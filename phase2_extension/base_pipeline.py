"""Separately frozen BAliBASE controlled-deletion experiment. No earlier writes."""
import sys
sys.dont_write_bytecode = True
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import io
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile
import time
import traceback
import uuid
import xml.etree.ElementTree as ET
import scoring as S
import storage as A

HERE = Path(__file__).resolve().parent
AA = set('ACDEFGHIKLMNPQRSTVWY')
ARRANGEMENTS = {'none':50,'partial':6,'high':1}
CONFIG = {'experiment':'phase2_balibase_v1','sequence_count':[4,60],'max_columns':2000,
          'deletion_length':12,'offsets':ARRANGEMENTS,'window':10,
          'caps':[0,2,10],'thread':1,'retree':2,'bootstrap_seed':3553002,
          'bootstrap_resamples':10000,'difficulty_threshold':0.95,
          'family_selection':'shared identifier or exact ungapped sequence connected components; first eligible case',
          'primary':'core-local normalized net gain cap2','timeout_seconds':600}

def csv_write(path, rows):
    if not rows: raise ValueError('No rows to write')
    out=io.StringIO(newline=''); w=csv.DictWriter(out,fieldnames=list(rows[0])); w.writeheader();w.writerows(rows)
    A.atomic(path,out.getvalue())

def fasta(rows):
    return ''.join('>'+n+'\n'+s+'\n' for n,s in rows.items())

def read_xml(data):
    root=ET.fromstring(data)
    rows={}
    for seq in root.findall('.//sequence'):
        name=seq.findtext('seq-name').strip()
        if name in rows: raise ValueError('Duplicate IDs')
        rows[name]=''.join(seq.findtext('seq-data').split()).upper().replace('.','-')
    if not rows or len({len(s) for s in rows.values()})!=1: raise ValueError('Unequal or empty rows')
    if any(set(s)-AA-{'-'} for s in rows.values()): raise ValueError('Non-standard protein alphabet')
    scores=[s for s in root.findall('.//column-score') if s.findtext('colsco-name')=='coreblock']
    if len(scores)!=1: raise ValueError('Missing or ambiguous core mask')
    raw=[int(s) for s in scores[0].findtext('colsco-data').split()]
    if len(raw)!=len(next(iter(rows.values()))) or set(raw)-{-1,0,1}: raise ValueError('Invalid core mask')
    mask=[int(v==1) for v in raw]
    return rows,mask

def choose_site(rows,mask):
    width=len(mask); best=None
    for a,b in itertools.combinations(sorted(rows),2):
        run=0
        for col in range(width):
            run=run+1 if mask[col] and rows[a][col]!='-' and rows[b][col]!='-' else 0
            if run>=82:
                # 0-based anchor; 82-column window begins 10 before it.
                anchor=col-71
                rank=(abs((anchor+30.5)-(width-1)/2),a,b,anchor)
                if best is None or rank<best: best=rank
    return None if best is None else {'focal':[best[1],best[2]],'anchor_0based':best[3]}

def identity(rows,mask):
    values=[]
    for a,b in itertools.combinations(rows,2):
        cols=[i for i,c in enumerate(mask) if c and rows[a][i]!='-' and rows[b][i]!='-']
        if cols: values.append(sum(rows[a][i]==rows[b][i] for i in cols)/len(cols))
    return {'mean_core_pairwise_identity':sum(values)/len(values) if values else None,
            'min_core_pairwise_identity':min(values) if values else None,
            'max_core_pairwise_identity':max(values) if values else None}

def prepare():
    manifest=HERE/'experiment_manifest.json'
    if manifest.exists():
        verify_frozen(); print('Frozen screening reused'); return
    archive=HERE/'downloads/BAliBASE_R1-5.tar.gz'
    families={}; metadata=[]; allrows={}
    with tarfile.open(archive) as tar:
        for item in tar:
            if not item.isfile(): continue
            dest=(HERE/'downloads'/item.name).resolve()
            if not dest.is_relative_to(HERE/'downloads'): raise ValueError('Unsafe archive path')
            A.atomic(dest,tar.extractfile(item).read())
        for path in sorted((HERE/'downloads/bb3_release').glob('RV*/BB*.xml')):
            if path.stem.startswith('BBS'): continue
            row={'case':path.stem,'set':path.parent.name,'source':str(path.relative_to(HERE)),
                 'eligible':False,'selected':False,'exclusion':'','n_sequences':None,'columns':None,
                 'core_columns':None,'mean_core_pairwise_identity':None,'min_core_pairwise_identity':None,
                 'max_core_pairwise_identity':None,'component':None,'focal1':None,'focal2':None,'anchor_1based':None}
            try:
                rows,mask=read_xml(path.read_bytes()); allrows[path.stem]=rows
                row.update(n_sequences=len(rows),columns=len(mask),core_columns=sum(mask),**identity(rows,mask))
                if not 4<=len(rows)<=60: raise ValueError('Sequence count outside 4–60')
                if len(mask)>2000: raise ValueError('More than 2000 columns')
                site=choose_site(rows,mask)
                if site is None: raise ValueError('No 82-column core gap-free focal interval')
                row.update(eligible=True,focal1=site['focal'][0],focal2=site['focal'][1],anchor_1based=site['anchor_0based']+1)
                families[path.stem]={'rows':rows,'core_mask':mask,**site,'source':row['source']}
            except ValueError as exc: row['exclusion']=str(exc)
            metadata.append(row)
    # Conservative connected components, including ineligible cases as bridges.
    parent={n:n for n in allrows}
    def find(x):
        while parent[x]!=x: x=parent[x]
        return x
    seen={}
    for n,rows in allrows.items():
        for name,seq in rows.items():
            for token in [('id',name),('seq',seq.replace('-',''))]:
                if token in seen:
                    a,b=find(n),find(seen[token]); parent[max(a,b)]=min(a,b)
                else: seen[token]=n
    groups={}
    for n in families: groups.setdefault(find(n),[]).append(n)
    selected=sorted(min(v) for v in groups.values())
    for row in metadata:
        if row['case'] in parent: row['component']=find(row['case'])
        row['selected']=row['case'] in selected
        if row['eligible'] and not row['selected']: row['exclusion']='Linked family component; first eligible representative used'
    csv_write(HERE/'family_screening.csv',metadata)
    A.json_write(HERE/'selected_families.json',{n:families[n] for n in selected})
    A.json_write(HERE/'config.json',CONFIG)
    source_files=sorted([*HERE.glob('*.py'),HERE/'PLAN.md',HERE/'config.json',HERE/'family_screening.csv',HERE/'selected_families.json'])
    frozen={str(p.relative_to(HERE)):A.sha(p) for p in source_files}
    files={str(p.relative_to(HERE)):A.sha(p) for p in (HERE/'mafft-win').rglob('*') if p.is_file()}
    provenance={'source_url':'https://www.lbgi.fr/balibase/BalibaseDownload/BAliBASE_R1-5.tar.gz',
                'index_url':'https://www.lbgi.fr/balibase/BalibaseDownload/',
                'version':'3.0, Dec 2004 (included README)', 'archive_sha256':A.sha(archive),
                'paper':'https://pubmed.ncbi.nlm.nih.gov/16044462/',
                'license':'No explicit data redistribution license found in archive README. Public download is not an explicit license. Cite BAliBASE; verify reuse permissions before redistribution/publication.',
                'full_length_cases':len(metadata),'truncated_xml':len(list((HERE/'downloads/bb3_release').glob('RV*/BBS*.xml'))),
                'eligible_cases':len(families),'selected_independent_components':len(selected),
                'linked_components_total':len({find(n) for n in allrows}),
                'mafft_copied_from':'Original official MAFFT Windows 7.526 installation; original files only read',
                'mafft_files_sha256':files,'formats':['XML with coreblock mask','MSF','TFA ungapped inputs'],
                'frozen_files':frozen,'frozen_before_mafft_utc_epoch':time.time()}
    A.json_write(manifest,provenance)
    print(json.dumps({k:v for k,v in provenance.items() if k not in ['frozen_files','mafft_files_sha256']},indent=2))

def verify_frozen():
    obj=A.read_json(HERE/'experiment_manifest.json')
    for path,h in obj['frozen_files'].items():
        if A.sha(HERE/path)!=h: raise RuntimeError('Frozen source/config changed: '+path)
    if A.sha(HERE/'downloads/BAliBASE_R1-5.tar.gz')!=obj['archive_sha256']: raise RuntimeError('Archive changed')

def dataset(family,arrangement):
    original=family['rows']; anchor=family['anchor_0based']; focal=family['focal']; width=len(family['core_mask'])
    deleted={focal[0]:set(range(anchor,anchor+12)),focal[1]:set(range(anchor+ARRANGEMENTS[arrangement],anchor+ARRANGEMENTS[arrangement]+12))}
    overlap=len(deleted[focal[0]]&deleted[focal[1]])
    assert overlap=={'none':0,'partial':6,'high':11}[arrangement]
    rows={}; ancestry={}; events=[]
    for n,seq in original.items():
        drop=deleted.get(n,set()); rows[n]=''.join('-' if i in drop else aa for i,aa in enumerate(seq))
        assert all(seq[i]!='-' for i in drop)
        ancestry[n]=[i+1 for i,aa in enumerate(seq) if aa!='-' and i not in drop]
        assert len(seq.replace('-',''))-len(rows[n].replace('-',''))==(12 if n in focal else 0)
        assert all(rows[n][c-1]==seq[c-1] for c in ancestry[n])
    for n,drop in deleted.items():
        assert len(drop)==12 and min(drop)>=10 and max(drop)<width-10
        events.append({'sequence_id':n,'start':min(drop)+1,'end':max(drop)+1,'reference_columns':sorted(i+1 for i in drop)})
    meta={'ancestry':ancestry,'events':events,'region_window':10,'core_mask':family['core_mask'],
          'shared_deleted_columns':overlap,'union_deleted_columns':24-overlap,
          'focal':focal,'original_column_count':width,'arrangement':arrangement}
    inputs={n:s.replace('-','') for n,s in rows.items()}
    S.validate_alignment(rows,inputs)
    return rows,inputs,meta

def truth_pairs(rows,meta):
    # Replacing non-core columns by all-gap columns would renumber residues.
    # Instead filter true pairs through immutable survivor-to-original-column maps.
    truth=S.pairs(rows)
    return {p for p in truth if meta['core_mask'][meta['ancestry'][p[0][0]][p[0][1]-1]-1]}

def prepare_case(name,fam,arr):
    root=HERE/'cases'/name/arr
    deps={'family':A.digest(A.canonical(fam)),'arrangement':arr,'config':A.digest(A.canonical(CONFIG))}
    if not A.complete(root/'dataset_done.json',deps):
        truth,inputs,meta=dataset(fam,arr)
        A.quarantine(root,[root/'input.fasta',root/'reference.fasta',root/'mapping.json',root/'dataset_done.json'])
        A.atomic(root/'input.fasta',fasta(inputs)); A.atomic(root/'reference.fasta',fasta(truth)); A.json_write(root/'mapping.json',meta)
        A.marker(root/'dataset_done.json',[root/'input.fasta',root/'reference.fasta',root/'mapping.json'],deps)
    return root

def align(root,cap,invoke=None):
    inp=S.parse_fasta((root/'input.fasta').read_text()); target=root/'alignments'/f'cap{cap}'
    target.mkdir(parents=True,exist_ok=True)
    deps={'input_sha256':A.sha(root/'input.fasta'),'flags':['--amino','--thread','1','--retree','2','--maxiterate',str(cap)],'mafft':'7.526'}
    validator=lambda:S.validate_alignment(S.parse_fasta((target/'aligned.fasta').read_text()),inp)
    if A.complete(target/'done.json',deps,validator): return False
    A.quarantine(target,[target/'aligned.fasta',target/'stderr.log',target/'done.json'])
    cmd=[str(HERE/'mafft-win/mafft.bat'),*deps['flags'],str(root/'input.fasta')]
    started=time.time()
    try:
        if invoke: out,err,code=invoke(cmd)
        else:
            env=os.environ.copy()
            for key in ['TMP','TEMP','TMPDIR','MAFFT_TMPDIR']: env[key]=str(HERE/'tmp')
            proc=subprocess.run(cmd,cwd=HERE,env=env,capture_output=True,timeout=CONFIG['timeout_seconds'])
            out,err,code=proc.stdout,proc.stderr,proc.returncode
        A.atomic(target/'stderr.log',err)
        if code: raise RuntimeError('MAFFT exit '+str(code))
        rows=S.parse_fasta(out.decode()); S.validate_alignment(rows,inp)
        A.atomic(target/'aligned.fasta',out)
        text=err.decode(errors='replace')
        A.marker(target/'done.json',[target/'aligned.fasta',target/'stderr.log'],deps,
                 {'command':cmd,'runtime_seconds':time.time()-started,'exit_code':code,
                  'convergence_reported':bool(re.search(r'converg',text,re.I)),
                  'diagnostic_lines':[line for line in text.splitlines() if re.search(r'converg|iterate|cycle|STEP|identical',line,re.I)]})
        return True
    except Exception:
        A.json_write(target/('failure_'+uuid.uuid4().hex+'.json'),{'command':cmd,'traceback':traceback.format_exc()})
        raise

def load_case(root):
    ref=S.parse_fasta((root/'reference.fasta').read_text()); meta=A.read_json(root/'mapping.json')
    truth=truth_pairs(ref,meta); boundary=S.boundary_residues(meta)
    return truth,boundary,meta

def score(root,cap,fail=False):
    target=root/'scores'; path=target/f'cap{cap}.json'; mark=target/f'cap{cap}_done.json'
    deps={'alignment':A.sha(root/'alignments'/f'cap{cap}'/'aligned.fasta'),'mapping':A.sha(root/'mapping.json'),
          'reference':A.sha(root/'reference.fasta'),'scorer':A.sha(HERE/'scoring.py'),'phase2_code':A.sha(HERE/'experiment.py')}
    if A.complete(mark,deps): return A.read_json(path)
    if fail: raise RuntimeError('Injected scoring failure after saved alignment')
    truth,boundary,meta=load_case(root)
    pred=S.pairs(S.parse_fasta((root/'alignments'/f'cap{cap}'/'aligned.fasta').read_text()))
    baseline=S.pairs(S.parse_fasta((root/'alignments/cap0/aligned.fasta').read_text()))
    # Restrict evaluable predictions to residues individually annotated in core.
    # Wrong cross-column pairings between those residues still count as false.
    core_res={(n,i) for n,anc in meta['ancestry'].items() for i,c in enumerate(anc,1) if meta['core_mask'][c-1]}
    pred_core={p for p in pred if p[0] in core_res and p[1] in core_res}
    base_core={p for p in baseline if p[0] in core_res and p[1] in core_res}
    acc=S.accuracy(pred_core,truth,boundary); comp=S.compare(base_core,pred_core,truth,boundary)
    affected=set(itertools.chain.from_iterable(pred^baseline))
    local_affected=affected&boundary
    done=A.read_json(root/'alignments'/f'cap{cap}'/'done.json')['details']
    result={'accuracy':acc,'comparison':comp,'all_predicted_pair_set_changed':pred!=baseline,
            'residues_with_changed_partners':len(affected),'boundary_residues_with_changed_partners':len(local_affected),
            'boundary_residues':len(boundary),'runtime_seconds':done['runtime_seconds'],
            'convergence_reported':done['convergence_reported'],'cap':cap,
            'baseline_core_local_missed':len(S.local_subset(truth-base_core,boundary)),
            'remaining_core_local_missed':acc['local']['missed_pairs'],
            'focal_neighborhoods':{n:{'reference_pairs':len(S.local_subset(truth,{r for r in boundary if r[0]==n})),
                'missed_pairs':len(S.local_subset(truth-pred_core,{r for r in boundary if r[0]==n}))} for n in meta['focal']}}
    A.json_write(path,result); A.marker(mark,[path],deps)
    return result

def baseline_table(selected):
    rows=[]
    for name in selected:
        for arr in ARRANGEMENTS:
            root=HERE/'cases'/name/arr; s=score(root,0); a=s['accuracy']['local']
            recall=a['sp_recall']
            if recall is None: raise ValueError('Undefined local reference denominator')
            rows.append({'family':name,'arrangement':arr,'local_sp_recall':recall,'local_missed_pairs':a['missed_pairs'],
                         'local_reference_pairs':a['reference_pairs'],'global_core_sp_recall':s['accuracy']['global']['sp_recall'],
                         'perfect':a['missed_pairs']==0,'difficulty':'perfect' if recall==1 else 'near-perfect' if recall>=.95 else 'difficult',
                         'perfect_neighborhoods':sum(v['missed_pairs']==0 for v in s['focal_neighborhoods'].values()),
                         'total_neighborhoods':len(s['focal_neighborhoods']),
                         **identity(selected[name]['rows'],selected[name]['core_mask'])})
    csv_write(HERE/'baseline_difficulty.csv',rows)
    A.json_write(HERE/'baseline_difficulty_done.json',{'sha256':A.sha(HERE/'baseline_difficulty.csv'),'utc_epoch':time.time(),
                 'all_cap0_finished_before_refinement':not any((HERE/'cases').glob('*/*/alignments/cap2/done.json'))})

def run(workers=4):
    verify_frozen()
    if not (HERE/'validation_results.json').exists() or not A.read_json(HERE/'validation_results.json')['passed']:
        raise RuntimeError('Pre-execution validation required')
    selected=A.read_json(HERE/'selected_families.json')
    with A.writer_lock(HERE):
        roots=[prepare_case(n,f,a) for n,f in selected.items() for a in ARRANGEMENTS]
        counts={'completed_now':0,'skipped':0,'failed':0}
        for cap in [0,2,10]:
            def work(root):
                try:
                    executed=align(root,cap); score(root,cap)
                    return {'case':str(root.relative_to(HERE)),'cap':cap,'state':'completed_now' if executed else 'skipped'}
                except Exception:
                    return {'case':str(root.relative_to(HERE)),'cap':cap,'state':'failed','traceback':traceback.format_exc()}
            with ThreadPoolExecutor(max_workers=workers) as pool:
                for result in pool.map(work,roots):
                    counts[result['state']]+=1
                    log=HERE/'logs'/f'{uuid.uuid4().hex}.json'; A.json_write(log,result)
                    print(json.dumps(result),flush=True)
            if counts['failed']: raise RuntimeError('Failures saved; resume same command after resolving cause')
            if cap==0 and not (HERE/'baseline_difficulty_done.json').exists(): baseline_table(selected)
        A.json_write(HERE/'last_run_status.json',counts)
        print(counts)

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('action',choices=['prepare','run']); parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    if args.action=='prepare': prepare()
    else: run(args.workers)
