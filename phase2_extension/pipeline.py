import sys
sys.dont_write_bytecode=True
from concurrent.futures import ThreadPoolExecutor
import argparse
import itertools
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess
import time
import traceback
import uuid
import base_pipeline as B
import scoring as S
import screen
import storage as A
H=Path(__file__).resolve().parent
C=A.read_json(H/'config.json')
B.HERE=H;B.CONFIG=C;B.ARRANGEMENTS=C['offsets']
AA=screen.AA

def dataset(family,arr):
    d=C['deletion_length'];anchor=family['anchor_0based'];focal=family['focal'];offset=C['offsets'][arr]
    deleted={focal[0]:set(range(anchor,anchor+d)),focal[1]:set(range(anchor+offset,anchor+offset+d))}
    original=family['rows'];rows={};ancestry={};events=[]
    for n,s in original.items():
        drop=deleted.get(n,set());assert all(s[i] in AA and family['core_mask'][i] for i in drop)
        rows[n]=''.join('-' if i in drop else aa for i,aa in enumerate(s))
        ancestry[n]=[i+1 for i,aa in enumerate(s) if aa!='-' and i not in drop]
        assert len(s.replace('-',''))-len(rows[n].replace('-',''))==(d if n in focal else 0)
        assert rows[n].replace('-','')==''.join(s[i-1] for i in ancestry[n])
    for n,drop in deleted.items():events.append({'sequence_id':n,'start':min(drop)+1,'end':max(drop)+1})
    meta={'ancestry':ancestry,'events':events,'region_window':C['window'],'core_mask':family['core_mask'],'focal':focal,
          'shared_deleted_columns':len(deleted[focal[0]]&deleted[focal[1]]),'union_deleted_columns':len(deleted[focal[0]]|deleted[focal[1]]),
          'arrangement':arr,'deletion_length':d,'offset':offset,'component':family['component']}
    assert meta['shared_deleted_columns']==max(0,d-offset)
    assert meta['union_deleted_columns']==2*d-meta['shared_deleted_columns']
    inputs={n:s.replace('-','') for n,s in rows.items()};S.validate_alignment(rows,inputs)
    return rows,inputs,meta
B.dataset=dataset

def read_case(root):
    rows=S.parse_fasta((root/'reference.fasta').read_text());meta=A.read_json(root/'mapping.json')
    evaluable={(n,i) for n,s in rows.items() for i,col in enumerate(meta['ancestry'][n],1) if meta['core_mask'][col-1] and s[col-1] in AA}
    truth={p for p in S.pairs(rows) if p[0] in evaluable and p[1] in evaluable}
    boundary=S.boundary_residues(meta)&evaluable
    return truth,boundary,evaluable,meta

def score(root,cap,fail=False):
    target=root/'scores';out=target/f'cap{cap}.json';mark=target/f'cap{cap}_done.json'
    deps={'reference':A.sha(root/'reference.fasta'),'mapping':A.sha(root/'mapping.json'),
          'baseline':A.sha(root/'alignments/cap0/aligned.fasta'),'alignment':A.sha(root/f'alignments/cap{cap}/aligned.fasta'),
          'scoring_source':A.sha(H/'scoring.py'),'pipeline_source':A.sha(H/'pipeline.py')}
    if A.complete(mark,deps):return A.read_json(out)
    if fail:raise RuntimeError('Injected failure after raw MAFFT completed')
    truth,boundary,evaluable,meta=read_case(root)
    pred=S.pairs(S.parse_fasta((root/f'alignments/cap{cap}/aligned.fasta').read_text()))
    baseline=pred if cap==0 else S.pairs(S.parse_fasta((root/'alignments/cap0/aligned.fasta').read_text()))
    pred_core={p for p in pred if p[0] in evaluable and p[1] in evaluable}
    base_core=pred_core if cap==0 else {p for p in baseline if p[0] in evaluable and p[1] in evaluable}
    acc=S.accuracy(pred_core,truth,boundary);comp=S.compare(base_core,pred_core,truth,boundary)
    affected=set(itertools.chain.from_iterable(pred^baseline));false=S.local_subset(pred_core-truth,boundary)
    distances=[abs(meta['ancestry'][a[0]][a[1]-1]-meta['ancestry'][b[0]][b[1]-1]) for a,b in false]
    details=A.read_json(root/f'alignments/cap{cap}/done.json')['details']
    result={'accuracy':acc,'comparison':comp,'residues_with_changed_partners':len(affected),
        'boundary_residues_with_changed_partners':len(affected&boundary),'boundary_residues':len(boundary),
        'all_pair_set_changed':pred!=baseline,'runtime_seconds':details['runtime_seconds'],'convergence_reported':details['convergence_reported'],
        'mean_false_local_partner_reference_displacement':statistics.mean(distances) if distances else None,
        'perfect_neighborhoods':sum(not S.local_subset(truth-pred_core,{r for r in boundary if r[0]==n}) for n in meta['focal']),
        'cap':cap,'baseline_missing_local':len(S.local_subset(truth-base_core,boundary))}
    A.json_write(out,result);A.marker(mark,[out],deps)
    return result

def identity(rows,mask,a=None,b=None):
    pairs=[(a,b)] if a is not None else itertools.combinations(rows,2);values=[]
    for a,b in pairs:
        cols=[i for i,m in enumerate(mask) if m and rows[a][i] in AA and rows[b][i] in AA]
        if cols:values.append(sum(rows[a][i]==rows[b][i] for i in cols)/len(cols))
    return statistics.mean(values) if values else None

def characteristics(name,fam):
    rows=fam['rows'];mask=fam['core_mask'];a=fam['anchor_0based'];d=C['deletion_length'];w=C['window']
    cols=range(a-w,a+2*d+w);conservation=[];entropy=[];gap_runs=[]
    for col in cols:
        values=[s[col] for s in rows.values() if s[col] in AA]
        if mask[col] and values:
            counts=[values.count(x)/len(values) for x in set(values)]
            conservation.append(max(counts));entropy.append(-sum(p*math.log(p) for p in counts)/math.log(20))
    for seq in rows.values():gap_runs.extend(len(m.group()) for m in re.finditer('-+',seq[a-w:a+2*d+w]))
    n=len(rows)*len(cols);local_symbols=[s[i] for s in rows.values() for i in cols]
    return {'family':name,'component':fam['component'],'reference_set':fam['reference_set'],'n_sequences':len(rows),
        'reference_columns':len(mask),'deletion_length':d,'window':w,'focal1':fam['focal'][0],'focal2':fam['focal'][1],
        'anchor_1based':a+1,'mean_core_pairwise_identity':identity(rows,mask),
        'focal_core_pairwise_identity':identity(rows,mask,*fam['focal']),
        'local_conservation':statistics.mean(conservation),'local_normalized_entropy':statistics.mean(entropy),
        'local_natural_gap_density':sum(x=='-' for x in local_symbols)/n,'local_gap_run_count':len(gap_runs),
        'local_gap_run_mean_length':statistics.mean(gap_runs) if gap_runs else 0,'local_gap_run_max_length':max(gap_runs,default=0),
        'local_ambiguous_fraction':sum(x not in AA and x!='-' for x in local_symbols)/n,
        'local_core_coverage':sum(mask[i] for i in cols)/len(cols),'family_core_coverage':sum(mask)/len(mask),
        'local_scorable_residues':sum(s[i] in AA and mask[i] for s in rows.values() for i in cols),
        'ambiguous_residues_family':sum(c not in AA and c!='-' for s in rows.values() for c in s),
        'region_start_1based':a-w+1,'region_end_1based':a+2*d+w}

def freeze():
    if (H/'experiment_manifest.json').exists():verify();print('Existing frozen manifest verified');return
    validation=A.read_json(H/'validation_latest.json');assert validation['passed']
    selected=A.read_json(H/'selected_families.json');assert len(selected)==C['primary_families']
    screen.csv_write(H/'pre_alignment_characteristics.csv',[characteristics(n,f) for n,f in selected.items()])
    A.json_write(H/'validation_before_execution.json',validation)
    src=[*H.glob('*.py'),*H.glob('*PLAN.md'),H/'config.json',H/'selected_families.json',H/'screening_sensitivity.csv',H/'screening_audit_218.csv',H/'screening_done.json',H/'pre_alignment_characteristics.csv']
    mafft=H/'mafft-win/mafft.bat';env=os.environ.copy()
    for key in ['TMP','TEMP','TMPDIR','MAFFT_TMPDIR']:env[key]=str(H/'tmp')
    proc=subprocess.run([str(mafft),'--version'],cwd=H,env=env,capture_output=True,timeout=60);assert proc.returncode==0
    version=(proc.stdout+proc.stderr).decode(errors='replace');assert '7.526' in version
    A.json_write(H/'mafft_version.json',{'command':[str(mafft),'--version'],'exit_code':proc.returncode,'output':version})
    A.json_write(H/'experiment_manifest.json',{'frozen_utc_epoch':time.time(),'config':C,'BAliBASE_sha256':screen.ARCHIVE_SHA,
        'source_url':'https://www.lbgi.fr/balibase/BalibaseDownload/BAliBASE_R1-5.tar.gz',
        'files':{p.relative_to(H).as_posix():A.sha(p) for p in src},
        'mafft_files':{p.relative_to(H).as_posix():A.sha(p) for p in (H/'mafft-win').rglob('*') if p.is_file()},
        'source_note':'base_pipeline.py and storage.py copied unchanged from validated Phase2. scoring.py extends accepted alphabet only; pipeline.py restricts evaluable residues explicitly. No original module imported or changed.',
        'license_note':'Original public archive README has no explicit redistribution license. Preserve attribution and verify reuse terms.'})
    print('Source/configuration/31 independent family selection frozen before MAFFT')

def verify():
    obj=A.read_json(H/'experiment_manifest.json')
    for p,digest in obj['files'].items():
        if A.sha(H/p)!=digest:raise RuntimeError('Frozen file changed: '+p)
    assert A.sha(H/'source/BAliBASE_R1-5.tar.gz')==obj['BAliBASE_sha256']

def baseline_table(selected):
    rows=[]
    for n in selected:
        for arr in C['offsets']:
            s=score(H/'cases'/n/arr,0);local=s['accuracy']['local'];assert local['reference_pairs']>0
            rows.append({'family':n,'arrangement':arr,'local_sp_recall':local['sp_recall'],'missed_pairs':local['missed_pairs'],
                'correct_pairs':local['correct_pairs'],'local_reference_pairs':local['reference_pairs'],'evaluable_false_pairs':local['false_pairs'],
                'global_core_sp_recall':s['accuracy']['global']['sp_recall'],'perfect_neighborhoods':s['perfect_neighborhoods'],
                'mean_false_partner_reference_displacement':s['mean_false_local_partner_reference_displacement'],
                'difficulty':'perfect' if local['missed_pairs']==0 else 'near-perfect' if local['sp_recall']>=.95 else 'difficult'})
    screen.csv_write(H/'baseline_characteristics.csv',rows)
    A.json_write(H/'baseline_characteristics_done.json',{'sha256':A.sha(H/'baseline_characteristics.csv'),'utc_epoch':time.time(),
        'before_refinement':not any((H/'cases').glob('*/*/alignments/cap2/done.json'))})

def run():
    verify();selected=A.read_json(H/'selected_families.json');validation=A.read_json(H/'validation_before_execution.json');assert validation['passed']
    counts={'completed_now':0,'skipped':0,'failed':0};invocation=uuid.uuid4().hex
    with A.writer_lock(H):
        roots=[B.prepare_case(n,f,arr) for n,f in selected.items() for arr in C['offsets']]
        for cap in C['caps']:
            def work(root):
                try:
                    called=B.align(root,cap);score(root,cap)
                    result={'family':root.parent.name,'arrangement':root.name,'cap':cap,'state':'completed_now' if called else 'skipped'}
                except Exception:
                    result={'family':root.parent.name,'arrangement':root.name,'cap':cap,'state':'failed','traceback':traceback.format_exc()}
                A.json_write(H/'logs'/invocation/(uuid.uuid4().hex+'.json'),result);return result
            with ThreadPoolExecutor(max_workers=C['workers']) as pool:
                for result in pool.map(work,roots):counts[result['state']]+=1;print(json.dumps(result),flush=True)
            if counts['failed']:raise RuntimeError('Failed stages logged and successful alignments preserved. Resume same command.')
            if cap==0 and not (H/'baseline_characteristics_done.json').exists():baseline_table(selected)
        A.json_write(H/'invocations'/f'{invocation}.json',counts);A.json_write(H/'last_run_status.json',counts);print(counts)

if __name__=='__main__':
    action=sys.argv[1]
    {'freeze':freeze,'run':run}[action]()
