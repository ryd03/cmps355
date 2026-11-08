"""Source-only screening; never imports or invokes MAFFT or previous analyses."""
import sys
sys.dont_write_bytecode=True
import csv
import io
import itertools
from pathlib import Path
import tarfile
import time
import xml.etree.ElementTree as X
import storage as A
H=Path(__file__).resolve().parent
AA=set('ACDEFGHIKLMNPQRSTVWY')
ARCHIVE_SHA='5378cf4e2e0741f788a0175d4fcb78c54387102b243e8de16c5e2cafed928190'
SOURCE=H.parent/'phase2_balibase_real_proteins/downloads/BAliBASE_R1-5.tar.gz'

def csv_write(path,rows):
    if not rows:raise ValueError('Empty table')
    out=io.StringIO(newline='');w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows);A.atomic(path,out.getvalue())

def load_sources():
    archive=H/'source/BAliBASE_R1-5.tar.gz'
    if not archive.exists():A.atomic(archive,SOURCE.read_bytes())
    assert A.sha(archive)==ARCHIVE_SHA
    result={}
    with tarfile.open(archive) as tar:
        for member in tar:
            if not member.isfile():continue
            name=Path(member.name).name
            if member.name.endswith('/README') or (name.endswith('.xml') and name.startswith('BB') and not name.startswith('BBS')):
                dest=(H/'source'/member.name).resolve();assert dest.is_relative_to(H/'source')
                if not dest.exists():A.atomic(dest,tar.extractfile(member).read())
            if not(name.endswith('.xml') and name.startswith('BB') and not name.startswith('BBS')):continue
            root=X.fromstring(tar.extractfile(member).read());rows={}
            for s in root.findall('.//sequence'):
                n=s.findtext('seq-name').strip();assert n not in rows
                rows[n]=''.join(s.findtext('seq-data').split()).upper().replace('.','-')
            masks=[c for c in root.findall('.//column-score') if c.findtext('colsco-name')=='coreblock'];assert len(masks)==1
            raw=list(map(int,masks[0].findtext('colsco-data').split()));assert set(raw)<={-1,0,1}
            assert len({len(s) for s in rows.values()})==1 and len(raw)==len(next(iter(rows.values())))
            result[Path(name).stem]={'rows':rows,'core_mask':[int(c==1) for c in raw],'reference_set':Path(member.name).parent.name,
                'xml_sha256':A.digest(tar.extractfile(member).read()),'source_xml':'source/'+member.name}
    assert len(result)==218
    return result

def components(families):
    parents={n:n for n in families};seen={}
    def find(x):
        while parents[x]!=x:x=parents[x]
        return x
    for n,f in sorted(families.items()):
        for seqid,seq in f['rows'].items():
            for token in [('id',seqid),('seq',seq.replace('-',''))]:
                if token in seen:
                    a,b=find(n),find(seen[token]);parents[max(a,b)]=min(a,b)
                else:seen[token]=n
    return {n:find(n) for n in families}

def prefix(values):
    p=[0]
    for v in values:p.append(p[-1]+v)
    return p

def site(f,d,w,reliability='targeted',symbols='retained',old=False):
    rows=f['rows'];mask=f['core_mask'];width=len(mask)
    if not 4<=len(rows)<=60 or width>2000:return None
    if any(set(s)-AA-set('XBZ-') for s in rows.values()):return None
    if symbols=='strict' and any(set(s)-AA-{'-'} for s in rows.values()):return None
    span=62 if old else 2*d
    if width<span+2*w:return None
    core=prefix(mask);std={n:prefix([int(c in AA) for c in s]) for n,s in rows.items()}
    scorable={n:prefix([int(c in AA and mask[i]) for i,c in enumerate(s)]) for n,s in rows.items()}
    names=sorted(rows);offsets=[50,6,1] if old else [d,d//2,1]
    candidates=sorted(range(w,width-span-w+1),key=lambda a:(abs(a+(span-1)/2-(width-1)/2),a))
    def okay(n,start,length):return std[n][start+length]-std[n][start]==length
    def flanks(n,start):
        p=scorable[n]
        return p[start]-p[start-w]>=3 and p[start+d+w]-p[start+d]>=3
    for anchor in candidates:
        if old or reliability=='continuous':
            if core[anchor+span+w]-core[anchor-w]!=span+2*w:continue
            eligible=[n for n in names if okay(n,anchor-w,span+2*w)]
            if len(eligible)>=2:return {'anchor_0based':anchor,'focal':eligible[:2]}
        else:
            if core[anchor+span]-core[anchor]!=span:continue
            first=[n for n in names if okay(n,anchor,d) and flanks(n,anchor)]
            second=[n for n in names if all(okay(n,anchor+offset,d) and flanks(n,anchor+offset) for offset in offsets)]
            for a in first:
                b=next((b for b in second if b!=a),None)
                if b:return {'anchor_0based':anchor,'focal':[a,b]}
    return None

def screen():
    if (H/'screening_done.json').exists():
        m=A.read_json(H/'screening_done.json')
        for p,digest in m['files'].items():assert A.sha(H/p)==digest
        print('Outcome-blind screening already complete; reused');return
    families=load_sources();comp=components(families);pilot_roots={comp['BB12007'],comp['BB12031']}
    protocols=[{'length':12,'window':10,'reliability':'original82','symbols':'strict'}]
    protocols +=[{'length':d,'window':w,'reliability':r,'symbols':s} for d in [12,10,8,6] for w in [10,5] for r in ['continuous','targeted'] for s in ['strict','retained']]
    sensitivity=[];sites={};audits=[]
    for p in protocols:
        key=f"d{p['length']}_w{p['window']}_{p['reliability']}_{p['symbols']}"; eligible={}
        for n,f in families.items():
            chosen=site(f,p['length'],p['window'],p['reliability'],p['symbols'],old=p['reliability']=='original82')
            if chosen:eligible[n]=chosen
        roots={comp[n] for n in eligible};new=roots-pilot_roots;sites[key]=eligible
        sensitivity.append({'protocol':key,**p,'eligible_cases':len(eligible),'independent_components':len(roots),
            'pilot_components_qualifying':len(roots&pilot_roots),'new_independent_components':len(new),
            'exact_overlaps':'0,6,11' if p['reliability']=='original82' else f"0,{p['length']//2},{p['length']-1}",
            'none_start_offset':50 if p['reliability']=='original82' else p['length']})
        print(sensitivity[-1],flush=True)
    primary=[r for r in sensitivity if r['reliability']=='targeted' and r['symbols']=='retained' and r['window']==10]
    qualifying=[r for r in primary if r['new_independent_components']>=20]
    choice=max(qualifying,key=lambda r:r['length']) if qualifying else max(primary,key=lambda r:(r['new_independent_components'],r['length']))
    eligible=sites[choice['protocol']];representatives={}
    for n in sorted(eligible):
        if comp[n] not in pilot_roots:representatives.setdefault(comp[n],n)
    selected={n:{**families[n],**eligible[n],'component':root} for root,n in sorted(representatives.items())}
    for n,f in families.items():
        rows=f['rows'];width=len(f['core_mask']);amb=sum(c not in AA and c!='-' for s in rows.values() for c in s)
        original=site(f,12,10,'continuous','strict',old=True)
        # Audit all criteria, not merely the first failed one.
        geometric_original=site(f,12,10,'continuous','retained',old=True)
        audits.append({'case':n,'component':comp[n],'n_sequences':len(rows),'columns':width,'core_columns':sum(f['core_mask']),
            'ambiguous_residues':amb,'original_excluded_for_ambiguity':amb>0,'outside_sequence_count_limit':not 4<=len(rows)<=60,
            'outside_column_limit':width>2000,'original82_eligible':original is not None,
            'original82_geometry_with_ambiguity_retained':geometric_original is not None,
            'primary_extension_site_valid':n in eligible,'pilot_component_excluded_from_primary':comp[n] in pilot_roots,
            'primary_selected':n in selected,'reference_set':f['reference_set'],
            'primary_exclusion':'original pilot component' if comp[n] in pilot_roots else 'no valid primary site/size' if n not in eligible else 'linked eligible case; one representative' if n not in selected else ''})
    config={'experiment':'phase2_extension_v1','deletion_length':choice['length'],'window':10,'minimum_core_residues_per_flank':3,
        'offsets':{'none':choice['length'],'partial':choice['length']//2,'high':1},'caps':[0,2,10],'thread':1,'retree':2,
        'bootstrap_resamples':10000,'bootstrap_seed':3553202,'difficulty_threshold':.95,'workers':2,'timeout_seconds':600,
        'protocol':choice['protocol'],'primary':'cap2 normalized core-local net gain','pilot_components_excluded':sorted(pilot_roots),
        'primary_families':len(selected),'one_site_per_component':True,'ambiguous_symbols':'X B Z retained; scoring excludes either ambiguous endpoint',
        'regression':'exploratory univariate baseline-adjusted descriptive model only if >=30 families; at most2 predictors plus intercept; no p-values',
        'selection_rule':'longest d among12,10,8,6 with >=20 new components; else most new components, ties longer'}
    csv_write(H/'screening_sensitivity.csv',sensitivity);csv_write(H/'screening_audit_218.csv',audits)
    A.json_write(H/'selected_families.json',selected);A.json_write(H/'config.json',config)
    A.json_write(H/'family_components.json',{'components':comp,'total_components':len(set(comp.values())),'pilot_roots':sorted(pilot_roots)})
    files=[H/'screening_sensitivity.csv',H/'screening_audit_218.csv',H/'selected_families.json',H/'config.json',H/'family_components.json',H/'SCREENING_PLAN.md',H/'screen.py']
    A.json_write(H/'screening_done.json',{'files':{p.name:A.sha(p) for p in files},'utc_epoch':time.time(),'chosen_protocol':choice,
        'mafft_executed_before_screening':False})
    print('PRIMARY CONFIG',config,flush=True)

if __name__=='__main__':screen()
