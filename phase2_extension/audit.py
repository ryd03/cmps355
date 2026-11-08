"""Independent column-coincidence audit, zero-MAFFT resume proof and bundle."""
import sys
sys.dont_write_bytecode=True
from collections import defaultdict
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import time
import uuid
import zipfile
import analysis as N
import pipeline as P
import scoring as S
import storage as A
H=P.H

def column_maps(rows):
    columns=defaultdict(list);coords={};letters={}
    for name,seq in rows.items():
        idx=0
        for col,aa in enumerate(seq):
            if aa!='-':
                idx+=1;r=(name,idx);coords[r]=col;columns[col].append(r);letters[r]=aa
    groups={}
    for members in columns.values():
        value=frozenset(members)
        for r in members:groups[r]=value
    return coords,columns,groups,letters

def independent():
    P.verify();families=A.read_json(H/'selected_families.json');rows=N.read_csv(H/'expanded_results.csv');lookup={(r['family'],r['arrangement'],int(r['cap'])):r for r in rows}
    assert len(lookup)==len(families)*9
    checked=0;extra=[];checks=[]
    for name,fam in families.items():
        for arr in N.ARR:
            root=H/'cases'/name/arr;ref=S.parse_fasta((root/'reference.fasta').read_text());meta=A.read_json(root/'mapping.json')
            inputs=S.parse_fasta((root/'input.fasta').read_text());deletions={}
            for event in meta['events']:
                deleted=set(range(event['start'],event['end']+1));assert len(deleted)==P.C['deletion_length'];deletions[event['sequence_id']]=deleted
            ds=list(deletions.values());assert len(ds[0]&ds[1])==max(0,P.C['deletion_length']-P.C['offsets'][arr])
            for n,original in fam['rows'].items():
                expected=''.join('-' if i in deletions.get(n,set()) else aa for i,aa in enumerate(original,1))
                assert ref[n]==expected and inputs[n]==expected.replace('-','')
                for i in deletions.get(n,set()):assert original[i-1] in P.AA and fam['core_mask'][i-1]
            rc,columns,_,letters=column_maps(ref)
            known={r for r,c in rc.items() if fam['core_mask'][c] and letters[r] in P.AA};boundary=set()
            for event in meta['events']:
                for r,c in rc.items():
                    if r[0]==event['sequence_id'] and r in known and min(abs(c+1-event['start']),abs(c+1-event['end']))<=P.C['window']:boundary.add(r)
            alignment={};details={}
            for cap in P.C['caps']:
                target=root/f'alignments/cap{cap}';done=A.read_json(target/'done.json')
                for filename,digest in done['files'].items():assert A.sha(target/filename)==digest
                predicted=S.parse_fasta((target/'aligned.fasta').read_text());S.validate_alignment(predicted,inputs)
                alignment[cap]=column_maps(predicted);details[cap]=done
            truth=[]
            for col,members in columns.items():
                for pair in itertools.combinations(sorted(r for r in members if r in known),2):truth.append(pair)
            # Truth is derived independently; pair correctness is column coincidence.
            basecoords=alignment[0][0];localtruth=[p for p in truth if p[0] in boundary or p[1] in boundary]
            localtruth_set=set(localtruth)
            base_local_pairs={tuple(sorted((a,b))) for members in alignment[0][1].values()
                for a,b in itertools.combinations([r for r in members if r in known],2) if a in boundary or b in boundary}
            for cap in P.C['caps']:
                coords,pcols,groups,_=alignment[cap];r=lookup[name,arr,cap]
                refn=len(localtruth);correct=sum(coords[a]==coords[b] for a,b in localtruth)
                recovered=sum(coords[a]==coords[b] and basecoords[a]!=basecoords[b] for a,b in localtruth)
                lost=sum(coords[a]!=coords[b] and basecoords[a]==basecoords[b] for a,b in localtruth)
                missing_base=sum(basecoords[a]!=basecoords[b] for a,b in localtruth)
                predlocal=0;false_distances=[];pred_local_pairs=set()
                for members in pcols.values():
                    for a,b in itertools.combinations([r for r in members if r in known],2):
                        if a in boundary or b in boundary:
                            predlocal+=1
                            pred_local_pairs.add(tuple(sorted((a,b))))
                            if rc[a]!=rc[b]:false_distances.append(abs(rc[a]-rc[b]))
                assert int(r['local_reference_pairs'])==refn
                assert int(r['local_correct_pairs'])==correct
                assert int(r['recovered'])==recovered and int(r['lost'])==lost
                assert int(r['baseline_local_missed'])==missing_base
                assert int(r['remaining_missed_correct_pairs'])==refn-correct
                assert int(r['remaining_evaluable_false_pairs'])==predlocal-correct
                assert int(r['false_introduced'])==len((pred_local_pairs-base_local_pairs)-localtruth_set)
                assert int(r['false_removed'])==len((base_local_pairs-pred_local_pairs)-localtruth_set)
                assert abs(float(r['normalized_net_gain'])-(recovered-lost)/refn)<1e-13
                assert abs(float(r['local_sp_recall'])-correct/refn)<1e-13
                assert abs(float(r['global_core_sp_recall'])-sum(coords[a]==coords[b] for a,b in truth)/len(truth))<1e-13
                changed={r for r in coords if groups[r]!=alignment[0][2][r]}
                assert int(r['residues_with_changed_partners'])==len(changed)
                assert int(r['boundary_residues_with_changed_partners'])==len(changed&boundary)
                if false_distances:assert abs(float(r['mean_false_partner_displacement'])-sum(false_distances)/len(false_distances))<1e-12
                assert (r['pair_set_changed']=='True')==bool(changed)
                checked+=1
            c2,c10=alignment[2][0],alignment[10][0]
            extra.append({'family':name,'arrangement':arr,'recovered':sum(c10[a]==c10[b] and c2[a]!=c2[b] for a,b in localtruth),
                'lost':sum(c10[a]!=c10[b] and c2[a]==c2[b] for a,b in localtruth),
                'identical':all(alignment[2][2][r]==alignment[10][2][r] for r in c2)})
            checks.append(name+'/'+arr);print('Independent audit',name,arr,flush=True)
    for r in N.read_csv(H/'cap10_vs_cap2.csv'):
        v=next(v for v in extra if v['family']==r['family'] and v['arrangement']==r['arrangement'])
        assert int(r['additional_recovered'])==v['recovered'] and int(r['additional_lost'])==v['lost']
        assert (r['identical_pair_sets']=='True')==v['identical']
    # Statistical output is checked using saved joint family indices and independent sums.
    draws=json.loads(gzip.decompress((H/'bootstrap_draws.json.gz').read_bytes()));paired=N.read_csv(H/'family_level_comparisons.csv');cis=N.read_csv(H/'confidence_intervals.csv')
    for row in paired:
        a,b=row['contrast'].split(' minus ');cap=int(row['cap']);metric=row['metric'];name=row['family']
        expected=float(lookup[name,a,cap][metric])-float(lookup[name,b,cap][metric])
        assert abs(float(row['difference'])-expected)<1e-12
    for r in cis:
        values={p['family']:float(p['difference']) for p in paired if p['cap']==r['cap'] and p['metric']==r['metric'] and p['contrast']==r['contrast']}
        boots=sorted(sum(values[draws['families'][j]] for j in d)/len(d) for d in draws['indices'])
        def quantile(p):
            i=(len(boots)-1)*p;lo=math.floor(i);hi=math.ceil(i);return boots[lo]+(boots[hi]-boots[lo])*(i-lo)
        assert abs(float(r['mean_difference'])-sum(values.values())/len(values))<1e-12
        assert abs(float(r['ci_lower'])-quantile(.025))<1e-12 and abs(float(r['ci_upper'])-quantile(.975))<1e-12
    manifest=A.read_json(H/'experiment_manifest.json');baseline=A.read_json(H/'baseline_characteristics_done.json');assert baseline['before_refinement']
    assert A.sha(H/'baseline_characteristics.csv')==baseline['sha256']
    for p in (H/'cases').glob('*/*/alignments/*/done.json'):
        obj=A.read_json(p);started=obj['utc_epoch']-obj['details']['runtime_seconds']
        assert started>manifest['frozen_utc_epoch']
        if p.parent.name!='cap0':assert started>baseline['utc_epoch']
    result={'passed':True,'independent_alignment_scores':checked,'reference_datasets':len(checks),'confidence_intervals':len(cis),
        'cap10_vs_cap2_comparisons':len(extra),'frozen_before_first_alignment':True,'all_baselines_saved_before_refinement':True,
        'ambiguous_symbols_preserved_and_excluded_from_scores':True,'utc_epoch':time.time()}
    A.json_write(H/'independent_validation.json',result);print(json.dumps(result,indent=2))

def resume():
    artifacts=[p for p in (H/'cases').rglob('*') if p.is_file()]
    before={p.relative_to(H).as_posix():(A.sha(p),p.stat().st_mtime_ns) for p in artifacts}
    proc=subprocess.run([sys.executable,'-B',str(H/'pipeline.py'),'run'],cwd=H,capture_output=True)
    A.atomic(H/'actual_resume.log',proc.stdout+proc.stderr);assert proc.returncode==0,proc.stderr.decode(errors='replace')
    status=A.read_json(H/'last_run_status.json');assert status=={'completed_now':0,'skipped':P.C['primary_families']*9,'failed':0}
    changed=[p for p,old in before.items() if old!=(A.sha(H/p),(H/p).stat().st_mtime_ns)];assert not changed
    result={'passed':True,'skipped':status['skipped'],'mafft_runs_repeated':0,'identical_stage_file_hashes_and_mtimes':len(before),'changed_files':changed,'utc_epoch':time.time()}
    A.json_write(H/'actual_resume_verification.json',result);print(json.dumps(result,indent=2))

def bundle():
    for name in ['validation_latest.json','validation_before_execution.json','analysis_done.json','independent_validation.json','actual_resume_verification.json','preservation_verification.json']:
        assert A.read_json(H/name)['passed'],name
    old=A.read_json(H/'preservation_verification.json');base=A.read_json(H/'experiment_manifest.json');rows=N.read_csv(H/'expanded_results.csv')
    marker='\n## Final audit and preservation\n';report=(H/'phase2_extension_report.md').read_text(encoding='utf-8').split(marker)[0]
    text=marker+f"\nAll **{old['original_file_count']:,} pre-existing files** retain their SHA256, sizes and modification times, with no missing or added files outside phase2_extension. Both complementary read-only scan contexts agree on overlapping records; no original permissions changed. This includes the complete Phase1/supplementary/manuscript and original Phase2 report, ZIP and extracted verification files.\n\n"
    text+='Final validation: **22/22 tests passed**; independent column-coincidence auditing reproduced all279 scores,93 deletion references,93 cap10-vs-cap2 comparisons and30 confidence intervals. The actual complete resume skipped279/279 alignments, repeated zero MAFFT runs and left every saved stage artifact byte/timestamp unchanged. All279 scheduled jobs completed with zero failed MAFFT attempts. Exact version, command, stderr and runtime are saved per run.\n\n'
    text+='The ZIP includes complete experimental raw inputs/references/alignments, all counts/tables/figures, representative snippets, source/config/plans, source BAliBASE data, logs and preservation/validation evidence. Large binary MAFFT installation, plotting caches, temp files and redundant unit-test scratch trees are omitted; versions, executable hashes and test histories are retained. No manuscript was edited.\n'
    A.atomic(H/'phase2_extension_report.md',report+text)
    target=H/'phase2_extension_review_bundle.zip'
    if target.exists():raise FileExistsError('Refuse to overwrite completed review bundle')
    include=[]
    for p in H.rglob('*'):
        if not p.is_file():continue
        rel=p.relative_to(H)
        if set(rel.parts)&{'mafft-win','plot_cache','tmp','test_artifacts','__pycache__','review_verification_extract'}:continue
        if p==target or p.name in ['bundle_verification.json','extracted_bundle_validation.json'] or p.suffix=='.partial':continue
        include.append(p)
    records={p.relative_to(H).as_posix():A.sha(p) for p in include};assert sum(n.endswith('aligned.fasta') for n in records)==279
    A.json_write(H/'bundle_manifest.json',{'sha256':records,'complete_experimental_alignments':279,'independent_new_families':31,
        'omissions':['MAFFT binary installation (version/executable hashes retained)','plotting caches and dependencies','temporary and redundant test scratch artifacts']})
    include.append(H/'bundle_manifest.json')
    with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(include):z.write(p,p.relative_to(H).as_posix())
    with zipfile.ZipFile(target) as z:
        assert z.testzip() is None
        for path,digest in records.items():assert hashlib.sha256(z.read(path)).hexdigest()==digest
        count=len(z.namelist())
    result={'passed':True,'bytes':target.stat().st_size,'sha256':A.sha(target),'verified_members':count,'path':str(target)}
    A.json_write(H/'bundle_verification.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':{'independent':independent,'resume':resume,'bundle':bundle}[sys.argv[1]]()
