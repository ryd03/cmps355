"""Prespecified family-paired extension analysis; prior experiments read only."""
import sys
sys.dont_write_bytecode=True
from collections import Counter
import csv
import gzip
import itertools
import json
import math
import os
from pathlib import Path
import random
import statistics as st
import time
import pipeline as P
import scoring as S
import screen
import storage as A
H=P.H;C=P.C;ARR=['none','partial','high']

def read_csv(path):
    with Path(path).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
def avg(v):return st.mean(v) if v else None
def q(v,p):
    v=sorted(v);i=(len(v)-1)*p;a=int(i);b=min(a+1,len(v)-1);return v[a]+(v[b]-v[a])*(i-a)
def category(baseline_missed,net):
    if net>0:return 'C_improved'
    if net<0:return 'D_worsened'
    return 'A_baseline_perfect' if baseline_missed==0 else 'B_errors_no_net_gain'

def collect():
    selected=A.read_json(H/'selected_families.json');rows=[];caps=[];jobs=[]
    baseline={(r['family'],r['arrangement']):r for r in read_csv(H/'baseline_characteristics.csv')}
    for name,fam in selected.items():
        for arr in ARR:
            root=H/'cases'/name/arr
            for cap in C['caps']:
                path=root/f'alignments/cap{cap}';done=A.read_json(path/'done.json')
                for f,digest in done['files'].items():assert A.sha(path/f)==digest
                s=P.score(root,cap);local=s['accuracy']['local'];comp=s['comparison']['local'];b=baseline[name,arr]
                assert local['reference_pairs']>0
                row={'family':name,'component':fam['component'],'arrangement':arr,'cap':cap,'deletion_length':C['deletion_length'],
                    'offset':C['offsets'][arr],'overlap':max(0,C['deletion_length']-C['offsets'][arr]),
                    'overlap_fraction':max(0,C['deletion_length']-C['offsets'][arr])/C['deletion_length'],
                    'union_deleted_columns':C['deletion_length']+min(C['offsets'][arr],C['deletion_length']),
                    'local_reference_pairs':local['reference_pairs'],'local_correct_pairs':local['correct_pairs'],
                    'remaining_missed_correct_pairs':local['missed_pairs'],'remaining_evaluable_false_pairs':local['false_pairs'],
                    'local_sp_recall':local['sp_recall'],'local_precision':local['precision'],'global_core_sp_recall':s['accuracy']['global']['sp_recall'],
                    'baseline_local_recall':float(b['local_sp_recall']),'baseline_local_missed':int(b['missed_pairs']),
                    'baseline_missing_fraction':1-float(b['local_sp_recall']),'baseline_false_pairs':int(b['evaluable_false_pairs']),
                    'baseline_perfect':int(b['missed_pairs'])==0,'baseline_difficulty':b['difficulty'],
                    'baseline_mean_false_partner_displacement':float(b['mean_false_partner_reference_displacement']) if b['mean_false_partner_reference_displacement'] else None,
                    'recovered':comp['recovered'],'lost':comp['lost'],'net_pairs':comp['recovered']-comp['lost'],
                    'normalized_net_gain':comp['net_gain'],'recovery_rate':comp['recovery_rate'],'loss_rate':comp['loss_rate'],
                    'false_introduced':comp['false_introduced'],'false_removed':comp['false_removed'],
                    'pair_set_changed':s['all_pair_set_changed'],'residues_with_changed_partners':s['residues_with_changed_partners'],
                    'boundary_residues_with_changed_partners':s['boundary_residues_with_changed_partners'],
                    'mean_false_partner_displacement':s['mean_false_local_partner_reference_displacement'],
                    'runtime_seconds':s['runtime_seconds'],'convergence_message':s['convergence_reported'],
                    'outcome_category':category(int(b['missed_pairs']),comp['recovered']-comp['lost'])}
                rows.append(row)
                jobs.append({'family':name,'arrangement':arr,'cap':cap,'state':'completed','alignment_sha256':A.sha(path/'aligned.fasta'),
                    'input_sha256':A.sha(root/'input.fasta'),'command':json.dumps(done['details']['command']),
                    'failure_attempts_preserved':len(list(path.glob('failure_*.json')))})
            # Separate cap10 vs cap2 stage; successful MAFFT never rerun for analysis.
            marker=root/'scores/cap10_vs_cap2_done.json';out=root/'scores/cap10_vs_cap2.json'
            deps={'cap2':A.sha(root/'alignments/cap2/aligned.fasta'),'cap10':A.sha(root/'alignments/cap10/aligned.fasta'),
                'mapping':A.sha(root/'mapping.json'),'analysis_code':A.sha(H/'analysis.py')}
            if not A.complete(marker,deps):
                truth,boundary,known,meta=P.read_case(root)
                p2=S.pairs(S.parse_fasta((root/'alignments/cap2/aligned.fasta').read_text()))
                p10=S.pairs(S.parse_fasta((root/'alignments/cap10/aligned.fasta').read_text()))
                p2c={p for p in p2 if p[0] in known and p[1] in known};p10c={p for p in p10 if p[0] in known and p[1] in known}
                comparison=S.compare(p2c,p10c,truth,boundary)
                value={'all_pair_sets_identical':p2==p10,'comparison':comparison}
                A.json_write(out,value);A.marker(marker,[out],deps)
            value=A.read_json(out);local=value['comparison']['local']
            caps.append({'family':name,'arrangement':arr,'identical_pair_sets':value['all_pair_sets_identical'],
                'additional_recovered':local['recovered'],'additional_lost':local['lost'],
                'net_pairs':local['recovered']-local['lost'],'normalized_cap10_minus_cap2':local['net_gain']})
            print('Analysis collected',name,arr,flush=True)
    assert len(rows)==C['primary_families']*9
    screen.csv_write(H/'expanded_results.csv',rows);screen.csv_write(H/'cap10_vs_cap2.csv',caps);screen.csv_write(H/'job_manifest.csv',jobs)
    return rows,caps

def paired(rows):
    names=sorted({r['family'] for r in rows});idx={(r['family'],r['arrangement'],r['cap']):r for r in rows}
    rng=random.Random(C['bootstrap_seed']);draws=[[rng.randrange(len(names)) for _ in names] for _ in range(C['bootstrap_resamples'])]
    A.atomic(H/'bootstrap_draws.json.gz',gzip.compress(A.canonical({'families':names,'seed':C['bootstrap_seed'],'indices':draws}),mtime=0))
    comparisons=[];estimates=[];loo=[]
    for cap in [2,10]:
        for metric in ['normalized_net_gain','recovered','lost','baseline_local_recall','baseline_local_missed']:
            for a,b in [('partial','none'),('high','none'),('partial','high')]:
                values=[idx[n,a,cap][metric]-idx[n,b,cap][metric] for n in names];boots=[sum(values[i] for i in draw)/len(names) for draw in draws]
                for n,value in zip(names,values):comparisons.append({'family':n,'cap':cap,'metric':metric,'contrast':a+' minus '+b,'difference':value})
                estimates.append({'cap':cap,'metric':metric,'contrast':a+' minus '+b,'n_families':len(names),
                    'mean_difference':avg(values),'ci_lower':q(boots,.025),'ci_upper':q(boots,.975),
                    'positive_families':sum(v>0 for v in values),'negative_families':sum(v<0 for v in values),'zero_families':sum(v==0 for v in values),
                    'primary':cap==2 and metric=='normalized_net_gain','bootstrap_resamples':len(draws)})
                if cap==2 and metric=='normalized_net_gain':
                    for i,n in enumerate(names):loo.append({'omitted_family':n,'contrast':a+' minus '+b,'remaining_families':len(names)-1,
                        'mean_difference':avg(values[:i]+values[i+1:]),'full_mean':avg(values)})
    screen.csv_write(H/'family_level_comparisons.csv',comparisons);screen.csv_write(H/'confidence_intervals.csv',estimates);screen.csv_write(H/'leave_one_family_out.csv',loo)
    return estimates

def groups(rows):
    output=[]
    for cap in C['caps']:
        for arr in ARR:
            for group in ['all','perfect','imperfect','near-perfect','difficult']:
                d=[r for r in rows if r['cap']==cap and r['arrangement']==arr and (group=='all' or
                    group=='perfect' and r['baseline_perfect'] or group=='imperfect' and not r['baseline_perfect'] or group==r['baseline_difficulty'])]
                output.append({'cap':cap,'arrangement':arr,'baseline_group':group,'n_families':len(d),
                    'mean_local_recall':avg([r['local_sp_recall'] for r in d]),'mean_global_core_recall':avg([r['global_core_sp_recall'] for r in d]),
                    'mean_baseline_recall':avg([r['baseline_local_recall'] for r in d]),'mean_normalized_net_gain':avg([r['normalized_net_gain'] for r in d]),
                    'recovered_sum':sum(r['recovered'] for r in d),'lost_sum':sum(r['lost'] for r in d),
                    'n_improved':sum(r['net_pairs']>0 for r in d),'n_worsened':sum(r['net_pairs']<0 for r in d),'n_unchanged':sum(r['net_pairs']==0 for r in d),
                    'changed_pair_sets':sum(r['pair_set_changed'] for r in d),'mean_runtime_seconds':avg([r['runtime_seconds'] for r in d])})
    screen.csv_write(H/'group_summary.csv',output);return output

def ranks(v):
    order=sorted(range(len(v)),key=lambda i:v[i]);result=[0.0]*len(v);j=0
    while j<len(v):
        k=j+1
        while k<len(v) and v[order[k]]==v[order[j]]:k+=1
        for z in order[j:k]:result[z]=(j+k-1)/2+1
        j=k
    return result
def pearson(a,b):
    ma,mb=avg(a),avg(b);x=[v-ma for v in a];y=[v-mb for v in b];den=math.sqrt(sum(v*v for v in x)*sum(v*v for v in y))
    return sum(v*w for v,w in zip(x,y))/den if den else None

def explanatory(rows):
    before=read_csv(H/'pre_alignment_characteristics.csv');features=[]
    numeric=['n_sequences','reference_columns','mean_core_pairwise_identity','focal_core_pairwise_identity','local_conservation',
        'local_normalized_entropy','local_natural_gap_density','local_gap_run_count','local_gap_run_mean_length','local_gap_run_max_length',
        'local_core_coverage','family_core_coverage','local_scorable_residues','local_ambiguous_fraction']
    for f in before:
        n=f['family'];d=[r for r in rows if r['family']==n and r['cap']==2]
        out={**f,**{k:float(f[k]) for k in numeric},'baseline_missing_fraction':avg([r['baseline_missing_fraction'] for r in d]),
            'baseline_mean_missed_pairs':avg([r['baseline_local_missed'] for r in d]),
            'mean_normalized_net_gain_cap2':avg([r['normalized_net_gain'] for r in d]),
            'mean_recovered_cap2':avg([r['recovered'] for r in d]),'mean_lost_cap2':avg([r['lost'] for r in d]),
            'n_improved_arrangements':sum(r['net_pairs']>0 for r in d),'n_worsened_arrangements':sum(r['net_pairs']<0 for r in d),
            'baseline_perfect_arrangements':sum(r['baseline_perfect'] for r in d),'baseline_mean_local_reference_pairs':avg([r['local_reference_pairs'] for r in d])}
        features.append(out)
    screen.csv_write(H/'family_characteristics.csv',features)
    correlations=[];outcomes=['mean_normalized_net_gain_cap2','mean_recovered_cap2','mean_lost_cap2']
    for covariate in ['baseline_missing_fraction',*numeric]:
        for outcome in outcomes:
            x=[r[covariate] for r in features];y=[r[outcome] for r in features]
            correlations.append({'feature':covariate,'outcome':outcome,'n_families':len(features),'spearman_rho':pearson(ranks(x),ranks(y)),
                'interpretation':'Exploratory descriptive association, no p-value or causal attribution'})
    screen.csv_write(H/'explanatory_correlations.csv',correlations)
    regression={'fitted':False,'reason':'Fewer than30 independent families','n_families':len(features)}
    if len(features)>=30:
        # Only this small prespecified exploratory model; no outcome-selected predictors.
        sys.path.insert(0,str(H.parent/'msa_refinement/.deps'))
        import numpy as np
        y=np.array([r['mean_normalized_net_gain_cap2'] for r in features]);base=np.array([[1,r['baseline_missing_fraction']] for r in features])
        full=np.array([[1,r['baseline_missing_fraction'],r['local_natural_gap_density']] for r in features])
        def fit(x):
            beta,_,rank,_=np.linalg.lstsq(x,y,rcond=None);res=y-x@beta;sst=float(sum((y-y.mean())**2));sse=float(sum(res**2))
            return {'coefficients':beta.tolist(),'rank':int(rank),'training_r_squared':1-sse/sst if sst else None,'residuals':res.tolist()}
        b,f=fit(base),fit(full)
        regression={'fitted':True,'n_families':len(features),'baseline_only':b,'baseline_plus_natural_gap_density':f,
            'coefficient_order':['intercept','baseline_missing_fraction','local_natural_gap_density'],
            'incremental_training_r_squared':f['training_r_squared']-b['training_r_squared'],
            'caution':'Descriptive training fit, not causal explained variation or validated prediction; baseline errors mechanically limit recovery; only31 families.'}
    A.json_write(H/'exploratory_regression.json',regression)
    return features,correlations,regression

def compare_previous(rows,caps):
    sources={'phase1':H.parent/'msa_refinement/runs/full_deletion_overlap_v2_0815b46e6b9208a2/full_summary.csv',
        'pilot':H.parent/'phase2_balibase_real_proteins/phase2_summary.csv'}
    previous={};provenance={}
    for k,source in sources.items():
        dest=H/'comparison_sources'/source.name
        if source.exists():
            if not dest.exists():A.atomic(dest,source.read_bytes())
            assert A.sha(dest)==A.sha(source)
        else:
            stored=A.read_json(H/'comparison_source_provenance.json')[k]
            assert A.sha(dest)==stored['sha256']
        previous[k]=read_csv(dest);provenance[k]={'original':str(source),'sha256':A.sha(dest),'local_copy':dest.relative_to(H).as_posix()}
    output=[]
    for cap in [2,10]:
        for arr in ARR:
            for design,div in [('Phase1','0.25'),('Phase1','0.40'),('Phase2 pilot','real'),('Phase2 extension','real')]:
                if design=='Phase1':
                    d=[r for r in previous['phase1'] if float(r['divergence'])==float(div) and int(r['iteration_cap'])==cap and r['condition']==arr]
                    missing=[int(r['baseline_local_missed']) for r in d];den=[int(r['local_reference_pairs']) for r in d]
                    rec=[int(r['local_recovered']) for r in d];lost=[int(r['local_lost']) for r in d];gain=[float(r['local_net_gain']) for r in d]
                    changed=[r['alignment_pair_set_changed']=='True' for r in d];recall=[1-m/n for m,n in zip(missing,den)]
                    agree=len(d);unit='evolutionary seed';spacing=50
                elif design=='Phase2 pilot':
                    d=[r for r in previous['pilot'] if int(r['cap'])==cap and r['arrangement']==arr]
                    missing=[int(r['baseline_local_missed']) for r in d];rec=[int(r['recovered']) for r in d];lost=[int(r['lost']) for r in d]
                    gain=[float(r['normalized_net_gain']) for r in d];recall=[float(r['baseline_local_recall']) for r in d]
                    changed=[r['pair_set_changed']=='True' for r in d];agree=sum(r['cap2_cap10_same_pairs']=='True' for r in d);unit='pilot family';spacing=50
                else:
                    d=[r for r in rows if r['cap']==cap and r['arrangement']==arr];missing=[r['baseline_local_missed'] for r in d]
                    rec=[r['recovered'] for r in d];lost=[r['lost'] for r in d];gain=[r['normalized_net_gain'] for r in d]
                    recall=[r['baseline_local_recall'] for r in d];changed=[r['pair_set_changed'] for r in d]
                    agree=sum(r['identical_pair_sets'] for r in caps if r['arrangement']==arr);unit='new family component';spacing=C['deletion_length']
                rsum,lsum=sum(rec),sum(lost)
                output.append({'design':design,'divergence':div,'arrangement':arr,'cap':cap,'n_units':len(d),'unit':unit,
                    'deletion_length':12,'none_offset':spacing,'baseline_mean_recall':avg(recall),'baseline_mean_missed_pairs':avg(missing),
                    'baseline_perfect_fraction':sum(m==0 for m in missing)/len(d),'recovered_sum':rsum,'lost_sum':lsum,
                    'aggregate_cancellation_fraction':2*min(rsum,lsum)/(rsum+lsum) if rsum+lsum else None,
                    'mean_normalized_net_gain':avg(gain),'changed_alignments':sum(changed),'cap2_cap10_identical_alignments':agree})
    screen.csv_write(H/'previous_experiments_comparison.csv',output);A.json_write(H/'comparison_source_provenance.json',provenance)
    return output

def examples(rows):
    selected=A.read_json(H/'selected_families.json');manifest=[]
    for label in ['A_baseline_perfect','B_errors_no_net_gain','C_improved','D_worsened']:
        data=sorted([r for r in rows if r['cap']==2 and r['outcome_category']==label],key=lambda r:(r['family'],r['arrangement']))
        if not data:
            manifest.append({'category':label,'count':0,'representative':None});continue
        r=data[0];root=H/'cases'/r['family']/r['arrangement'];meta=A.read_json(root/'mapping.json')
        boundary=S.boundary_residues(meta);dest=H/'examples'/label;dest.mkdir(parents=True,exist_ok=True)
        mapping=[]
        for title,path in [('reference',root/'reference.fasta'),('cap0',root/'alignments/cap0/aligned.fasta'),('cap2',root/'alignments/cap2/aligned.fasta'),('cap10',root/'alignments/cap10/aligned.fasta')]:
            aligned=S.parse_fasta(path.read_text());coords={};counter={n:0 for n in aligned}
            for n,seq in aligned.items():
                for col,aa in enumerate(seq):
                    if aa!='-':counter[n]+=1;coords[n,counter[n]]=col
            columns=[c for residue,c in coords.items() if residue in boundary];lo=max(0,min(columns)-3);hi=min(len(next(iter(aligned.values()))),max(columns)+4)
            A.atomic(dest/(title+'_snippet.fasta'),P.B.fasta({n:s[lo:hi] for n,s in aligned.items()}))
            for residue in sorted(boundary):
                n,i=residue;mapping.append({'alignment':title,'sequence_id':n,'ungapped_survivor_index':i,'original_reference_column':meta['ancestry'][n][i-1],
                    'alignment_column_1based':coords[residue]+1,'snippet_column_1based':coords[residue]-lo+1,'amino_acid':aligned[n][coords[residue]]})
        screen.csv_write(dest/'focal_residue_mapping.csv',mapping)
        A.json_write(dest/'case.json',r)
        manifest.append({'category':label,'count':len(data),'representative':r['family']+'/'+r['arrangement'],'recovered':r['recovered'],'lost':r['lost'],
            'baseline_missed':r['baseline_local_missed'],'caveat':'First lexicographic example, not a frequency estimate; snippet coordinates differ between alignments'})
    A.json_write(H/'examples/manifest.json',manifest);return manifest

def main():
    P.verify();rows,caps=collect();estimates=paired(rows);summary=groups(rows);features,corr,reg=explanatory(rows);comparison=compare_previous(rows,caps);representatives=examples(rows)
    import figures_report
    figures_report.generate(rows,caps,estimates,summary,features,corr,reg,comparison,representatives)
    A.json_write(H/'analysis_done.json',{'passed':True,'rows':len(rows),'families':len(features),'confidence_intervals':len(estimates),'utc_epoch':time.time()})
    print('All extension analysis outputs complete',len(rows),'runs',flush=True)
if __name__=='__main__':main()
