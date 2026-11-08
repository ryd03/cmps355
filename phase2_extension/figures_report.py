"""Transparent descriptive plots and report; figures do not drive eligibility."""
import sys
sys.dont_write_bytecode=True
from collections import Counter
import os
import statistics as st
from pathlib import Path
import analysis as N
import storage as A
H=N.H;ARR=N.ARR

def plots(rows,caps,estimates,summary,features,correlations,regression,comparison):
    for key in ['MPLCONFIGDIR','XDG_CACHE_HOME','TMP','TEMP','TMPDIR']:
        p=H/('plot_cache' if key in ['MPLCONFIGDIR','XDG_CACHE_HOME'] else 'tmp');p.mkdir(exist_ok=True);os.environ[key]=str(p)
    sys.path.insert(0,str(H.parent/'msa_refinement/.deps'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    colors=['#0072B2','#D55E00','#009E73'];dest=H/'figures';dest.mkdir(exist_ok=True)
    def save(fig,name):
        fig.savefig(dest/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(dest/(name+'.svg'),bbox_inches='tight');plt.close(fig)
    def jitter(i):return ((i*7)%17-8)/65
    def points(ax,data,field,scale=1):
        for j,a in enumerate(ARR):
            d=[r for r in data if r['arrangement']==a]
            ax.scatter([j+jitter(i) for i in range(len(d))],[scale*r[field] for r in d],c=colors[j],s=22,alpha=.7)
            ax.plot([j-.22,j+.22],[scale*st.mean(r[field] for r in d)]*2,color='black',lw=2)
        ax.set_xticks(range(3),ARR)
    sensitivity=N.read_csv(H/'screening_sensitivity.csv')
    fig,ax=plt.subplots(figsize=(8,4));ds=[12,10,8,6]
    for k,(mode,color) in enumerate([('continuous','#999999'),('targeted','#0072B2')]):
        v=[next(int(r['new_independent_components']) for r in sensitivity if int(r['length'])==d and int(r['window'])==10 and r['reliability']==mode and r['symbols']=='retained') for d in ds]
        bars=ax.bar([i+(k-.5)*.35 for i in range(4)],v,width=.35,label=mode+' core rule',color=color)
        ax.bar_label(bars,padding=3)
    ax.set_xticks(range(4),ds);ax.set_xlabel('Candidate deletion length (amino acids)');ax.set_ylabel('New independent family components');ax.legend()
    ax.set_title('Source-only screening: window 10; X/B/Z retained\nPrimary keeps length 12 (31 new families); original 82 rule: 0 new');fig.tight_layout();save(fig,'screening_sensitivity')
    base=[r for r in rows if r['cap']==0]
    fig,axes=plt.subplots(1,2,figsize=(9,4));points(axes[0],base,'baseline_local_recall');points(axes[1],base,'baseline_local_missed')
    axes[0].set_ylabel('Cap-0 local core SP recall');axes[0].axhline(.95,color='.5',ls=':');axes[1].set_ylabel('Missing correct local pairs at cap0')
    fig.suptitle('Complete baseline distribution: 31 families × 3 arrangements');fig.tight_layout();save(fig,'baseline_difficulty')
    fig,axes=plt.subplots(1,2,figsize=(10,4),sharey=True)
    for ax,cap in zip(axes,[2,10]):
        d=[r for r in rows if r['cap']==cap]
        for j,a in enumerate(ARR):
            v=[r for r in d if r['arrangement']==a]
            ax.scatter([j-.12+jitter(i)/2 for i in range(len(v))],[r['recovered'] for r in v],c=colors[j],s=18)
            ax.scatter([j+.12+jitter(i)/2 for i in range(len(v))],[-r['lost'] for r in v],c=colors[j],marker='v',s=20)
        ax.set_xticks(range(3),ARR);ax.axhline(0,color='.5');ax.set_title(f'Cap {cap} vs cap0')
    axes[0].set_ylabel('Correct local pairs recovered (+) / lost (−)');fig.tight_layout();save(fig,'recovered_lost')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    d=[r for r in rows if r['cap']==2];points(axes[0],d,'normalized_net_gain',100);axes[0].axhline(0,color='.5');axes[0].set_ylabel('Cap2 local net gain (percentage points)')
    primary=[r for r in estimates if r['primary']]
    paired=N.read_csv(H/'family_level_comparisons.csv')
    for j,r in enumerate(primary):
        values=[100*float(p['difference']) for p in paired if p['cap']=='2' and p['metric']=='normalized_net_gain' and p['contrast']==r['contrast']]
        axes[1].scatter(values,[j+jitter(i) for i in range(len(values))],c='.6',s=16)
        axes[1].plot([100*r['ci_lower'],100*r['ci_upper']],[j,j],color='#0072B2',lw=3);axes[1].scatter([100*r['mean_difference']],[j],c='#D55E00',marker='D',zorder=3)
    axes[1].set_yticks(range(3),[r['contrast'] for r in primary]);axes[1].axvline(0,color='.4',ls=':');axes[1].set_xlabel('Paired difference (percentage points)')
    axes[1].set_title('31 paired family values; mean and marginal 95% CI');fig.tight_layout();save(fig,'matched_gains')
    labels=['A_baseline_perfect','B_errors_no_net_gain','C_improved','D_worsened'];palette=['#999999','#E69F00','#009E73','#CC79A7']
    fig,ax=plt.subplots(figsize=(8,4));bottom=[0,0,0]
    for label,col in zip(labels,palette):
        counts=[sum(r['outcome_category']==label and r['arrangement']==a for r in d) for a in ARR]
        ax.bar(ARR,counts,bottom=bottom,label=label.replace('_',' '),color=col);bottom=[a+b for a,b in zip(bottom,counts)]
    ax.set_ylabel('Family datasets (cap2 vs cap0)');ax.legend(fontsize=8,bbox_to_anchor=(1.02,1));fig.tight_layout();save(fig,'outcome_categories')
    xfields=[('baseline_missing_fraction','Baseline missing fraction'),('n_sequences','Protein sequences'),('mean_core_pairwise_identity','Family core identity'),
        ('focal_core_pairwise_identity','Focal pair core identity'),('local_conservation','Local conservation'),('local_natural_gap_density','Local natural gap density'),
        ('local_gap_run_count','Local natural gap runs'),('local_core_coverage','Local core coverage'),('local_normalized_entropy','Local normalized entropy')]
    fig,axes=plt.subplots(3,3,figsize=(11,9),sharey=True)
    for ax,(field,label) in zip(axes.flat,xfields):
        dot=ax.scatter([r[field] for r in features],[100*r['mean_normalized_net_gain_cap2'] for r in features],
            c=[r['baseline_missing_fraction'] for r in features],cmap='viridis',vmin=0,vmax=1,s=28)
        ax.axhline(0,color='.7',lw=.7);ax.set_xlabel(label)
        rho=next(r['spearman_rho'] for r in correlations if r['feature']==field and r['outcome']=='mean_normalized_net_gain_cap2')
        ax.set_title('Spearman ρ='+('undefined' if rho is None else f'{rho:.2f}'),fontsize=9)
    for ax in axes[:,0]:ax.set_ylabel('Family-mean cap2 net gain (pp)')
    fig.suptitle('Exploratory family-level associations; colour = baseline missing fraction');fig.tight_layout(rect=(0,0,1,.96));save(fig,'family_explanatory_characteristics')
    fig,axes=plt.subplots(1,2,figsize=(10,4));points(axes[0],caps,'normalized_cap10_minus_cap2',100);axes[0].axhline(0,color='.5');axes[0].set_ylabel('Cap10 minus cap2 local recall (pp)')
    axes[1].scatter([r['boundary_residues_with_changed_partners'] for r in d],[100*r['normalized_net_gain'] for r in d],c=[colors[ARR.index(r['arrangement'])] for r in d],s=25)
    axes[1].axhline(0,color='.5');axes[1].set_xlabel('Boundary residues with changed partners');axes[1].set_ylabel('Cap2 local net gain (pp)');fig.tight_layout();save(fig,'caps_and_partner_changes')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for design,div in [('Phase1','0.25'),('Phase1','0.40'),('Phase2 pilot','real'),('Phase2 extension','real')]:
        data=[next(r for r in comparison if r['design']==design and r['divergence']==div and r['cap']==2 and r['arrangement']==a) for a in ARR]
        axes[0].plot(range(3),[r['baseline_mean_recall'] for r in data],marker='o',label=design+' '+div)
        axes[1].plot(range(3),[100*r['mean_normalized_net_gain'] for r in data],marker='o')
    for ax in axes:ax.set_xticks(range(3),ARR)
    axes[0].set_ylabel('Mean baseline local recall');axes[1].set_ylabel('Mean cap2 local net gain (pp)');axes[0].legend(fontsize=8)
    fig.suptitle('Separate protocols: no pooling; extension none-offset12 vs previous50');fig.tight_layout();save(fig,'previous_experiments')

def report(rows,caps,estimates,summary,features,correlations,regression,comparison,representatives):
    d=[r for r in rows if r['cap']==2];base=[r for r in rows if r['cap']==0];screening=N.read_csv(H/'screening_sensitivity.csv')
    n=len(features);primary=[r for r in estimates if r['primary']]
    lines=['# Phase 2 extension — controlled deletions across real-protein families','',
      f'Completed **{n} new independent family components**, {len(d)} matched datasets and **{len(rows)} MAFFT alignments**. Every selected family is retained. Original pilot families BB12007/BB12031 and their shared-protein components are excluded from primary estimates. The extension tests a changed screening/spacing protocol, not a replacement for the original pilot.','',
      '## Outcome-blind protocol and screening','',
      'Before MAFFT, SCREENING_PLAN.md specified a finite sensitivity grid and a deterministic selection rule: retain window10 and targeted reliability, then choose the longest deletion length attaining20 new family components (otherwise the most components, ties longer). This selected **length12**, with31 new components; lengths10/8/6 and window5 were screened only, not run as additional factors. No eligibility rule used baseline or refined alignment outcomes.','',
      'The verified original [BAliBASE3.0 distribution](https://www.lbgi.fr/balibase/BalibaseDownload/) is reused read only, SHA256 `'+N.P.screen.ARCHIVE_SHA+'`. Full-length BB references are used; BBS truncated copies are not extra families. All218 full-length cases are audited in screening_audit_218.csv. All cases, including exclusions, enter the shared-ID/exact-sequence connected-component graph. This operational clustering reduces reuse of the same proteins but does not guarantee all components are unrelated. One eligible case/site per new component is selected deterministically.','',
      'Changes: the pilot required an82-column uninterrupted core-and-gap-free span with no-overlap offset50. The extension requires reliable deletion columns and at least3 reliable standard residues on each side of every focal event within10 columns, permitting noncore/gapped flanks where scoring remains reliable. No-overlap offset is now12, partial6 and high1. Overlap0/6/11 (0%,50%,91.67%); union24/18/13. Thus spacing and site eligibility changed, even though length12 is retained; do not pool the protocols. Practical limits remain4–60 sequences and<=2000 columns; no cropping/downsampling is used.','',
      'X/B/Z remain unchanged in inputs and outputs: no invented amino acids or silent replacements. [MAFFT documentation](https://mafft.cbrc.jp/alignment/software/anysymbol.html) accepts those symbols and treats them as unknown for scoring. Either ambiguous endpoint is excluded from reference denominators and evaluable predicted pairs; focal deletions and required core flanks must contain standard amino acids. Curated coreblock=1 is the reliability criterion. Other columns are not asserted to be biologically false.','',
      '| Length | Continuous core: new families | Targeted core: new families |','|---:|---:|---:|']
    for length in [12,10,8,6]:
        vals=[next(r['new_independent_components'] for r in screening if int(r['length'])==length and r['window']=='10' and r['symbols']=='retained' and r['reliability']==mode) for mode in ['continuous','targeted']]
        lines.append(f'| {length} | {vals[0]} | {vals[1]} |')
    lines+=['','For length12/window10, targeted screening yields23 new components when any family with ambiguity is excluded, versus31 when ambiguities are retained and excluded only from scoring. Total eligible components are33, including the two pilot components; **31 are genuinely new**. The complete33-protocol grid is in screening_sensitivity.csv. This design prefers the longer already-studied deletion length rather than maximizing observations with shorter events.','',
      '## Measurement and statistical units','',
      'Survivors retain immutable IDs (sequence ID, ungapped index) and original reference columns. Exactly12 residues are removed in each of two focal proteins. Correct local pairs are core-standard reference pairs involving an extant focal residue within10 reference columns inclusive of its OWN deletion edges. Correct pairs recovered/lost are compared against cap0. Normalized net gain=(recovered−lost)/true local-pair count; undefined denominators/rates remain null. Remaining evaluable false pairs and missed true pairs are separate quantities. One residue can change many partners; partner-change counts are complementary, not a count of evolutionary displacements.','',
      'Primary estimates compare partial−none, high−none and partial−high at cap2 using all31 matched family components. Paired bootstrap:10000 joint family resamples, seed3553202, percentile marginal95% intervals with linear interpolation. Pairs and the three arrangements are not independent replicates. All primary/secondary contrasts, individual-family differences and leave-one-family-out estimates are saved. No multiplicity-adjusted discovery or independent causal-overlap claim is made. Union, spacing, deleted identities and local context change together.','',
      '## Baseline difficulty, before refinement','',
      f"All93 cap0 jobs and baseline_characteristics.csv were saved before any cap2/cap10 jobs. Baseline local recall spans {min(r['baseline_local_recall'] for r in base):.5f}–{max(r['baseline_local_recall'] for r in base):.5f}; {sum(r['baseline_perfect'] for r in base)}/93 datasets are already perfect. Secondary groups were fixed beforehand: perfect, near-perfect [0.95,1), difficult <0.95, plus all-imperfect. Group_summary.csv retains empty groups and all cases; primary results never select only successful refinements.",'',
      '| Arrangement | Perfect baselines | Mean baseline recall | Mean missing correct pairs |','|---|---:|---:|---:|']
    for a in ARR:
        b=[r for r in base if r['arrangement']==a];lines.append(f"| {a} | {sum(r['baseline_perfect'] for r in b)}/{n} | {N.avg([r['baseline_local_recall'] for r in b]):.5f} | {N.avg([r['baseline_local_missed'] for r in b]):.2f} |")
    lines+=['','## Refinement outcomes','',
      '| Cap | Arrangement | Recovered | Lost | Mean net gain (pp) | Improved / worsened / unchanged | Changed alignments |','|---:|---|---:|---:|---:|---|---:|']
    for r in summary:
        if r['baseline_group']=='all' and r['cap']:
            lines.append(f"| {r['cap']} | {r['arrangement']} | {r['recovered_sum']} | {r['lost_sum']} | {100*r['mean_normalized_net_gain']:.5f} | {r['n_improved']} / {r['n_worsened']} / {r['n_unchanged']} | {r['changed_pair_sets']}/{n} |")
    lines+=['','Counts are descriptive totals, not independent residue-pair sample sizes. Expanded_results.csv supplies every dataset/cap, exact denominator, local/global core SP recall, recovered/lost, false pairs, recovery/loss rates, remaining errors, runtime and partner changes.','',
      '| Primary contrast | Mean (pp) | Marginal95% interval (pp) | Positive / negative / zero family differences |','|---|---:|---|---|']
    for r in primary:lines.append(f"| {r['contrast']} | {100*r['mean_difference']:.5f} | [{100*r['ci_lower']:.5f}, {100*r['ci_upper']:.5f}] | {r['positive_families']} / {r['negative_families']} / {r['zero_families']} |")
    consistency=all(r['ci_lower']>0 or r['ci_upper']<0 for r in primary)
    lines+=['', 'The overlap contrasts '+('have intervals away from zero in this sample, but remain protocol-specific descriptive evidence rather than independent causal effects.' if consistency else 'do not all have intervals away from zero; differences are not uniformly consistent across the expanded families.'),
      ' Family-level signs and leave-one-out estimates show whether the mean is shared broadly or driven by a few cases; no formal significance declaration follows from these intervals.','',
      '## Why some families gain and others do not','',
      'These are measured associations and illustrative cases, not proven optimization mechanisms. A perfect baseline cannot recover additional correct local pairs; an imperfect baseline can still remain unchanged, recover/loss-cancel, or worsen. A positive net result requires more corrected relationships than destroyed correct ones; global alignment changes alone do not establish local correction.','',
      '| Category | Datasets at cap2 | Representative |','|---|---:|---|']
    for r in representatives:lines.append(f"| {r['category']} | {r['count']} | {r['representative'] or 'Absent; no example manufactured'} |")
    lines+=['','Observed details of these illustrations:']
    for example in representatives:
        if not example['representative']:continue
        family,arr=example['representative'].split('/')
        r=next(r for r in d if r['family']==family and r['arrangement']==arr)
        f=next(f for f in features if f['family']==family)
        lines+=['',f"- **{example['category']}: {family}/{arr}** — baseline recall {r['baseline_local_recall']:.5f}, missing{r['baseline_local_missed']}; recovered{r['recovered']}, lost{r['lost']}, remaining missed{r['remaining_missed_correct_pairs']}. Changed boundary partner sets for{r['boundary_residues_with_changed_partners']} residues. Family mean core identity{f['mean_core_pairwise_identity']:.4f}, local conservation{f['local_conservation']:.4f}, natural-gap fraction{f['local_natural_gap_density']:.4f}. The saved reference and cap0/2/10 snippets show the actual placements; these observed feature values alone do not establish why the internal search retained or changed them."]
    lines+=['','Each representative directory contains before/reference/cap2/cap10 amino-acid snippets, all focal survivor IDs/original columns/predicted columns and exact counts. First lexicographic cases are illustrations, not chosen to exaggerate effects. Snippets have different predicted column coordinates and are not independent samples. CategoryB can contain recovered/lost cancellation; inspect raw counts.','',
      'Original-reference characteristics were saved before alignment: family/focal core identity, sequence count, local majority-residue conservation and normalized Shannon entropy, natural gap fraction and clipped run lengths, ambiguity, core coverage and scorable information. The common local region is44 reference columns spanning both possible deletions plus10-column flanks. Family_characteristics.csv adds cap0 difficulty and cap2 outcome averages across arrangements; every family has equal weight. Baseline false-partner displacement is mean absolute original reference-column separation of wrongly paired evaluable local residues, not a biological ambiguity probability.','',
      '| Feature | Spearman with net gain | With recovered pairs | With lost pairs |','|---|---:|---:|---:|']
    for f in ['baseline_missing_fraction','n_sequences','mean_core_pairwise_identity','focal_core_pairwise_identity','local_conservation','local_natural_gap_density','local_gap_run_count','local_core_coverage']:
        cells=[]
        for outcome in ['mean_normalized_net_gain_cap2','mean_recovered_cap2','mean_lost_cap2']:
            rho=next(r['spearman_rho'] for r in correlations if r['feature']==f and r['outcome']==outcome)
            cells.append('undefined' if rho is None else f'{rho:.5f}')
        lines.append('| '+f+' | '+' | '.join(cells)+' |')
    lines+=['','In this sample, family-mean gains have a positive rank association with baseline missing fraction (ρ≈0.35) and negative associations with local conservation (ρ≈−0.47) and family-wide core identity (ρ≈−0.30). These describe where larger gains were observed, not deterministic explanations. Sequence count and focal-pair identity show little rank association with net gain. Gross pair losses show a positive exploratory association with family-wide core coverage (ρ≈0.56) and a negative association with sequence count (ρ≈−0.41), while local conservation/gap density have little association with gross loss. Raw pair counts depend on available/scorable relationships, so these correlations must not be interpreted as per-residue damage risks. All prespecified correlations, including weak and opposing results, are retained; no selection of only favorable associations is used.']
    if regression['fitted']:
        b=regression['baseline_only'];f=regression['baseline_plus_natural_gap_density']
        lines+=['',f"The prespecified small exploratory family-level linear model uses mean baseline missing fraction and original local natural-gap density only (31 families, intercept plus2 predictors). Baseline-only training R²={b['training_r_squared']:.5f}; adding gap density gives R²={f['training_r_squared']:.5f}, increment={regression['incremental_training_r_squared']:.5f}. Coefficients (fraction gain units): {f['coefficients']}. Rank={f['rank']}. These are descriptive in-sample fits, not causal attribution, validated prediction or a decomposition of biological mechanisms. Baseline error mechanically bounds possible recovery; correlated features and influential cases can affect apparent associations. No p-values are calculated and no covariates are added after outcomes."]
    lines+=['','Low identity, gap-rich regions and lower conservation can plausibly create multiple competing placements, but these summaries cannot prove that MAFFT chose an internally optimal yet reference-incorrect alignment. Persistent errors can reflect remaining objective/search limitations or imperfect reference homology. The data establish which pairs changed and what features covary with those changes; explaining the optimization mechanism would require a separate investigation.','',
      '## Does cap10 add benefit?','',
      f"Cap2 and cap10 pair sets are identical in {sum(r['identical_pair_sets'] for r in caps)}/{len(caps)} datasets. Relative to cap2, cap10 improves local recall in {sum(r['net_pairs']>0 for r in caps)}, worsens it in {sum(r['net_pairs']<0 for r in caps)}, and leaves it unchanged in {sum(r['net_pairs']==0 for r in caps)}. It additionally recovers {sum(r['additional_recovered'] for r in caps)} and loses {sum(r['additional_lost'] for r in caps)} correct local pairs. Cap10_vs_cap2.csv retains every case; more iterations need not monotonically improve curated-reference accuracy. MAFFT convergence diagnostics describe internal sub-alignments and may coexist with reached-cap messages, so they do not guarantee cap2/cap10 equality.",'',
      '## Comparison with earlier experiments','',
      'Saved previous results are copied byte-for-byte into comparison_sources; no earlier pipeline is invoked. Keep all designs separate: Phase1 has300 simulated datasets from50 seed blocks at each divergence; the original pilot has2 families/6 datasets; the extension has31 NEW family components/93 datasets. No-overlap spacing50 versus12 and changed site reliability/ambiguous-endpoint treatment limit direct comparisons despite length12 throughout. Simulations have exact ancestral truth; real proteins have curated core-supported references.','',
      '| Study | Divergence | Arrangement | Perfect baseline fraction | Mean baseline recall | Mean cap2 gain (pp) | Recovered / lost | Cancellation fraction |','|---|---|---|---:|---:|---:|---|---:|']
    for r in comparison:
        if r['cap']==2:lines.append(f"| {r['design']} | {r['divergence']} | {r['arrangement']} | {r['baseline_perfect_fraction']:.4f} | {r['baseline_mean_recall']:.5f} | {100*r['mean_normalized_net_gain']:.5f} | {r['recovered_sum']} / {r['lost_sum']} | {r['aggregate_cancellation_fraction']} |")
    lines+=['','Aggregate cancellation=2min(total recovered,total lost)/(total recovered+total lost), null if both zero. Totals are not comparable independent sample sizes. Original pilot cap2 recovery107 was the sum17+17+73 across BB12031 arrangements; BB12007 recovered0. Phase1 cap2/cap10 were identical in all300 datasets, pilot0/6, extension as reported above. Neither the pilot pattern nor agreement with simulations was an eligibility target.','',
      'The pilot pattern did not persist as a uniform ordering: the extension mean baseline recalls are about79.0% none,84.3% partial and85.0% high, whereas high overlap was hardest in both pilot families. None is tied for the lowest baseline recall in17 extension families, partial in8 and high in10 (ties count in each arrangement, so these totals exceed31). The extension mean cap2 gain is highest for none, rather than high as in the pilot. Different families, screened sites, unknown-endpoint treatment and the explicit spacing change prevent attributing that reversal to real-protein biology or overlap alone. Phase1/2 differences remain descriptive findings rather than evidence that one protocol estimates a pure causal effect.','',
      '## Reliability, limitations and reproduction','',
      'Twenty-two automated tests passed before any alignment, including survivor identity, ambiguity, exact overlap/length/core filtering, hand-worked scoring, marker/hash corruption, successful-output reuse after scoring failure, and a real process crash after two completed alignments. Immutable validation histories, a pre-execution validation record, source/config/site freeze, per-stage hashes, stderr, command/version/runtime, failure logs and resume records are retained. MAFFT7.526 is a local copy of the verified original distribution; only --maxiterate changes. All jobs use --amino --thread1 --retree2. Runtime includes launcher overhead and two-job contention. Independent post-run audits and the full actual resume check are reported separately and appended on completion.','',
      'The31 components are selected benchmark families, not a random sample of natural protein evolution. Shared-protein clustering can miss distant relationships. Reliable core windows favor conserved regions; ambiguous endpoints are omitted, not inferred. Scorable denominator and sequence count vary by family/arrangement. Reference reliability is imperfect; natural gaps are not proof of ancestral deletions. Local conservation/natural gaps are correlated features and descriptive regressions do not isolate mechanisms. Baseline-stratified subsets are secondary and may contain different families across arrangements. Statistical intervals are marginal with no multiple-testing discovery claim. Overlap remains confounded with spacing, context and union-gap geometry. No current manuscript or original experiment is modified.','',
      'See README.md for reproducibility, resume, validation and packaging. All new files remain in phase2_extension/. The original BAliBASE archive README contains no explicit redistribution license; preserve attribution and verify reuse terms. The report makes no novelty or independent causal-effect claim.']
    A.atomic(H/'phase2_extension_report.md','\n'.join(lines)+'\n')

def generate(rows,caps,estimates,summary,features,corr,reg,comparison,representatives):
    plots(rows,caps,estimates,summary,features,corr,reg,comparison)
    report(rows,caps,estimates,summary,features,corr,reg,comparison,representatives)
