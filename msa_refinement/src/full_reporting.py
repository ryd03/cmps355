"""Verified Phase 1 inventory, reports and figures; never writes Phase 0 artifacts."""
from collections import Counter
import csv
import io
import json
from pathlib import Path
import statistics
from . import pipeline, simulation, scoring, runner
from .storage import atomic, canonical, complete, digest, json_write, read_json, sha
from .reporting import csv_text
from .full_pipeline import batch_label, PROJECT


def inventory(root):
    root=Path(root)
    config=read_json(root/"config.json")
    manifest=read_json(root/"experiment_manifest.json")
    cfg_hash=digest(canonical(config))
    jobs=[]; rows=[]
    skips=Counter()
    for path in (root/"execution_events").glob("*.json"):
        obj=read_json(path)
        if obj["state"]=="skipped":
            skips[obj["stage"].replace("\\","/")]+=1
    for p in config["divergence_levels"]:
        batch=batch_label(p)
        for seed in config["simulator"]["replicate_seeds"]:
            block=root/"blocks"/f"{batch}_seed{seed}"
            block_ok=complete(block/"block_done.json",{"seed":seed,"p":p,"config":cfg_hash,"simulator":sha(simulation.__file__)})
            for cond in config["conditions"]:
                case=pipeline.case_path(root,batch,seed,cond)
                sim_dep={"config":cfg_hash,"seed":seed,"condition":cond,"p":p,"simulator":sha(simulation.__file__)}
                sim_ok=block_ok and complete(case/"simulation_done.json",sim_dep)
                meta=read_json(case/"ancestry.json") if sim_ok else None
                inputs=scoring.parse_fasta((case/"input.fasta").read_text()) if sim_ok else None
                h={name:sha(case/name) for name in ["input.fasta","truth_alignment.fasta","ancestry.json"]} if sim_ok else {}
                valid={}; outputs={}; scores={}; run_details={}
                for cap in [0,2,10]:
                    directory=case/"alignments"/f"maxiter_{cap:03d}"
                    output=directory/"aligned.fasta"
                    dep={"input":h.get("input.fasta"),"flags":["--amino","--thread","1","--retree","2","--maxiterate",str(cap)],
                         "mafft_version":manifest["mafft_version"],"runner":sha(runner.__file__)}
                    def validate():
                        scoring.validate_alignment(scoring.parse_fasta(output.read_text()),inputs)
                    valid[cap]=sim_ok and complete(directory/"run_done.json",dep,validate)
                    if valid[cap]:
                        outputs[cap]=sha(output)
                        run_details[cap]=read_json(directory/"run_done.json")["details"]
                        score_path=case/"scores"/f"maxiter_{cap:03d}.json"
                        score_dep={"aligned":outputs[cap],"truth":h["truth_alignment.fasta"],"ancestry":h["ancestry.json"],
                                   "scoring_version":scoring.SCORING_VERSION,"scoring_code":sha(scoring.__file__)}
                        if complete(score_path.with_name(score_path.name+".done.json"),score_dep):
                            scores[cap]=read_json(score_path)
                for cap in [0,2,10]:
                    directory=case/"alignments"/f"maxiter_{cap:03d}"
                    done=directory/"run_done.json"
                    state="completed" if valid[cap] else "failed" if (directory/"run_done.json.failed.json").exists() else "interrupted_or_invalid" if directory.exists() else "pending"
                    comp=None
                    if cap and cap in scores and 0 in scores:
                        comp_path=case/"comparisons"/f"against_maxiter_{cap:03d}.json"
                        dep={**scores[cap]["dependencies"],"baseline":outputs[0]}
                        if complete(comp_path.with_name(comp_path.name+".done.json"),dep):
                            comp=read_json(comp_path)
                    jobs.append({"batch":batch,"divergence":p,"seed":seed,"condition":cond,"iteration_cap":cap,
                                 "case":case.name,"alignment_state":state,"score_complete":cap in scores,
                                 "comparison_complete":True if cap==0 else comp is not None,
                                 "historical_alignment_failures":len(list((directory/"failures").glob("*.json"))),
                                 "historical_scoring_failures":len(list((case/"scores/failures").glob("*.json"))),
                                 "resume_skip_count":skips[str(done.relative_to(root)).replace("\\","/")],
                                 "path":directory.relative_to(root).as_posix()})
                    if cap not in scores:
                        continue
                    score=scores[cap]
                    r={"divergence":p,"batch":batch,"seed":seed,"condition":cond,"iteration_cap":cap,
                       "expected_root_tip_difference":meta["expected_root_to_tip_difference"],
                       "realized_overlap":meta["realized_overlap"],"union_gap_columns":meta["union_gap_columns"],
                       "deletion1_start":meta["events"][0]["start"],"deletion2_start":meta["events"][1]["start"],
                       "mean_pairwise_identity":meta["mean_pairwise_identity_shared_extant_ancestral_sites"],
                       "mean_root_tip_endpoint_differences":statistics.fmean(meta["root_tip_endpoint_differences"].values()),
                       "sum_branch_endpoint_changes":sum(e["endpoint_changes"] for e in meta["edges"]),
                       "surviving_residues":meta["total_surviving_residues"],"aligned_length":score["aligned_length"],
                       "runtime_seconds":run_details[cap]["runtime_seconds"],"input_sha256":h["input.fasta"],"aligned_sha256":outputs[cap],
                       "baseline_local_missed":scores.get(0,{}).get("local",{}).get("missed_pairs"),
                       "baseline_global_missed":scores.get(0,{}).get("global",{}).get("missed_pairs")}
                    for region in ["local","far","global"]:
                        r.update({region+"_"+key:value for key,value in score[region].items()})
                        if comp:
                            r.update({region+"_"+key:value for key,value in comp[region].items()})
                    if comp:
                        r["alignment_pair_set_changed"]=comp["alignment_pair_set_changed"]
                    rows.append(r)
    return jobs,rows


def cell_summaries(rows):
    cells=[]
    for p in [.25,.4]:
        for cond in ["none","partial","high"]:
            for cap in [0,2,10]:
                selected=[r for r in rows if r["divergence"]==p and r["condition"]==cond and r["iteration_cap"]==cap]
                if not selected:
                    continue
                row={"divergence":p,"condition":cond,"union_gap_columns":selected[0]["union_gap_columns"],"iteration_cap":cap,"n_seed_blocks":len(selected),
                     "baseline_local_perfect":sum(r["baseline_local_missed"]==0 for r in selected),
                     "baseline_global_perfect":sum(r["baseline_global_missed"]==0 for r in selected),
                     "mean_pairwise_identity":statistics.fmean(r["mean_pairwise_identity"] for r in selected)}
                for metric in ["local_sp_recall","global_sp_recall","local_precision","global_precision","runtime_seconds"]:
                    values=[r[metric] for r in selected if r[metric] is not None]
                    row["mean_"+metric]=statistics.fmean(values) if values else None
                    row["median_"+metric]=statistics.median(values) if values else None
                if cap and all("local_net_gain" in r for r in selected):
                    for key in ["local_recovered","local_lost","local_false_introduced","local_false_removed","global_recovered","global_lost"]:
                        row["total_"+key]=sum(r[key] for r in selected)
                        row["mean_"+key]=statistics.fmean(r[key] for r in selected)
                    row["mean_net_local_gain"]=statistics.fmean(r["local_net_gain"] for r in selected)
                    row["alignment_changed_datasets"]=sum(r["alignment_pair_set_changed"] for r in selected)
                    row["local_correct_pair_changed_datasets"]=sum(r["local_recovered"]+r["local_lost"]>0 for r in selected)
                    row["net_gain_positive_datasets"]=sum(r["local_net_gain"]>0 for r in selected)
                    row["net_gain_negative_datasets"]=sum(r["local_net_gain"]<0 for r in selected)
                    row["net_gain_zero_datasets"]=sum(r["local_net_gain"]==0 for r in selected)
                cells.append(row)
    return cells


def read_intervals(root):
    result=[]
    with (root/"confidence_intervals.csv").open(newline="",encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            for key in ["estimate","ci_lower","ci_upper","confidence_level"]:
                row[key]=float(row[key])
            for key in ["iteration_cap","n_seed_blocks","bootstrap_replicates","bootstrap_seed"]:
                row[key]=int(row[key])
            result.append(row)
    return result


def figures(root,rows,cells,intervals):
    import os
    os.environ.setdefault("MPLCONFIGDIR",str(PROJECT/"runs/phase1_plot_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors={"none":"#0072B2","partial":"#D55E00","high":"#009E73"}
    directory=root/"figures"
    directory.mkdir(exist_ok=True)
    def save(fig,name):
        fig.tight_layout()
        for ext in ["png","svg"]:
            buffer=io.BytesIO(); fig.savefig(buffer,format=ext,dpi=160)
            atomic(directory/(name+"."+ext),buffer.getvalue())
        plt.close(fig)
    for cap in [2,10]:
        fig,axes=plt.subplots(1,2,figsize=(12,4.5),sharey=True)
        for ax,p in zip(axes,[.25,.4]):
            for seed in range(2001,2051):
                subset=[next(r for r in rows if r["divergence"]==p and r["seed"]==seed and r["condition"]==c and r["iteration_cap"]==cap) for c in ["none","partial","high"]]
                ax.plot(range(3),[r["local_net_gain"] for r in subset],color="#888888",alpha=.22,lw=.7)
            for x,cond in enumerate(["none","partial","high"]):
                ci=next(r for r in intervals if r["family"]=="mean_net_gain" and float(r["divergence"])==p and r["iteration_cap"]==cap and r["contrast"]==cond)
                ax.plot(x,ci["estimate"],"o",color=colors[cond],ms=8)
                ax.vlines(x,ci["ci_lower"],ci["ci_upper"],color=colors[cond],lw=3)
            ax.axhline(0,color="black",lw=.8)
            ax.set_xticks(range(3),["none\nunion=24","partial\nunion=18","high\nunion=13"])
            ax.set(title=f"Expected root-to-tip difference {p:.2f}",ylabel="Normalized focal-local net gain")
        fig.suptitle(f"Cap {cap} versus cap 0: 50 paired seed blocks; means and marginal 95% CIs")
        save(fig,f"net_local_cap{cap}")
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for ax,cap in zip(axes,[2,10]):
        selected=[r for r in intervals if r["family"]=="geometry_contrast" and r["iteration_cap"]==cap]
        for y,r in enumerate(selected):
            ax.hlines(y,r["ci_lower"],r["ci_upper"],color="#0072B2",lw=2)
            ax.plot(r["estimate"],y,"o",color="#0072B2")
        ax.set_yticks(range(len(selected)),[f"p={float(r['divergence']):.2f}: {r['contrast'].replace('_',' ')}" for r in selected])
        ax.axvline(0,color="black",lw=.8)
        ax.set(title="Primary cap 2" if cap==2 else "Secondary cap 10",xlabel="Mean paired difference in normalized local net gain")
        ax.invert_yaxis()
    fig.suptitle("Seed-paired deletion geometry contrasts; marginal 95% bootstrap intervals")
    save(fig,"paired_geometry_intervals")
    for metric,ylabel,name in [("total_local_recovered","Correct pairs recovered","recovered_pairs"),("total_local_lost","Correct pairs lost","lost_pairs"),
                               ("alignment_changed_datasets","Datasets with changed pair sets / 50","alignment_changes"),
                               ("mean_local_sp_recall","Mean local SP recall","local_sp_recall"),
                               ("mean_global_sp_recall","Mean global SP recall","global_sp_recall"),
                               ("median_runtime_seconds","Median wall-clock seconds (4 concurrent blocks)","runtime")]:
        fig,axes=plt.subplots(1,2,figsize=(11,4.5))
        caps=[0,2,10] if "sp_recall" in metric or "runtime" in metric else [2,10]
        for ax,p in zip(axes,[.25,.4]):
            for i,cap in enumerate(caps):
                vals=[next(c for c in cells if c["divergence"]==p and c["condition"]==cond and c["iteration_cap"]==cap)[metric] for cond in ["none","partial","high"]]
                width=.24
                ax.bar([x+(i-(len(caps)-1)/2)*width for x in range(3)],vals,width,label=f"cap {cap}")
            ax.set_xticks(range(3),["none (24)","partial (18)","high (13)"])
            ax.set(title=f"Expected root-to-tip difference {p:.2f}",ylabel=ylabel,xlabel="Geometry (union gap columns)")
            ax.legend()
        fig.suptitle(ylabel)
        save(fig,name)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),sharey=True)
    for ax,p in zip(axes,[.25,.4]):
        s=[c for c in cells if c["divergence"]==p and c["iteration_cap"]==0]
        ax.bar([i-.17 for i in range(3)],[c["baseline_local_perfect"] for c in s],.34,label="Local-perfect")
        ax.bar([i+.17 for i in range(3)],[c["baseline_global_perfect"] for c in s],.34,label="Global-perfect")
        ax.set_xticks(range(3),["none (24)","partial (18)","high (13)"])
        ax.set(title=f"Expected root-to-tip difference {p:.2f}",ylabel="Baseline alignments already perfect / 50",ylim=(0,50))
        ax.legend()
    fig.suptitle("Ceiling effects at cap 0; separate local and global definitions")
    save(fig,"ceiling_effects")


def report_text(root,config,jobs,rows,cells,intervals):
    manifest=read_json(root/"experiment_manifest.json")
    text=["# Phase 1 full research report","",
          "This report contains only the 50 new evolutionary seed blocks 2001-2050. Both divergence levels use these same matched seed IDs. "
          "The independent experimental unit is the seed block; amino-acid pairs are not independent samples. Phase 0 observations are excluded.","",
          "## Frozen design and execution","",
          f"Experiment `{root.name}`; configuration SHA256 `{manifest['frozen_config_file_sha256']}`; frozen UTC epoch {manifest['frozen_utc_epoch']}. "
          f"Python {manifest['python']}, {manifest['platform']}, MAFFT `{manifest['mafft_version']}`.","",
          "The full 12-tip fixed-tree v2 simulation, uniform 240-residue root, symmetric 20-state endpoint CTMC, "
          "independent terminal length-12 deletions in protein1/protein7, anchor jitter, offsets 50/6/1, ancestry-based local window and scoring are unchanged. "
          "p=.25 and p=.40 are expected root-to-tip differences, not pairwise divergence. Every branch is evolved. "
          "Each seed/divergence realization is generated once and shared across conditions. There is no adaptive calibration in Phase 1.","",
          "Every dataset is run with `--amino --thread 1 --retree 2 --maxiterate 0|2|10`. "
          "Four seed blocks may be orchestrated concurrently, each with one active single-threaded MAFFT call. "
          "Elapsed times include launcher overhead and possible resource contention; they do not establish a pure algorithmic speed effect. "
          "Iteration caps do not prove actual iteration counts. No internal WSP measurement is claimed.","",
          f"Planned 300 datasets / 900 alignments. Validated alignments: {sum(j['alignment_state']=='completed' for j in jobs)}; "
          f"currently failed: {sum(j['alignment_state']=='failed' for j in jobs)}; "
          f"uncompleted: {sum(j['alignment_state']!='completed' for j in jobs)}; score rows: {len(rows)}; "
          f"historical alignment failures: {sum(j['historical_alignment_failures'] for j in jobs)}. "
          "All planned jobs, failures, comparison completion and skip counts are in `full_job_manifest.csv`; invocation and raw stderr logs are retained.","",
          "## Primary analysis: cap 2 versus cap 0","",
          "For each seed and condition, normalize `(recovered_local - lost_local)` by that condition's true-local reference-pair denominator, "
          "then compute partial-minus-none and high-minus-none contrasts. The estimand is the mean of 50 within-seed differences. "
          "Bootstrap 10,000 samples of 50 seed IDs with replacement, jointly preserving all conditions, divergences and caps; "
          "fixed RNG seed 3552050. Intervals are percentile 95% intervals with linear quantile interpolation. "
          "These are marginal intervals for four primary contrasts, not simultaneous familywise intervals. No p-values or significance labels are used.","",
          "| Divergence | Contrast | Mean paired difference | 95% CI lower | 95% CI upper | Seed blocks |",
          "|---|---|---:|---:|---:|---:|"]
    for r in intervals:
        if r["analysis_role"]=="primary":
            text.append(f"| {r['divergence']} | {r['contrast']} | {r['estimate']:.8f} | {r['ci_lower']:.8f} | {r['ci_upper']:.8f} | {r['n_seed_blocks']} |")
    text.append("")
    for r in intervals:
        if r["analysis_role"]=="primary":
            location="spans zero" if r["ci_lower"]<=0<=r["ci_upper"] else "lies above zero" if r["ci_lower"]>0 else "lies below zero"
            text.append(f"At p={r['divergence']}, the {r['contrast']} estimate is {r['estimate']:.6f} "
                        f"({100*r['estimate']:.3f} percentage points); its marginal interval {location}. "
                        "This describes uncertainty in the paired simulation contrast, not an independent causal effect of overlap.")
            text.append("")
    text.extend(["","A positive contrast means refinement's normalized correct-pair net change was greater for the indicated geometry than for no overlap. "
                 "It does not imply that refinement improved every alignment, that the absolute gains were large, or that overlap is independently causal.","",
                 "## Individual gains, losses, accuracy and ceiling effects","",
                 "`full_summary.csv` has one row per dataset/mode; `paired_comparisons.csv` preserves each seed's two primary contrasts and separate recovered/lost counts and denominators. "
                 "Raw counts below are summed descriptively over 50 datasets within a cell, never used as independent statistical samples. "
                 "The same pair may be counted at both caps; do not sum across caps as unique events.","",
                 "| p | Geometry (union) | cap | changed /50 | recovered local | lost local | mean net local | mean local SP | mean global SP | median seconds | baseline local-perfect | baseline global-perfect |",
                 "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"])
    for c in cells:
        def fmt(key):
            v=c.get(key)
            return "" if v is None else f"{v:.6f}" if isinstance(v,float) else str(v)
        keys=["divergence","condition","iteration_cap","alignment_changed_datasets","total_local_recovered","total_local_lost","mean_net_local_gain","mean_local_sp_recall","mean_global_sp_recall","median_runtime_seconds","baseline_local_perfect","baseline_global_perfect"]
        vals=[fmt(k) for k in keys]; vals[1]+=f" ({c['union_gap_columns']})"
        text.append("| "+" | ".join(vals)+" |")
    text.extend(["","A local-perfect baseline has zero missed focal-local reference pairs. A globally perfect baseline has zero missed global reference pairs. "
                 "They are separate ceiling definitions; local-perfect does not imply global-perfect. Undefined rates remain null (blank in CSV), not zero. "
                 "Global SP recall is secondary; its high value can mask local displacement errors.","",
                 "## Secondary cap 10 and divergence comparisons","",
                 "| Family | p / difference | cap | Contrast | Estimate | 95% CI lower | 95% CI upper |",
                 "|---|---|---:|---|---:|---:|---:|"])
    for r in intervals:
        if r["family"]!="mean_net_gain" and r["analysis_role"]!="primary":
            text.append(f"| {r['family']} | {r['divergence']} | {r['iteration_cap']} | {r['contrast']} | {r['estimate']:.8f} | {r['ci_lower']:.8f} | {r['ci_upper']:.8f} |")
    text.extend(["","Divergence differences are paired 0.40 minus 0.25 comparisons by seed, retaining the shared root/anchor and unchanged named RNG streams. "
                 "Geometry interaction rows compare the geometry contrast at .40 with that at .25. "
                 "All secondary intervals are descriptive marginal intervals without multiplicity adjustment. "
                 "`confidence_intervals.csv` also contains per-condition mean net-gain intervals; `divergence_paired_comparisons.csv` lists every seed-level divergence difference.","",
                 "## Validation, recovery, preservation and reproduction","",
                 "The unmodified v2 scientific and seven-scenario recovery suite plus Phase 1-specific schedule, normalization, bootstrap, incomplete-data, concurrency and crash-restart tests ran before execution. "
                 "See `validation_before_execution.json` and the bundled phase1_validation logs/JUnit files. "
                 "The separate final audit regenerates all datasets and recomputes all scores and comparisons, verifies exact flags and the frozen config, checks unique rows, and reproduces the CI analysis. "
                 "See `full_artifact_audit.json` when present. `full_resume_verification.json` records the completed-run resume check. "
                 "`phase0_preservation.json` verifies unchanged bytes and timestamps against the saved Phase 0 baseline.","",
                 "All stage publications use the original atomic-write/hash-marker implementation. Successful raw MAFFT outputs are separate from scoring. "
                 "An interrupted unfinished call reruns; successful calls are reused. Failed jobs retain stderr, partial stdout and traceback, and can be explicitly retried. "
                 "Statistical estimates require all 900 validated outputs and all comparisons, so no failed seed is silently dropped.","",
                 "```", "python -m src.full_cli validate", "python -m src.full_cli freeze --config config/full.json --workdir ./runs",
                 "python -m src.full_cli run --config config/full.json --workdir ./runs",
                 "python -m src.full_cli status --workdir ./runs",
                 "python -m src.full_cli run --config config/full.json --workdir ./runs",
                 "python -m src.full_cli audit --workdir ./runs", "python -m src.full_cli report --workdir ./runs",
                 "python -m src.full_cli bundle --workdir ./runs", "```","",
                 "Detailed installation, resume and archive instructions are in `PHASE1_REPRODUCE.md`. Exact per-call commands, versions, seeds and input/output hashes are stored alongside every alignment.","",
                 "## Limitations and interpretation","",
                 "The treatment is a combined deletion arrangement. Overlap is confounded with deletion spacing, the union of affected reference columns (24/18/13), "
                 "local sequence context, and which site identities were removed. An independent causal effect of overlap is not identified. "
                 "Confidence intervals quantify variation across 50 evolutionary seed blocks in this simulation, not uncertainty across natural protein families. "
                 "There is one fixed tree, a uniformly distributed amino-acid ancestor and a simple equal-rate substitution model, "
                 "two independent terminal deletion events, and no insertions. MAFFT guide-tree differences may mediate alignment changes. "
                 "No mechanism or novelty is proved, no unadjusted multiple-testing significance claim is made, and runtime comparisons are limited by the execution environment.","",
                 "The primary cap, contrasts, seed range, divergence levels and CI method were fixed before the first Phase 1 alignment. "
                 "The pilot is a separate completed experiment and contributes no observations to these estimates. "
                 "Phase 1 is the stopping point; no insertion experiments, added algorithms or extra seed batches are authorized.",""])
    text.extend(["## Actual verification evidence",""])
    for label,path in [("Before execution",root/"validation_before_execution.json"),("Latest validation suite",PROJECT/"runs/phase1_validation/latest.json")]:
        if path.exists():
            evidence=read_json(path)
            text.append(f"{label}: exit code {evidence['exit_code']}; {evidence.get('summary','see validation log')}.")
            text.append("")
    if (root/"full_artifact_audit.json").exists():
        evidence=read_json(root/"full_artifact_audit.json")
        text.extend([f"Independent audit passed: {evidence['datasets_regenerated']} datasets regenerated, "
                     f"{evidence['scores_recomputed']} scores and {evidence['comparisons_recomputed']} comparisons recomputed, "
                     f"{evidence['statistical_estimates_reproduced']} statistical estimates reproduced. "
                     f"Cap 2 and cap 10 had identical residue-pair sets in {evidence['cap2_cap10_identical_pair_sets']}/300 datasets.",""])
    if (root/"full_resume_verification.json").exists():
        evidence=read_json(root/"full_resume_verification.json")
        text.extend([f"Actual completed-run resume passed: {evidence['alignment_calls_made']} new alignment calls; "
                     f"{evidence['alignments_skipped']} alignments skipped; {evidence['artifacts_checked']} stage artifacts retained identical bytes and timestamps; summary CSV unchanged.",""])
    if (root/"phase0_preservation.json").exists():
        evidence=read_json(root/"phase0_preservation.json")
        text.extend([f"Phase 0 preservation: passed={evidence['passed']}; {evidence['phase0_files_checked']} original files checked against their saved hashes and timestamps.",""])
    return "\n".join(text)


def build(root):
    root=Path(root)
    config=read_json(root/"config.json")
    jobs,rows=inventory(root)
    atomic(root/"full_summary.csv",csv_text(rows))
    atomic(root/"full_job_manifest.csv",csv_text(jobs))
    json_write(root/"full_job_manifest.json",{"planned":len(jobs),"jobs":jobs,"completed":sum(j["alignment_state"]=="completed" for j in jobs),
               "failed":sum(j["alignment_state"]=="failed" for j in jobs),"score_complete":sum(j["score_complete"] for j in jobs),
               "comparison_complete":sum(j["comparison_complete"] for j in jobs)})
    ready=(len(jobs)==900 and len(rows)==900 and all(j["alignment_state"]=="completed" and j["score_complete"] and j["comparison_complete"] for j in jobs))
    if not ready:
        atomic(root/"full_report.md","# Phase 1 execution incomplete\n\nStatistical estimates are withheld until all 900 jobs and comparisons validate. See full_job_manifest.csv for every planned/failed/pending job.\n")
        return jobs,rows
    from . import full_analysis
    products=[root/name for name in ["paired_comparisons.csv","divergence_paired_comparisons.csv","confidence_intervals.csv","bootstrap_resampling.json.gz","bootstrap_estimates.csv.gz","analysis_done.json"]]
    deps={"summary":sha(root/"full_summary.csv"),"analysis_settings":config["analysis"],"analysis_source":sha(full_analysis.__file__)}
    pipeline.stage(root,root,"analysis_products_done.json",deps,products,lambda: full_analysis.analyze(rows,config,root) and {},
                   lambda: read_json(root/"analysis_done.json")["seed_blocks"],True,Counter())
    intervals=read_intervals(root)
    cells=cell_summaries(rows)
    atomic(root/"cell_summary.csv",csv_text(cells))
    figures(root,rows,cells,intervals)
    atomic(root/"full_report.md",report_text(root,config,jobs,rows,cells,intervals))
    return jobs,rows
