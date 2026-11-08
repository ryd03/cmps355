"""Read-only exploratory supplement: existing Phase 1 counts, distributions and LOO."""
import argparse
from collections import Counter
import csv
from fractions import Fraction
import hashlib
import io
import json
import math
import os
from pathlib import Path
import statistics
import sys
import time
import uuid

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parent
PROJECT = WORKSPACE / "msa_refinement"
DEFAULT_SOURCE = PROJECT / "runs/full_deletion_overlap_v2_0815b46e6b9208a2"
CONDITIONS = ["none", "partial", "high"]
LEVELS = [.25, .40]
SEEDS = list(range(2001, 2051))
PROCESS_ATTEMPTS = []
WRITE_VIOLATIONS = []
WRITE_FLAGS = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND


def inside_output(path):
    if isinstance(path, int):
        return True
    return Path(os.fsdecode(path)).resolve().is_relative_to(HERE)


def readonly_guard(event, args):
    if event in {"subprocess.Popen", "os.system", "os.posix_spawn", "os.spawn"}:
        PROCESS_ATTEMPTS.append({"event": event, "arguments": repr(args)})
        raise PermissionError("External process execution is forbidden in this supplementary analysis")
    targets = []
    if event == "open" and len(args) >= 3 and (args[2] or 0) & WRITE_FLAGS:
        targets = [args[0]]
    elif event in {"os.mkdir", "os.remove", "os.rmdir", "os.chmod", "os.utime", "os.truncate"}:
        targets = [args[0]]
    elif event in {"os.rename", "os.link", "os.symlink"}:
        targets = list(args[:2])
    for path in targets:
        if not inside_output(path):
            WRITE_VIOLATIONS.append({"event": event, "path": str(path)})
            raise PermissionError("Original files are read-only: " + str(path))


def install_guard():
    sys.addaudithook(readonly_guard)
    for name in ["cache/matplotlib", "tmp", "figures", "tables", "logs"]:
        (HERE / name).mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(HERE / "cache/matplotlib")
    os.environ["XDG_CACHE_HOME"] = str(HERE / "cache")
    for key in ["TMP", "TEMP", "TMPDIR"]:
        os.environ[key] = str(HERE / "tmp")
    if (PROJECT / ".deps").exists():
        sys.path.insert(0, str(PROJECT / ".deps"))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def output_path(name):
    path = (HERE / name).resolve()
    if not path.is_relative_to(HERE):
        raise PermissionError("All new files must be inside supplementary_analysis_phase1")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def atomic(name, data):
    path = output_path(name)
    temp = path.with_name(uuid.uuid4().hex[:16] + ".partial")
    with temp.open("wb") as handle:
        handle.write(data.encode("utf-8") if isinstance(data, str) else data)
        handle.flush(); os.fsync(handle.fileno())
    os.replace(temp, path)


def save_json(name, data):
    atomic(name, json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def csv_text(rows):
    if not rows:
        return ""
    fields = list(dict.fromkeys(k for r in rows for k in r if not k.startswith("_")))
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fields, lineterminator="\n", extrasaction="ignore")
    writer.writeheader(); writer.writerows(rows)
    return buffer.getvalue()


def read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def close(a, b):
    assert math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12), (a, b)


def quantile(values, probability):
    values = sorted(values)
    if not values:
        raise ValueError("Empty distribution")
    rank = (len(values)-1)*probability
    lo, hi = math.floor(rank), math.ceil(rank)
    return float(values[lo]) + (float(values[hi])-float(values[lo]))*(rank-lo)


def mean_exact(values):
    return sum(values, Fraction(0)) / len(values)


def sign(value):
    return (value > 0) - (value < 0)


def concentration(values):
    if any(v < 0 for v in values):
        raise ValueError("Concentration requires nonnegative contributions")
    ranked = sorted(values, reverse=True)
    total = sum(ranked)
    result = {"total":total, "nonzero_datasets":sum(v>0 for v in values)}
    for k in [1,5,10]:
        result[f"top_{k}_pairs"] = sum(ranked[:k])
        result[f"top_{k}_share"] = sum(ranked[:k])/total if total else None
    for pct in [50,80]:
        cumulative = 0
        n = None
        if total:
            for i, value in enumerate(ranked, 1):
                cumulative += value
                if 100*cumulative >= pct*total:
                    n = i; break
        result[f"datasets_for_{pct}_percent"] = n
    return result


def require_baseline():
    path = HERE / "preservation_before.json"
    if not path.exists():
        raise ValueError("A complete preservation baseline must be recorded before analysis")
    baseline = read_json(path)
    assert baseline["complete"] and not baseline["errors"]
    return baseline


def check_source_file(path, baseline):
    path = Path(path).resolve()
    name = path.relative_to(WORKSPACE).as_posix()
    expected = baseline["files"].get(name)
    assert expected is not None, "Input not covered by preservation baseline: " + name
    assert sha(path) == expected["sha256"], "Original input differs from baseline: " + name
    return {"path":name, "sha256":expected["sha256"], "size":expected["size"]}


def marker_check(path, data, baseline):
    check_source_file(path, baseline)
    marker = read_json(path)
    assert marker["dependencies"] == read_json(data)["dependencies"]
    assert marker["files"][Path(data).name] == sha(data)


def load_data(source, baseline):
    source = Path(source).resolve()
    assert not source.is_relative_to(HERE), "Original Phase 1 data must remain outside the supplementary directory"
    named = ["config.json", "experiment_manifest.json", "full_summary.csv", "paired_comparisons.csv",
             "confidence_intervals.csv", "cell_summary.csv", "full_artifact_audit.json"]
    inputs = [check_source_file(source/name, baseline) for name in named]
    config = read_json(source/"config.json")
    assert config["phase"] == "phase1"
    assert config["simulator"]["replicate_seeds"] == SEEDS and config["divergence_levels"] == LEVELS
    assert config["conditions"] == CONDITIONS
    audit = read_json(source/"full_artifact_audit.json")
    assert audit["passed"] and audit["alignment_jobs_verified"] == 900
    raw = read_csv(source/"full_summary.csv")
    indexed = {(float(r["divergence"]),int(r["seed"]),r["condition"],int(r["iteration_cap"])):r for r in raw}
    expected = {(p,s,c,k) for p in LEVELS for s in SEEDS for c in CONDITIONS for k in [0,2,10]}
    assert len(raw) == len(indexed) == 900 and set(indexed) == expected
    rows = []; stage_markers_checked = 0
    for p in LEVELS:
        for seed in SEEDS:
            for condition in CONDITIONS:
                b, r = indexed[p,seed,condition,0], indexed[p,seed,condition,2]
                case = source/"cases"/f"{r['batch']}_seed{seed}_{condition}"
                bpath, rpath = case/"scores/maxiter_000.json", case/"scores/maxiter_002.json"
                cpath = case/"comparisons/against_maxiter_002.json"
                for path in [bpath,rpath,cpath]:
                    check_source_file(path, baseline)
                    marker_check(path.with_name(path.name+".done.json"),path,baseline)
                    stage_markers_checked += 1
                bs,rs,comparison = read_json(bpath),read_json(rpath),read_json(cpath)
                local = comparison["local"]
                recovered,lost,denom = int(r["local_recovered"]),int(r["local_lost"]),int(r["local_reference_denominator"])
                assert (recovered,lost,denom) == (local["recovered"],local["lost"],local["reference_denominator"])
                assert bs["local"]["reference_pairs"] == rs["local"]["reference_pairs"] == denom > 0
                assert rs["local"]["correct_pairs"] - bs["local"]["correct_pairs"] == recovered-lost
                assert rs["local"]["false_pairs"] - bs["local"]["false_pairs"] == local["false_introduced"]-local["false_removed"]
                for score,csv_row in [(bs,b),(rs,r)]:
                    for region in ["local","far","global"]:
                        for metric in ["correct_pairs","false_pairs","missed_pairs","reference_pairs","predicted_pairs"]:
                            assert score[region][metric] == int(csv_row[region+"_"+metric])
                        close(score[region]["sp_recall"], csv_row[region+"_sp_recall"])
                    for metric in ["correct_pairs","false_pairs","missed_pairs","reference_pairs","predicted_pairs"]:
                        assert score["global"][metric] == score["local"][metric]+score["far"][metric]
                delta = Fraction(recovered-lost,denom)
                close(delta,r["local_net_gain"])
                close(delta,rs["local"]["sp_recall"]-bs["local"]["sp_recall"])
                close(local["net_gain"],delta)
                close(Fraction(bs["local"]["correct_pairs"],denom),bs["local"]["sp_recall"])
                close(Fraction(rs["local"]["correct_pairs"],denom),rs["local"]["sp_recall"])
                # Link saved score dependencies to the baseline hashes of raw inputs.
                for cap,score in [(0,bs),(2,rs)]:
                    for dep,name in [("aligned",f"alignments/maxiter_{cap:03d}/aligned.fasta"),("truth","truth_alignment.fasta"),("ancestry","ancestry.json")]:
                        rel=(case/name).relative_to(WORKSPACE).as_posix()
                        assert score["dependencies"][dep] == baseline["files"][rel]["sha256"]
                net = recovered-lost
                rows.append({"divergence":p,"seed":seed,"condition":condition,"iteration_cap":2,
                    "comparison":"cap_2_minus_cap_0","n_local_reference_pairs":denom,
                    "local_accuracy_cap0":bs["local"]["sp_recall"],"local_accuracy_cap2":rs["local"]["sp_recall"],
                    "local_accuracy_change":float(delta),"local_accuracy_change_percentage_points":100*float(delta),
                    "recovered_local":recovered,"lost_local":lost,"net_correct_local":net,
                    "positive_net_pairs":max(net,0),"negative_net_pairs":max(-net,0),
                    "classification":"improved" if net>0 else "worsened" if net<0 else "unchanged",
                    "baseline_local_recall_ceiling":bs["local"]["missed_pairs"]==0,
                    "alignment_pair_set_changed":comparison["alignment_pair_set_changed"],
                    "union_gap_columns":int(r["union_gap_columns"]),"_exact_change":delta})
    save_json("source_manifest.json",{"exploratory_post_hoc":True,"experiment_directory":str(source),"inputs":inputs,
        "dataset_rows":len(rows),"source_summary_rows":len(raw),"validated_score_and_comparison_markers":stage_markers_checked,
        "original_scoring_definition":"Focal-boundary reference-pair SP recall; integer correct-pair net divided by original true-local pair denominator",
        "original_data_copied_into_bundle":False})
    return rows,indexed


def descriptive_tables(rows):
    frequency,distribution,shares,ranks = [],[],[],[]
    for p in LEVELS:
        for cond in CONDITIONS:
            group = [r for r in rows if r["divergence"]==p and r["condition"]==cond]
            assert len(group)==50 and {r["seed"] for r in group}==set(SEEDS)
            counts=Counter(r["classification"] for r in group)
            changes=[r["_exact_change"] for r in group]
            f={"divergence":p,"condition":cond,"independent_seed_blocks":len(group),"union_gap_columns":group[0]["union_gap_columns"]}
            for label in ["improved","worsened","unchanged"]:
                f[label+"_datasets"]=counts[label]; f[label+"_percent"]=100*counts[label]/len(group)
            f.update({"mean_local_accuracy_change":float(mean_exact(changes)),"median_local_accuracy_change":statistics.median(changes),
                      "mean_change_percentage_points":100*float(mean_exact(changes)),
                      "median_change_percentage_points":100*float(statistics.median(changes)),
                      "alignment_changed_but_net_accuracy_unchanged":sum(r["alignment_pair_set_changed"] and r["net_correct_local"]==0 for r in group)})
            f["median_local_accuracy_change"]=float(f["median_local_accuracy_change"])
            frequency.append(f)
            for metric in ["recovered_local","lost_local","net_correct_local","local_accuracy_change"]:
                values=[r[metric] for r in group]
                q1,q3=quantile(values,.25),quantile(values,.75)
                distribution.append({"divergence":p,"condition":cond,"metric":metric,"n_seed_blocks":len(group),
                    "total":sum(values),"mean":statistics.fmean(values),"median":statistics.median(values),
                    "q1":q1,"q3":q3,"iqr":q3-q1,"min":min(values),"max":max(values),"zero_datasets":sum(v==0 for v in values)})
            for metric in ["recovered_local","lost_local","positive_net_pairs","negative_net_pairs"]:
                ordered=sorted(group,key=lambda r:(-r[metric],r["seed"]))
                values=[r[metric] for r in ordered]
                info=concentration(values)
                shares.append({"divergence":p,"condition":cond,"metric":metric,"n_seed_blocks":len(group),**info,
                    "top_1_seed":ordered[0]["seed"] if info["total"] else None})
                cumulative=0
                for rank,r in enumerate(ordered,1):
                    cumulative+=r[metric]
                    ranks.append({"divergence":p,"condition":cond,"metric":metric,"rank":rank,"seed":r["seed"],
                        "contribution_pairs":r[metric],"cumulative_pairs":cumulative,
                        "cumulative_share":cumulative/info["total"] if info["total"] else None})
    return frequency,distribution,shares,ranks


def leave_one_out(rows,source):
    indexed={(r["divergence"],r["seed"],r["condition"]):r for r in rows}
    original_paired={(float(r["divergence"]),int(r["seed"]),r["contrast"]):r for r in read_csv(source/"paired_comparisons.csv") if int(r["iteration_cap"])==2}
    original_cis={(float(r["divergence"]),r["contrast"]):r for r in read_csv(source/"confidence_intervals.csv") if r["analysis_role"]=="primary"}
    assert len(original_paired)==200 and len(original_cis)==4
    pairs,loo,summary=[],[],[]
    for p in LEVELS:
        for cond in ["partial","high"]:
            label=cond+"_minus_none"
            vector={s:indexed[p,s,cond]["_exact_change"]-indexed[p,s,"none"]["_exact_change"] for s in SEEDS}
            full=mean_exact(list(vector.values()))
            reference=original_cis[p,label]
            close(full,reference["estimate"])
            for s in SEEDS:
                original=original_paired[p,s,label]
                close(vector[s],original["paired_difference"])
                pair={"divergence":p,"seed":s,"contrast":label,"none_gain":float(indexed[p,s,"none"]["_exact_change"]),
                    "treatment_gain":float(indexed[p,s,cond]["_exact_change"]),"paired_difference":float(vector[s]),
                    "original_paired_difference":original["paired_difference"]}
                for prefix,case in [("none",indexed[p,s,"none"]),("treatment",indexed[p,s,cond])]:
                    for metric in ["recovered_local","lost_local","net_correct_local","n_local_reference_pairs"]:
                        pair[prefix+"_"+metric]=case[metric]
                pairs.append(pair)
                kept=[value for key,value in vector.items() if key!=s]
                estimate=mean_exact(kept)
                shift=estimate-full
                assert shift==(full-vector[s])/49
                loo.append({"divergence":p,"contrast":label,"omitted_seed":s,"retained_seed_blocks":49,
                    "original_full_estimate":reference["estimate"],"recomputed_full_estimate":float(full),
                    "omitted_seed_contrast":float(vector[s]),"leave_one_out_estimate":float(estimate),
                    "change_from_full_estimate":float(shift),"change_percentage_points":100*float(shift),
                    "absolute_change":abs(float(shift)),"direction_changed":sign(estimate)!=sign(full),
                    "original_ci_lower_reference_only":reference["ci_lower"],"original_ci_upper_reference_only":reference["ci_upper"],
                    "loo_ci_recomputed":False,"_exact_estimate":estimate})
            selected=[r for r in loo if r["divergence"]==p and r["contrast"]==label]
            influential=max(selected,key=lambda r:r["absolute_change"])
            summary.append({"divergence":p,"contrast":label,"full_seed_blocks":50,"loo_seed_blocks":49,
                "original_full_estimate":reference["estimate"],"original_ci_lower_reference_only":reference["ci_lower"],
                "original_ci_upper_reference_only":reference["ci_upper"],
                "min_leave_one_out_estimate":min(r["leave_one_out_estimate"] for r in selected),
                "max_leave_one_out_estimate":max(r["leave_one_out_estimate"] for r in selected),
                "max_absolute_change":influential["absolute_change"],"max_absolute_change_percentage_points":100*influential["absolute_change"],
                "most_influential_seed":influential["omitted_seed"],
                "direction_changes_out_of_50":sum(r["direction_changed"] for r in selected),
                "loo_ci_recomputed":False})
    return pairs,loo,summary


def unit_and_data_validation(rows,frequency,distribution,shares,ranks,loo,source):
    # Hand-worked tests, including null totals and exact pairing influence.
    close(quantile([0,1,2,3],.25),.75); close(quantile([0,1,2,3],.75),2.25)
    assert concentration([10,5,5,0])["top_1_share"]==.5
    assert concentration([10,5,5,0])["datasets_for_80_percent"]==3
    assert concentration([0,0])["top_1_share"] is None
    assert concentration([0,0])["datasets_for_50_percent"] is None
    try:
        concentration([-1,2]); raise AssertionError("Negative contribution accepted")
    except ValueError:
        pass
    toy=[Fraction(-2),Fraction(1),Fraction(4)]
    assert mean_exact(toy)==1
    assert [mean_exact([v for j,v in enumerate(toy) if i!=j]) for i in range(3)]==[Fraction(5,2),Fraction(1),Fraction(-1,2)]
    assert len(rows)==300 and len(loo)==200 and len(ranks)==1200
    # Independent array-library checks of all descriptive quantities.
    import numpy as np
    for f in frequency:
        group=[r for r in rows if r["divergence"]==f["divergence"] and r["condition"]==f["condition"]]
        net=np.array([r["net_correct_local"] for r in group])
        delta=np.array([r["local_accuracy_change"] for r in group])
        for label,condition in [("improved",net>0),("worsened",net<0),("unchanged",net==0)]:
            assert f[label+"_datasets"]==int(np.count_nonzero(condition))
            close(f[label+"_percent"],100*np.mean(condition))
        assert sum(f[label+"_datasets"] for label in ["improved","worsened","unchanged"])==50
        close(f["mean_local_accuracy_change"],np.mean(delta)); close(f["median_local_accuracy_change"],np.median(delta))
    for d in distribution:
        vals=np.array([r[d["metric"]] for r in rows if r["divergence"]==d["divergence"] and r["condition"]==d["condition"]])
        q1,q3=np.quantile(vals,[.25,.75],method="linear")
        close(d["q1"],q1); close(d["q3"],q3); close(d["iqr"],q3-q1)
        close(d["mean"],np.mean(vals)); close(d["median"],np.median(vals))
    for c in shares:
        vals=sorted([r[c["metric"]] for r in rows if r["divergence"]==c["divergence"] and r["condition"]==c["condition"]],reverse=True)
        assert c["total"]==sum(vals)
        if c["total"]:
            close(c["top_5_share"],sum(vals[:5])/sum(vals))
            selected=[r for r in ranks if r["divergence"]==c["divergence"] and r["condition"]==c["condition"] and r["metric"]==c["metric"]]
            close(selected[-1]["cumulative_share"],1)
    original_cells=[r for r in read_csv(source/"cell_summary.csv") if r["iteration_cap"]=="2"]
    for c in original_cells:
        group=[r for r in rows if r["divergence"]==float(c["divergence"]) and r["condition"]==c["condition"]]
        assert sum(r["recovered_local"] for r in group)==int(c["total_local_recovered"])
        assert sum(r["lost_local"] for r in group)==int(c["total_local_lost"])
        close(mean_exact([r["_exact_change"] for r in group]),c["mean_net_local_gain"])
    assert not PROCESS_ATTEMPTS and not WRITE_VIOLATIONS
    assert not any(name in sys.modules for name in ["src.runner","src.simulation","src.pipeline","src.full_pipeline"])
    return {"passed":True,"dataset_comparisons_validated":300,"original_score_jsons_checked":600,
        "original_comparison_jsons_checked":300,"completion_markers_checked":900,"descriptive_groups":6,
        "distribution_summaries_checked":24,"concentration_summaries_checked":24,
        "original_primary_seed_contrasts_matched":200,"original_primary_estimates_matched":4,
        "leave_one_seed_out_estimates_checked":200,"quartiles_independently_checked_with_numpy":True,
        "mafft_calls_executed":0,"all_external_process_execution_blocked":True,"external_process_attempts":PROCESS_ATTEMPTS,
        "original_write_attempts_blocked":WRITE_VIOLATIONS,"original_simulator_or_runner_imported":False,
        "new_simulations_or_reference_alignments_generated":0,"original_confidence_intervals_recomputed_or_modified":False,
        "numpy_version":np.__version__,"python_version":sys.version,"exploratory_post_hoc":True}


def make_figures(rows,frequency,ranks,loo,loo_summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":300})
    colors={"none":"#0072B2","partial":"#D55E00","high":"#009E73"}
    positions={s:(i-24.5)/600 for i,s in enumerate(SEEDS)}
    def save(fig,name):
        fig.tight_layout(rect=(0,0,1,.95))
        for ext in ["png","svg"]:
            buffer=io.BytesIO(); fig.savefig(buffer,format=ext,dpi=300)
            atomic(f"figures/{name}.{ext}",buffer.getvalue())
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),sharey=True)
    for ax,p in zip(axes,LEVELS):
        bottom=[0.,0.,0.]
        for label,color in [("improved","#0072B2"),("worsened","#D55E00"),("unchanged","#B8B8B8")]:
            vals=[next(f for f in frequency if f["divergence"]==p and f["condition"]==c)[label+"_percent"] for c in CONDITIONS]
            ax.bar(range(3),vals,bottom=bottom,label=label.capitalize(),color=color,width=.65)
            for i,v in enumerate(vals):
                if v>=8:
                    ax.text(i,bottom[i]+v/2,f"{v:.0f}%\n({round(v/2)}/50)",ha="center",va="center",fontsize=9)
            bottom=[a+b for a,b in zip(bottom,vals)]
        ax.set(title=f"Expected root-to-tip difference {p:.2f}; n=50 per group",ylabel="Datasets (%)",ylim=(0,100))
        ax.set_xticks(range(3),["none (24)","partial (18)","high (13)"])
        ax.set_xlabel("Deletion geometry (union gap columns)")
    axes[1].legend(loc="upper center",bbox_to_anchor=(.5,-.16),ncol=3,frameon=False)
    fig.suptitle("Exploratory/post-hoc: how often cap 2 changes focal-local accuracy versus cap 0")
    save(fig,"accuracy_outcome_frequencies")
    def points_box(ax,p,metric,scale=1):
        groups=[[scale*r[metric] for r in rows if r["divergence"]==p and r["condition"]==c] for c in CONDITIONS]
        ax.boxplot(groups,positions=range(3),widths=.4,showfliers=False,patch_artist=True,
                   boxprops={"facecolor":"#E8E8E8"},medianprops={"color":"black","linewidth":1.4})
        for x,c in enumerate(CONDITIONS):
            group=[r for r in rows if r["divergence"]==p and r["condition"]==c]
            ax.scatter([x+positions[r["seed"]] for r in group],[scale*r[metric] for r in group],s=15,alpha=.6,color=colors[c],zorder=3)
            ax.plot(x,statistics.fmean(scale*r[metric] for r in group),"D",mfc="white",mec="black",ms=5,zorder=5)
        ax.set_xticks(range(3),CONDITIONS); ax.axhline(0,color="black",lw=.7,alpha=.5)
        ax.grid(axis="y",alpha=.15)
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),sharey=True)
    for ax,p in zip(axes,LEVELS):
        points_box(ax,p,"local_accuracy_change",100)
        ax.set(title=f"Root-to-tip difference {p:.2f}",ylabel="Local SP recall change (percentage points)")
    fig.suptitle("Exploratory/post-hoc: 50 observations per group; boxes show quartiles, diamonds means")
    save(fig,"accuracy_change_distributions")
    fig,axes=plt.subplots(2,3,figsize=(13,8),sharex=True)
    for i,p in enumerate(LEVELS):
        for j,(metric,label) in enumerate([("recovered_local","Correct local pairs recovered"),("lost_local","Correct local pairs lost"),("net_correct_local","Net correct local pairs")]):
            points_box(axes[i,j],p,metric)
            axes[i,j].set(title=f"p={p:.2f}: {label}",ylabel="Pairs per dataset")
    fig.suptitle("Exploratory/post-hoc: recovered, lost and net pairs; every dataset and outlier shown")
    save(fig,"pair_count_distributions")
    fig,axes=plt.subplots(2,2,figsize=(11,8),sharex=True,sharey=True)
    for i,p in enumerate(LEVELS):
        for j,(metric,label) in enumerate([("recovered_local","Gross recovered pairs"),("lost_local","Gross lost pairs")]):
            ax=axes[i,j]
            for c in CONDITIONS:
                selected=[r for r in ranks if r["divergence"]==p and r["condition"]==c and r["metric"]==metric]
                if selected[-1]["cumulative_share"] is not None:
                    ax.plot([0]+[r["rank"] for r in selected],[0]+[100*r["cumulative_share"] for r in selected],label=c,color=colors[c],lw=1.8)
            ax.plot([0,50],[0,100],"--",color="#999999",lw=1,label="Equal contributions" if i==j==0 else None)
            ax.axvline(5,color="#777777",ls=":",lw=.8); ax.axvline(10,color="#777777",ls=":",lw=.8)
            ax.set(title=f"p={p:.2f}: {label}",xlabel="Top-ranked datasets (of 50)",ylabel="Share of pairs (%)",xlim=(0,50),ylim=(0,102))
    axes[0,0].legend(fontsize=8,frameon=False)
    fig.suptitle("Exploratory/post-hoc: concentration of gross recoveries and losses (ranked separately)")
    save(fig,"concentration_curves")
    fig,axes=plt.subplots(1,2,figsize=(11,5),sharey=True)
    for ax,p in zip(axes,LEVELS):
        for seed in SEEDS:
            selected=[next(r for r in rows if r["divergence"]==p and r["seed"]==seed and r["condition"]==c) for c in CONDITIONS]
            ax.plot([x+positions[seed] for x in range(3)],[r["local_accuracy_change"] for r in selected],color="#777777",alpha=.25,lw=.7,marker=".",ms=3)
        for x,c in enumerate(CONDITIONS):
            vals=[r["local_accuracy_change"] for r in rows if r["divergence"]==p and r["condition"]==c]
            ax.plot(x,statistics.fmean(vals),"o",color=colors[c],ms=7,zorder=5)
            ax.plot(x,statistics.median(vals),"D",mfc="white",mec="black",ms=5,zorder=6)
        ax.axhline(0,color="black",lw=.8)
        ax.set_xticks(range(3),["none\nunion=24","partial\nunion=18","high\nunion=13"])
        ax.set(title=f"Root-to-tip difference {p:.2f}",ylabel="Normalized local net refinement gain")
    handles=[Line2D([0],[0],color="#777777",lw=1,label="One matched seed block"),Line2D([0],[0],marker="o",color="#0072B2",ls="",label="Mean"),
             Line2D([0],[0],marker="D",mfc="white",mec="black",color="black",ls="",label="Median")]
    axes[1].legend(handles=handles,fontsize=8,frameon=False)
    fig.suptitle("Exploratory/post-hoc: all 50 seed blocks remain paired across deletion geometries")
    save(fig,"paired_seed_net_gain")
    for metric,name,ylabel in [("leave_one_out_estimate","leave_one_seed_out_estimates","Mean paired contrast (n=49)"),
                               ("change_percentage_points","leave_one_seed_out_influence","LOO shift (percentage points)")]:
        fig,axes=plt.subplots(2,2,figsize=(12,8))
        for i,p in enumerate(LEVELS):
            for j,c in enumerate(["partial","high"]):
                ax=axes[i,j]; label=c+"_minus_none"
                selected=[r for r in loo if r["divergence"]==p and r["contrast"]==label]
                ref=next(r for r in loo_summary if r["divergence"]==p and r["contrast"]==label)
                if metric=="leave_one_out_estimate":
                    ax.axhspan(float(ref["original_ci_lower_reference_only"]),float(ref["original_ci_upper_reference_only"]),color="#DADADA",alpha=.45,label="Original n=50 CI (reference only)")
                    ax.axhline(float(ref["original_full_estimate"]),color="#0072B2",ls="--",lw=1,label="Original n=50 mean")
                    ax.plot([r["omitted_seed"] for r in selected],[r[metric] for r in selected],".",color="#333333",ms=6)
                else:
                    for r in selected:
                        color="#0072B2" if r[metric]>=0 else "#D55E00"
                        ax.vlines(r["omitted_seed"],0,r[metric],color=color,alpha=.7,lw=1)
                        ax.plot(r["omitted_seed"],r[metric],".",color=color,ms=5)
                    influential=max(selected,key=lambda r:r["absolute_change"])
                    ax.annotate(str(influential["omitted_seed"]),(influential["omitted_seed"],influential[metric]),xytext=(4,4),textcoords="offset points",fontsize=8)
                ax.axhline(0,color="black",lw=.7)
                ax.set(title=f"p={p:.2f}: {c} minus none",xlabel="Omitted evolutionary seed",ylabel=ylabel)
                ax.set_xticks([2001,2010,2020,2030,2040,2050]); ax.tick_params(axis="x",labelrotation=30)
        if metric=="leave_one_out_estimate":
            axes[0,0].legend(fontsize=7,frameon=False)
        fig.suptitle("Exploratory/post-hoc: leave-one-seed-out sensitivity; 49 retained blocks; no new CIs")
        save(fig,name)
    return matplotlib.__version__


def write_report(rows,frequency,distribution,shares,loo_summary,validation,source):
    text=["# Phase 1 supplementary analysis — exploratory/post-hoc","",
        "This supplement describes existing cap-2 versus cap-0 results only. It does not replace the original prespecified analysis, "
        "confidence intervals, reports or manuscript figures. No MAFFT calls, simulations, reference regeneration or original-pipeline commands were executed.","",
        f"Read-only source: `{source}`. There are six groups of 50 independent seed blocks each (2001–2050): three deletion geometries at each expected root-to-tip difference .25/.40. "
        "The same seeds are matched across conditions and divergences; the 300 datasets are not 300 independent blocks. These are simulated protein sets, not observed natural protein families.","",
        "## 1. How often does refinement help?","",
        "Local accuracy is the original focal-boundary reference-pair SP recall. Its exact change equals `(recovered_local-lost_local)/true_local_pairs`. "
        "An integer positive/negative/zero numerator classifies improvement/worsening/no change. Unchanged local accuracy does not mean that every aligned pair stayed unchanged. "
        "All percentages use n=50 within the indicated group. Changes below are percentage points, not relative percentage improvements.","",
        "| p | Geometry (union gap columns) | Improved | Worsened | Unchanged | Mean change (pp) | Median change (pp) |",
        "|---|---|---:|---:|---:|---:|---:|"]
    for r in frequency:
        fmt=lambda label:f"{r[label+'_datasets']}/50 ({r[label+'_percent']:.0f}%)"
        text.append(f"| {r['divergence']:.2f} | {r['condition']} ({r['union_gap_columns']}) | {fmt('improved')} | {fmt('worsened')} | {fmt('unchanged')} | {r['mean_change_percentage_points']:.4f} | {r['median_change_percentage_points']:.4f} |")
    text.extend(["","Zero changes, positive means and zero medians can coexist: counts and medians provide a different descriptive view from means. "
        "See `tables/improvement_frequencies.csv` and `tables/per_dataset_results.csv` for full precision and all individual observations.","",
        "## 2. Distribution and concentration of recovered and lost pairs","",
        "Counts include all 50 datasets, including zeros. Q1/Q3 are linearly interpolated empirical quartiles; IQR=Q3-Q1. "
        "Net counts are recovered minus lost, before normalization. Pair counts are descriptive measurements, not independent experimental samples.","",
        "| p | Geometry | Metric | Total | Median | Q1 | Q3 | IQR | Maximum |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|"])
    for r in distribution:
        if r["metric"]!="local_accuracy_change":
            text.append(f"| {r['divergence']:.2f} | {r['condition']} | {r['metric']} | {r['total']} | {r['median']:g} | {r['q1']:g} | {r['q3']:g} | {r['iqr']:g} | {r['max']} |")
    text.extend(["","Concentration ranks datasets separately for each nonnegative measurement. Top 5 means 10% of the 50 datasets. "
        "Gross recoveries can coexist with losses in the same dataset; they are not automatically net benefits. "
        "Signed net totals are not used as concentration denominators because positive/negative cancellation can produce misleading shares. "
        "Additional positive-net and negative-net concentration results are in the CSV. Zero-total shares are undefined.","",
        "| p | Geometry | Metric | Nonzero contributors | Top 1 share | Top 5 share | Top 10 share | Datasets for 80% |",
        "|---|---|---|---:|---:|---:|---:|---:|"])
    def pct(v):
        return "undefined" if v is None else f"{100*v:.1f}%"
    for r in shares:
        if r["metric"] in ["recovered_local","lost_local"]:
            text.append(f"| {r['divergence']:.2f} | {r['condition']} | {r['metric']} | {r['nonzero_datasets']}/50 | {pct(r['top_1_share'])} | {pct(r['top_5_share'])} | {pct(r['top_10_share'])} | {r['datasets_for_80_percent']} |")
    gross=[r for r in shares if r["metric"] in ["recovered_local","lost_local"] and r["total"]]
    text.extend(["",f"Observed top-five shares of gross recovered/lost pairs range from {100*min(r['top_5_share'] for r in gross):.1f}% "
        f"to {100*max(r['top_5_share'] for r in gross):.1f}%. This directly quantifies concentration in a minority of simulated datasets; "
        "it does not show that particular natural protein-family types cause refinement to help or harm. "
        "`tables/concentration_rankings.csv` identifies every contributor and cumulative share without hiding zero-contribution datasets.","",
        "## 3. Matched per-seed comparisons","",
        "`figures/paired_seed_net_gain.png` joins none, partial and high geometry for every seed in separate divergence panels. "
        "Small deterministic x offsets reveal overlapping observations and do not alter outcomes or use/change any random seed. "
        "Colored circles are means; hollow diamonds are medians. All matching is retained in `tables/per_seed_paired_contrasts.csv`.","",
        "## 4. Leave-one-seed-out (LOO) sensitivity","",
        "One evolutionary seed block is removed at a time from all matched conditions/divergences. Each estimate uses the remaining 49 blocks. "
        "The four original cap-2 contrasts are recomputed as means of paired, already-normalized condition differences. "
        "Original means and original n=50 confidence intervals are read-only references. No LOO confidence intervals or p-values are computed. "
        "A point estimate changing sign is sensitivity of its direction, not a new significance result.","",
        "| p | Contrast | Original mean | LOO minimum | LOO maximum | Largest absolute shift (pp) | Omitted seed for largest shift | Direction changes /50 |",
        "|---|---|---:|---:|---:|---:|---:|---:|"])
    for r in loo_summary:
        text.append(f"| {r['divergence']:.2f} | {r['contrast']} | {float(r['original_full_estimate']):.8f} | {r['min_leave_one_out_estimate']:.8f} | {r['max_leave_one_out_estimate']:.8f} | {r['max_absolute_change_percentage_points']:.4f} | {r['most_influential_seed']} | {r['direction_changes_out_of_50']}/50 |")
    for r in loo_summary:
        direction="changes direction after at least one omission" if r["direction_changes_out_of_50"] else "retains its direction for every single-seed omission"
        text.extend(["",f"At p={r['divergence']:.2f}, {r['contrast']} {direction}; the largest one-block shift is "
            f"{r['max_absolute_change_percentage_points']:.4f} percentage points (seed {r['most_influential_seed']}). "
            "Use the size and sign of these shifts to assess point-estimate influence, alongside the original uncertainty."])
    text.extend(["","All four original primary 95% bootstrap confidence intervals included zero. These original intervals remain unchanged. "
        "LOO point estimates alone cannot establish that a previously uncertain contrast becomes supported or unsupported after deletion, "
        "nor prove that results would be robust to different simulation models or protein families. The 200 overlapping LOO estimates are dependent diagnostics, not 200 independent samples.","",
        "## Supported results, exploratory observations and limits","",
        "Supported descriptive results are the exact improvement/worsening counts, empirical distributions and contributor shares, "
        "and the arithmetic influence of omitting one observed seed block. Their calculations were checked against the original per-mode scores, "
        "comparison JSONs, completion markers, original paired contrasts and original cell means. All individual observations remain visible.","",
        "The concentration and LOO interpretations are exploratory/post-hoc. They assess how the existing sample is distributed; "
        "they do not supply new confirmatory evidence, explain the biological cause of an outlier, or justify excluding an influential dataset. "
        "No seeds are dropped from the original estimates. No original experimental configuration, scoring rule or statistical result is changed.","",
        "Deletion overlap remains confounded with spacing, affected-site identities and union gap columns (24/18/13). "
        "An independent causal effect of overlap cannot be identified. One fixed tree, a uniform random amino-acid root, equal-rate substitutions "
        "and only two independent terminal deletions restrict generalization. Data from different conditions and divergences share seed blocks; "
        "residue pairs and all 300 condition-specific datasets must not be treated as independent replicates.","",
        "## Validation and strict preservation","",
        f"Calculation validation passed: {validation['dataset_comparisons_validated']} dataset comparisons, "
        f"{validation['original_score_jsons_checked']} saved scores, {validation['original_comparison_jsons_checked']} saved comparisons, "
        f"{validation['completion_markers_checked']} completion markers, 200 original seed-level contrasts and 200 LOO estimates. "
        "Independent NumPy checks agree on means, medians and quartiles; hand-worked tests cover concentration, zero totals and LOO pairing.","",
        "The supplementary script enforces a write boundary around this directory and blocks external process launches. "
        "Python bytecode writes are disabled; plotting caches and temporary files are confined here. "
        "Executed MAFFT calls: **0**. New simulations/reference alignments: **0**. Original simulation/MAFFT pipeline modules imported: **0**.","",
        "The full before/after SHA256 and timestamp checks are in `preservation_before.json` and `preservation_verification.json`. "
        "They cover all original workspace files, not only the reports and alignments. See `validation_results.json`, `logs/`, "
        "`source_manifest.json` and `README.md` for validation evidence, source paths and reproduction instructions. "
        "The final packaging step appends the actual preservation outcome after verification; it does not modify an original report.","",
        "## Figures","",
        "All new figures are explicitly exploratory/post-hoc and exported as 300-DPI PNG plus SVG. "
        "Boxes use Q1/median/Q3 and 1.5-IQR whiskers, with all observations including outliers displayed separately. "
        "No original manuscript figure is modified.",""])
    for name in ["accuracy_outcome_frequencies","accuracy_change_distributions","pair_count_distributions","concentration_curves","paired_seed_net_gain","leave_one_seed_out_estimates","leave_one_seed_out_influence"]:
        text.append(f"- [{name}](figures/{name}.png)")
    atomic("supplementary_analysis.md","\n".join(text)+"\n")


def main():
    parser=argparse.ArgumentParser(description="Exploratory supplement using existing results only")
    parser.add_argument("--source",type=Path,default=DEFAULT_SOURCE)
    args=parser.parse_args()
    install_guard()
    baseline=require_baseline()
    source=args.source.resolve()
    started=time.time()
    rows,_=load_data(source,baseline)
    frequency,distribution,shares,ranks=descriptive_tables(rows)
    pairs,loo,loo_summary=leave_one_out(rows,source)
    validation=unit_and_data_validation(rows,frequency,distribution,shares,ranks,loo,source)
    outputs={"per_dataset_results":rows,"improvement_frequencies":frequency,"distribution_summaries":distribution,
        "concentration_summaries":shares,"concentration_rankings":ranks,"per_seed_paired_contrasts":pairs,
        "leave_one_seed_out_results":loo,"leave_one_seed_out_summary":loo_summary}
    for name,data in outputs.items():
        atomic("tables/"+name+".csv",csv_text(data))
    validation["matplotlib_version"]=make_figures(rows,frequency,ranks,loo,loo_summary)
    assert not PROCESS_ATTEMPTS and not WRITE_VIOLATIONS
    validation["external_process_attempts"]=PROCESS_ATTEMPTS
    validation["original_write_attempts_blocked"]=WRITE_VIOLATIONS
    validation["elapsed_seconds"]=time.time()-started
    save_json("validation_results.json",validation)
    write_report(rows,frequency,distribution,shares,loo_summary,validation,source)
    log={"completed":True,"exploratory_post_hoc":True,"utc_epoch":time.time(),"source":str(source),
         "datasets":len(rows),"independent_seed_blocks_per_group":50,"mafft_calls_executed":0,
         "external_process_attempts":PROCESS_ATTEMPTS,"original_write_attempts":WRITE_VIOLATIONS,
         "original_pipeline_imported":False,"new_script_sha256":sha(__file__),"validation_passed":True}
    save_json("logs/analysis_execution.json",log)
    print(json.dumps(log,indent=2))


if __name__=="__main__":
    main()
