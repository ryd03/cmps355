"""Rebuild reports from verified artifacts, never append summary rows."""
import csv
from collections import Counter
import io
import json
from pathlib import Path
import platform
import statistics
import zipfile
from .storage import atomic, canonical, complete, digest, json_write, read_json, sha, writer_lock
from . import scoring, simulation, runner


def verified(path):
    try:
        marker = read_json(path)
        return complete(path, marker["dependencies"])
    except (OSError, ValueError, KeyError):
        return False


def inventory(root):
    config = read_json(root / "config.json")
    manifest = read_json(root / "experiment_manifest.json")
    jobs, rows = [], []
    skips = Counter()
    for log in (root / "execution_events").glob("*.json"):
        entry = read_json(log)
        if entry["state"] == "skipped":
            skips[entry["stage"].replace("\\", "/")] += 1
    batches = ["initial"]
    if (root / "harder_decision.json").exists() or any((root / "cases").glob("harder_*")):
        batches.append("harder")
    for batch in batches:
        for seed in config["simulator"]["replicate_seeds"]:
            for cond in config["simulator"]["second_deletion_start_delta_by_condition"]:
                case = root / "cases" / f"{batch}_seed{seed}_{cond}"
                sim_ok = verified(case / "simulation_done.json")
                meta = read_json(case / "ancestry.json") if sim_ok else None
                valid = {}
                for cap in config["mafft"]["maxiterate"]:
                    directory = case / "alignments" / f"maxiter_{cap:03d}"
                    done = directory / "run_done.json"
                    ok = sim_ok and verified(done)
                    if ok:
                        try:
                            d = read_json(done)["dependencies"]
                            assert d["input"] == sha(case / "input.fasta")
                            assert d["mafft_version"] == manifest["identity"]["mafft_version"]
                            scoring.validate_alignment(scoring.parse_fasta((directory / "aligned.fasta").read_text()),
                                                       scoring.parse_fasta((case / "input.fasta").read_text()))
                        except (OSError, ValueError, AssertionError, KeyError):
                            ok = False
                    valid[cap] = ok
                    state = "completed" if ok else ("failed" if (directory / "run_done.json.failed.json").exists()
                                                     else "interrupted_or_invalid" if directory.exists() else "pending")
                    scorefile = case / "scores" / f"maxiter_{cap:03d}.json"
                    score_ok = ok and verified(scorefile.with_name(scorefile.name + ".done.json"))
                    if score_ok:
                        dep = read_json(scorefile)["dependencies"]
                        score_ok = (dep["aligned"] == sha(directory / "aligned.fasta") and
                                    dep["truth"] == sha(case / "truth_alignment.fasta") and
                                    dep["ancestry"] == sha(case / "ancestry.json") and
                                    dep["scoring_code"] == sha(scoring.__file__))
                    jobs.append({"case": case.name, "batch": batch, "seed": seed, "condition": cond,
                                 "iteration_cap": cap, "alignment_state": state, "score_complete": score_ok,
                                 "verified_resume_skip_count": skips[str(done.relative_to(root)).replace("\\", "/")],
                                 "historical_failures": len(list((directory / "failures").glob("*.json"))),
                                 "path": str(directory.relative_to(root))})
                    if not score_ok:
                        continue
                    score = read_json(scorefile)
                    run = read_json(done)["details"]
                    baseline = case / "scores/maxiter_000.json"
                    row = {"batch": batch, "seed": seed, "condition": cond, "iteration_cap": cap,
                           "expected_root_tip_difference": meta["expected_root_to_tip_difference"],
                           "realized_overlap": meta["realized_overlap"], "union_gap_columns": meta["union_gap_columns"],
                           "deletion1_start": meta["events"][0]["start"], "deletion2_start": meta["events"][1]["start"],
                           "mean_pairwise_identity": meta["mean_pairwise_identity_shared_extant_ancestral_sites"],
                           "mean_root_tip_endpoint_differences": statistics.mean(meta["root_tip_endpoint_differences"].values()),
                           "sum_branch_endpoint_changes": sum(e["endpoint_changes"] for e in meta["edges"]),
                           "surviving_residues": meta["total_surviving_residues"],
                           "aligned_length": score["aligned_length"], "runtime_seconds": run["runtime_seconds"],
                           "input_sha256": sha(case / "input.fasta"), "aligned_sha256": sha(directory / "aligned.fasta")}
                    if verified(baseline.with_name(baseline.name + ".done.json")):
                        row["baseline_local_missed"] = read_json(baseline)["local"]["missed_pairs"]
                    for region in ["local", "far", "global"]:
                        row.update({region + "_" + k: v for k, v in score[region].items()})
                    compfile = case / "comparisons" / f"against_maxiter_{cap:03d}.json"
                    if cap and valid.get(0) and verified(compfile.with_name(compfile.name + ".done.json")):
                        comp = read_json(compfile)
                        dep = comp["dependencies"]
                        if dep["aligned"] == row["aligned_sha256"] and dep["baseline"] == sha(case / "alignments/maxiter_000/aligned.fasta") and dep["scoring_code"] == sha(scoring.__file__):
                            for region in ["local", "far", "global"]:
                                row.update({region + "_" + k: v for k, v in comp[region].items()})
                            row["alignment_pair_set_changed"] = comp["alignment_pair_set_changed"]
                    rows.append(row)
    return jobs, rows


def csv_text(rows):
    if not rows:
        return ""
    fields = list(dict.fromkeys(k for row in rows for k in row))
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def plots(root, rows):
    if not rows:
        return
    import os
    os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    directory = root / "plots"
    directory.mkdir(exist_ok=True)
    colors = {1001: "#0072B2", 1002: "#D55E00", 1003: "#009E73", 1004: "#CC79A7"}
    conditions = ["none", "partial", "high"]
    batches = list(dict.fromkeys(r["batch"] for r in rows))
    for metric, title, ylabel in [
        ("local_net_gain", "Focal-boundary net correct-pair change", "(recovered - lost) / local reference pairs"),
        ("local_recovered", "Recovered focal-local reference pairs", "Recovered pairs"),
        ("local_lost", "Lost focal-local reference pairs", "Lost pairs"),
        ("global_sp_recall", "Global reference-pair SP recall", "SP recall"),
        ("baseline_local_missed", "Baseline focal-local missed pairs", "Missed reference pairs"),
        ("runtime_seconds", "MAFFT wall-clock runtime", "Seconds (Windows overhead included)")]:
        caps = [0] if metric == "baseline_local_missed" else [0, 2, 10] if metric in ["global_sp_recall", "runtime_seconds"] else [2, 10]
        fig, axs = plt.subplots(len(batches), len(caps), figsize=(5 * len(caps), 3.8 * len(batches)), squeeze=False)
        for i, batch in enumerate(batches):
            for j, cap in enumerate(caps):
                ax = axs[i][j]
                for seed, color in colors.items():
                    selected = [next((r for r in rows if r["batch"] == batch and r["iteration_cap"] == cap and
                                      r["seed"] == seed and r["condition"] == cond), {}) for cond in conditions]
                    ys = [r.get(metric) for r in selected]
                    ax.plot(range(3), [float("nan") if y is None else y for y in ys], "o-", color=color, label=str(seed), alpha=.8)
                ax.set_xticks(range(3), ["none\nunion=24", "partial\nunion=18", "high\nunion=13"])
                ax.set(title=f"{batch}; iteration cap {cap}", ylabel=ylabel)
                ax.grid(alpha=.2)
                if metric == "local_net_gain":
                    ax.axhline(0, color="black", linewidth=.7)
                ax.legend(title="Block seed", fontsize=8)
        fig.suptitle(title + " (all four matched blocks)")
        fig.tight_layout()
        for ext in ["png", "svg"]:
            buffer = io.BytesIO()
            fig.savefig(buffer, format=ext, dpi=160)
            atomic(directory / (metric + "." + ext), buffer.getvalue())
        plt.close(fig)


def build(root):
    root = Path(root)
    jobs, rows = inventory(root)
    atomic(root / "pilot_summary.csv", csv_text(rows))
    json_write(root / "job_manifest.json", {"jobs": jobs, "completed": sum(j["alignment_state"] == "completed" for j in jobs),
               "failed": sum(j["alignment_state"] == "failed" for j in jobs),
               "pending_or_invalid": sum(j["alignment_state"] not in ["completed", "failed"] for j in jobs),
               "skip_history": "See execution_events and invocations; skips reuse completed jobs, not missing observations."})
    atomic(root / "job_manifest.csv", csv_text(jobs))
    plots(root, rows)
    config = read_json(root / "config.json")
    manifest = read_json(root / "experiment_manifest.json")
    decisions = [read_json(f) for f in [root / "initial_decision.json", root / "harder_decision.json"] if f.exists()]
    text = ["# Phase 0 pilot report", "",
            "This is a descriptive feasibility pilot with four independent seed blocks (1001–1004). "
            "Conditions and iteration caps are paired within each block. Residue pairs are not independent replicates. "
            "No significance test or novelty claim is made. Phase 0 is the stopping point.", "",
            "## Execution and provenance", "",
            f"Experiment: `{root.name}`. Python {platform.python_version()}; {manifest['platform']}. "
            f"MAFFT version output: `{manifest['identity']['mafft_version']}`.", "",
            "Commands (from project root):", "", "```",
            "python -m src.cli validate",
            "python -m src.cli pilot --config config/pilot.json --workdir ./runs",
            "python -m src.cli status --workdir ./runs",
            "python -m src.cli pilot --config config/pilot.json --workdir ./runs",
            "python -m src.cli report --workdir ./runs", "```", "",
            "Every MAFFT call uses `--amino --thread 1 --retree 2 --maxiterate 0|2|10` on the same input within a case. "
            "Exact executable argv, input/output hashes, elapsed seconds, stderr, and version are in each `run_done.json` and `stderr.log`. "
            "The iteration cap does not establish the number of cycles actually performed. No internal WSP measurement is claimed.", "",
            f"Validated completed MAFFT jobs: **{sum(j['alignment_state'] == 'completed' for j in jobs)}**; "
            f"currently failed: **{sum(j['alignment_state'] == 'failed' for j in jobs)}**; "
            f"pending/interrupted/invalid: **{sum(j['alignment_state'] not in ['completed','failed'] for j in jobs)}**. "
            f"Historical alignment failures: {sum(j['historical_failures'] for j in jobs)}. "
            "`job_manifest.csv` lists all expected jobs. `invocations/` and `execution_events/` retain completed/failed/skipped stage counts per invocation.", "",
            "## Prespecified calibration decision", ""]
    for d in decisions:
        text.extend(["```json", json.dumps(d, indent=2), "```", ""])
    if decisions and decisions[0].get("complete"):
        text.append("The harder calibration was triggered and retained as a separate batch." if decisions[0]["low_signal"] else
                    "The harder calibration was not triggered. Only the predefined initial 36-job batch was executed.")
        text.append("")
    text.append("Harder calibration ran only if the initial batch met the numerical rule: fewer than 2/12 cases with any local correct-pair gain/loss under either cap, OR at least 10/12 baselines with no local missed correct pairs. A maximum of one harder batch is allowed.")
    text.extend(["", "### Invocation history", "", "| Invocation | batch | newly completed | failed | completed work skipped | failed work skipped |", "|---|---|---:|---:|---:|---:|"])
    for invocation_file in sorted((root / "invocations").glob("*.json")):
        for batch, details in read_json(invocation_file)["batches"].items():
            c = details["alignments"]
            text.append(f"| {invocation_file.stem} | {batch} | {c.get('completed',0)} | {c.get('failed',0)} | {c.get('skipped',0)} | {c.get('failed_skipped',0)} |")
    text.extend(["", "## Observations", "",
                 "Counts below sum descriptive events over cases and caps: the same recovered pair may be counted once for cap 2 and once for cap 10. These are not unique biological events or independent samples.", "",
                 "| Batch | cap | pair-set changes /12 | recovered local | lost local | false added local | false removed local | mean net local | mean global SP |",
                 "|---|---:|---:|---:|---:|---:|---:|---:|---:|"])
    for batch in list(dict.fromkeys(r["batch"] for r in rows)):
        for cap in [2, 10]:
            selected = [r for r in rows if r["batch"] == batch and r["iteration_cap"] == cap and "local_recovered" in r]
            if not selected:
                continue
            sums = [sum(r[k] for r in selected) for k in ["alignment_pair_set_changed", "local_recovered", "local_lost", "local_false_introduced", "local_false_removed"]]
            net = statistics.mean(r["local_net_gain"] for r in selected if r["local_net_gain"] is not None)
            sp = statistics.mean(r["global_sp_recall"] for r in selected)
            text.append(f"| {batch} | {cap} | {sums[0]}/{len(selected)} | {sums[1]} | {sums[2]} | {sums[3]} | {sums[4]} | {net:.6f} | {sp:.6f} |")
    text.extend(["", "### Geometry summaries (each of four blocks retained in figures)", "",
                 "| Batch | geometry (union) | cap | recovered [min,max] | lost [min,max] | mean net local |",
                 "|---|---|---:|---:|---:|---:|"])
    for batch in list(dict.fromkeys(r["batch"] for r in rows)):
        for cond in ["none", "partial", "high"]:
            for cap in [2, 10]:
                s = [r for r in rows if r["batch"] == batch and r["condition"] == cond and r["iteration_cap"] == cap and "local_recovered" in r]
                if s:
                    text.append(f"| {batch} | {cond} ({s[0]['union_gap_columns']}) | {cap} | "
                                f"{sum(r['local_recovered'] for r in s)} [{min(r['local_recovered'] for r in s)},{max(r['local_recovered'] for r in s)}] | "
                                f"{sum(r['local_lost'] for r in s)} [{min(r['local_lost'] for r in s)},{max(r['local_lost'] for r in s)}] | "
                                f"{statistics.mean(r['local_net_gain'] for r in s):.6f} |")
    text.extend(["", "### Every case and mode", "",
                 "Blank comparison cells indicate baseline or unavailable comparisons. Full precision, near/far/global decompositions, raw counts, denominators, null rates, hashes and times are in `pilot_summary.csv`.", "",
                 "| Batch | seed | geometry | cap | identity | baseline local missed | recovered | lost | true local denominator | net local | global SP | seconds |",
                 "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"])
    for r in rows:
        def fmt(key):
            v = r.get(key)
            return "" if v is None else f"{v:.6f}" if isinstance(v, float) else str(v)
        text.append("| " + " | ".join(fmt(k) for k in ["batch", "seed", "condition", "iteration_cap", "mean_pairwise_identity", "baseline_local_missed", "local_recovered", "local_lost", "local_reference_pairs", "local_net_gain", "global_sp_recall", "runtime_seconds"]) + " |")
    validation = Path(__file__).resolve().parents[1] / "validation_results/latest.json"
    text.extend(["", "## Validation and interruption recovery", ""])
    if validation.exists():
        evidence = read_json(validation)
        text.append(f"Latest automated suite exit code: {evidence['exit_code']}. {evidence.get('summary','')} "
                    "See `validation_results/` for complete pytest output, JUnit test names and failures, source hashes, and timing.")
    else:
        text.append("Validation evidence unavailable; results must not be accepted without the validation suite.")
    text.extend(["", "Recovery uses OS writer locks released on process death, same-directory unique partial files with flush/fsync and atomic replace, and final hash-verified markers. "
                 "The automated fake-MAFFT tests kill a real child pipeline process after two completed alignments, count subsequent invocations, inject a post-alignment scoring failure, "
                 "test corrupt/partial outputs, rerun idempotence, changed configurations, failure/timeout retry, and exclusive locking. "
                 "A mid-alignment interruption reruns that unfinished call. Completed raw outputs remain reusable if only scoring changes; scoring changes have their own source hash and quarantined prior artifacts.", "",
                 "## Interpretation and limits", ""])
    for batch in list(dict.fromkeys(r["batch"] for r in rows)):
        for cap in [2, 10]:
            selected = [r for r in rows if r["batch"] == batch and r["iteration_cap"] == cap and "local_net_gain" in r]
            if selected:
                seeds = sorted({r["seed"] for r in selected if r["local_recovered"] + r["local_lost"] > 0})
                gains = sum(r["local_net_gain"] > 0 for r in selected)
                losses = sum(r["local_net_gain"] < 0 for r in selected)
                means = {cond: statistics.mean(r["local_net_gain"] for r in selected if r["condition"] == cond)
                         for cond in ["none", "partial", "high"] if any(r["condition"] == cond for r in selected)}
                text.append(f"{batch}, cap {cap}: {gains} cases had positive net local change, {losses} had negative net change, "
                            f"and {len(selected)-gains-losses} had zero net change. Correct-pair changes occurred in seed blocks {seeds}. "
                            "Mean normalized net local changes by geometry were " + ", ".join(f"{c}={v:.6f}" for c,v in means.items()) +
                            ". This is descriptive variation across combined geometries; the small pilot cannot establish a reliable independent overlap effect.")
                text.append("")
    if decisions and decisions[-1].get("complete"):
        if decisions[-1]["low_signal"]:
            text.append("The final permitted batch remains low-signal under the locked rule. There is insufficient demonstrated dynamic range to launch the full study without design review. Stop here; do not search for favorable settings.")
        else:
            text.append("The executed pilot passes the prespecified dynamic-range screen. It supports the feasibility of studying local gains and losses, but four blocks cannot establish an overlap effect or justify a confirmatory sample size. Full-study design, power/precision planning, and explicit authorization are still required.")
    text.extend(["", "Geometry-specific observations above describe whether recovery/loss varied with the combined deletion arrangement. "
                 "Overlap is intrinsically confounded with union gap-column count (24/18/13), deletion spacing, removed site identities, and local sequence context. "
                 "No independent causal effect of overlap can be identified. Realized mean identity uses pairwise shared extant ancestral sites, not MAFFT output; substitutions recorded are endpoint differences, not unobserved numbers of CTMC jumps.", "",
                 "The root is a uniform 20-amino-acid random protein. The symmetric equal-rate CTMC is a pilot simplification, not a realistic protein model. "
                 "Every edge is evolved, using total root-to-tip scale `-19/20 log(1-20p/19)` for p=.25 (and .40 only if triggered). "
                 "These are expected root-to-tip endpoint differences, not specified pairwise identities. One fixed 12-tip tree, two independent terminal deletions, and no insertions limit generalization. "
                 "MAFFT infers its own guide tree, which may change between geometries. Runtime includes Windows launcher overhead and should not be generalized to Linux.", "",
                 "Implementation assumptions: deletion boundary distances use inclusive deleted start/end coordinates, so extant residues start one coordinate away; "
                 "the window is inclusive at distance 10. `Random` streams use the full big-endian SHA256 of canonical JSON `[seed, stream_name]`; names are `root`, `edge:<zero-based child-index path>`, `anchor`. "
                 "A complete substitution-only realization is generated once per seed/batch, checkpointed in `blocks/`, and reused for all three conditions. Tip root-to-tip distances are unit length. No scientific design or metric was changed.", "",
                 "Official MAFFT documentation describes FFT-NS-i as refining an FFT-NS-2 starting alignment using WSP: "
                 "https://mafft.ddbj.nig.ac.jp/alignment/software/algorithms/ . This is a possible mechanism, not a measured explanation for these results. "
                 "Related work already studies indel/alignment effects; Ezawa (2016) used mammalian DNA and E-INS methods, distinct from this protein FFT-NS pilot. No literature priority is established here.", "",
                 "All experimental datasets, truth and ancestry, MAFFT outputs, logs, seeds and frozen configuration are retained under `cases/`. "
                 "Any installation or validation corrections are recorded in `IMPLEMENTATION_NOTES.md` and validation logs. No full experiment has been started.", ""])
    atomic(root / "pilot_report.md", "\n".join(text))
    return jobs, rows


def bundle(project, workdir):
    project, workdir = Path(project), Path(workdir)
    target = project / "pilot_review_bundle.zip"
    partial = target.with_suffix(".zip.partial")
    files = []
    for directory in ["src", "tests", "config", "data", "validation_results", "tools"]:
        candidates = (project / directory).glob("*") if directory == "validation_results" else (project / directory).rglob("*")
        for path in candidates:
            if path.is_file() and "__pycache__" not in path.parts and (directory != "tools" or path.name in ["installation.json", "setup_mafft.py"]):
                files.append(path)
    files += [f for f in project.glob("*.md")]
    files += [f for f in project.glob("*.txt")]
    for path in workdir.rglob("*"):
        if path.is_file() and path.name != ".writer.lock" and not path.name.endswith(".partial"):
            files.append(path)
    with zipfile.ZipFile(partial, "w", zipfile.ZIP_DEFLATED) as archive:
        hashes = {}
        for path in sorted(set(files)):
            name = path.relative_to(project).as_posix()
            archive.write(path, name)
            hashes[name] = sha(path)
        manifests = list(workdir.glob("*/experiment_manifest.json"))
        if len(manifests) == 1:
            root = manifests[0].parent
            for path in [root / "pilot_report.md", root / "pilot_summary.csv", root / "job_manifest.csv"]:
                if path.exists():
                    archive.write(path,path.name)
                    hashes[path.name] = sha(path)
        archive.writestr("bundle_sha256_manifest.json", json.dumps(hashes, indent=2, sort_keys=True))
    with zipfile.ZipFile(partial) as archive:
        bad = archive.testzip()
        if bad:
            raise ValueError("Corrupt ZIP member " + bad)
        for name, expected in hashes.items():
            assert digest(archive.read(name)) == expected
    import os
    os.replace(partial, target)
    json_write(project / "bundle_verification.json", {"zip_sha256": sha(target), "members_verified": len(hashes),
               "zip_crc_and_all_member_sha256_passed": True})
    return target
