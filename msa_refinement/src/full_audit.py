"""Independent Phase 1 content audit and Phase 0 preservation verification."""
from pathlib import Path
from . import scoring, simulation
from .full_pipeline import PROJECT, check_design, batch_label, source_hashes
from .full_reporting import inventory
from .reporting import csv_text
from .storage import canonical, digest, json_write, read_json, sha


def preserve_phase0(root):
    baseline=PROJECT/"runs/phase1_provenance/phase0_before.json"
    if not baseline.exists():
        raise ValueError("Missing Phase 0 preservation baseline")
    evidence=Path(root)/"phase0_preservation.json"
    if evidence.exists():
        previous=read_json(evidence)
        if (previous.get("passed") and previous.get("workspace_project_path") and
            previous["workspace_project_path"]!=str(PROJECT.resolve()) and previous.get("baseline_sha256")==sha(baseline)):
            # The ZIP intentionally excludes Phase 0 raw files. Keep the original proof,
            # without pretending to check unavailable files/timestamps after relocation.
            return {**previous,"check_scope":"verified_original_workspace_evidence_only","current_workspace_checked":False}
    records=read_json(baseline)["files"]
    changed=[]
    for name,expected in records.items():
        path=PROJECT/name
        if not path.exists() or sha(path)!=expected["sha256"] or path.stat().st_mtime_ns!=expected["mtime_ns"]:
            changed.append(name)
    result={"passed":not changed,"phase0_files_checked":len(records),"bytes_and_mtimes_unchanged":not changed,"changed_files":changed,
            "baseline_sha256":sha(baseline),"workspace_project_path":str(PROJECT.resolve()),
            "check_scope":"all_original_workspace_files","current_workspace_checked":True}
    json_write(Path(root)/"phase0_preservation.json",result)
    if changed:
        raise ValueError("Phase 0 changed: "+str(changed))
    return result


def audit(root):
    root=Path(root)
    config=read_json(root/"config.json")
    check_design(config)
    manifest=read_json(root/"experiment_manifest.json")
    assert sha(root/"config.json")==manifest["frozen_config_file_sha256"]
    assert "FAKE" not in manifest["mafft_version"].upper()
    jobs,rows=inventory(root)
    assert len(jobs)==900 and len(rows)==900
    assert all(j["alignment_state"]=="completed" and j["score_complete"] and j["comparison_complete"] for j in jobs)
    assert len({(r["divergence"],r["seed"],r["condition"],r["iteration_cap"]) for r in rows})==900
    assert (root/"full_summary.csv").read_text(encoding="utf-8")==csv_text(rows)
    scores_checked=0; comparisons_checked=0; cap2_cap10_identical=0; case_count=0
    for p in config["divergence_levels"]:
        batch=batch_label(p)
        for seed in config["simulator"]["replicate_seeds"]:
            shared=simulation.realization(config,seed,p)
            block=read_json(root/"blocks"/f"{batch}_seed{seed}"/"realization.json")
            assert shared==block
            for cond in config["conditions"]:
                case=root/"cases"/f"{batch}_seed{seed}_{cond}"
                inp,truth,meta=simulation.simulate(config,seed,cond,p,shared)
                assert scoring.parse_fasta((case/"input.fasta").read_text())==inp
                assert scoring.parse_fasta((case/"truth_alignment.fasta").read_text())==truth
                assert read_json(case/"ancestry.json")==meta
                t=scoring.pairs(truth); boundary=scoring.boundary_residues(meta); modes={}
                for cap in [0,2,10]:
                    directory=case/"alignments"/f"maxiter_{cap:03d}"
                    alignment=scoring.parse_fasta((directory/"aligned.fasta").read_text())
                    scoring.validate_alignment(alignment,inp)
                    modes[cap]=scoring.pairs(alignment)
                    saved=read_json(case/"scores"/f"maxiter_{cap:03d}.json")
                    recomputed=scoring.accuracy(modes[cap],t,boundary)
                    assert all(saved[k]==v for k,v in recomputed.items())
                    done=read_json(directory/"run_done.json")
                    assert done["utc_epoch"]>=manifest["frozen_utc_epoch"]
                    assert done["details"]["command"][-8:-1]==["--amino","--thread","1","--retree","2","--maxiterate",str(cap)]
                    for key in ["reference_pairs","predicted_pairs","correct_pairs","false_pairs","missed_pairs"]:
                        assert saved["global"][key]==saved["local"][key]+saved["far"][key]
                    scores_checked+=1
                for cap in [2,10]:
                    saved=read_json(case/"comparisons"/f"against_maxiter_{cap:03d}.json")
                    actual=scoring.compare(modes[0],modes[cap],t,boundary)
                    assert all(saved[k]==v for k,v in actual.items())
                    comparisons_checked+=1
                cap2_cap10_identical+=modes[2]==modes[10]
                case_count+=1
            if seed%10==0:
                print(f"Audit: p={p}, seed {seed}; {scores_checked} scores independently recomputed",flush=True)
    from .full_analysis import analyze
    paired,intervals=analyze(rows,config)
    assert (root/"paired_comparisons.csv").read_text()==csv_text(paired)
    assert (root/"confidence_intervals.csv").read_text()==csv_text(intervals)
    validation=read_json(PROJECT/"runs/phase1_validation/latest.json")
    assert validation["exit_code"]==0 and validation["source_sha256"]==source_hashes()
    phase0=preserve_phase0(root)
    result={"passed":True,"datasets_regenerated":case_count,"alignment_jobs_verified":len(jobs),"scores_recomputed":scores_checked,
            "comparisons_recomputed":comparisons_checked,"cap2_cap10_identical_pair_sets":cap2_cap10_identical,
            "statistical_estimates_reproduced":len(intervals),"paired_comparisons_reproduced":len(paired),"seed_blocks_per_level":50,
            "pilot_blocks_in_estimates":0,"frozen_before_first_alignment":True,"phase0_preservation_passed":phase0["passed"],
            "phase0_preservation_check_scope":phase0["check_scope"],
            "source_sha256":source_hashes()}
    json_write(root/"full_artifact_audit.json",result)
    return result
