"""Phase 1 orchestration using the unchanged validated Phase 0 stage functions."""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import os
import platform
import sys
import time
from . import pipeline, simulation, scoring, runner
from .storage import atomic, canonical, digest, json_write, read_json, sha, writer_lock

PROJECT = Path(__file__).resolve().parents[1]
CONDITIONS = ["none", "partial", "high"]
CORE_FILES = ["simulation.py", "scoring.py", "runner.py", "pipeline.py", "storage.py", "reporting.py"]


def source_hashes():
    return {p.relative_to(PROJECT).as_posix(): sha(p) for folder in ["src", "tests"]
            for p in sorted((PROJECT / folder).glob("*.py"))}


def check_design(config, test_design=False):
    if config["phase"] != "phase1":
        raise ValueError("Phase 1 configuration required")
    original = read_json(PROJECT / "config/pilot.json")
    for key, value in original["simulator"].items():
        if key not in ["replicate_seeds", "expected_root_to_tip_substitution_probability"]:
            if config["simulator"].get(key) != value:
                raise ValueError("Changed v2 simulation field: " + key)
    if config["mafft"] != original["mafft"]:
        raise ValueError("Phase 1 must retain the exact v2 MAFFT flags")
    if config["conditions"] != CONDITIONS:
        raise ValueError("Exactly the three original conditions are required")
    seeds = config["simulator"]["replicate_seeds"]
    levels = config["divergence_levels"]
    if len(set(seeds)) != len(seeds) or any(s not in range(2001,2051) for s in seeds):
        raise ValueError("Only unique new seeds 2001-2050 are admissible")
    if not test_design and (seeds != list(range(2001,2051)) or levels != [.25,.40]):
        raise ValueError("The full experiment requires all 50 seeds and both divergence levels")
    if config["analysis"]["primary_iteration_cap"] != 2 or config["analysis"]["primary_contrasts"] != ["partial_minus_none","high_minus_none"]:
        raise ValueError("Frozen primary analysis mismatch")
    if not test_design:
        if config["expected_mafft_runs"] != 900 or config["expected_datasets"] != 300:
            raise ValueError("Full design must contain 900 alignments and 300 datasets")
    for name, expected in config["provenance"]["unchanged_core_files"].items():
        if sha(PROJECT / "src" / name) != expected:
            raise ValueError("Frozen core source was changed: " + name)
    if sha(PROJECT / "config/pilot.json") != config["provenance"]["phase0_config_sha256"]:
        raise ValueError("Original pilot config changed")


def batch_label(p):
    return f"div{round(p*100):03d}"


def experiment_identity(config, version):
    _, key = pipeline.identity(config, version)
    # Statistical/reporting changes do not invalidate valid raw alignments.
    key["phase1_orchestrator_sha256"] = sha(__file__)
    return config["experiment_id"] + "_" + digest(canonical(key))[:16], key


def prepare(config_path, workdir, test_design=False):
    config = read_json(config_path)
    check_design(config, test_design)
    prefix = runner.executable_prefix()
    version = runner.version(prefix)
    if not test_design and ("FAKE" in version.upper() or any(k.startswith(("PILOT_TEST_", "FULL_TEST_", "FAKE_")) for k in os.environ)):
        raise ValueError("Test wrappers/fault injections cannot be used for real Phase 1 experiments")
    name, identity = experiment_identity(config, version)
    root = Path(workdir).resolve() / name
    with writer_lock(root):
        manifest = root / "experiment_manifest.json"
        if manifest.exists():
            saved = read_json(manifest)
            if saved["identity"] != identity or sha(root / "config.json") != saved["frozen_config_file_sha256"]:
                raise ValueError("Frozen experiment identity/config mismatch; refusing stale reuse")
            if sha(root / "analysis_plan.md") != saved["analysis_plan_sha256"]:
                raise ValueError("Frozen analysis plan changed")
        else:
            # Save exact original bytes as well as a canonical scientific identity.
            atomic(root / "config.json", Path(config_path).read_bytes())
            atomic(root / "analysis_plan.md", (PROJECT / "PHASE1_DESIGN.md").read_bytes())
            for folder in ["src", "tests", "config"]:
                for path in (PROJECT / folder).glob("*"):
                    if path.is_file():
                        atomic(root / "frozen_source" / folder / path.name, path.read_bytes())
            json_write(manifest, {"experiment_id": name, "phase": "phase1", "identity": identity,
                       "frozen_config_file_sha256": sha(root / "config.json"),
                       "analysis_plan_sha256": sha(root / "analysis_plan.md"),
                       "frozen_utc_epoch": time.time(), "python": platform.python_version(),
                       "platform": platform.platform(), "mafft_version": version, "mafft_prefix": prefix,
                       "source_sha256": source_hashes(), "planned_alignments": len(config["divergence_levels"])*len(config["simulator"]["replicate_seeds"])*9})
            pipeline.event(root, "freeze", "completed", frozen_config_file_sha256=sha(root / "config.json"))
        validation = PROJECT / "runs/phase1_validation/latest.json"
        if validation.exists() and not (root / "validation_before_execution.json").exists():
            atomic(root / "validation_before_execution.json", validation.read_bytes())
    return root, config, prefix, version


def execute(config_path, workdir, retry=False, test_design=False, finalize=True):
    root, config, prefix, version = prepare(config_path,workdir,test_design)
    with writer_lock(root):
        pipeline.event(root,"phase1","invoked",argv=sys.argv,retry_failed=retry)
        # Retry is an operational parameter, never part of data dependencies/config hashes.
        tasks = [(p,seed) for seed in config["simulator"]["replicate_seeds"] for p in config["divergence_levels"]]
        results=[]
        # Thread workers orchestrate subprocesses; each block has disjoint artifacts and counters.
        with ThreadPoolExecutor(max_workers=config["run"]["concurrent_seed_blocks"]) as executor:
            futures={executor.submit(run_block,root,config,p,seed,prefix,version,retry):(p,seed) for p,seed in tasks}
            for future in as_completed(futures):
                p,seed=futures[future]
                try:
                    record=future.result()
                    results.append(record)
                    json_write(root / "block_invocations" / (str(time.time_ns())+f"_{seed}_{batch_label(p)}.json"),record)
                    total=sum(r["alignments"].get("completed",0) for r in results)
                    skipped=sum(r["alignments"].get("skipped",0) for r in results)
                    failed=sum(r["alignments"].get("failed",0)+r["alignments"].get("failed_skipped",0) for r in results)
                    print(f"Phase 1 blocks {len(results)}/{len(tasks)}; {batch_label(p)} seed {seed}; alignment stages newly complete={total}, skipped={skipped}, failed={failed}",flush=True)
                except Exception as exc:
                    pipeline.failure(root,root / "orchestrator_failures",f"{batch_label(p)}_seed{seed}",exc)
        summary={"utc_epoch":time.time(),"blocks":results,"retry_failed":retry}
        json_write(root / "last_invocation.json",summary)
        json_write(root / "invocations" / (str(time.time_ns())+".json"),summary)
        if finalize:
            from .full_reporting import build
            build(root)
    return root


def run_block(root, config, p, seed, prefix, version, retry):
    # Each worker shares only immutable configuration; artifact paths and counters are disjoint.
    simulated, alignments, scores = Counter(),Counter(),Counter()
    batch=batch_label(p)
    directory=root/"blocks"/f"{batch}_seed{seed}"
    path=directory/"realization.json"
    deps={"seed":seed,"p":p,"config":digest(canonical(config)),"simulator":sha(simulation.__file__)}
    def generate():
        json_write(path,simulation.realization(config,seed,p))
    if pipeline.stage(root,directory,"block_done.json",deps,[path],generate,lambda: read_json(path)["substitution_only_tips"],retry,simulated):
        shared=read_json(path)
        for cond in config["conditions"]:
            case=pipeline.case_path(root,batch,seed,cond)
            if pipeline.dataset(root,case,config,seed,cond,p,retry,simulated,shared):
                for cap in config["mafft"]["maxiterate"]:
                    if pipeline.alignment(root,case,config,cap,prefix,version,retry,alignments):
                        if os.environ.get("FULL_TEST_PAUSE_AFTER_ALIGN") and alignments["completed"]==int(os.environ["FULL_TEST_PAUSE_AFTER_ALIGN"]):
                            atomic(root/"test_pause_ready","ready")
                            while True:
                                time.sleep(.1)
                        pipeline.score_case(root,case,cap,retry,scores)
    return {"batch":batch,"seed":seed,"simulation":dict(simulated),"alignments":dict(alignments),"scores":dict(scores)}
