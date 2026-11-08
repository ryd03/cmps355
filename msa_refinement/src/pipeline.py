"""Separate verified simulation, alignment, scoring, and comparison checkpoints."""
from collections import Counter
import copy
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
import uuid
from . import simulation, scoring, runner
from .storage import (atomic, canonical, complete, digest, json_write, marker, quarantine,
                      read_json, sha, writer_lock)


def code_hash(module):
    return sha(module.__file__)


def identity(config, mafft_version):
    key = {"config": config, "simulation_code": code_hash(simulation),
           "runner_code": code_hash(runner), "runner_version": runner.RUNNER_VERSION,
           "simulator_version": simulation.SIMULATOR_VERSION,
           "mafft_version": mafft_version, "python": platform.python_version()}
    return config["experiment_id"] + "_" + digest(canonical(key))[:16], key


def event(root, stage, state, **details):
    # Individual immutable events are authoritative; JSONL is only a convenience.
    obj = {"stage": stage, "state": state, "utc_epoch": time.time(), **details}
    json_write(root / "execution_events" / (str(time.time_ns()) + "_" + uuid.uuid4().hex + ".json"), obj)
    with (root / "events.jsonl").open("a", encoding="utf-8") as out:
        out.write(canonical(obj).decode() + "\n")


def failure(root, directory, stage, exc):
    details = {"stage": stage, "type": type(exc).__name__, "error": str(exc),
               "traceback": traceback.format_exc(), "utc_epoch": time.time(),
               "returncode": getattr(exc, "returncode", None), "command": getattr(exc, "cmd", None),
               "preserved_partial_artifacts": [p.name for p in directory.glob("*.partial")]}
    json_write(directory / "failures" / (str(time.time_ns()) + ".json"), details)
    json_write(directory / "failure.json", details)
    event(root, stage, "failed", **{k: v for k, v in details.items() if k != "stage"})


def stage(root, directory, name, dependencies, files, action, validator, retry, counts):
    directory.mkdir(parents=True, exist_ok=True)
    done = directory / name
    stage_id = str(done.relative_to(root))
    if complete(done, dependencies, validator):
        counts["skipped"] += 1
        event(root, stage_id, "skipped")
        return True
    if (directory / (name + ".failed.json")).exists() and not retry:
        counts["failed_skipped"] += 1
        event(root, stage_id, "failed_skipped", reason="use --retry-failed")
        return False
    interrupted = done.exists() or any(directory.glob("*.partial")) or (directory / (name + ".running.json")).exists()
    quarantine(directory, [done, *files])
    json_write(directory / (name + ".running.json"), {"pid": os.getpid(), "utc_epoch": time.time()})
    event(root, stage_id, "started", interrupted_or_invalid=interrupted)
    try:
        details = action() or {}
        validator()
        marker(done, files, dependencies, details)
        for p in [directory / (name + ".running.json"), directory / (name + ".failed.json")]:
            p.unlink(missing_ok=True)
        counts["completed"] += 1
        event(root, stage_id, "completed", **details)
        return True
    except Exception as exc:
        failure(root, directory, stage_id, exc)
        json_write(directory / (name + ".failed.json"), {"error": str(exc), "type": type(exc).__name__})
        counts["failed"] += 1
        return False


def case_path(root, batch, seed, condition):
    return root / "cases" / f"{batch}_seed{seed}_{condition}"


def dataset(root, case, config, seed, condition, p, retry, counts, shared=None):
    files = [case / f for f in ("input.fasta", "truth_alignment.fasta", "ancestry.json")]
    deps = {"config": digest(canonical(config)), "seed": seed, "condition": condition, "p": p,
            "simulator": code_hash(simulation)}

    def validate():
        inp, truth = [scoring.parse_fasta(f.read_text()) for f in files[:2]]
        simulation.validate_dataset(inp, truth, read_json(files[2]))

    def generate():
        inp, truth, meta = simulation.simulate(config, seed, condition, p, shared)
        atomic(files[0], simulation.fasta(inp))
        atomic(files[1], simulation.fasta(truth))
        json_write(files[2], meta)

    return stage(root, case, "simulation_done.json", deps, files, generate, validate, retry, counts)


def alignment(root, case, config, cap, prefix, mafft_version, retry, counts):
    directory = case / "alignments" / f"maxiter_{cap:03d}"
    output, stderr = directory / "aligned.fasta", directory / "stderr.log"
    source = case / "input.fasta"
    argv = runner.command(prefix, source, cap, config)
    deps = {"input": sha(source), "flags": argv[len(prefix):-1], "mafft_version": mafft_version,
            "runner": code_hash(runner)}

    def validate():
        scoring.validate_alignment(scoring.parse_fasta(output.read_text()), scoring.parse_fasta(source.read_text()))

    def execute():
        tmp = output.with_name("output_" + uuid.uuid4().hex[:12] + ".partial")
        logtmp = stderr.with_name("stderr_" + uuid.uuid4().hex[:12] + ".partial")
        before = time.perf_counter()
        with tmp.open("wb") as out, logtmp.open("wb") as err:
            process = subprocess.Popen(argv, stdout=out, stderr=err)
            try:
                code = process.wait(timeout=config["run"]["timeout_seconds_per_mafft"])
            except subprocess.TimeoutExpired:
                if os.name == "nt":
                    try:
                        subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"], capture_output=True, timeout=5)
                    except subprocess.TimeoutExpired:
                        pass
                    if process.poll() is None:
                        process.kill()
                else:
                    process.kill()
                process.wait()
                raise
            finally:
                out.flush(); err.flush()
                os.fsync(out.fileno()); os.fsync(err.fileno())
        os.replace(logtmp, stderr)
        if code:
            raise subprocess.CalledProcessError(code, argv)
        scoring.validate_alignment(scoring.parse_fasta(tmp.read_text()), scoring.parse_fasta(source.read_text()))
        os.replace(tmp, output)
        return {"command": argv, "mafft_version": mafft_version, "runtime_seconds": time.perf_counter() - before,
                "exit_code": code, "input_sha256": sha(source), "aligned_sha256": sha(output)}

    return stage(root, directory, "run_done.json", deps, [output, stderr], execute, validate, retry, counts)


def score_case(root, case, cap, retry, counts):
    score_dir = case / "scores"
    output = score_dir / f"maxiter_{cap:03d}.json"
    align = case / "alignments" / f"maxiter_{cap:03d}" / "aligned.fasta"
    truth_file, ancestry = case / "truth_alignment.fasta", case / "ancestry.json"
    deps = {"aligned": sha(align), "truth": sha(truth_file), "ancestry": sha(ancestry),
            "scoring_version": scoring.SCORING_VERSION, "scoring_code": code_hash(scoring)}

    def execute():
        # Explicit fault injection used only by automated recovery tests.
        if os.environ.get("PILOT_TEST_SCORE_FAIL") == "1":
            raise RuntimeError("Injected scoring failure after completed alignment")
        pred = scoring.pairs(scoring.parse_fasta(align.read_text()))
        truth = scoring.pairs(scoring.parse_fasta(truth_file.read_text()))
        result = scoring.accuracy(pred, truth, scoring.boundary_residues(read_json(ancestry)))
        result.update({"aligned_length": len(next(iter(scoring.parse_fasta(align.read_text()).values()))),
                       "dependencies": deps})
        json_write(output, result)

    ok = stage(root, score_dir, output.name + ".done.json", deps, [output], execute,
               lambda: read_json(output)["global"], retry, counts)
    if not ok or cap == 0:
        return ok
    baseline = case / "alignments" / "maxiter_000" / "aligned.fasta"
    baseline_marker = baseline.parent / "run_done.json"
    from .reporting import verified
    if not verified(baseline_marker) or read_json(baseline_marker)["dependencies"]["input"] != sha(case / "input.fasta"):
        return False
    comp_dir = case / "comparisons"
    comp = comp_dir / f"against_maxiter_{cap:03d}.json"
    comp_deps = {**deps, "baseline": sha(baseline)}

    def comparison():
        b, r, t = [scoring.pairs(scoring.parse_fasta(f.read_text())) for f in [baseline, align, truth_file]]
        json_write(comp, {**scoring.compare(b, r, t, scoring.boundary_residues(read_json(ancestry))),
                          "dependencies": comp_deps})

    return stage(root, comp_dir, comp.name + ".done.json", comp_deps, [comp], comparison,
                 lambda: read_json(comp)["local"], retry, counts)


def evaluate_signal(batch, observations, config):
    corrections, ceilings = 0, 0
    for baseline_missed, changes in observations:
        corrections += any(recovered + lost > 0 for recovered, lost in changes)
        ceilings += baseline_missed == 0
    preset = config["predeclared_calibration"]
    a = preset["trigger_if_fewer_than_n_of_12_cases_with_any_local_true_pair_correction_or_loss"]
    b = preset["trigger_if_at_least_n_of_12_baselines_have_no_local_missed_reference_pairs"]
    return {"batch": batch, "complete": True, "datasets": len(observations), "datasets_with_local_correct_pair_change": corrections,
            "zero_baseline_missed_datasets": ceilings, "change_threshold_strict_less_than": a,
            "ceiling_threshold_at_least": b, "low_signal": corrections < a or ceilings >= b}


def decision(root, batch, config):
    from .reporting import inventory
    _, rows = inventory(root)
    observations = []
    for seed in config["simulator"]["replicate_seeds"]:
        for cond in config["simulator"]["second_deletion_start_delta_by_condition"]:
            selected = {r["iteration_cap"]: r for r in rows if r["batch"] == batch and r["seed"] == seed and r["condition"] == cond}
            if set(selected) != {0, 2, 10} or any("local_recovered" not in selected[k] for k in [2, 10]):
                return {"batch": batch, "complete": False, "low_signal": None}
            observations.append((selected[0]["local_missed_pairs"], [(selected[k]["local_recovered"], selected[k]["local_lost"]) for k in [2,10]]))
    return evaluate_signal(batch, observations, config)


def run_batch(root, batch, p, config, prefix, mafft_version, retry, counts):
    align_count = Counter()
    simulated = Counter()
    for seed in config["simulator"]["replicate_seeds"]:
        block = root / "blocks" / f"{batch}_seed{seed}"
        block_file = block / "realization.json"
        block_deps = {"seed": seed, "p": p, "config": digest(canonical(config)), "simulator": code_hash(simulation)}
        def generate_block():
            json_write(block_file, simulation.realization(config, seed, p))
        if not stage(root, block, "block_done.json", block_deps, [block_file], generate_block,
                     lambda: read_json(block_file)["substitution_only_tips"], retry, simulated):
            continue
        shared = read_json(block_file)
        for cond in config["simulator"]["second_deletion_start_delta_by_condition"]:
            case = case_path(root, batch, seed, cond)
            if not dataset(root, case, config, seed, cond, p, retry, simulated, shared):
                continue
            for cap in config["mafft"]["maxiterate"]:
                if alignment(root, case, config, cap, prefix, mafft_version, retry, align_count):
                    if os.environ.get("PILOT_TEST_PAUSE_AFTER_ALIGN"):
                        n = int(os.environ["PILOT_TEST_PAUSE_AFTER_ALIGN"])
                        if align_count["completed"] == n:
                            atomic(root / "test_pause_ready", "ready")
                            while True:
                                time.sleep(0.1)
                    score_case(root, case, cap, retry, counts)
            print(f"{batch} seed={seed} {cond}: alignment stages {dict(align_count)}", flush=True)
    return {"alignments": dict(align_count), "simulation": dict(simulated)}


def pilot(config_path, workdir, retry=False):
    config = read_json(config_path)
    if config["phase"] != "phase0":
        raise ValueError("Only Phase 0 is implemented")
    prefix = runner.executable_prefix()
    mafft_version = runner.version(prefix)
    experiment_id, key = identity(config, mafft_version)
    root = Path(workdir).resolve() / experiment_id
    with writer_lock(root):
        manifest_path = root / "experiment_manifest.json"
        if manifest_path.exists():
            if read_json(manifest_path)["identity"] != key:
                raise ValueError("Namespace collision: refusing stale reuse")
        else:
            json_write(root / "config.json", config)
            json_write(manifest_path, {"experiment_id": experiment_id, "identity": key,
                       "config_sha256": digest(canonical(config)), "mafft_prefix": prefix,
                       "platform": platform.platform(), "python_executable": sys.executable,
                       "created_utc_epoch": time.time()})
        # Record source snapshot, scoring is intentionally outside experiment identity.
        for source in Path(__file__).parent.glob("*.py"):
            snapshot = root / "source_snapshots" / sha(source)[:16] / source.name
            if not snapshot.exists():
                atomic(snapshot, source.read_bytes())
        event(root, "pilot", "invoked", argv=sys.argv, retry_failed=retry)
        counts = Counter()
        batches = {}
        p = config["simulator"]["expected_root_to_tip_substitution_probability"]
        batches["initial"] = run_batch(root, "initial", p, config, prefix, mafft_version, retry, counts)
        first = decision(root, "initial", config)
        json_write(root / "initial_decision.json", first)
        event(root, "initial_decision", "evaluated", **first)
        if first["complete"] and first["low_signal"]:
            harder = config["predeclared_calibration"]["expected_root_to_tip_substitution_probability"]
            batches["harder"] = run_batch(root, "harder", harder, config, prefix, mafft_version, retry, counts)
            second = decision(root, "harder", config)
            json_write(root / "harder_decision.json", second)
            event(root, "harder_decision", "evaluated", **second)
        invocation = {"batches": batches, "score_stages": dict(counts), "utc_epoch": time.time()}
        json_write(root / "last_invocation.json", invocation)
        json_write(root / "invocations" / (str(time.time_ns()) + ".json"), invocation)
        from .reporting import build
        build(root)
    return root
