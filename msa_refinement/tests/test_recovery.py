import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import pytest
from src import pipeline
from src.storage import read_json, json_write, sha, writer_lock

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def setup(tmp_path, monkeypatch):
    from src import reporting
    monkeypatch.setattr(reporting, "plots", lambda *args: None)
    config = read_json(ROOT / "config/pilot.json")
    config["simulator"]["replicate_seeds"] = [1001]
    config["simulator"]["second_deletion_start_delta_by_condition"] = {"none":50}
    config["predeclared_calibration"]["trigger_if_fewer_than_n_of_12_cases_with_any_local_true_pair_correction_or_loss"] = 0
    config["predeclared_calibration"]["trigger_if_at_least_n_of_12_baselines_have_no_local_missed_reference_pairs"] = 100
    path = tmp_path / "config.json"
    json_write(path, config)
    count = tmp_path / "invocations.jsonl"
    monkeypatch.setenv("PILOT_MAFFT_COMMAND", json.dumps([sys.executable,str(ROOT / "tests/fake_mafft.py")]))
    monkeypatch.setenv("FAKE_COUNT", str(count))
    def run(retry=False):
        return pipeline.pilot(path, tmp_path / "runs", retry)
    return path, count, run


def calls(count):
    return count.read_text().count("\n") if count.exists() else 0


def raw(root):
    return sorted(root.glob("cases/*/alignments/*/aligned.fasta"))


def snapshot(paths):
    return {str(p): (p.read_bytes(),p.stat().st_mtime_ns) for p in paths}


def test_1_kill_after_two_completed_jobs_restart_counts(setup, monkeypatch):
    path, count, run = setup
    monkeypatch.setenv("PILOT_TEST_PAUSE_AFTER_ALIGN", "2")
    # Run the actual pipeline in a killable OS child. Fake MAFFT is a subprocess too.
    code = "from src.pipeline import pilot; import sys; pilot(sys.argv[1],sys.argv[2])"
    child = subprocess.Popen([sys.executable,"-c",code,str(path),str(path.parent/"runs")], cwd=ROOT,
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    deadline = time.monotonic()+40
    try:
        while not list((path.parent / "runs").glob("*/test_pause_ready")) and child.poll() is None and time.monotonic()<deadline:
            time.sleep(.05)
        assert list((path.parent / "runs").glob("*/test_pause_ready")), child.communicate(timeout=1) if child.poll() is not None else "no pause marker"
        child.kill()
        child.communicate(timeout=10)
    finally:
        if child.poll() is None:
            child.kill(); child.communicate()
    assert calls(count)==2
    root = next((path.parent / "runs").iterdir())
    before = snapshot(raw(root))
    monkeypatch.delenv("PILOT_TEST_PAUSE_AFTER_ALIGN")
    run()
    assert calls(count)==3
    assert all(snapshot([Path(p)])[p]==v for p,v in before.items())
    assert len(raw(root))==3


def test_2_scoring_failure_preserves_raw_bytes_and_timestamps(setup, monkeypatch):
    path,count,run = setup
    monkeypatch.setenv("PILOT_TEST_SCORE_FAIL","1")
    root=run()
    before=snapshot(raw(root))
    assert calls(count)==3
    assert list(root.glob("cases/*/scores/failures/*.json"))
    monkeypatch.delenv("PILOT_TEST_SCORE_FAIL")
    run(retry=True)
    assert calls(count)==3 and snapshot(raw(root))==before
    assert len(list(root.glob("cases/*/scores/maxiter_???.json")))==3


def test_3_truncated_partial_only_affected_stage_reruns(setup):
    path,count,run=setup
    root=run()
    files=raw(root)
    unaffected=snapshot(files[1:])
    first=files[0]
    (first.parent/"run_done.json").unlink()
    first.unlink()
    (first.parent/"aligned.fasta.dead.partial").write_text(">protein1\nA")
    run()
    assert calls(count)==4 and snapshot(files[1:])==unaffected
    assert list(first.parent.glob("quarantine/*/*.partial"))


@pytest.mark.parametrize("damage", ["marker", "hash"])
def test_4_damaged_marker_or_mismatched_hash(setup,damage):
    path,count,run=setup
    root=run()
    first=raw(root)[0]
    if damage=="marker":
        (first.parent/"run_done.json").write_text("{broken")
    else:
        first.write_text(first.read_text()+"\n")
    run()
    assert calls(count)==4


def test_5_second_run_zero_mafft_identical_summary_no_duplicate_rows(setup):
    path,count,run=setup
    root=run()
    before=(root/"pilot_summary.csv").read_bytes()
    raws=snapshot(raw(root))
    run()
    assert calls(count)==3 and (root/"pilot_summary.csv").read_bytes()==before
    assert len(before.decode().splitlines())==4
    assert snapshot(raw(root))==raws
    assert read_json(root/"last_invocation.json")["batches"]["initial"]["alignments"]=={"skipped":3}


def test_6_changed_config_isolated_and_original_untouched(setup):
    path,count,run=setup
    root=run()
    before=snapshot(raw(root))
    config=read_json(path)
    config["simulator"]["replicate_seeds"]=[1002]
    json_write(path,config)
    other=run()
    assert root!=other and calls(count)==6 and snapshot(raw(root))==before


@pytest.mark.parametrize("mode", ["fail","timeout","invalid"])
def test_7_failures_timeouts_invalid_keep_logs_and_retry(setup,monkeypatch,mode):
    path,count,run=setup
    config=read_json(path)
    config["run"]["timeout_seconds_per_mafft"]=2 if mode=="timeout" else 10
    json_write(path,config)
    monkeypatch.setenv("FAKE_MODE",mode)
    root=run()
    assert not raw(root)
    assert calls(count)==3
    failures=list(root.glob("cases/*/alignments/*/failures/*.json"))
    assert len(failures)==3
    assert all(read_json(p)["traceback"] and read_json(p)["type"] for p in failures)
    assert list(root.glob("cases/*/alignments/*/*.partial"))
    logs=list(root.glob("cases/*/alignments/*/stderr.log"))+list(root.glob("cases/*/alignments/*/stderr_*.partial"))
    assert len(logs)==3
    run() # ordinary resume does not retry failed jobs
    assert calls(count)==3
    monkeypatch.delenv("FAKE_MODE")
    run(retry=True)
    assert calls(count)==6 and len(raw(root))==3


def test_exclusive_writer_lock_and_release_after_exit(tmp_path):
    with writer_lock(tmp_path):
        with pytest.raises(RuntimeError,match="writer lock"):
            with writer_lock(tmp_path):
                pass
    with writer_lock(tmp_path):
        pass


def test_scoring_only_change_reuses_raw_alignments(setup,monkeypatch):
    path,count,run=setup
    root=run()
    before=snapshot(raw(root))
    monkeypatch.setattr(pipeline.scoring,"SCORING_VERSION","test-new-scoring-version")
    run()
    assert calls(count)==3 and snapshot(raw(root))==before
    assert all(read_json(p)["dependencies"]["scoring_version"]=="test-new-scoring-version" for p in root.glob("cases/*/scores/maxiter_???.json"))


def test_numeric_trigger_thresholds(tmp_path):
    config=read_json(ROOT/"config/pilot.json")
    for changed,ceilings,expected in [(1,0,True),(2,9,False),(2,10,True),(12,0,False),(0,12,True)]:
        observations = [(0 if i<ceilings else 1, [(int(i<changed),0),(0,0)]) for i in range(12)]
        assert pipeline.evaluate_signal("initial",observations,config)["low_signal"]==expected
