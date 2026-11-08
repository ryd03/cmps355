import copy
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time
import pytest
from src import full_pipeline as full, full_analysis as analysis, full_reporting, simulation
from src.storage import read_json,json_write,sha

ROOT=Path(__file__).resolve().parents[1]


def test_full_frozen_design_exact_900_and_original_science():
    config=read_json(ROOT/"config/full.json")
    full.check_design(config)
    assert config["simulator"]["replicate_seeds"]==list(range(2001,2051))
    schedule={(p,s,c,k) for p in config["divergence_levels"] for s in config["simulator"]["replicate_seeds"] for c in config["conditions"] for k in config["mafft"]["maxiterate"]}
    assert len(schedule)==900
    assert {s for p,s,c,k in schedule}.isdisjoint({1001,1002,1003,1004})
    for key in ["topology_newick_unit_scale","deletion_length_per_event","focal_1_deletion_anchor_1based","second_deletion_start_delta_by_condition","region_window_ancestral_positions"]:
        changed=copy.deepcopy(config)
        changed["simulator"][key]="invalid"
        with pytest.raises(ValueError,match="Changed v2"):
            full.check_design(changed)


def test_bootstrap_quantiles_manual_and_constant():
    # All four size-two resamples of blocks [1,3] have means 1,2,2,3.
    result,draws=analysis.interval([1,3],[[0,0],[0,1],[1,0],[1,1]])
    assert draws==[1,2,2,3]
    assert result["estimate"]==2
    assert result["ci_lower"]==pytest.approx(1.075)
    assert result["ci_upper"]==pytest.approx(2.925)
    indices=analysis.resample_indices(50,100,77)
    constant,_=analysis.interval([.125]*50,indices)
    assert constant["estimate"]==constant["ci_lower"]==constant["ci_upper"]==.125
    assert indices==analysis.resample_indices(50,100,77)
    assert len(indices)==100 and all(len(row)==50 and all(0<=i<50 for i in row) for row in indices)
    with pytest.raises(ValueError):
        analysis.interval([None],[ [0] ])


def synthetic_rows():
    result=[]
    for p in [.25,.4]:
        for seed in range(2001,2051):
            for cond,denom,gain in [("none",440,11),("partial",420,11),("high",420,22)]:
                for cap in [0,2,10]:
                    multiplier=2 if p==.4 else 1
                    result.append({"divergence":p,"seed":seed,"condition":cond,"iteration_cap":cap,
                        "local_net_gain":gain*multiplier/denom,"local_recovered":gain*multiplier,"local_lost":0,
                        "local_reference_denominator":denom,"local_baseline_missed_denominator":100,"local_baseline_correct_denominator":denom-100})
    return result


def test_normalize_each_condition_then_pair_no_pilot():
    config=read_json(ROOT/"config/full.json")
    rows=synthetic_rows()
    indexed,records=analysis.paired_records(list(reversed(rows)),config)
    assert len(records)==400
    r=next(r for r in records if r["divergence"]==.25 and r["iteration_cap"]==2 and r["contrast"]=="partial_minus_none")
    assert r["paired_difference"]==pytest.approx(11/420-11/440)
    assert r["paired_difference"]!=0 # equal raw gains are unequal normalized gains
    config["analysis"]["bootstrap_replicates"]=100
    paired,intervals=analysis.analyze(rows,config)
    primary=[r for r in intervals if r["analysis_role"]=="primary"]
    assert len(primary)==4 and all(r["n_seed_blocks"]==50 for r in intervals)
    assert all(r["estimate"]==pytest.approx(r["ci_lower"])==pytest.approx(r["ci_upper"]) for r in primary)
    for wrong in [rows[:-1],rows+[rows[0]],rows+[dict(rows[0],seed=1001)]]:
        with pytest.raises(ValueError):
            analysis.paired_records(wrong,config)


def test_bootstrap_preserves_pairing_instead_of_independent_arms():
    indices=analysis.resample_indices(4,100,355)
    background=[-1000,20,500,4000]
    treatment=[v+.5 for v in background]
    differences=[t-b for t,b in zip(treatment,background)]
    estimate,_=analysis.interval(differences,indices)
    assert estimate["ci_lower"]==estimate["ci_upper"]==.5


@pytest.mark.parametrize("p",[.25,.4])
def test_all_new_seed_block_geometry_and_shared_histories(p):
    config=read_json(ROOT/"config/full.json")
    for seed in range(2001,2051):
        shared=simulation.realization(config,seed,p)
        metas=[]
        for cond,overlap,union in [("none",0,24),("partial",6,18),("high",11,13)]:
            inp,truth,meta=simulation.simulate(config,seed,cond,p,shared)
            assert meta["realized_overlap"]==overlap and meta["union_gap_columns"]==union
            assert all(truth[s].replace('-','')==seq for s,seq in inp.items())
            metas.append(meta)
        assert all(m["substitution_only_tips"]==shared["substitution_only_tips"] for m in metas)
        assert len({m["events"][0]["start"] for m in metas})==1


@pytest.fixture
def full_setup(tmp_path,monkeypatch):
    config=read_json(ROOT/"config/full.json")
    config["experiment_id"]="testfull"
    config["simulator"]["replicate_seeds"]=[2001]
    config["divergence_levels"]=[.25]
    config["run"]["concurrent_seed_blocks"]=1
    path=tmp_path/"config.json"; json_write(path,config)
    count=tmp_path/"calls.jsonl"
    monkeypatch.setenv("PILOT_MAFFT_COMMAND",json.dumps([sys.executable,str(ROOT/"tests/fake_mafft.py")]))
    monkeypatch.setenv("FAKE_COUNT",str(count))
    def run(retry=False):
        return full.execute(path,tmp_path/"runs",retry=retry,test_design=True,finalize=False)
    return path,count,run


def count_calls(path):
    return len(path.read_text().splitlines()) if path.exists() else 0


def raw_snapshot(root):
    return {p.relative_to(root).as_posix():(sha(p),p.stat().st_mtime_ns) for p in root.glob("cases/*/alignments/*/aligned.fasta")}


def test_phase1_scoring_failure_raw_reuse_and_incomplete_analysis_gate(full_setup,monkeypatch):
    path,count,run=full_setup
    monkeypatch.setenv("PILOT_TEST_SCORE_FAIL","1")
    root=run(); before=raw_snapshot(root)
    assert count_calls(count)==9 and len(before)==9
    monkeypatch.delenv("PILOT_TEST_SCORE_FAIL")
    run(retry=True)
    assert count_calls(count)==9 and raw_snapshot(root)==before
    jobs,rows=full_reporting.build(root)
    assert len(jobs)==9 and len(rows)==9
    assert not (root/"confidence_intervals.csv").exists()
    assert all(j["score_complete"] and j["comparison_complete"] for j in jobs)
    raw=(root/"full_summary.csv").read_bytes()
    run(); full_reporting.build(root)
    assert raw_snapshot(root)==before and count_calls(count)==9
    assert (root/"full_summary.csv").read_bytes()==raw


def test_phase1_concurrent_blocks_disjoint_checkpoints(full_setup):
    path,count,run=full_setup
    config=read_json(path)
    config["simulator"]["replicate_seeds"]=[2001,2002]
    config["run"]["concurrent_seed_blocks"]=2
    json_write(path,config)
    root=run()
    jobs,rows=full_reporting.inventory(root)
    assert len(jobs)==18 and len(rows)==18 and count_calls(count)==18
    assert all(j["score_complete"] and j["comparison_complete"] for j in jobs)
    before=raw_snapshot(root)
    run()
    assert count_calls(count)==18 and raw_snapshot(root)==before


def test_phase1_actual_process_kill_and_resume_two_successes(full_setup,monkeypatch):
    path,count,run=full_setup
    monkeypatch.setenv("FULL_TEST_PAUSE_AFTER_ALIGN","2")
    code="from src.full_pipeline import execute; import sys; execute(sys.argv[1],sys.argv[2],test_design=True,finalize=False)"
    process=subprocess.Popen([sys.executable,"-c",code,str(path),str(path.parent/"runs")],cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    deadline=time.monotonic()+40
    try:
        while not list((path.parent/"runs").glob("*/test_pause_ready")) and process.poll() is None and time.monotonic()<deadline:
            time.sleep(.05)
        assert list((path.parent/"runs").glob("*/test_pause_ready")),process.communicate(timeout=1) if process.poll() is not None else "pause not reached"
        process.kill(); process.communicate(timeout=10)
    finally:
        if process.poll() is None:
            process.kill(); process.communicate()
    root=next((path.parent/"runs").glob("testfull_*"))
    before=raw_snapshot(root)
    assert count_calls(count)==2 and len(before)==2
    monkeypatch.delenv("FULL_TEST_PAUSE_AFTER_ALIGN")
    run()
    assert count_calls(count)==9
    after=raw_snapshot(root)
    assert all(after[key]==value for key,value in before.items())


def test_phase1_comparison_truth_hash_damage_rejected(full_setup):
    path,count,run=full_setup
    root=run()
    comp=next(root.glob("cases/*/comparisons/against_maxiter_002.json"))
    marker=comp.with_name(comp.name+".done.json")
    marker.write_text('{"files":{}}')
    jobs,rows=full_reporting.inventory(root)
    assert sum(not j["comparison_complete"] for j in jobs)==1
    before=raw_snapshot(root)
    run()
    assert raw_snapshot(root)==before and count_calls(count)==9
    assert all(j["comparison_complete"] for j in full_reporting.inventory(root)[0])


def test_phase1_status_available_under_writer_lock(full_setup,monkeypatch,capsys):
    from src import full_cli
    from src.storage import writer_lock
    path,count,run=full_setup
    root=run()
    monkeypatch.setattr(full_cli,"roots",lambda workdir: iter([root]))
    monkeypatch.setattr(sys,"argv",["full_cli","status","--workdir",str(path.parent/"runs")])
    with writer_lock(root):
        assert full_cli.main()==0
    assert "9 completed" in capsys.readouterr().out


def test_relocated_archive_uses_original_preservation_evidence(tmp_path,monkeypatch):
    from src import full_audit
    monkeypatch.setattr(full_audit,"PROJECT",tmp_path)
    baseline=tmp_path/"runs/phase1_provenance/phase0_before.json"
    json_write(baseline,{"files":{"unavailable_original_pilot.zip":{"sha256":"original","mtime_ns":1}}})
    root=tmp_path/"new"
    json_write(root/"phase0_preservation.json",{"passed":True,"workspace_project_path":"original_workspace",
               "baseline_sha256":sha(baseline),"phase0_files_checked":1})
    result=full_audit.preserve_phase0(root)
    assert result["passed"] and result["current_workspace_checked"] is False
    assert result["check_scope"]=="verified_original_workspace_evidence_only"


def test_original_workspace_preservation_detects_changed_file(tmp_path,monkeypatch):
    from src import full_audit
    monkeypatch.setattr(full_audit,"PROJECT",tmp_path)
    path=tmp_path/"original.txt"; path.write_text("original")
    json_write(tmp_path/"runs/phase1_provenance/phase0_before.json",{"files":{"original.txt":{"sha256":sha(path),"mtime_ns":path.stat().st_mtime_ns}}})
    root=tmp_path/"new"
    assert full_audit.preserve_phase0(root)["passed"]
    path.write_text("changed")
    with pytest.raises(ValueError,match="Phase 0 changed"):
        full_audit.preserve_phase0(root)


def test_full_bundle_relative_root_zip_roundtrip(tmp_path,monkeypatch):
    """Archive integration fixture only; these are not experimental alignments."""
    import zipfile
    from src import full_bundle
    from src.storage import digest
    monkeypatch.setattr(full_bundle,"PROJECT",tmp_path)
    monkeypatch.chdir(tmp_path)
    root=Path("runs/full_test_archive")
    for name in ["full_artifact_audit.json","full_resume_verification.json","phase0_preservation.json"]:
        json_write(root/name,{"passed":True,"test_fixture_only":True})
    for name in ["full_report.md","full_summary.csv","paired_comparisons.csv","confidence_intervals.csv","cell_summary.csv","full_job_manifest.csv"]:
        path=root/name; path.parent.mkdir(parents=True,exist_ok=True); path.write_text("test fixture only\n")
    (tmp_path/"tools").mkdir()
    for name in ["setup_mafft.py","installation.json"]:
        (tmp_path/"tools"/name).write_text("test fixture only\n")
    (tmp_path/"PHASE1_REPRODUCE.md").write_text("# Test reproduction README\n")
    json_write(tmp_path/"runs/phase1_provenance/phase0_before.json",{"test_fixture_only":True})
    for i in range(900):
        path=root/"cases"/f"test{i:03d}"/"aligned.fasta"
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(">test_fixture_only\nA\n")
    target=full_bundle.bundle(root)
    with zipfile.ZipFile(target) as archive:
        assert archive.testzip() is None
        manifest=json.loads(archive.read("bundle_sha256_manifest.json"))
        assert len([name for name in manifest if name.endswith("/aligned.fasta")])==900
        assert archive.read("README.md").startswith(b"# Test reproduction")
        assert all(digest(archive.read(name))==expected for name,expected in manifest.items())
