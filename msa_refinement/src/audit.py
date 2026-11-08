"""Independent final artifact verification; recomputes scores from raw alignments."""
from pathlib import Path
from . import scoring, simulation
from .pipeline import decision
from .reporting import inventory, verified, csv_text
from .storage import json_write, read_json, sha


def audit(root):
    root = Path(root)
    config = read_json(root / "config.json")
    jobs, rows = inventory(root)
    assert len(jobs) in [36,72], "Pilot must contain 36 jobs or a triggered 72 jobs"
    assert all(j["alignment_state"] == "completed" and j["score_complete"] for j in jobs)
    assert len(rows) == len(jobs)
    assert len({(r["batch"],r["seed"],r["condition"],r["iteration_cap"]) for r in rows}) == len(rows)
    assert (root / "pilot_summary.csv").read_text() == csv_text(rows)
    assert "FAKE" not in read_json(root / "experiment_manifest.json")["identity"]["mafft_version"]
    n=0
    for case in sorted((root / "cases").iterdir()):
        meta = read_json(case / "ancestry.json")
        inp = scoring.parse_fasta((case / "input.fasta").read_text())
        truth_rows = scoring.parse_fasta((case / "truth_alignment.fasta").read_text())
        simulation.validate_dataset(inp,truth_rows,meta)
        target_p = config["predeclared_calibration"]["expected_root_to_tip_substitution_probability"] if case.name.startswith("harder_") else config["simulator"]["expected_root_to_tip_substitution_probability"]
        expected = simulation.simulate(config, meta["seed"], meta["condition"], target_p)
        assert inp == expected[0] and truth_rows == expected[1] and meta == expected[2]
        t = scoring.pairs(truth_rows)
        boundary = scoring.boundary_residues(meta)
        modes = {}
        for cap in [0,2,10]:
            d = case / "alignments" / f"maxiter_{cap:03d}"
            rows_pred = scoring.parse_fasta((d / "aligned.fasta").read_text())
            scoring.validate_alignment(rows_pred,inp)
            modes[cap] = scoring.pairs(rows_pred)
            score = read_json(case / "scores" / f"maxiter_{cap:03d}.json")
            actual = scoring.accuracy(modes[cap],t,boundary)
            assert all(score[k]==v for k,v in actual.items())
            argv = read_json(d / "run_done.json")["details"]["command"]
            assert argv[-8:-1] == ["--amino","--thread","1","--retree","2","--maxiterate",str(cap)]
            for k in ["reference_pairs","predicted_pairs","correct_pairs","false_pairs","missed_pairs"]:
                assert score["global"][k] == score["local"][k] + score["far"][k]
        for cap in [2,10]:
            comparison = read_json(case / "comparisons" / f"against_maxiter_{cap:03d}.json")
            actual = scoring.compare(modes[0],modes[cap],t,boundary)
            assert all(comparison[k] == v for k,v in actual.items())
        n+=1
    for batch in ["initial","harder"]:
        path = root / (batch + "_decision.json")
        if path.exists():
            assert read_json(path) == decision(root,batch,config)
    initial = read_json(root / "initial_decision.json")
    assert (root / "harder_decision.json").exists() == initial["low_signal"]
    if len(jobs)==72:
        assert n==24
    else:
        assert n==12
    result = {"passed":True,"datasets_verified":n,"alignment_jobs_verified":len(jobs),
              "scores_recomputed":len(jobs),"comparisons_recomputed":2*n,
              "summary_unique_rows":len(rows),"configuration_and_seeds_reproduced":True,
              "trigger_verified":True,"source_sha256":{p.name:sha(p) for p in Path(__file__).parent.glob("*.py")}}
    json_write(root / "artifact_audit.json",result)
    return result
