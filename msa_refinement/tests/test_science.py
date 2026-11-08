import copy
import itertools
import math
from pathlib import Path
import random
import pytest
from src import scoring as s, simulation as sim
from src.pipeline import identity
from src.storage import canonical, digest, read_json

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def config():
    return read_json(ROOT / "config/pilot.json")


def test_supplied_toy_exact():
    t = s.pairs(s.parse_fasta((ROOT / "data/toy/reference.fasta").read_text()))
    b = s.pairs(s.parse_fasta((ROOT / "data/toy/incorrect.fasta").read_text()))
    assert (len(t), len(b), len(t & b), len(b-t), len(t-b)) == (29, 30, 28, 2, 1)
    c = s.compare(b, t, t, set())["global"]
    assert (c["recovered"], c["lost"], c["false_removed"], c["false_introduced"]) == (1, 0, 2, 0)


def test_reordering_and_all_gap_columns():
    rows = s.parse_fasta((ROOT / "data/toy/reference.fasta").read_text())
    altered = {k: "--" + v[:4] + "-" + v[4:] + "--" for k,v in reversed(list(rows.items()))}
    assert s.pairs(rows) == s.pairs(altered)


@pytest.mark.parametrize("reverse", [False, True])
def test_two_hand_worked_examples(reverse):
    # All letters identical: only ancestral coordinates define truth.
    # True columns: A1/B1/C1; B2/C2; A2/B3/C3 => 7 pairs.
    # Wrong: A1/B1/C1; A2/B2/C2; B3/C3 => 7 pairs, TP5 FP2 FN2.
    truth = s.pairs({"protein1": "A-A", "b": "AAA", "c": "AAA"})
    wrong = s.pairs({"protein1": "AA-", "b": "AAA", "c": "AAA"})
    boundary = {("protein1", 2)}
    b, r = (truth, wrong) if reverse else (wrong, truth)
    c = s.compare(b, r, truth, boundary)
    assert len(truth) == len(wrong) == 7
    assert len(truth & wrong) == 5
    assert c["local"]["reference_denominator"] == 2
    assert c["local"]["recovered"] == (0 if reverse else 2)
    assert c["local"]["lost"] == (2 if reverse else 0)
    assert c["local"]["false_introduced"] == (2 if reverse else 0)
    assert c["local"]["false_removed"] == (0 if reverse else 2)
    assert c["local"]["net_gain"] == (-1 if reverse else 1)
    assert c["local"]["recovery_rate"] == (None if reverse else 1)
    assert c["local"]["loss_rate"] == (1 if reverse else None)
    for key in ["recovered", "lost", "false_introduced", "false_removed", "reference_denominator", "baseline_missed_denominator", "baseline_correct_denominator"]:
        assert c["global"][key] == c["local"][key] + c["far"][key]
    a = s.accuracy(wrong, truth, boundary)
    for key in ["reference_pairs", "predicted_pairs", "correct_pairs", "false_pairs", "missed_pairs"]:
        assert a["global"][key] == a["local"][key] + a["far"][key]
    assert s.compare(b, r, truth, set())["local"]["net_gain"] is None
    assert s.accuracy(set(), set(), set())["global"]["precision"] is None


def test_boundary_uses_own_event_and_true_ancestry():
    meta = {"events": [{"sequence_id":"protein1", "start":20, "end":31}, {"sequence_id":"protein7", "start":70, "end":81}],
            "ancestry": {"protein1":[1,9,10,19,32,41,42,70], "protein7":[20,59,60,69,82,91,92]}, "region_window":10}
    boundary = s.boundary_residues(meta)
    assert boundary == {("protein1",3),("protein1",4),("protein1",5),("protein1",6),
                        ("protein7",3),("protein7",4),("protein7",5),("protein7",6)}
    assert not s.local_subset({(("background",1),("background2",1))}, boundary)
    assert s.local_subset({(("background",1),("protein1",3))}, boundary)


def test_tree_and_ctmc_expectations(config):
    tree = sim.parse_tree(config["simulator"]["topology_newick_unit_scale"])
    sim.validate_tree(tree)
    assert len(sim.tips(tree)) == 12
    for p in [.25,.40]:
        scale = sim.scale_for(p)
        assert sim.probability(scale) == pytest.approx(p)
        # Semigroup composition for the symmetric 20-state transition kernel.
        x, y = sim.probability(.25*scale), sim.probability(.75*scale)
        assert x+y-(20/19)*x*y == pytest.approx(p)
        realized = sim.realization(config, 1001, p)
        assert len(realized["edges"]) == 22
        assert all(e["endpoint_changes"] > 0 for e in realized["edges"])
        assert all(e["unit_length"] > 0 for e in realized["edges"])
        assert all(sum(a!=b for a,b in zip(realized["ancestor"],seq)) > 0 for seq in realized["substitution_only_tips"].values())


@pytest.mark.parametrize("p", [.25, .40])
@pytest.mark.parametrize("seed", [1001,1002,1003,1004])
def test_matched_simulation_geometry_ancestry(config, seed, p):
    datasets = [sim.simulate(config, seed, cond, p) for cond in ["none","partial","high"]]
    for key in ["ancestor","substitution_only_tips","edges"]:
        assert all(d[2][key] == datasets[0][2][key] for d in datasets)
    for (inp, truth, meta), overlap, union in zip(datasets, [0,6,11], [24,18,13]):
        assert meta["realized_overlap"] == overlap
        assert meta["union_gap_columns"] == union
        assert sum(any(row[c]=="-" for row in truth.values()) for c in range(240)) == union
        assert meta["total_surviving_residues"] == 2856
        assert 0 < meta["mean_pairwise_identity_shared_extant_ancestral_sites"] < 1
        assert all(len(seq)==(228 if name in ["protein1","protein7"] else 240) for name,seq in inp.items())
        assert all(e["branch_type"]=="terminal" and 1 < e["start"] < e["end"] < 240 and e["end"]-e["start"]+1==12 for e in meta["events"])
        assert len({e["id"] for e in meta["events"]}) == 2
        s.validate_alignment(truth, inp)
        # Independent coordinate dictionary reconstructs all true pairs.
        coordinates = {}
        for name, ancestors in meta["ancestry"].items():
            for i, coord in enumerate(ancestors, 1):
                coordinates.setdefault(coord, []).append((name,i))
        manual = {tuple(sorted(pair)) for residues in coordinates.values() for pair in itertools.combinations(residues,2)}
        assert manual == s.pairs(truth)
        assert s.accuracy(manual, manual, s.boundary_residues(meta))["local"]["reference_pairs"] > 0
    again = sim.simulate(config,seed,"none",p)
    assert digest(sim.fasta(again[0]).encode()) == digest(sim.fasta(datasets[0][0]).encode())


@pytest.mark.parametrize("change", ["seed", "count", "length", "flags"])
def test_config_isolates_namespaces(config, change):
    changed = copy.deepcopy(config)
    if change=="seed": changed["simulator"]["replicate_seeds"][0] += 100
    if change=="count": changed["simulator"]["num_sequences"] = 13
    if change=="length": changed["simulator"]["deletion_length_per_event"] = 13
    if change=="flags": changed["mafft"]["retree"] = 1
    assert identity(config,"test")[0] != identity(changed,"test")[0]
    assert identity(config,"test")[0] == identity(copy.deepcopy(config),"test")[0]


@pytest.mark.parametrize("text", [">a\nAA\n>a\nAA\n", ">a\nAA\n", ">a\nAA\n>x\nAA\n", ">a\nAC\n>b\nAA\n", ">a\nAA-\n>b\nAA\n"])
def test_invalid_mafft_output_rejected(text):
    with pytest.raises(ValueError):
        s.validate_alignment(s.parse_fasta(text), {"a":"AA","b":"AA"})


def test_correct_loss_bookkeeping_random_shifts():
    rng = random.Random(77)
    t = s.pairs({"a":"AAAA", "b":"AAAA", "c":"AAAA"})
    for _ in range(40):
        rows = [{name: "-"*(n:=rng.randrange(3))+"AAAA"+"-"*(2-n) for name in ["a","b","c"]} for _ in range(2)]
        b,r = map(s.pairs,rows)
        c = s.compare(b,r,t,{("a",2)})
        assert c["global"]["recovered"]-c["global"]["lost"] == len(r&t)-len(b&t)
