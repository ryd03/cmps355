"""Explicit symmetric 20-state endpoint CTMC on every fixed-tree edge."""
import itertools
import math
import random
import re
from .storage import canonical, digest

SIMULATOR_VERSION = "fixed-tree-ctmc-v2.1"


def stream(seed, name):
    # Full SHA256 big-endian integer; stable named streams, no built-in hash().
    return random.Random(int(digest(canonical([seed, name])), 16))


def parse_tree(newick):
    tokens = re.findall(r"[(),:;]|[^(),:;\s]+", newick)
    position = 0

    def node():
        nonlocal position
        children = []
        label = None
        if tokens[position] == "(":
            position += 1
            children.append(node())
            while tokens[position] == ",":
                position += 1
                children.append(node())
            assert tokens[position] == ")"
            position += 1
        else:
            label = tokens[position]
            position += 1
        length = 0.0
        if tokens[position] == ":":
            position += 1
            length = float(tokens[position])
            position += 1
        return {"label": label, "length": length, "children": children}

    root = node()
    assert tokens[position:] == [";"]
    return root


def tips(node, distance=0.0):
    distance += node["length"]
    if not node["children"]:
        return [(node["label"], distance)]
    return list(itertools.chain.from_iterable(tips(c, distance) for c in node["children"]))


def validate_tree(tree):
    leaves = tips(tree)
    assert len(leaves) == 12
    assert {s for s, _ in leaves} == {f"protein{i}" for i in range(1, 13)}
    assert all(abs(d - 1) < 1e-10 for _, d in leaves)
    assert len(tree["children"]) == 2
    groups = [{s for s, _ in tips(c)} for c in tree["children"]]
    assert any("protein1" in g and "protein7" not in g for g in groups)


def probability(t):
    return 19 / 20 * (-math.expm1(-20 * t / 19))


def scale_for(p):
    return -19 / 20 * math.log1p(-20 * p / 19)


def realization(config, seed, p):
    sim = config["simulator"]
    tree = parse_tree(sim["topology_newick_unit_scale"])
    validate_tree(tree)
    alphabet = sim["amino_acid_alphabet"]
    assert len(set(alphabet)) == 20
    root_rng = stream(seed, "root")
    ancestor = "".join(root_rng.choice(alphabet) for _ in range(sim["ancestor_length"]))
    states, edges = {}, []
    scale = scale_for(p)

    def evolve(node, parent, path):
        rng = stream(seed, "edge:" + path)
        prob = probability(node["length"] * scale)
        child = "".join(rng.choice(alphabet.replace(aa, "")) if rng.random() < prob else aa
                        for aa in parent)
        edges.append({"path": path, "unit_length": node["length"], "p_endpoint_change": prob,
                      "endpoint_changes": sum(a != b for a, b in zip(parent, child)),
                      "state": child, "terminal_label": node["label"]})
        if node["children"]:
            for i, sub in enumerate(node["children"]):
                evolve(sub, child, path + "." + str(i))
        else:
            states[node["label"]] = child

    for i, sub in enumerate(tree["children"]):
        evolve(sub, ancestor, str(i))
    return {"ancestor": ancestor, "substitution_only_tips": states, "edges": edges,
            "root_to_tip_scale": scale, "expected_root_to_tip_difference": probability(scale),
            "rng_derivation": "Random(int(SHA256(canonical_JSON([seed, stream_name])),16)); root, edge:<child-index-path>, anchor",
            "seed": seed}


def simulate(config, seed, condition, p, shared=None):
    sim = config["simulator"]
    shared = realization(config, seed, p) if shared is None else shared
    assert shared["seed"] == seed
    assert abs(shared["expected_root_to_tip_difference"] - p) < 1e-12
    length = sim["ancestor_length"]
    start = sim["focal_1_deletion_anchor_1based"] + stream(seed, "anchor").randint(*sim["anchor_jitter_range"])
    starts = [start, start + sim["second_deletion_start_delta_by_condition"][condition]]
    events = [{"id": f"terminal_deletion_{name}", "sequence_id": name, "branch_type": "terminal",
               "start": s, "end": s + sim["deletion_length_per_event"] - 1}
              for name, s in zip(sim["focal_lineages"], starts)]
    deleted = {e["sequence_id"]: set(range(e["start"], e["end"] + 1)) for e in events}
    assert all(1 < e["start"] <= e["end"] < length for e in events)
    sets = list(deleted.values())
    overlap, union = len(sets[0] & sets[1]), len(sets[0] | sets[1])
    assert overlap == sim["expected_shared_deleted_positions"][condition]
    assert union == sim["expected_union_gap_columns"][condition]
    truth, sequences, ancestry = {}, {}, {}
    for name, seq in shared["substitution_only_tips"].items():
        truth[name] = "".join("-" if i in deleted.get(name, ()) else aa for i, aa in enumerate(seq, 1))
        sequences[name] = truth[name].replace("-", "")
        ancestry[name] = [i for i in range(1, length + 1) if i not in deleted.get(name, ())]
    identities = []
    for a, b in itertools.combinations(truth, 2):
        comparable = [(x, y) for x, y in zip(truth[a], truth[b]) if x != "-" and y != "-"]
        identities.append(sum(x == y for x, y in comparable) / len(comparable))
    metadata = {**shared, "condition": condition, "events": events, "ancestry": ancestry,
                "realized_overlap": overlap, "union_gap_columns": union,
                "region_window": sim["region_window_ancestral_positions"],
                "mean_pairwise_identity_shared_extant_ancestral_sites": sum(identities) / len(identities),
                "pairwise_identity_denominator": "each sequence pair's shared extant ancestral sites; mean of 66 identities",
                "root_tip_endpoint_differences": {s: sum(a != b for a, b in zip(shared["ancestor"], seq))
                                                   for s, seq in shared["substitution_only_tips"].items()},
                "total_surviving_residues": sum(map(len, sequences.values()))}
    validate_dataset(sequences, truth, metadata)
    return sequences, truth, metadata


def validate_dataset(sequences, truth, metadata):
    assert set(sequences) == set(truth) == set(metadata["ancestry"])
    assert all(truth[s].replace("-", "") == seq for s, seq in sequences.items())
    for name, seq in truth.items():
        assert metadata["ancestry"][name] == [i for i, aa in enumerate(seq, 1) if aa != "-"]
    assert len({len(seq) for seq in truth.values()}) == 1


def fasta(rows):
    return "".join(f">{name}\n{seq}\n" for name, seq in rows.items())
