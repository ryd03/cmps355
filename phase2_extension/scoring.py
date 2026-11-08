"""Pairs identify residues by sequence ID and ungapped 1-based index, never letter."""
from itertools import combinations

SCORING_VERSION = "ancestral-focal-boundary-v2.1"


def parse_fasta(text):
    rows, name = {}, None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in rows:
                raise ValueError("Duplicate FASTA ID: " + name)
            rows[name] = ""
        else:
            if name is None:
                raise ValueError("Sequence before FASTA header")
            rows[name] += line.upper()
    if not rows or any(not seq for seq in rows.values()):
        raise ValueError("Empty FASTA")
    if any(set(seq) - set("ACDEFGHIKLMNPQRSTVWYXBZ-") for seq in rows.values()):
        raise ValueError("Invalid protein alphabet")
    return rows


def validate_alignment(rows, inputs):
    if set(rows) != set(inputs):
        raise ValueError("Missing or unknown sequence IDs")
    if len({len(s) for s in rows.values()}) != 1:
        raise ValueError("Unequal aligned row lengths")
    if any(rows[s].replace("-", "") != inputs[s] for s in rows):
        raise ValueError("Changed ungapped sequence")


def pairs(rows):
    if len({len(s) for s in rows.values()}) != 1:
        raise ValueError("Unequal aligned lengths")
    ordered = sorted(rows)
    indexes = dict.fromkeys(ordered, 0)
    result = set()
    for column in zip(*(rows[s] for s in ordered)):
        residues = []
        for name, aa in zip(ordered, column):
            if aa != "-":
                indexes[name] += 1
                residues.append((name, indexes[name]))
        result.update(combinations(residues, 2))
    return result


def boundary_residues(metadata):
    result = set()
    for event in metadata["events"]:
        name = event["sequence_id"]
        for idx, ancestor in enumerate(metadata["ancestry"][name], 1):
            if min(abs(ancestor - event["start"]), abs(ancestor - event["end"])) <= metadata["region_window"]:
                result.add((name, idx))
    return result


def local_subset(pair_set, boundary):
    return {p for p in pair_set if p[0] in boundary or p[1] in boundary}


def ratio(a, b):
    return a / b if b else None


def accuracy(predicted, truth, boundary):
    result = {}
    for label, select in (("global", lambda s: s),
                          ("local", lambda s: local_subset(s, boundary)),
                          ("far", lambda s: s - local_subset(s, boundary))):
        pred, ref = select(predicted), select(truth)
        tp = len(pred & ref)
        result[label] = {"reference_pairs": len(ref), "predicted_pairs": len(pred), "correct_pairs": tp,
                         "false_pairs": len(pred - ref), "missed_pairs": len(ref - pred),
                         "sp_recall": ratio(tp, len(ref)), "precision": ratio(tp, len(pred))}
    return result


def compare(baseline, refined, truth, boundary):
    result = {}
    for label, select in (("global", lambda s: s),
                          ("local", lambda s: local_subset(s, boundary)),
                          ("far", lambda s: s - local_subset(s, boundary))):
        b, r, t = select(baseline), select(refined), select(truth)
        recovered, lost = (r - b) & t, (b - r) & t
        assert not recovered & lost
        assert len(r & t) - len(b & t) == len(recovered) - len(lost)
        result[label] = {"recovered": len(recovered), "lost": len(lost),
                         "false_introduced": len((r - b) - t), "false_removed": len((b - r) - t),
                         "reference_denominator": len(t), "baseline_missed_denominator": len(t - b),
                         "baseline_correct_denominator": len(b & t),
                         "net_gain": ratio(len(recovered) - len(lost), len(t)),
                         "recovery_rate": ratio(len(recovered), len(t - b)),
                         "loss_rate": ratio(len(lost), len(b & t))}
    result["alignment_pair_set_changed"] = baseline != refined
    return result
