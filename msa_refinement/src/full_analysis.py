"""Paired seed-block analysis. No pilot data and no residue-pair pseudoreplication."""
import gzip
import math
import random
import statistics
from .storage import atomic, canonical, digest, json_write
from .reporting import csv_text


def quantile(values,q):
    if not values:
        raise ValueError("Cannot take a quantile of empty values")
    ordered=sorted(values)
    rank=(len(ordered)-1)*q
    lo=math.floor(rank); hi=math.ceil(rank)
    return ordered[lo]+(ordered[hi]-ordered[lo])*(rank-lo)


def resample_indices(n,replicates,seed):
    if n<=0 or replicates<1:
        raise ValueError("Positive block and bootstrap counts are required")
    rng=random.Random(seed)
    return [[rng.randrange(n) for _ in range(n)] for _ in range(replicates)]


def interval(values,indices,confidence=.95):
    if not values or any(v is None or not math.isfinite(v) for v in values):
        raise ValueError("Undefined/nonfinite block outcome; do not silently exclude blocks")
    draws=[statistics.fmean(values[i] for i in row) for row in indices]
    alpha=(1-confidence)/2
    return {"estimate":statistics.fmean(values),"ci_lower":quantile(draws,alpha),
            "ci_upper":quantile(draws,1-alpha),"n_seed_blocks":len(values)},draws


def paired_records(rows,config):
    seeds=config["simulator"]["replicate_seeds"]
    levels=config["divergence_levels"]
    if seeds != list(range(2001,2051)) or levels != [.25,.4]:
        raise ValueError("Full statistical estimates require exactly seeds 2001-2050 and both levels")
    expected={(p,seed,cond,cap) for p in levels for seed in seeds for cond in config["conditions"] for cap in [0,2,10]}
    indexed={}
    for r in rows:
        key=(r["divergence"],r["seed"],r["condition"],r["iteration_cap"])
        if key in indexed or key not in expected:
            raise ValueError("Duplicate or out-of-design observation (pilot seeds are excluded)")
        indexed[key]=r
    if set(indexed)!=expected:
        raise ValueError("Incomplete experiment; no confidence intervals may be estimated")
    records=[]
    for p in levels:
        for cap in [2,10]:
            for seed in seeds:
                none=indexed[p,seed,"none",cap]
                for condition in ["partial","high"]:
                    treatment=indexed[p,seed,condition,cap]
                    if any(r.get("local_net_gain") is None for r in [none,treatment]):
                        raise ValueError("Missing normalized net outcome")
                    row={"divergence":p,"seed":seed,"iteration_cap":cap,"contrast":condition+"_minus_none",
                         "analysis_role":"primary" if cap==2 else "secondary",
                         "none_net_local_gain":none["local_net_gain"],"treatment_net_local_gain":treatment["local_net_gain"],
                         "paired_difference":treatment["local_net_gain"]-none["local_net_gain"]}
                    for label,case in [("none",none),("treatment",treatment)]:
                        for metric in ["local_recovered","local_lost","local_reference_denominator","local_baseline_missed_denominator","local_baseline_correct_denominator"]:
                            row[label+"_"+metric]=case[metric]
                    records.append(row)
    return indexed,records


def analyze(rows,config,root=None):
    indexed,paired=paired_records(rows,config)
    seeds=config["simulator"]["replicate_seeds"]
    settings=config["analysis"]
    indices=resample_indices(len(seeds),settings["bootstrap_replicates"],settings["bootstrap_seed"])
    vectors=[]
    for cap in [2,10]:
        for p in config["divergence_levels"]:
            for condition in ["partial","high"]:
                values=[indexed[p,s,condition,cap]["local_net_gain"]-indexed[p,s,"none",cap]["local_net_gain"] for s in seeds]
                vectors.append(({"family":"geometry_contrast","divergence":p,"iteration_cap":cap,"contrast":condition+"_minus_none",
                                 "analysis_role":"primary" if cap==2 else "secondary"},values))
            for cond in config["conditions"]:
                vectors.append(({"family":"mean_net_gain","divergence":p,"iteration_cap":cap,"contrast":cond,
                                 "analysis_role":"secondary_descriptive"},[indexed[p,s,cond,cap]["local_net_gain"] for s in seeds]))
        for cond in config["conditions"]:
            values=[indexed[.4,s,cond,cap]["local_net_gain"]-indexed[.25,s,cond,cap]["local_net_gain"] for s in seeds]
            vectors.append(({"family":"divergence_difference","divergence":"0.40_minus_0.25","iteration_cap":cap,
                             "contrast":cond,"analysis_role":"secondary"},values))
        for cond in ["partial","high"]:
            values=[(indexed[.4,s,cond,cap]["local_net_gain"]-indexed[.4,s,"none",cap]["local_net_gain"])-
                    (indexed[.25,s,cond,cap]["local_net_gain"]-indexed[.25,s,"none",cap]["local_net_gain"]) for s in seeds]
            vectors.append(({"family":"divergence_geometry_interaction","divergence":"0.40_minus_0.25","iteration_cap":cap,
                             "contrast":cond+"_minus_none","analysis_role":"secondary"},values))
    results=[]; distributions=[]
    for i,(label,values) in enumerate(vectors):
        estimate,draws=interval(values,indices,settings["confidence_level"])
        key=f"estimate_{i:02d}"
        results.append({"estimate_id":key,**label,**estimate,"bootstrap_replicates":len(indices),
                        "bootstrap_seed":settings["bootstrap_seed"],"ci_method":settings["interval_method"],
                        "confidence_level":settings["confidence_level"],"coverage":"marginal_not_simultaneous"})
        distributions.append(draws)
    divergence_pairs=[]
    for s in seeds:
        for cap in [2,10]:
            for cond in config["conditions"]:
                lo,hi=indexed[.25,s,cond,cap],indexed[.4,s,cond,cap]
                divergence_pairs.append({"seed":s,"condition":cond,"iteration_cap":cap,
                    "net_gain_div025":lo["local_net_gain"],"net_gain_div040":hi["local_net_gain"],
                    "paired_divergence_difference":hi["local_net_gain"]-lo["local_net_gain"],
                    "recovered_div025":lo["local_recovered"],"lost_div025":lo["local_lost"],
                    "recovered_div040":hi["local_recovered"],"lost_div040":hi["local_lost"]})
    if root is not None:
        atomic(root/"paired_comparisons.csv",csv_text(paired))
        atomic(root/"divergence_paired_comparisons.csv",csv_text(divergence_pairs))
        atomic(root/"confidence_intervals.csv",csv_text(results))
        data=canonical({"seeds":seeds,"indices_zero_based":indices,"bootstrap_seed":settings["bootstrap_seed"]})
        atomic(root/"bootstrap_resampling.json.gz",gzip.compress(data,mtime=0))
        draw_rows=[{"resample":j+1,**{results[i]["estimate_id"]:distributions[i][j] for i in range(len(results))}} for j in range(len(indices))]
        atomic(root/"bootstrap_estimates.csv.gz",gzip.compress(csv_text(draw_rows).encode(),mtime=0))
        json_write(root/"analysis_done.json",{"seed_blocks":len(seeds),"pilot_blocks_used":0,"estimate_count":len(results),
                   "bootstrap_resampling_sha256":digest(data),"settings":settings})
    return paired,results
