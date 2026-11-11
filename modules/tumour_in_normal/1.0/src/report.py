#!/usr/bin/env python
# pairs.tsv: one row per pair, tool and variant (deTiN arm, TINC copy-number mode), with a call.
# normals.tsv: one row per normal and tool, from the tool's primary variant (deTiN joint if run, else ssnv; TINC's
# first copy-number mode). A normal serving several tumours takes the most severe call among its pairs, in the
# order contaminated, low_level, clean, not_assessable; the estimate, interval and site count are those of the pair
# that gave the call, the highest estimate on a tie.
import pandas as pd

opts = snakemake.params.opts
SEVERITY = ["contaminated", "low_level", "clean", "not_assessable"]
COLUMNS = ["normal_id", "tumour_id", "seq_type", "genome_build", "tool", "variant", "call", "estimate",
           "ci_low", "ci_high", "n_sites", "reason"]


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def detin_rows(path):
    rows = []
    for _, r in pd.read_csv(path, sep="\t", dtype=str).iterrows():
        tin, low, high = num(r.TiN), num(r.TiN_CI_low), num(r.TiN_CI_high)
        if "insuffcient data" in str(r.messages):
            call = "not_assessable"
        elif low > opts["detin_call"]["ci_low_above"]:
            call = "contaminated"
        elif tin > 0:
            call = "low_level"
        else:
            call = "clean"
        rows.append([r.normal_id, r.tumour_id, r.seq_type, r.genome_build, "detin", r.arm, call, tin, low, high,
                     num(r.ssnv_sites), r.messages])
    return rows


def tinc_rows(path):
    rows = []
    for _, r in pd.read_csv(path, sep="\t", dtype=str).iterrows():
        rule = opts["tinc_call"]
        if r.status != "ok" or num(r.get("TIT")) < rule["purity_at_least"]:
            call = "not_assessable"
            reason = r.reason if r.status != "ok" else "tumour purity below " + str(rule["purity_at_least"])
        else:
            tin = num(r.get(rule["estimate"]))
            call = "contaminated" if tin > rule["above"] else ("low_level" if tin > 0 else "clean")
            reason = "NA"
        rows.append([r.normal_id, r.tumour_id, r.seq_type, r.genome_build, "tinc", r.cna_mode, call,
                     num(r.get(rule["estimate"])), float("nan"), float("nan"), num(r.get("n_mutations_total")),
                     reason])
    return rows


rows = []
for path in snakemake.input.detin:
    rows += detin_rows(path)
for path in snakemake.input.tinc:
    rows += tinc_rows(path)
for pair in snakemake.params.not_run:
    rows.append([pair["normal_id"], pair["tumour_id"], pair["seq_type"], pair["genome_build"], pair["tool"],
                 "NA", "not_assessable", float("nan"), float("nan"), float("nan"), float("nan"), pair["reason"]])
pairs = pd.DataFrame(rows, columns=COLUMNS)
pairs.to_csv(snakemake.output.pairs, sep="\t", index=False)

primary = []
for (tool, normal, seq_type, build), g in pairs.groupby(["tool", "normal_id", "seq_type", "genome_build"]):
    if tool == "detin":
        order = ["joint", "ssnv", "NA"]
    else:
        order = list(opts["tinc_cna_modes"]) + ["NA"]
    g = g[g.variant == next(v for v in order if v in set(g.variant))]
    g = g.assign(rank=g.call.map(SEVERITY.index)).sort_values(["rank", "estimate"], ascending=[True, False])
    best = g.iloc[0]
    primary.append(dict(best.drop("rank"), n_pairs=len(g), pairs=",".join(g.tumour_id + "--" + g.normal_id)))
normals = pd.DataFrame(primary, columns=COLUMNS + ["n_pairs", "pairs"]).rename(columns={"tumour_id": "from_tumour"})
normals.to_csv(snakemake.output.normals, sep="\t", index=False)
