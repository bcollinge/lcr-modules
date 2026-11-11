#!/usr/bin/env python
# One row per deTiN run of a pair: the written estimate, interval and recovered-SNV count, beside each model's own
# estimate and the messages deTiN printed, both taken from its stdout.
import re

patterns = {
    "ssnv_TiN": r"SSNV based TiN estimate converged: TiN = (\S+) based on",
    "ssnv_sites": r"SSNV based TiN estimate converged: TiN = \S+ based on (\S+) sites",
    "ascna_TiN": r"aSCNA based TiN estimate(?: from selected TiN cluster :| TiN =) +(\S+)",
    "ascna_clusters": r"detected (\d+) clusters",
    "joint_TiN": r"joint TiN estimate = (\S+)",
}
messages = [
    "reselected cluster based on SSNVs",
    "SSNV based TiN estimate exceed 0.3 using only aSCNA based estimate",
    "No aSCNAs only using SSNV based model",
    "insuffcient data to generate TiN estimate.",
    "insufficient evidence to justify TiN > 0",
    "Only using aSCNA data",
]

columns = ["tumour_id", "normal_id", "seq_type", "genome_build", "arm", "TiN", "TiN_CI_low", "TiN_CI_high",
           "n_SSNVs_added"] + list(patterns) + ["messages"]
w = snakemake.wildcards
with open(snakemake.output.tsv, "w") as out:
    out.write("\t".join(columns) + "\n")
    for arm, tin, ci, added, stdout in zip(snakemake.params.arms, snakemake.input.tin, snakemake.input.ci,
                                           snakemake.input.added, snakemake.input.stdout):
        low, high = [x.strip() for x in open(ci).read().split(" - ")]
        log = open(stdout).read()
        found = []
        for key, pattern in patterns.items():
            hits = re.findall(pattern, log)
            found.append(hits[-1] if hits else "NA")
        fired = [m for m in messages if m in log]
        out.write("\t".join([w.tumour_id, w.normal_id, w.seq_type, w.genome_build, arm, open(tin).read().strip(),
                             low, high, open(added).read().strip()] + found + ["; ".join(fired) or "NA"]) + "\n")
