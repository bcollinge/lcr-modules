#!/usr/bin/env python
# Sorted unique contig and position of a call-stats table, read the way deTiN reads it: "#" lines skipped,
# and two leading rows skipped when the header is not found on the first line.
import pandas as pd

fields = ["contig", "position"]
try:
    table = pd.read_csv(snakemake.input.callstats, sep="\t", index_col=False, comment="#", usecols=fields,
                        dtype={"contig": str})
except (ValueError, LookupError):
    table = pd.read_csv(snakemake.input.callstats, sep="\t", index_col=False, comment="#", skiprows=2,
                        usecols=fields, dtype={"contig": str})
table = table.drop_duplicates().sort_values(["contig", "position"])
table.to_csv(snakemake.output.positions, sep="\t", header=False, index=False)
