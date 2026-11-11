#!/usr/bin/env python
# Writes deTiN's population pickle in the layout of deTiN_utilities.build_exac_pickle: a dict keyed
# "<contig>_<position>" holding ref_allele, alt_allele, AC/AF, for records whose summed AF is at least min_af.
import pickle
import sys

sites_tsv, min_af, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
sites = {}
with open(sites_tsv) as handle:
    for line in handle:
        contig, pos, ref, alt, af = line.rstrip("\n").split("\t")
        if af in (".", ""):
            continue
        if sum(float(x) for x in af.split(",") if x != ".") >= min_af:
            sites[contig + "_" + pos] = {"ref_allele": ref, "alt_allele": alt, "AF": af}
with open(out, "wb") as handle:
    pickle.dump(sites, handle, protocol = pickle.HIGHEST_PROTOCOL)
