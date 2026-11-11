#!/usr/bin/env python
# deTiN call-stats table from a FilterMutectCalls VCF of a matched pair. Biallelic SNVs only, counts from AD.
# PASS is KEEP. A record whose filters all lie in candidate_filters is REJECT with failure_reasons
# "alt_allele_in_normal", which deTiN's candidate rule accepts; any other record is REJECT with its FILTER
# string (";" as ","), which the rule does not accept.
import gzip

candidate_filters = set(snakemake.params.candidate_filters)

columns = ["contig", "position", "ref_allele", "alt_allele", "tumor_name", "normal_name", "t_alt_count",
           "t_ref_count", "n_alt_count", "n_ref_count", "failure_reasons", "judgement"]

tumour = normal = None
with gzip.open(snakemake.input.vcf, "rt") as vcf, open(snakemake.output.callstats, "w") as out:
    out.write("\t".join(columns) + "\n")
    for line in vcf:
        if line.startswith("##tumor_sample="):
            tumour = line.rstrip("\n").split("=", 1)[1]
        elif line.startswith("##normal_sample="):
            normal = line.rstrip("\n").split("=", 1)[1]
        elif line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            assert tumour in samples and normal in samples, "VCF header lacks tumor_sample or normal_sample"
            t_col, n_col = samples.index(tumour), samples.index(normal)
        elif not line.startswith("#"):
            f = line.rstrip("\n").split("\t")
            ref, alt = f[3], f[4]
            if len(ref) != 1 or len(alt) != 1 or alt == "*":
                continue
            keys = f[8].split(":")
            t_ad = dict(zip(keys, f[9 + t_col].split(":")))["AD"].split(",")
            n_ad = dict(zip(keys, f[9 + n_col].split(":")))["AD"].split(",")
            filters = f[6].split(";")
            if f[6] == "PASS":
                reasons, judgement = "", "KEEP"
            elif set(filters) <= candidate_filters:
                reasons, judgement = "alt_allele_in_normal", "REJECT"
            else:
                reasons, judgement = ",".join(filters), "REJECT"
            out.write("\t".join([f[0], f[1], ref, alt, tumour, normal, t_ad[1], t_ad[0], n_ad[1], n_ad[0],
                                 reasons, judgement]) + "\n")
