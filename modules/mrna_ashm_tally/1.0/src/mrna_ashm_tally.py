#!/usr/bin/env python
"""Per-base counts over a set of regions in one RNA alignment.

Every position of every region is reported, whether or not a read covers it, so
that zero depth is a value rather than a missing row. Positions are emitted with
the reference base and the count of each observed base; calling thresholds are
applied downstream, not here.

  mrna_ashm_tally.py --bam X.cram --index X.cram.crai --reference genome.fa \
      --regions regions.tsv --sample SAMPLE --out SAMPLE.tally.tsv

regions: a headed TSV with columns region, chrom, start, end (1-based inclusive).
Contigs are matched with or without a "chr" prefix, whichever the alignment uses.
"""
import argparse
import csv
import os
import sys

import pysam

BASES = ("A", "C", "G", "T")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    # given explicitly, because the index may not sit beside the alignment or
    # carry the extension its format implies
    ap.add_argument("--index", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-mapq", type=int, default=1)
    ap.add_argument("--min-baseq", type=int, default=10)
    ap.add_argument("--max-depth", type=int, default=100000)
    args = ap.parse_args()

    with open(args.regions) as fh:
        regions = list(csv.DictReader(fh, delimiter="\t"))
    if not regions:
        sys.exit(f"no regions in {args.regions}")

    af = pysam.AlignmentFile(args.bam, reference_filename=args.reference,
                             index_filename=args.index)
    fa = pysam.FastaFile(args.reference)
    prefixed = any(c.startswith("chr") for c in af.references)

    # written to a temp path and renamed on success, so an interrupted run
    # cannot leave a partial file that looks complete
    tmp = f"{args.out}.{os.getpid()}.partial"
    n_fail = 0

    with open(tmp, "w") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "region", "chrom", "pos", "ref", "depth", *BASES])

        for r in regions:
            bare = r["chrom"][3:] if r["chrom"].startswith("chr") else r["chrom"]
            chrom = ("chr" + bare) if prefixed else bare
            start, end = int(r["start"]), int(r["end"])
            refseq = fa.fetch(chrom, start - 1, end).upper()
            counts = {p: [0, 0, 0, 0] for p in range(start, end + 1)}

            # an unreadable container raises during iteration rather than at
            # pileup(); the region is then reported at zero depth instead of
            # aborting the sample
            try:
                for col in af.pileup(chrom, start - 1, end, truncate=True,
                                     stepper="samtools",
                                     min_mapping_quality=args.min_mapq,
                                     min_base_quality=args.min_baseq,
                                     ignore_overlaps=True,
                                     max_depth=args.max_depth,
                                     multiple_iterators=False):
                    c = counts[col.reference_pos + 1]
                    for b in col.get_query_sequences(add_indels=False):
                        b = b.upper()
                        if b in BASES:
                            c[BASES.index(b)] += 1
            except (ValueError, OSError):
                n_fail += 1
                counts = {p: [0, 0, 0, 0] for p in range(start, end + 1)}

            for p in range(start, end + 1):
                c = counts[p]
                w.writerow([args.sample, r["region"], bare, p,
                            refseq[p - start], sum(c), *c])

    af.close()
    fa.close()
    os.replace(tmp, args.out)

    msg = f"{args.sample}: {len(regions)} regions"
    if n_fail:
        msg += f" | {n_fail} unreadable"
    print(msg, file=sys.stderr)


if __name__ == "__main__":
    main()
