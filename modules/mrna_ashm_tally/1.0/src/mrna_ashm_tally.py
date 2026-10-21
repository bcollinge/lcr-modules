#!/usr/bin/env python
"""Per-base counts over a set of regions in one RNA alignment.

Every position of every region is reported, whether or not a read covers it, so
that zero depth is a value rather than a missing row. Positions are emitted with
the reference base and the count of each observed base; calling thresholds are
applied downstream, not here.

  mrna_ashm_tally.py --bam X.cram --index X.cram.crai --reference genome.fa \
      --regions regions.tsv --sample SAMPLE --out SAMPLE.tally.tsv.gz

regions: a headed TSV with columns region, chrom, start, end (1-based inclusive),
sorted and non-overlapping. The output is bgzip-compressed and tabix-indexed
(OUT.tbi) on chrom and pos.
Contigs are matched with or without a "chr" prefix, whichever the alignment uses.

Consecutive regions on a chromosome that share a name (a comma-separated list)
are read in one pileup spanning them, so a spliced read covering several exons
of a gene is processed once rather than once per exon. Only region positions
are reported.
"""
import argparse
import csv
import os
import sys

import pysam

BASES = ("A", "C", "G", "T")


def spans(regions):
    """Group consecutive regions on one chromosome that share a name."""
    groups = []
    for r in regions:
        names = set(r["region"].split(","))
        if groups and groups[-1][-1]["chrom"] == r["chrom"] and names & groups[-1][-1]["_names"]:
            groups[-1].append(r)
        else:
            groups.append([r])
        r["_names"] = names
    return groups


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

        for group in spans(regions):
            bare = group[0]["chrom"][3:] if group[0]["chrom"].startswith("chr") else group[0]["chrom"]
            chrom = ("chr" + bare) if prefixed else bare
            first, last = int(group[0]["start"]), int(group[-1]["end"])
            counts = {p: [0, 0, 0, 0] for r in group
                      for p in range(int(r["start"]), int(r["end"]) + 1)}

            # an unreadable container raises during iteration rather than at
            # pileup(); the span is then reported at zero depth instead of
            # aborting the sample
            try:
                for col in af.pileup(chrom, first - 1, last, truncate=True,
                                     stepper="samtools",
                                     min_mapping_quality=args.min_mapq,
                                     min_base_quality=args.min_baseq,
                                     ignore_overlaps=True,
                                     max_depth=args.max_depth,
                                     multiple_iterators=False):
                    c = counts.get(col.reference_pos + 1)
                    if c is None:
                        continue
                    for b in col.get_query_sequences(add_indels=False):
                        b = b.upper()
                        if b in BASES:
                            c[BASES.index(b)] += 1
            except (ValueError, OSError):
                n_fail += len(group)
                counts = {p: [0, 0, 0, 0] for p in counts}

            for r in group:
                start, end = int(r["start"]), int(r["end"])
                refseq = fa.fetch(chrom, start - 1, end).upper()
                for p in range(start, end + 1):
                    c = counts[p]
                    w.writerow([args.sample, r["region"], bare, p,
                                refseq[p - start], sum(c), *c])

    af.close()
    fa.close()
    pysam.tabix_index(tmp, force=True, seq_col=2, start_col=3, end_col=3,
                      line_skip=1)
    os.replace(tmp + ".gz.tbi", args.out + ".tbi")
    os.replace(tmp + ".gz", args.out)

    msg = f"{args.sample}: {len(regions)} regions"
    if n_fail:
        msg += f" | {n_fail} unreadable"
    print(msg, file=sys.stderr)


if __name__ == "__main__":
    main()
