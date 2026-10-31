#!/usr/bin/env python3
# maf_rows_to_vcf.py MAF FASTA OUT_VCF
# Writes the data rows of a MAF as VCF records in the fasta's contig naming, ID = 1-based data-row ordinal, for an
# exact-allele lookup. SNP/DNP/TNP/ONP keep Start_Position; insertions and deletions get the anchor base from the fasta.
# Rows on contigs absent from the fasta, or with alleles that cannot be written, are left out (counted on stderr).
import sys

maf, fasta, out = sys.argv[1:4]
fai = {}
order = []
for line in open(fasta + ".fai"):
    name, length, offset, linebases, linewidth = line.split("\t")[:5]
    fai[name] = (int(length), int(offset), int(linebases), int(linewidth))
    order.append(name)

fh = open(fasta, "rb")
def base(contig, pos):
    length, offset, linebases, linewidth = fai[contig]
    q = pos - 1
    fh.seek(offset + (q // linebases) * linewidth + q % linebases)
    return fh.read(1).decode().upper()

def contig_of(chrom):
    if chrom in fai:
        return chrom
    alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
    return alt if alt in fai else None

records, skipped, n = [], 0, 0
with open(maf) as h:
    cols = None
    for line in h:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if cols is None:
            cols = {c: i for i, c in enumerate(f)}
            continue
        n += 1
        contig = contig_of(f[cols["Chromosome"]])
        start, ref, alt, vt = int(f[cols["Start_Position"]]), f[cols["Reference_Allele"]], f[cols["Tumor_Seq_Allele2"]], f[cols["Variant_Type"]]
        if contig is None or not set(ref + alt) <= set("ACGTN-"):
            skipped += 1
            continue
        if vt in ("SNP", "DNP", "TNP", "ONP") and len(ref) == len(alt) and "-" not in ref + alt:
            records.append((contig, start, ref, alt, n))
        elif vt == "INS" and ref == "-" and "-" not in alt:
            b = base(contig, start)
            records.append((contig, start, b, b + alt, n))
        elif vt == "DEL" and alt == "-" and "-" not in ref and start > 1:
            b = base(contig, start - 1)
            records.append((contig, start - 1, b + ref, b, n))
        else:
            skipped += 1

rank = {c: i for i, c in enumerate(order)}
records.sort(key=lambda r: (rank[r[0]], r[1]))
with open(out, "w") as o:
    o.write("##fileformat=VCFv4.2\n")
    for c in order:
        o.write(f"##contig=<ID={c},length={fai[c][0]}>\n")
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for c, p, r, a, i in records:
        o.write(f"{c}\t{p}\t{i}\t{r}\t{a}\t.\t.\t.\n")
print(f"{n} rows, {len(records)} written, {skipped} left out", file=sys.stderr)
