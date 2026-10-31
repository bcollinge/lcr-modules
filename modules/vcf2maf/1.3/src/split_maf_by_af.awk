# awk -v c=CUTOFF -v d=DROPPED -f split_maf_by_af.awk AF_TSV MAF > KEPT
# AF_TSV: MAF data-row ordinal and gnomAD AF ("," between several matching records; "." when absent). Rows whose
# largest AF is above c go to DROPPED, all other rows to stdout. Comment and header lines go to both.
BEGIN { FS = OFS = "\t" }
FILENAME == ARGV[1] { m = 0; k = split($2, v, ","); for (i = 1; i <= k; i++) if (v[i] + 0 > m) m = v[i] + 0; af[$1] = m; next }
/^#/ { print; print > d; next }
!h { print; print > d; h = 1; next }
{ n++; if ((n in af) && af[n] > c) print > d; else print }
