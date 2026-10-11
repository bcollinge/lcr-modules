#!/usr/bin/env bash
# gnomad_lookup.sh IN FASTA RELEASE TAG ATOMIZE TMPDIR
# Writes IN to stdout with INFO/TAG (Number=A): per ALT, the release INFO/AF matched by exact allele after splitting and
# left-aligning a copy of the record (ATOMIZE=1: MNVs also split into SNVs, max AF over bases); 0 when absent.
# The original records are written unchanged apart from INFO/TAG.
set -euo pipefail
IN=$1; FASTA=$2; RES=$3; TAG=$4; ATOMIZE=$5; T=$6
mkdir -p "$T"; trap 'rm -rf "$T"' EXIT
A=""; [ "$ATOMIZE" = 1 ] && A="-a"
# Copy: ID = ordinal of the original record, GLK_ALTIDX = its ALT index, carried through the split
bcftools view -G "$IN" |
awk 'BEGIN{FS=OFS="\t"} /^##/{print; next}
  /^#CHROM/{print "##INFO=<ID=GLK_ALTIDX,Number=A,Type=Integer,Description=\"ALT index in the original record\">"; print; next}
  {n++; k=split($5,a,","); s="1"; for(i=2;i<=k;i++) s=s","i; $3=n; $8="GLK_ALTIDX="s; print}' |
bcftools norm -m -any $A -f "$FASTA" -c w 2> "$T/norm.log" |
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/GLK_ALTIDX\n' > "$T/copy.tsv"
cut -f1,2 "$T/copy.tsv" | sort -u > "$T/sites.tsv"
if [ -s "$T/sites.tsv" ]; then
    bcftools query -R "$T/sites.tsv" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$RES" > "$T/res.tsv"
else
    : > "$T/res.tsv"
fi
# ordinal, ALT index, AF (exact allele; one release record may serve several copy records)
awk 'BEGIN{FS=OFS="\t"} NR==FNR {k=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4; v=($5=="."?0:$5+0); if (!(k in af) || v>af[k]) af[k]=v; next}
  $6!="." {k=$1 SUBSEP $2 SUBSEP $4 SUBSEP $5; print $3, $6, ((k in af) ? af[k] : 0)}' "$T/res.tsv" "$T/copy.tsv" |
sort -t "$(printf '\t')" -k1,1n -k2,2n > "$T/af.tsv"
bcftools view "$IN" | awk -v tag="$TAG" -v f="$T/af.tsv" '
  BEGIN{FS=OFS="\t"; more=(getline line < f)>0}
  index($0,"##INFO=<ID=" tag ",")==1 {next}
  /^##/{print; next}
  /^#CHROM/{print "##INFO=<ID=" tag ",Number=A,Type=Float,Description=\"gnomAD AF, exact allele after split and left-align of the call; 0 when absent\">"; print; next}
  {n++; delete af
   while (more) {split(line,x,"\t"); if (x[1]+0!=n) break; v=x[3]+0; if (!(x[2] in af) || v>af[x[2]]) af[x[2]]=v; more=(getline line < f)>0}
   k=split($5,a,","); s=""; for(i=1;i<=k;i++) s=s (i>1?",":"") ((i"") in af ? af[i""] : 0)
   m=split($8,kv,";"); info=""; for(j=1;j<=m;j++) if (kv[j]!="." && kv[j] !~ ("^" tag "(=|$)")) info=info (info==""?"":";") kv[j]
   $8=(info==""?"":info";") tag "=" s; print}'
