#!/usr/bin/env bash
# One samtools bedcov pass of a BAM/CRAM over the kit panel.
# Usage: fingerprint.sh <bam> <panel.bed> <ref.fa> <out.tsv>
# Panel contigs are bare; they get a chr prefix when the BAM header uses one.
# Output: panel columns plus summed base depth (MAPQ >= 20, flags 0xF04 excluded).
set -euo pipefail

bam=$1 panel=$2 ref=$3 out=$4
tmp=$(mktemp "$out.panel.XXXX")
trap 'rm -f "$tmp"' EXIT

sn=$(samtools view -H --reference "$ref" "$bam" | awk '$1 == "@SQ" && !n++ {sub(/^SN:/, "", $2); print $2}')
if [[ $sn == chr* ]]; then sed 's/^/chr/' "$panel" > "$tmp"; else cp "$panel" "$tmp"; fi

samtools bedcov -Q 20 -G 0xF04 -j --reference "$ref" "$tmp" "$bam" > "$out.part"
mv "$out.part" "$out"
