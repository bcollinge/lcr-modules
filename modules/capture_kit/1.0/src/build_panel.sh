#!/usr/bin/env bash
# Build the kit-discriminating panel for one build.
# Usage: build_panel.sh <kits.tsv> <build> <annotation_dir> <outdir>
#   kits.tsv: kit, build, path (.bb or BED); build: grch37 | hg38
#   annotation_dir: UCSC genomicSuperDups.<ucsc>.txt.gz, gap.<ucsc>.txt.gz,
#     <ucsc>.chrom.sizes and, for hg38, centromeres.hg38.txt.gz
# Needs bedtools and bigBedToBed on PATH.
# Output: <outdir>/panel.bed (BED4, set label), panel.sizes.tsv, kits/<kit>.bed
# Sets: A--not--B = A minus (B slopped by GAP), segdups removed, regions >= MIN_REGION;
#   core = intersection of all kits; floor = 10 kb windows clear of every kit by 5 kb,
#   of gaps, centromeres and segdups. Each set thinned evenly to PANEL_SIZE regions.
set -euo pipefail

kits_tsv=$1 build=$2 ann=$3 out=$4
GAP=${GAP:-1000}
MIN_REGION=${MIN_REGION:-100}
PANEL_SIZE=${PANEL_SIZE:-2000}
FLOOR_SLOP=${FLOOR_SLOP:-5000}
FLOOR_WIN=${FLOOR_WIN:-10000}
export LC_ALL=C

case $build in
  grch37) ucsc=hg19 ;;
  hg38) ucsc=hg38 ;;
  *) echo "unknown build $build" >&2; exit 1 ;;
esac

mkdir -p "$out/kits" "$out/sets"
tmp=$(mktemp -d "$out/tmp.XXXX")
trap 'rm -rf "$tmp"' EXIT

# bare autosome names, sorted
clean() { grep -vE '^(track|browser|#)' | cut -f1-3 | sed 's/^chr//' | awk -F'\t' '$1 ~ /^[0-9]+$/' | sort -k1,1 -k2,2n; }
ucsc_table() { gzip -dc "$1" | cut -f2-4 | clean; }

genome=$tmp/genome.txt
sed 's/^chr//' "$ann/$ucsc.chrom.sizes" | awk -F'\t' '$1 ~ /^[0-9]+$/' | sort -k1,1 > "$genome"

ucsc_table "$ann/genomicSuperDups.$ucsc.txt.gz" | bedtools merge > "$tmp/segdups.bed"
{
  ucsc_table "$ann/gap.$ucsc.txt.gz"
  [[ -f $ann/centromeres.$ucsc.txt.gz ]] && ucsc_table "$ann/centromeres.$ucsc.txt.gz"
  cat "$tmp/segdups.bed"
} | sort -k1,1 -k2,2n | bedtools merge > "$tmp/mask.bed"

kits=()
while IFS=$'\t' read -r kit b path; do
  [[ $kit == kit || $b != "$build" ]] && continue
  if [[ $path == *.bb ]]; then bigBedToBed "$path" /dev/stdout; else cat "$path"; fi \
    | clean | bedtools merge > "$out/kits/$kit.bed"
  kits+=("$kit")
done < "$kits_tsv"
(( ${#kits[@]} >= 2 )) || { echo "fewer than two kits for $build" >&2; exit 1; }

# keep regions >= MIN_REGION after removing segdups
finish() { bedtools subtract -a - -b "$tmp/segdups.bed" | awk -v m="$MIN_REGION" '$3 - $2 >= m'; }

for a in "${kits[@]}"; do
  for b in "${kits[@]}"; do
    [[ $a == "$b" ]] && continue
    bedtools slop -i "$out/kits/$b.bed" -g "$genome" -b "$GAP" \
      | bedtools subtract -a "$out/kits/$a.bed" -b - | finish > "$out/sets/$a--not--$b.bed"
  done
done

cp "$out/kits/${kits[0]}.bed" "$tmp/core.bed"
for k in "${kits[@]:1}"; do
  bedtools intersect -a "$tmp/core.bed" -b "$out/kits/$k.bed" > "$tmp/core.next"
  mv "$tmp/core.next" "$tmp/core.bed"
done
finish < "$tmp/core.bed" > "$out/sets/core.bed"

cat "$out"/kits/*.bed | sort -k1,1 -k2,2n \
  | bedtools slop -i - -g "$genome" -b "$FLOOR_SLOP" \
  | cat - "$tmp/mask.bed" | sort -k1,1 -k2,2n | bedtools merge \
  | bedtools complement -i - -g "$genome" \
  | bedtools makewindows -b - -w "$FLOOR_WIN" \
  | awk -v w="$FLOOR_WIN" '$3 - $2 == w' > "$out/sets/floor.bed"

printf 'set\tn_all\tbp_all\tn\tbp\n' > "$out/panel.sizes.tsv"
: > "$tmp/panel.bed"
for f in "$out"/sets/*.bed; do
  s=$(basename "$f" .bed)
  n=$(wc -l < "$f" | tr -d ' ')
  step=$(( n > PANEL_SIZE ? (n + PANEL_SIZE - 1) / PANEL_SIZE : 1 ))
  awk -v s="$s" -v k="$step" 'BEGIN{OFS="\t"} (NR - 1) % k == 0 {print $1, $2, $3, s}' "$f" > "$tmp/thin.bed"
  awk -v s="$s" -v n="$n" -v f="$f" 'BEGIN{OFS="\t"} {c++; bp += $3 - $2}
    END {while ((getline l < f) > 0) {split(l, x, "\t"); ba += x[3] - x[2]}; print s, n, ba + 0, c + 0, bp + 0}' \
    "$tmp/thin.bed" >> "$out/panel.sizes.tsv"
  cat "$tmp/thin.bed" >> "$tmp/panel.bed"
done
sort -k1,1 -k2,2n "$tmp/panel.bed" > "$out/panel.bed"
