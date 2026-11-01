#!/usr/bin/env Rscript
# Rascal purity and ploidy fit of a capture-mode CBS profile (gatk4_cnv_cbs.R output on targets) scaled to its median MedianRatio.
# Usage: gatk4_cnv_fit_capture.R <ratio.tsv> <out purity.tsv> <out solutions.tsv> <run> <contigs, comma-separated> <segment weight>
#        <distance function> <penploidy> <penalty> <ambiguous distance> <ambiguous purity>
# As gatk4_cnv_fit.R, except segments are runs of equal MedianRatio within a contig over targets (Start..End), weighted by
# "span" (first target start to last target end, bp) or "targets" (target count).
suppressPackageStartupMessages({library(data.table); library(rascal)})
a = commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 11, a[6] %in% c("span", "targets"))
contigs = strsplit(a[5], ",", fixed = TRUE)[[1]]
penploidy = as.numeric(a[8]); penalty = as.numeric(a[9])

w = fread(a[1], colClasses = list(character = "Chromosome"))
stopifnot(nrow(w) > 0, !anyNA(w$MedianRatio), "End" %in% names(w))
m = median(w$MedianRatio)
sg = copy_number_segments(data.frame(sample = "s", chromosome = match(w$Chromosome, contigs), start = w$Start,
                                     end = w$End, segmented = w$MedianRatio / m))
wt = if (a[6] == "span") sg$end - sg$start + 1 else sg$bin_count
s = as.data.table(find_best_fit_solutions(sg$copy_number, wt, distance_function = a[7]))
s[, penalised_distance := distance * (1 + abs(ploidy - 2))^penploidy / cellularity^penalty]
if (nrow(s)) {
  i = order(s$penalised_distance, -s$cellularity)[1]
  b = order(s$distance, -s$cellularity)[1]
  out = data.table(run = a[4], purity = s$cellularity[i], ploidy = s$ploidy[i],
                   ambiguous_fit = any(s$distance <= (1 + as.numeric(a[10])) * s$distance[b] & abs(s$cellularity - s$cellularity[b]) > as.numeric(a[11])),
                   solutions = nrow(s))
} else {
  out = data.table(run = a[4], purity = NA_real_, ploidy = NA_real_, ambiguous_fit = NA, solutions = 0L)
}
fwrite(out, a[2], sep = "\t", na = "NA")
fwrite(s, a[3], sep = "\t", na = "NA")
cat("segments", nrow(sg), "weight", a[6], "solutions", nrow(s), "purity", out$purity, "ploidy", out$ploidy, "\n", file = stderr())
