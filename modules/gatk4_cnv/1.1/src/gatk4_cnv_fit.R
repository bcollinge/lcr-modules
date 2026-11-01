#!/usr/bin/env Rscript
# Rascal purity and ploidy fit of a CBS profile (gatk4_cnv_cbs.R output) scaled to its median MedianRatio.
# Usage: gatk4_cnv_fit.R <ratio.tsv> <out purity.tsv> <out solutions.tsv> <run> <contigs, comma-separated> <bin length>
#        <distance function> <penploidy> <penalty> <ambiguous distance> <ambiguous purity>
# Segments are runs of equal MedianRatio within a contig; each bin counts as bin length bp, so segment weights are bin counts.
# Solution of record: lowest distance x (1 + |ploidy - 2|)^penploidy / cellularity^penalty, ties to the higher cellularity.
# ambiguous_fit: another solution within (1 + ambiguous distance) x the lowest raw distance has a cellularity more than
# ambiguous purity from the lowest-distance solution's.
# Out: purity.tsv (run, purity, ploidy, ambiguous_fit, solutions; purity NA when Rascal returns no solution);
# solutions.tsv (every Rascal solution with its penalised distance).
suppressPackageStartupMessages({library(data.table); library(rascal)})
a = commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 11)
contigs = strsplit(a[5], ",", fixed = TRUE)[[1]]
bin_length = as.numeric(a[6])
penploidy = as.numeric(a[8]); penalty = as.numeric(a[9])

w = fread(a[1], colClasses = list(character = "Chromosome"))
stopifnot(nrow(w) > 0, !anyNA(w$MedianRatio))
m = median(w$MedianRatio)
sg = copy_number_segments(data.frame(sample = "s", chromosome = match(w$Chromosome, contigs), start = w$Start,
                                     end = w$Start + bin_length - 1, segmented = w$MedianRatio / m))
s = as.data.table(find_best_fit_solutions(sg$copy_number, sg$weight, distance_function = a[7]))
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
cat("segments", nrow(sg), "solutions", nrow(s), "purity", out$purity, "ploidy", out$ploidy, "\n", file = stderr())
