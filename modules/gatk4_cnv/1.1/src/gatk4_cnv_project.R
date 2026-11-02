#!/usr/bin/env Rscript
# Capture-mode profile (gatk4_cnv_correct.R output on targets) on the target grid (gatk4_cnv_grid.R output).
# Usage: gatk4_cnv_project.R <ratio.tsv> <grid.tsv> <out grid.tsv> <min targets per bin>
# Out, one row per grid bin: Chromosome, Start (bin_start), End (bin_end), n_targets (kit), n_called (targets with a ratio,
# by midpoint), Ratio (median of their Ratio), MedianRatio (of the segment holding most of them, ties to the first), CopyNumber
# (round half up of 2 x MedianRatio); Ratio, MedianRatio and CopyNumber NA when n_called < min targets per bin.
suppressPackageStartupMessages(library(data.table))
a = commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 4)
min_t = as.numeric(a[4])

w = fread(a[1], colClasses = list(character = "Chromosome"))
g = fread(a[2], colClasses = list(character = "Chromosome"))[, `:=`(bin_start = as.integer(bin_start), bin_end = as.integer(bin_end))]
stopifnot(nrow(w) > 0, all(c("Start", "End") %in% names(w)))
bin = g$bin_end[1] - g$bin_start[1]
w[, seg := rleid(Chromosome, MedianRatio)][, bin_start := ((Start - 1 + End) %/% 2) %/% bin * bin]
p = w[, .(n_called = .N, Ratio = median(Ratio), MedianRatio = MedianRatio[match(names(which.max(table(factor(seg, levels = unique(seg))))), seg)]),
      by = .(Chromosome, bin_start)]
x = merge(g, p, by = c("Chromosome", "bin_start"), all.x = TRUE, sort = FALSE)
x[is.na(n_called), n_called := 0L]
x[n_called < min_t, `:=`(Ratio = NA_real_, MedianRatio = NA_real_)]
x[, CopyNumber := floor(2 * MedianRatio + 0.5)]
x = x[order(match(Chromosome, unique(g$Chromosome)), bin_start)]
fwrite(x[, .(Chromosome, Start = bin_start, End = bin_end, n_targets, n_called, Ratio, MedianRatio, CopyNumber)], a[3], sep = "\t", na = "NA")
cat("bins", nrow(x), "called", x[!is.na(MedianRatio), .N], "targets", nrow(w), "\n", file = stderr())
