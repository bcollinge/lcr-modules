#!/usr/bin/env Rscript
# Fixed grid over each contig of a target interval_list, with the targets per bin.
# Usage: gatk4_cnv_grid.R <targets.interval_list> <out grid.tsv> <contigs, comma-separated> <bin bp>
# Out, one row per bin of the listed contigs (lengths from the @SQ header): Chromosome, bin_start (0-based, from 0),
# bin_end (exclusive, capped at the contig length), n_targets (targets whose midpoint lies in the bin).
suppressPackageStartupMessages(library(data.table))
a = commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 4)
contigs = strsplit(a[3], ",", fixed = TRUE)[[1]]
bin = as.integer(a[4])

hdr = grep("^@SQ", readLines(a[1]), value = TRUE)
len = data.table(Chromosome = sub(".*\tSN:([^\t]+).*", "\\1", hdr), len = as.integer(sub(".*\tLN:([0-9]+).*", "\\1", hdr)))[Chromosome %in% contigs]
stopifnot(setequal(len$Chromosome, contigs))
t = fread(cmd = paste("grep -v '^@'", shQuote(a[1])), header = FALSE, select = 1:3, col.names = c("Chromosome", "start", "end"),
          colClasses = list(character = 1L))[Chromosome %in% contigs]
stopifnot(nrow(t) > 0)
t[, bin_start := ((start - 1 + end) %/% 2) %/% bin * bin]
g = len[, .(bin_start = seq(0L, len - 1L, by = bin)), by = .(Chromosome, len)][, bin_end := pmin(bin_start + bin, len)][, len := NULL]
g = merge(g, t[, .(n_targets = .N), by = .(Chromosome, bin_start)], by = c("Chromosome", "bin_start"), all.x = TRUE)
g[is.na(n_targets), n_targets := 0L]
g = g[order(match(Chromosome, contigs), bin_start), .(Chromosome, bin_start, bin_end, n_targets)]
fwrite(g, a[2], sep = "\t")
cat("bins", nrow(g), "with targets", g[n_targets > 0, .N], "targets", nrow(t), "\n", file = stderr())
