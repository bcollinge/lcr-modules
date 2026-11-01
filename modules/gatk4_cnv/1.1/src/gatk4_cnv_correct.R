#!/usr/bin/env Rscript
# Purity correction of a CBS profile (gatk4_cnv_cbs.R output) with the Rascal purity (gatk4_cnv_fit.R output), the lcr seg
# and the profile plot.
# Usage: gatk4_cnv_correct.R <ratio.tsv> <purity.tsv> <out ratio.tsv> <out seg> <out png> <seg ID> <contigs, comma-separated>
#        <bin length> <floor>
# Ratio and MedianRatio become max((r / m - (1 - p)) / p, 0), m the median MedianRatio, p the purity (6 decimals); CopyNumber is
# the round half up of 2 x the corrected MedianRatio. Without a purity the profile is only divided by m (not purity corrected).
# Seg: one row per run of equal corrected MedianRatio within a contig: ID, chrom, start, end, LOH_flag (NA), log.ratio
# (log2 of the corrected MedianRatio floored at floor).
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
a = commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 9)
contigs = strsplit(a[7], ",", fixed = TRUE)[[1]]
bin_length = as.numeric(a[8]); floor_cr = as.numeric(a[9])

w = fread(a[1], colClasses = list(character = "Chromosome"))
p = fread(a[2])$purity
stopifnot(nrow(w) > 0, length(p) == 1)
m = median(w$MedianRatio)
corr = function(r) if (is.na(p)) round(r / m, 6) else round(pmax((r / m - (1 - p)) / p, 0), 6)
x = copy(w)[, `:=`(Ratio = corr(Ratio), MedianRatio = corr(MedianRatio))]
x[, CopyNumber := floor(2 * MedianRatio + 0.5)]
fwrite(x, a[3], sep = "\t")

if (!"End" %in% names(x)) x[, End := Start + bin_length - 1]
seg = x[, .(chrom = Chromosome[1], start = Start[1], end = End[.N], value = MedianRatio[1]), by = .(k = rleid(Chromosome, MedianRatio))]
fwrite(seg[, .(ID = a[6], chrom, start, end, LOH_flag = NA, log.ratio = round(log2(pmax(value, floor_cr)), 6))], a[4], sep = "\t", na = "NA")

lab = if (is.na(p)) "copy number, not purity corrected" else sprintf("copy number, purity corrected (purity %.2f)", p)
d = rbind(w[, .(chr = factor(Chromosome, levels = contigs), Start, y = log2(pmax(Ratio, 2^-3)), seg = log2(pmax(MedianRatio, 2^-3)), panel = "log ratio")],
          x[, .(chr = factor(Chromosome, levels = contigs), Start, y = pmin(2 * Ratio, 6), seg = pmin(2 * MedianRatio, 6), panel = lab)])
d[, panel := factor(panel, levels = c("log ratio", lab))]
g = ggplot(d, aes(Start / 1e6)) +
  geom_point(aes(y = y), size = 0.15, colour = "grey60") +
  geom_step(aes(y = seg), colour = "#b2182b", linewidth = 0.5) +
  facet_grid(panel ~ chr, scales = "free", space = "free_x", switch = "x") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = unit(0.5, "mm"),
        strip.background = element_rect(fill = "white", colour = NA), panel.grid.minor = element_blank())
ggsave(a[5], g, width = 12, height = 4.5, dpi = 120, bg = "white")
cat("bins", nrow(x), "segments", nrow(seg), "purity", p, "\n", file = stderr())
