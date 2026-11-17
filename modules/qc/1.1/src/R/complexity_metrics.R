#!/usr/bin/env Rscript
# One row of coverage evenness, GC bias and complexity metrics for a sample.
# Base level (from the mosdepth coverage distribution: region.dist over the targets for capture, global.dist over the
#   genome otherwise):
#   cv = sd / mean of per-base coverage; pct_ge_0.2x_mean = share of bases at >= 0.2 x the mean.
# Region level (per target for capture, per bin otherwise; regions with mean > 0 in the dedup run):
#   region_cv = sd / mean of region means; pct_regions_0.2_5x_median = share within 0.2-5 x the median region.
# Bin dispersion (non-capture): fragments per autosomal bin (mean depth x bin / read length, halved for paired reads),
#   robust variance / median within GC deciles of bins with < 1% N, averaged over deciles; 1 is counting noise alone.
# Complexity: preseq expected distinct at the seq_type's fixed read count (NA beyond the curve); total reads and
#   the read count the fit was observed at.
# Each base / region metric is given for the dedup and withdup coverage.

suppressPackageStartupMessages(library(data.table))

i <- snakemake@input
w <- snakemake@wildcards
profile <- fread(i$profile)

base_level <- function(path, capture) {
  d <- fread(path, header = FALSE, col.names = c("chrom", "cov", "frac"))[chrom == "total"]
  setorder(d, cov)
  p <- d$frac - c(d$frac[-1], 0)
  m <- sum(d$cov * p)
  s <- sqrt(sum(p * (d$cov - m)^2))
  list(mean = m, cv = if (m > 0) s / m else NA_real_,
       pct_ge_0.2x_mean = if (m > 0) d[cov >= ceiling(0.2 * m), max(c(frac, 0))] else NA_real_)
}

region_level <- function(path) {
  r <- fread(cmd = paste("gzip -dc", shQuote(path)), header = FALSE)
  x <- r[[ncol(r)]]
  x <- x[x > 0]
  if (length(x) < 10) return(list(region_cv = NA_real_, pct_regions = NA_real_))
  med <- median(x)
  list(region_cv = sd(x) / mean(x), pct_regions = mean(x >= 0.2 * med & x <= 5 * med))
}

dispersion <- function(path, gc_path) {
  r <- fread(cmd = paste("gzip -dc", shQuote(path)), header = FALSE, col.names = c("chrom", "start", "end", "mean"))
  g <- fread(gc_path)
  r <- merge(r, g, by = c("chrom", "start", "end"))
  r <- r[grepl("^(chr)?[0-9]+$", chrom) & n_share < 0.01]
  len <- profile$mean_read_length
  if (!nrow(r) || is.na(len)) return(NA_real_)
  frag <- r$mean * (r$end - r$start) / len / (if (isTRUE(profile$paired_share > 0.5)) 2 else 1)
  r[, frag := frag]
  r[, dec := cut(gc, unique(quantile(gc, 0:10 / 10)), include.lowest = TRUE)]
  v <- r[, .(d = if (median(frag) > 0) mad(frag)^2 / median(frag) else NA_real_), by = dec]$d
  mean(v, na.rm = TRUE)
}

capture <- w$seq_type == "capture"
out <- data.table(sample_id = w$sample_id, seq_type = w$seq_type, genome_build = w$genome_build,
                  reads_profiled = profile$reads, mean_read_length = profile$mean_read_length,
                  paired_share = profile$paired_share)
for (k in c("dedup", "withdup")) {
  b <- base_level(i[[paste0(k, "_dist")]], capture)
  r <- region_level(i[[paste0(k, "_regions")]])
  out[, paste0(c("mean_cov_", "cv_", "pct_ge_0.2x_mean_", "region_cv_", "pct_regions_0.2_5x_median_"), k) :=
        list(b$mean, b$cv, b$pct_ge_0.2x_mean, r$region_cv, r$pct_regions)]
  out[, paste0("bin_dispersion_", k) := if (capture) NA_real_ else dispersion(i[[paste0(k, "_regions")]], i$bin_gc)]
}

gc <- fread(cmd = paste("grep -v '^#'", shQuote(i$gc_summary)), fill = TRUE)
out[, `:=`(gc_at_dropout = gc$AT_DROPOUT[1], gc_gc_dropout = gc$GC_DROPOUT[1])]

pc <- fread(i$preseq)
n <- as.numeric(snakemake@params$fixed_reads)
out[, preseq_fixed_reads := n]
out[, preseq_expected_distinct := if (nrow(pc) && n >= min(pc$TOTAL_READS) && n <= max(pc$TOTAL_READS))
      approx(pc$TOTAL_READS, pc$EXPECTED_DISTINCT, n)$y else NA_real_]
out[, preseq_fit := nrow(pc) > 0]

fwrite(out, snakemake@output$tsv, sep = "\t", na = "NA")
