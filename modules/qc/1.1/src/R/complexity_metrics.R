#!/usr/bin/env Rscript
# One row of complexity and GC bias metrics for a sample.
# Complexity: preseq expected distinct fragments (reads if unpaired; paired_share says which) at the seq_type's fixed
#   count, linear between 1M steps (NA beyond the curve), and the method behind it (c_curve: observed, interpolated;
#   lc_extrap: extrapolated past the fragments sequenced).
#   In paired mode preseq merges two mates into one fragment only when they are on the same chromosome and span at most
#   its seg_len (5,000 bp by default); otherwise each mate counts as its own unit, without a message. Libraries with
#   many inter-chromosomal or long-insert pairs therefore count some fragments twice.
# GC bias (gc_bias_seq_types only, NA otherwise): AT and GC dropout from CollectGcBiasMetrics over the whole genome.

suppressPackageStartupMessages(library(data.table))

i <- snakemake@input
w <- snakemake@wildcards
profile <- fread(i$profile)
out <- data.table(sample_id = w$sample_id, seq_type = w$seq_type, genome_build = w$genome_build,
                  reads_profiled = profile$reads, mean_read_length = profile$mean_read_length,
                  paired_share = profile$paired_share)

if (length(i$gc_summary)) {
  gc <- fread(cmd = paste("grep -v '^#'", shQuote(i$gc_summary)), fill = TRUE)
  out[, `:=`(gc_at_dropout = gc$AT_DROPOUT[1], gc_gc_dropout = gc$GC_DROPOUT[1])]
} else {
  out[, `:=`(gc_at_dropout = NA_real_, gc_gc_dropout = NA_real_)]
}

pc <- fread(i$preseq)
n <- as.numeric(snakemake@params$fixed_fragments)
hit <- nrow(pc) && n >= min(pc$TOTAL_FRAGMENTS) && n <= max(pc$TOTAL_FRAGMENTS)
out[, preseq_fixed_fragments := n]
out[, preseq_distinct_at_fixed := if (hit) approx(pc$TOTAL_FRAGMENTS, pc$EXPECTED_DISTINCT, n)$y else NA_real_]
out[, preseq_method := if (hit) pc$METHOD[which.min(abs(pc$TOTAL_FRAGMENTS - n))] else NA_character_]
out[, preseq_fragments_observed := if (nrow(pc)) max(pc[METHOD == "c_curve", TOTAL_FRAGMENTS]) else NA_real_]

fwrite(out, snakemake@output$tsv, sep = "\t", na = "NA")
