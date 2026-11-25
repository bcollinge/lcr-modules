#!/usr/bin/env Rscript
# TINC on one pair: the pair's SNVs from a two-sample VCF (AD counts), Battenberg copy number in CNAqc columns,
# then TINC::autofit with its defaults. A pair TINC cannot fit gets a not_assessable row carrying the reason.
suppressPackageStartupMessages({
  library(TINC)
  library(dplyr)
  library(tidyr) # TINC 0.1.0 calls separate() without importing tidyr
})

log_con <- file(snakemake@log[["stderr"]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

p <- snakemake@params

# chr1-chr22, FILTER PASS, single-base substitutions, ref and alt counts from AD, tumour VAF strictly between 0 and
# 1 (as TINC's own VCF parser).
vcf_sites <- function(path) {
  lines <- readLines(path)
  hdr <- lines[startsWith(lines, "#")]
  meta <- function(key) sub(paste0("^##", key, "="), "", grep(paste0("^##", key, "="), hdr, value = TRUE))
  samples <- strsplit(tail(hdr, 1), "\t")[[1]][-(1:9)]
  t_col <- 9 + match(meta("tumor_sample"), samples)
  n_col <- 9 + match(meta("normal_sample"), samples)
  body <- lines[!startsWith(lines, "#")]
  if (length(body) == 0) stop("no SNVs in the VCF")
  f <- do.call(rbind, strsplit(body, "\t"))
  chr <- ifelse(startsWith(f[, 1], "chr"), f[, 1], paste0("chr", f[, 1]))
  keep <- chr %in% paste0("chr", 1:22) & nchar(f[, 4]) == 1 & nchar(f[, 5]) == 1 & f[, 7] == "PASS"
  ad <- function(col) {
    i <- vapply(strsplit(f[keep, 9], ":"), function(x) match("AD", x), integer(1))
    v <- mapply(function(s, j) strsplit(s, ":")[[1]][j], f[keep, col], i)
    do.call(rbind, lapply(strsplit(v, ","), as.numeric))
  }
  n <- ad(n_col)
  t <- ad(t_col)
  pos <- as.numeric(f[keep, 2])
  tibble::tibble(chr = chr[keep], from = pos - 1, to = pos, ref = f[keep, 4], alt = f[keep, 5],
                 n_ref_count = n[, 1], n_alt_count = n[, 2], t_ref_count = t[, 1], t_alt_count = t[, 2]) %>%
    dplyr::filter(t_alt_count / (t_ref_count + t_alt_count) > 0, t_alt_count / (t_ref_count + t_alt_count) < 1)
}

# Battenberg subclones.txt to CNAqc segments. clonal: segments whose first state has fraction 1. majority: the
# state with the larger fraction. state1: the first state whatever its fraction. Rows with no state called
# (fraction 0, no second state) are dropped or kept as written (1:1) by `filler`.
battenberg_segments <- function(path, mode, filler) {
  bb <- read.delim(path, check.names = FALSE)
  frac1 <- as.numeric(bb$frac1_A)
  frac2 <- suppressWarnings(as.numeric(bb$frac2_A))
  frac2[is.na(frac2)] <- 0
  unfilled <- frac1 == 0 & frac2 == 0
  use2 <- mode == "majority" & frac2 > frac1
  seg <- tibble::tibble(
    chr = ifelse(startsWith(as.character(bb$chr), "chr"), as.character(bb$chr), paste0("chr", bb$chr)),
    from = as.numeric(bb$startpos),
    to = as.numeric(bb$endpos),
    Major = as.integer(ifelse(use2, bb$nMaj2_A, bb$nMaj1_A)),
    minor = as.integer(ifelse(use2, bb$nMin2_A, bb$nMin1_A)),
    keep = if (mode == "clonal") frac1 == 1 else !(unfilled & filler == "drop")
  )
  seg %>% dplyr::filter(keep, !is.na(Major), !is.na(minor)) %>% dplyr::select(-keep)
}

row <- data.frame(tumour_id = snakemake@wildcards[["tumour_id"]], normal_id = snakemake@wildcards[["normal_id"]],
                  seq_type = snakemake@wildcards[["seq_type"]], genome_build = snakemake@wildcards[["genome_build"]],
                  pair_status = snakemake@wildcards[["pair_status"]], cna_mode = p[["cna_mode"]],
                  seed = as.integer(p[["seed"]]), status = "ok", reason = NA_character_, n_input_snvs = NA_integer_)
set.seed(as.integer(p[["seed"]]))
fit <- tryCatch({
  mutations <- vcf_sites(snakemake@input[["vcf"]])
  row$n_input_snvs <- nrow(mutations)
  cna <- battenberg_segments(snakemake@input[["subclones"]], p[["cna_mode"]], p[["cna_filler"]])
  TINC::autofit(input = mutations, cna = cna, FAST = isTRUE(p[["fast"]]))
}, error = function(e) {
  row$status <<- "not_assessable"
  row$reason <<- gsub("[\t\n]+", " ", conditionMessage(e))
  NULL
})

if (!is.null(fit)) {
  result <- TINC::to_string(fit)
  result$TIN_class <- TINC:::classification_normal(fit$TIN)[["QC"]]
  result$TIT_class <- TINC:::classification_tumour(fit$TIT)[["QC"]]
  row <- cbind(row, result)
  row$call <- if (row$TIT < p[["call_purity"]]) "not_assessable" else
    if (row[[p[["call_estimate"]]]] > p[["call_above"]]) "contaminated" else "clean"
} else {
  row$call <- "not_assessable"
}
write.table(row, snakemake@output[["tsv"]], sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(fit, snakemake@output[["rds"]])

pdf(snakemake@output[["pdf"]], width = 16, height = 18, bg = "white")
if (!is.null(fit)) print(TINC::plot_full_page_report(fit)) else { plot.new(); text(0.5, 0.5, row$reason) }
invisible(dev.off())
