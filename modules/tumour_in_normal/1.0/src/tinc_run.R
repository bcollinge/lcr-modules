#!/usr/bin/env Rscript
# TINC on one pair: Strelka 2 PASS SNVs through TINC::load_VCF_Strelka, Battenberg copy number in CNAqc columns
# when given, then TINC::autofit with its defaults. A pair with a missing input or that TINC cannot fit gets a
# not_assessable row carrying the reason.
suppressPackageStartupMessages({
  library(TINC)
  library(dplyr)
  library(tidyr) # load_VCF_Strelka in TINC 0.1.0 calls separate() without importing tidyr
})

log_con <- file(snakemake@log[["stderr"]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

vcf <- snakemake@params[["vcf"]]
cna_file <- snakemake@params[["cna"]]
cna_mode <- snakemake@wildcards[["cna_mode"]]
filler <- snakemake@params[["filler"]]
seed <- as.integer(snakemake@params[["seed"]])

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
                  cna_mode = cna_mode, seed = seed, status = "ok", reason = NA_character_)
fit <- NULL
missing_input <- c(vcf = !nzchar(vcf) || !file.exists(vcf),
                   cna = cna_mode != "none" && (!nzchar(cna_file) || !file.exists(cna_file)))
if (any(missing_input)) {
  row$status <- "not_assessable"
  row$reason <- paste("no input:", paste(names(missing_input)[missing_input], collapse = ", "))
} else {
  # load_VCF_Strelka keeps chr1-chr22 only; unprefixed contigs are prefixed so the same autosomes are read.
  lines <- readLines(gzfile(vcf))
  body <- lines[!startsWith(lines, "#")]
  body <- sub("^([0-9]+)\t", "chr\\1\t", body)
  tmp <- tempfile(fileext = ".vcf")
  writeLines(body, tmp)
  cna <- if (cna_mode == "none") NULL else battenberg_segments(cna_file, cna_mode, filler)
  set.seed(seed)
  fit <- tryCatch({
    mutations <- TINC::load_VCF_Strelka(tmp)
    TINC::autofit(input = mutations, cna = cna, FAST = isTRUE(snakemake@params[["fast"]]))
  }, error = function(e) {
    row$status <<- "not_assessable"
    row$reason <<- gsub("[\t\n]+", " ", conditionMessage(e))
    NULL
  })
}

if (!is.null(fit)) {
  result <- TINC::to_string(fit)
  result$TIN_class <- TINC:::classification_normal(fit$TIN)[["QC"]]
  result$TIT_class <- TINC:::classification_tumour(fit$TIT)[["QC"]]
  row <- cbind(row, result)
}
write.table(row, snakemake@output[["tsv"]], sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(fit, snakemake@output[["rds"]])

pdf(snakemake@output[["pdf"]], width = 16, height = 18)
if (!is.null(fit)) print(TINC::plot_full_page_report(fit)) else { plot.new(); text(0.5, 0.5, row$reason) }
invisible(dev.off())
