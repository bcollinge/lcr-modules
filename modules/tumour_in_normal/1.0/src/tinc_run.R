#!/usr/bin/env Rscript
# TINC on one pair: Strelka 2 PASS SNVs through TINC::load_VCF_Strelka, or a Mutect2 site table (mutect2_sites), Battenberg copy number in CNAqc columns
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
vcf_format <- snakemake@params[["vcf_format"]]

# Matched Mutect2 call to TINC's input table: chr1-chr22, records PASS or filtered only as germline and/or
# normal_artifact, biallelic SNVs, ref and alt counts from AD, tumour VAF strictly between 0 and 1 (as TINC's parser).
mutect2_sites <- function(lines) {
  hdr <- lines[startsWith(lines, "#")]
  meta <- function(key) sub(paste0("^##", key, "="), "", grep(paste0("^##", key, "="), hdr, value = TRUE))
  samples <- strsplit(tail(hdr, 1), "\t")[[1]][-(1:9)]
  t_col <- 9 + match(meta("tumor_sample"), samples)
  n_col <- 9 + match(meta("normal_sample"), samples)
  f <- do.call(rbind, strsplit(lines[!startsWith(lines, "#")], "\t"))
  chr <- ifelse(startsWith(f[, 1], "chr"), f[, 1], paste0("chr", f[, 1]))
  filters <- strsplit(f[, 7], ";")
  keep <- chr %in% paste0("chr", 1:22) & nchar(f[, 4]) == 1 & nchar(f[, 5]) == 1 &
    vapply(filters, function(x) all(x == "PASS") || all(x %in% c("germline", "normal_artifact")), logical(1))
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
  if (vcf_format == "strelka") {
    body <- lines[!startsWith(lines, "#")]
    body <- sub("^([0-9]+)\t", "chr\\1\t", body)
    tmp <- tempfile(fileext = ".vcf")
    writeLines(body, tmp)
  }
  cna <- if (cna_mode == "none") NULL else battenberg_segments(cna_file, cna_mode, filler)
  set.seed(seed)
  fit <- tryCatch({
    mutations <- if (vcf_format == "strelka") TINC::load_VCF_Strelka(tmp) else mutect2_sites(lines)
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
