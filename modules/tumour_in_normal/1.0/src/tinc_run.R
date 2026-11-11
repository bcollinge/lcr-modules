#!/usr/bin/env Rscript
# TINC on one pair: Strelka 2 PASS SNVs through TINC::load_VCF_Strelka, Battenberg clone-A copy number
# in CNAqc columns when given, then TINC::autofit with its defaults.
suppressPackageStartupMessages({
  library(TINC)
  library(dplyr)
})

log_con <- file(snakemake@log[["stderr"]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

vcf <- snakemake@input[["vcf"]]
cna_file <- snakemake@input[["cna"]]

# load_VCF_Strelka keeps chr1-chr22 only; unprefixed contigs are prefixed so the same autosomes are read.
lines <- readLines(gzcon(file(vcf, "rb")))
body <- lines[!startsWith(lines, "#")]
body <- sub("^([0-9]+)\t", "chr\\1\t", body)
tmp <- tempfile(fileext = ".vcf")
writeLines(body, tmp)
mutations <- TINC::load_VCF_Strelka(tmp)

cna <- NULL
if (length(cna_file) > 0 && nzchar(cna_file)) {
  bb <- read.delim(cna_file, check.names = FALSE)
  cna <- tibble::tibble(
    chr = ifelse(startsWith(as.character(bb$chr), "chr"), as.character(bb$chr), paste0("chr", bb$chr)),
    from = as.numeric(bb$startpos),
    to = as.numeric(bb$endpos),
    Major = as.integer(bb$nMaj1_A),
    minor = as.integer(bb$nMin1_A)
  ) %>% dplyr::filter(!is.na(Major), !is.na(minor))
}

fit <- TINC::autofit(input = mutations, cna = cna, FAST = isTRUE(snakemake@params[["fast"]]))

result <- TINC::to_string(fit)
result$TIN_class <- TINC:::classification_normal(fit$TIN)[["QC"]]
result$TIT_class <- TINC:::classification_tumour(fit$TIT)[["QC"]]
result$tumour_id <- snakemake@wildcards[["tumour_id"]]
result$normal_id <- snakemake@wildcards[["normal_id"]]
result$cna_used <- !is.null(cna)
write.table(result, snakemake@output[["tsv"]], sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(fit, snakemake@output[["rds"]])

pdf(snakemake@output[["pdf"]], width = 16, height = 18)
print(TINC::plot_full_page_report(fit))
invisible(dev.off())
