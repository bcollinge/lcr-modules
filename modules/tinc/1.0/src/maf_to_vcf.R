#!/usr/bin/env Rscript
# A tumour-normal MAF to the two-sample VCF tinc_run.R reads: single-base substitutions only, AD = ref,alt counts
# of each sample, FILTER PASS, sample columns named by the ##tumor_sample / ##normal_sample header lines.
log_con <- file(snakemake@log[["stderr"]], open = "wt")
sink(log_con, type = "message")

tumour_id <- snakemake@wildcards[["tumour_id"]]
normal_id <- snakemake@wildcards[["normal_id"]]

maf <- read.delim(snakemake@input[["maf"]], comment.char = "#", quote = "", colClasses = "character",
                  na.strings = character(0))
needed <- c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type",
            "t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count")
missing <- setdiff(needed, names(maf))
if (length(missing) > 0) stop("MAF lacks columns: ", paste(missing, collapse = ", "))

snv <- maf[maf$Variant_Type == "SNP" & nchar(maf$Reference_Allele) == 1 & nchar(maf$Tumor_Seq_Allele2) == 1, ]
counts <- c("t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count")
has_counts <- Reduce(`&`, lapply(snv[counts], function(x) grepl("^[0-9]+$", x)))
message(nrow(maf), " MAF rows; ", nrow(snv), " single-base substitutions; ", sum(!has_counts),
        " of them without integer counts, dropped")
snv <- snv[has_counts, ]

header <- c("##fileformat=VCFv4.2",
            paste0("##tumor_sample=", tumour_id),
            paste0("##normal_sample=", normal_id),
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Ref and alt read counts from the MAF">',
            paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", normal_id, tumour_id),
                  collapse = "\t"))
body <- if (nrow(snv) == 0) character(0) else paste(
  snv$Chromosome, snv$Start_Position, ".", snv$Reference_Allele, snv$Tumor_Seq_Allele2, ".", "PASS", ".", "AD",
  paste(snv$n_ref_count, snv$n_alt_count, sep = ","), paste(snv$t_ref_count, snv$t_alt_count, sep = ","), sep = "\t")
writeLines(c(header, body), snakemake@output[["vcf"]])
