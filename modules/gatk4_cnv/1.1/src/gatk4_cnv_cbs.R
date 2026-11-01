#!/usr/bin/env Rscript
# CBS of a GATK4 denoised copy ratio profile.
# Usage: gatk4_cnv_cbs.R <denoisedCR.tsv> <out ratio.tsv> <contigs, comma-separated> <seed key> <seed base> <floor> <alpha> <undo SD>
# Out, one row per bin of the listed contigs in contig order: Chromosome, Start, Ratio (denoised copy ratio),
# MedianRatio (median copy ratio of the bin's CBS segment, 6 decimals), CopyNumber (round half up of 2 x MedianRatio), End.
suppressPackageStartupMessages({library(data.table); library(DNAcopy)})
a = commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 8)
contigs = strsplit(a[3], ",", fixed = TRUE)[[1]]
seed = as.numeric(a[5]) + {h = 0; for (b in utf8ToInt(a[4])) h = (h * 31 + b) %% 1000003; h}
floor_cr = as.numeric(a[6])

w = fread(cmd = paste("grep -v '^@'", shQuote(a[1])), colClasses = list(character = "CONTIG"))
w = w[CONTIG %in% contigs][order(match(CONTIG, contigs), START)]
stopifnot(nrow(w) > 0, all(is.finite(w$LOG2_COPY_RATIO)))
lr = log2(pmax(2^w$LOG2_COPY_RATIO, floor_cr))
cna = CNA(lr, match(w$CONTIG, contigs), w$START, data.type = "logratio", sampleid = "s", presorted = TRUE)
set.seed(seed)
s = segment(cna, alpha = as.numeric(a[7]), undo.splits = "sdundo", undo.SD = as.numeric(a[8]), verbose = 0)
lvl = rep(s$output$seg.mean, s$output$num.mark)
stopifnot(length(lvl) == length(lr))
mr = round(2^ave(lr, rleid(w$CONTIG, lvl), FUN = median), 6)
fwrite(w[, .(Chromosome = CONTIG, Start = START, Ratio = 2^LOG2_COPY_RATIO, MedianRatio = mr, CopyNumber = floor(2 * mr + 0.5), End = END)],
       a[2], sep = "\t")
cat("bins", nrow(w), "segments", nrow(s$output), "seed", seed, "\n", file = stderr())
