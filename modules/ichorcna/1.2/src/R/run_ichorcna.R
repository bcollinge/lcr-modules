#!/usr/bin/env Rscript
# Command-line entry to ichorCNA::run_ichorCNA. Wig, centromere and extdata
# arguments may be bare file names, resolved in the package extdata. Vector
# arguments are R expressions as strings, e.g. --ploidy "c(2,3)".

suppressPackageStartupMessages({
  library(optparse)
  library(ichorCNA)
})

chr <- function(name, default = NULL) make_option(paste0("--", name), type = "character", default = default)
num <- function(name, default) make_option(paste0("--", name), type = "double", default = default)
lgl <- function(name, default) make_option(paste0("--", name), type = "logical", default = default)

option_list <- list(
  chr("id"), chr("WIG"), chr("normalPanel"), chr("gcWig"), chr("mapWig"), chr("centromere"),
  chr("exons.bed"), chr("sex"), chr("outDir", "./"),
  chr("normal", "0.5"), chr("ploidy", "2"), chr("scStates", "NULL"),
  chr("chrs", "c(1:22, \"X\")"), chr("chrTrain", "c(1:22)"), chr("chrNormalize", "c(1:22)"),
  chr("genomeStyle", "NCBI"), chr("genomeBuild", "hg19"),
  chr("plotFileType", "pdf"), chr("plotYLim", "c(-2,2)"),
  num("maxCN", 7), num("minMapScore", 0.9), num("txnE", 0.9999999), num("txnStrength", 1e7),
  num("maxFracGenomeSubclone", 0.5), num("maxFracCNASubclone", 0.7),
  num("fracReadsInChrYForMale", 0.001), num("cores", 1),
  lgl("estimateNormal", TRUE), lgl("estimatePloidy", TRUE), lgl("estimateScPrevalence", TRUE),
  lgl("includeHOMD", FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))

extdata <- function(f) {
  if (is.null(f) || file.exists(f)) return(f)
  p <- system.file("extdata", f, package = "ichorCNA")
  if (p == "") stop("not a file and not in ichorCNA extdata: ", f)
  p
}

run <- function() run_ichorCNA(
  id = opt$id,
  tumor_wig = opt$WIG,
  normal_panel = opt$normalPanel,
  gcWig = extdata(opt$gcWig),
  mapWig = extdata(opt$mapWig),
  centromere = extdata(opt$centromere),
  exons.bed = opt$exons.bed,
  sex = opt$sex,
  normal = opt$normal,
  ploidy = opt$ploidy,
  scStates = opt$scStates,
  maxCN = opt$maxCN,
  estimateNormal = opt$estimateNormal,
  estimatePloidy = opt$estimatePloidy,
  estimateScPrevalence = opt$estimateScPrevalence,
  includeHOMD = opt$includeHOMD,
  chrs = opt$chrs,
  chrTrain = opt$chrTrain,
  chrNormalize = opt$chrNormalize,
  genomeStyle = opt$genomeStyle,
  genomeBuild = opt$genomeBuild,
  minMapScore = opt$minMapScore,
  txnE = opt$txnE,
  txnStrength = opt$txnStrength,
  maxFracGenomeSubclone = opt$maxFracGenomeSubclone,
  maxFracCNASubclone = opt$maxFracCNASubclone,
  fracReadsInChrYForMale = opt$fracReadsInChrYForMale,
  plotFileType = opt$plotFileType,
  plotYLim = opt$plotYLim,
  outDir = opt$outDir,
  cores = opt$cores
)

# Upstream bug (GavinHaLab/ichorCNA#29): with more than 10,000 bins and no repTimeWig, plotSolutions fails
# in plotFitCompareByChr after every output is written. Only that error, with all outputs present, passes.
outputs <- file.path(opt$outDir, paste0(opt$id, c(".cna.seg", ".seg", ".seg.txt", ".params.txt", ".correctedDepth.txt")))
tryCatch(run(), error = function(e) {
  if (grepl("values(x)[[y]]", conditionMessage(e), fixed = TRUE) && all(file.exists(outputs))) {
    message("ignored post-output plotting error (ichorCNA#29): ", conditionMessage(e))
  } else {
    stop(e)
  }
})
