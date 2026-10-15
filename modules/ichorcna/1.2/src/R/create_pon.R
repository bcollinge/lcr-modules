#!/usr/bin/env Rscript
# Command-line entry to ichorCNA::createPanelOfNormals. Wig, centromere and
# extdata arguments may be bare file names, resolved in the package extdata.

suppressPackageStartupMessages({
  library(optparse)
  library(ichorCNA)
})

option_list <- list(
  make_option("--filelist", type = "character", help = "File listing one normal wig per line. Required."),
  make_option("--outfile", type = "character", help = "Output prefix; writes <outfile>_<method>.rds. Required."),
  make_option("--gcWig", type = "character", help = "GC wig at the bin size. Required."),
  make_option("--mapWig", type = "character", default = NULL, help = "Mappability wig at the bin size."),
  make_option("--centromere", type = "character", help = "Centromere table. Required."),
  make_option("--chrs", type = "character", default = "c(1:22, \"X\")"),
  make_option("--chrNormalize", type = "character", default = "c(1:22)"),
  make_option("--genomeStyle", type = "character", default = "NCBI"),
  make_option("--genomeBuild", type = "character", default = "hg19"),
  make_option("--minMapScore", type = "double", default = 0.0),
  make_option("--fracReadsInChrYForMale", type = "double", default = 0.001),
  make_option("--method", type = "character", default = "median")
)
opt <- parse_args(OptionParser(option_list = option_list))

extdata <- function(f) {
  if (is.null(f) || file.exists(f)) return(f)
  p <- system.file("extdata", f, package = "ichorCNA")
  if (p == "") stop("not a file and not in ichorCNA extdata: ", f)
  p
}

createPanelOfNormals(
  gcWig = extdata(opt$gcWig),
  mapWig = extdata(opt$mapWig),
  filelist = opt$filelist,
  outfile = opt$outfile,
  centromere = extdata(opt$centromere),
  chrs = opt$chrs,
  chrNormalize = opt$chrNormalize,
  genomeStyle = opt$genomeStyle,
  genomeBuild = opt$genomeBuild,
  minMapScore = opt$minMapScore,
  fracReadsInChrYForMale = opt$fracReadsInChrYForMale,
  method = opt$method
)
