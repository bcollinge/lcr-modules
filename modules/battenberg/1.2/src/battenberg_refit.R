library(Battenberg)
library(optparse)
# Refit a finished Battenberg run at a preset rho/psi, following
# inst/example/battenberg_snp6_refitting.R with the WGS settings of
# battenberg_wgs_hg38.R. Inputs are read from the run directory; every
# output is written to the refit directory.
option_list = list(
  make_option(c("-t", "--tumourname"), type="character", default=NULL, help="Samplename of the tumour"),
  make_option(c("-i", "--input_dir"), type="character", default=NULL, help="Directory of the finished Battenberg run"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Directory where refit output will be written"),
  make_option(c("--rho"), type="numeric", default=NULL, help="Preset rho"),
  make_option(c("--psi"), type="numeric", default=NULL, help="Preset psi"),
  make_option(c("--sex"), type="character", default=NULL, help="Sex of the sample"),
  make_option(c("--reference"), type="character", default=NULL, help="Path to reference directory"),
  make_option(c("--chr_prefixed_genome"), type="logical", default=FALSE, action="store_true", help="Genome has chr prefixes")
)
opt = parse_args(OptionParser(option_list=option_list))

TUMOURNAME = opt$tumourname
IN_DIR = normalizePath(opt$input_dir)
IS.MALE = opt$sex=="male" | opt$sex=="Male"
CHR_PREFIXED = opt$chr_prefixed_genome
IMPUTEINFOFILE = paste0(opt$reference, "/impute_info.txt")

# Same values as battenberg_wgs_hg38.R
PLATFORM_GAMMA = 1
CLONALITY_DIST_METRIC = 0
ASCAT_DIST_METRIC = 1
MIN_PLOIDY = 1.6
MAX_PLOIDY = 4.8
MIN_RHO = 0.1
MIN_GOODNESS_OF_FIT = 0.63
BALANCED_THRESHOLD = 0.51
CALC_SEG_BAF_OPTION = 3

in_file = function(suffix) file.path(IN_DIR, paste0(TUMOURNAME, suffix))
bafsegmented_file = in_file(".BAFsegmented.txt")
baf_file = in_file("_mutantBAF.tab")
logr_file = in_file("_mutantLogR_gcCorrected.tab")
allelecounts_file = in_file("_alleleCounts.tab")

setwd(opt$output)
chrom_names = get.chrom.names(IMPUTEINFOFILE, IS.MALE)

fit.copy.number(samplename=TUMOURNAME,
                outputfile.prefix=paste0(TUMOURNAME, "_"),
                inputfile.baf.segmented=bafsegmented_file,
                inputfile.baf=baf_file,
                inputfile.logr=logr_file,
                dist_choice=CLONALITY_DIST_METRIC,
                ascat_dist_choice=ASCAT_DIST_METRIC,
                min.ploidy=MIN_PLOIDY,
                max.ploidy=MAX_PLOIDY,
                min.rho=MIN_RHO,
                min.goodness=MIN_GOODNESS_OF_FIT,
                uninformative_BAF_threshold=BALANCED_THRESHOLD,
                gamma_param=PLATFORM_GAMMA,
                use_preset_rho_psi=T,
                preset_rho=opt$rho,
                preset_psi=opt$psi,
                read_depth=30)

fig_prefix = if (CHR_PREFIXED) "_subclones_" else "_subclones_chr"
callSubclones(sample.name=TUMOURNAME,
              baf.segmented.file=bafsegmented_file,
              logr.file=logr_file,
              rho.psi.file=paste0(TUMOURNAME, "_rho_and_psi.txt"),
              output.file=paste0(TUMOURNAME, "_subclones.txt"),
              output.figures.prefix=paste0(TUMOURNAME, fig_prefix),
              output.gw.figures.prefix=paste0(TUMOURNAME, "_BattenbergProfile"),
              masking_output_file=paste0(TUMOURNAME, "_segment_masking_details.txt"),
              prior_breakpoints_file=NULL,
              chr_names=chrom_names,
              gamma=PLATFORM_GAMMA,
              segmentation.gamma=NA,
              siglevel=0.05,
              maxdist=0.01,
              noperms=1000,
              calc_seg_baf_option=CALC_SEG_BAF_OPTION)

make_posthoc_plots(samplename=TUMOURNAME,
                   logr_file=logr_file,
                   subclones_file=paste0(TUMOURNAME, "_subclones.txt"),
                   rho_psi_file=paste0(TUMOURNAME, "_rho_and_psi.txt"),
                   bafsegmented_file=bafsegmented_file,
                   logrsegmented_file=paste0(TUMOURNAME, ".logRsegmented.txt"),
                   allelecounts_file=allelecounts_file)

cnfit_to_refit_suggestions(samplename=TUMOURNAME,
                           subclones_file=paste0(TUMOURNAME, "_subclones.txt"),
                           rho_psi_file=paste0(TUMOURNAME, "_rho_and_psi.txt"),
                           gamma_param=PLATFORM_GAMMA)
