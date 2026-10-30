# gatk4_cnv

## Purpose

`gatk4_cnv` is a Level 2 module for copy number from shallow (low-pass) whole genomes without a matched normal. Read
counts in fixed-length bins are denoised against a panel of normals with the GATK4 somatic CNV tools
(`CreateReadCountPanelOfNormals`, `DenoiseReadCounts`), segmented with CBS (DNAcopy), fitted for purity and ploidy with
Rascal and purity corrected. No allelic counts are used.

## Rules

| Rule | Command | Output |
|---|---|---|
| `_gatk4_cnv_input_bam` | symlink BAM/CRAM and index | `00-inputs/bam/{seq_type}--{genome_build}/{sample_id}.{bam,cram}` |
| `_gatk4_cnv_preprocess_intervals` | `PreprocessIntervals -L <contig>... --bin-length --padding` | `01-intervals/{genome_build}/preprocessed.interval_list` |
| `_gatk4_cnv_annotate_intervals` | `AnnotateIntervals` (GC content) | `01-intervals/{genome_build}/annotated.tsv` |
| `_gatk4_cnv_collect_read_counts` | `CollectReadCounts --format HDF5` | `02-counts/{seq_type}--{genome_build}/{sample_id}.counts.hdf5` |
| `_gatk4_cnv_create_pon` | `CreateReadCountPanelOfNormals --annotated-intervals` | `03-pon/{genome_build}/{panel}.pon.hdf5` |
| `_gatk4_cnv_denoise` | `DenoiseReadCounts` | `04-denoised/{seq_type}--{genome_build}/{sample_id}--{panel}.denoisedCR.tsv` |
| `_gatk4_cnv_cbs` | `src/gatk4_cnv_cbs.R` | `05-cbs/.../{sample_id}--{panel}.ratio.tsv` |
| `_gatk4_cnv_fit` | `src/gatk4_cnv_fit.R` | `06-fit/.../{sample_id}--{panel}.{purity,solutions}.tsv` |
| `_gatk4_cnv_correct` | `src/gatk4_cnv_correct.R` | `07-corrected/.../{sample_id}--{panel}.{ratio.tsv,seg,png}` |
| `_gatk4_cnv_output_*` | symlinks | `99-outputs/{seg,ratio,purity,plot}/{seq_type}--{genome_build}/`, `99-outputs/pon/{genome_build}/` |

- **CBS**: copy ratios are floored at `options.cbs.floor` before log2; `segment(alpha, undo.splits = "sdundo", undo.SD)`
  with `set.seed` from `options.cbs.seed_key`, so reruns give the same segments. `MedianRatio` is the median copy ratio
  of the bin's segment.
- **Fit**: Rascal `find_best_fit_solutions` on the segments of the profile divided by its median `MedianRatio`. The
  solution of record has the lowest `distance x (1 + |ploidy - 2|)^penploidy / cellularity^penalty` (ACE's penalties);
  `ambiguous_fit` is TRUE when a solution within `ambiguous_distance` of the lowest distance has a purity more than
  `ambiguous_purity` away from the lowest-distance solution's.
- **Correction**: `r' = max((r / m - (1 - p)) / p, 0)`, `m` the median `MedianRatio`, `p` the purity. Without a fit the
  profile is divided by `m` only.
- **Seg**: `ID, chrom, start, end, LOH_flag (NA), log.ratio`, one row per run of equal corrected `MedianRatio`;
  `log.ratio` is log2 of the corrected value floored at `options.cbs.floor`. No blacklist is applied.

## Configuration

- `panels`: panel name -> normal sample_ids. `denoise`: sample_id -> panel name. A panel holding normals of the
  sample's patient is rebuilt without them as `<panel>_minus_<patient_id>`.
- `conda_envs.rascal`: R with DNAcopy, data.table, ggplot2 and Rascal 0.7.0, which is not on conda; the env's
  post-deploy script installs it from GitHub. An existing environment can be named by its prefix path instead.
