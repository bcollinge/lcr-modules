# Changelog

All notable changes to the `qc` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1] - 2026-11-15

- Adds per-sample contamination: GATK 4.1.8.1 GetPileupSummaries at common SNPs (`contamination_sites`, keyed by
  genome build and capture space) and CalculateContamination on the sample alone, for the seq_types in
  `contamination_seq_types`. Pileup and contamination tables are under `99-outputs/contamination/`, so a matched-normal
  CalculateContamination can be run downstream from the two pileups.
- Adds coverage evenness, GC bias and complexity per sample (`complexity_seq_types`): mosdepth with and without
  duplicates (base-level CV and share of bases at >= 0.2x the mean; per-target or per-bin CV and share of regions
  within 0.2-5x the median; for non-capture samples a GC-stratified bin dispersion index, 1 = counting noise),
  CollectGcBiasMetrics on every seq_type, and preseq lc_extrap reported at a fixed read count per seq_type
  (`preseq_fixed_reads`). One row per sample under `99-outputs/complexity/`.
- The per-sample metrics table also carries LibrarySize (HS_LIBRARY_SIZE), ProportionExcludedDuplicates, Fold80BasePenalty,
  ATDropout and GCDropout from the Picard metrics already collected (NA where the tool does not report one).
- All other 1.0 rules and outputs are unchanged.

## [1.0] - 2022-05-18

This release was authored by Kostia Dreval.

<!-- TODO: Explain each important module design decision below. -->

- This module replaces the current `picard_qc` module.
- Verion 1 uses same bed file for both baits and target regions. This can be separated in the future versions if this functiionality will be of interest.
- The module generates supplementary table-ready tab-deliminated file with QC metrics to be reported as a good practise described by MIRAGE.
- The `MeanCorrectedCoverage` means that the filters for unmapped, low-quality, duplicated, secondary, and supplemental reads has been applied before calculating coverage. The reported value is effective coverage.
- Additional flags for each tool can be specified via config.
