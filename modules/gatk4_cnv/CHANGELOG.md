# Changelog

All notable changes to the `gatk4_cnv` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-10-30

This release was authored by Brett Collinge.

- Shallow whole-genome copy number from read depth: GATK 4.1.8.1 `PreprocessIntervals`, `AnnotateIntervals` (GC),
  `CollectReadCounts`, `CreateReadCountPanelOfNormals` and `DenoiseReadCounts`, then CBS (DNAcopy) of the denoised
  log2 copy ratio, a Rascal purity and ploidy fit and the purity-corrected profile.
- A sample is denoised against its configured panel without the normals of its own patient.
- The Rascal solution of record is chosen with the ACE ploidy and cellularity penalties on Rascal's distance; fits with
  a near-tie runner-up of different purity are flagged as `ambiguous_fit`.
- Outputs per sample: the lcr seg, the purity-corrected bin profile, the purity table and a profile plot.
