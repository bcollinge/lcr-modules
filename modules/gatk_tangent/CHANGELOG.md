# Changelog

All notable changes to the `gatk_tangent` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-10-17

This release was authored by Brett Collinge.

- Runs GATK 4.1.8.1 `PreprocessIntervals` (fixed-length bins over configured contigs), `AnnotateIntervals` (GC),
  `CollectReadCounts`, `CreateReadCountPanelOfNormals`, `DenoiseReadCounts` and, optionally, `ModelSegments` on the
  denoised copy ratios alone.
- `options.denoise` lists `DenoiseReadCounts` runs by label (for example all eigensamples and a fixed count), each
  written under its own label.
- Panels and denoising runs are listed in the config (`panels`, `denoise`), so per-sample panels (patient exclusion,
  leave-one-out) are plain config.
- Stops before running if a panel spans genome builds or holds a sample of the denoised sample's patient.
- Reuses the `gatk_cnv` 1.0 GATK env and its pure-Java BLAS workaround.
