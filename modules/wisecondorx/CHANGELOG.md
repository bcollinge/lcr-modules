# Changelog

All notable changes to the `wisecondorx` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-10-15

This release was authored by Brett Collinge.

- Runs WisecondorX 1.2.9 (bioconda) in three steps: `convert` (5 kb bins, BAM or CRAM), `newref` (100 kb) on
  config-declared reference sets, and `predict` per sample against the reference sets named for it.
- Reference sets and predictions are listed in the config (`references`, `predict`), so one sample can be
  scored against several sets and per-sample sets (matched-normal exclusion, leave-one-out) are plain config.
- Stops before running if a set has fewer than 10 members (the `newref` minimum), spans genome builds, or holds a
  sample of the predicted sample's patient.
- WisecondorX 1.3.2, the latest bioconda build at release, drops chromosome 10 in `convert` (its chromosome regex
  skips `10`; fixed upstream in 1.3.3, not yet on bioconda), so 1.2.9 is pinned.
- `predict` runs with `--seed 1` so segmentation is reproducible.
