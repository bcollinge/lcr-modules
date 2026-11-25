# Changelog

All notable changes to the `tinc` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-11-24

This release was authored by Brett Collinge.

- TINC 0.1.0 (Caravagna lab, bioconda r-tinc) on matched tumour-normal pairs: tumour-in-normal (TIN) and
  tumour-in-tumour (TIT) estimates from a tumour-normal MAF that carries the normal's read counts, with Battenberg
  allele-specific copy number.
- The MAF's single-base substitutions are written as a two-sample VCF (AD from the MAF's ref and alt counts) and
  read by the runner; TINC::autofit runs with its defaults under a fixed seed.
- One row per pair with TINC's estimates, its own QC classes and a call (contaminated / clean / not_assessable)
  from `options.call`; a pair TINC cannot fit gets a not_assessable row with the reason. Merged per seq_type and
  genome build in `99-outputs/merged/`.
