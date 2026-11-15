# Changelog

All notable changes to the `somalier` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-11-15

Draft by Ryan Morin (2020), finished by Brett Collinge.

- `somalier extract` per sample at the somalier sites for its genome build, then one `somalier relate` over every
  digest in the run. Site files for GRCh37, hs37d5 and hg38 list the same sites in the same order, so samples aligned
  to different builds are compared directly.
- Outputs `99-outputs/relate.pairs.tsv` (relatedness, IBS0, IBS2, sites compared per pair) and
  `99-outputs/relate.samples.tsv` (per-sample depth, heterozygosity, X heterozygous sites, Y depth).
- Thresholds for calling pairs the same individual, and any sex call, are left to the caller.
