# Changelog

All notable changes to the `tumour_in_normal` module will be documented in this file.

## [1.0] - 2026-11-11

Initial release. Runs two published tumour-in-normal estimators on matched pairs:

- deTiN (Taylor-Weiner et al., Nat Methods 2019), fetched at commit d296e11 and run from its own directory. Its SNV
  table is built from a Mutect2 FilterMutectCalls VCF (run with `--genotype-germline-sites`): PASS records are
  KEEP, records filtered only by the configured germline filters are written as `alt_allele_in_normal`. With copy
  number it runs twice per pair, the joint estimate and `--only_ascnas`; without, SNVs only. Each model's own
  estimate and the messages deTiN prints are kept beside its written estimate. Runs on hg19-coordinate builds only.
- TINC (Caravagna lab, Nat Commun 2024), from bioconda r-tinc 0.1.0: Strelka 2 PASS SNVs through
  `TINC::load_VCF_Strelka`, with Battenberg copy number in CNAqc columns per configured state, then
  `TINC::autofit` with a fixed seed. A pair with no input, or that TINC cannot fit, is reported not_assessable.
- A report of pairs and of normals, with the call vocabulary contaminated, low_level, clean, not_assessable.
