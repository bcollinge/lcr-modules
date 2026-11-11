# Changelog

All notable changes to the `tumour_in_normal` module will be documented in this file.

## [1.0] - 2026-11-11

Initial release. Runs two published tumour-in-normal estimators on matched pairs:

- deTiN (Taylor-Weiner et al., Nat Methods 2018), fetched at commit d296e11 and run from its own directory,
  once per available arm (SNVs only, copy number only, both). Inputs are a MuTect call-stats table, a population
  VCF, and allelic segments with tumour and normal het counts; each is optional.
- TINC (Caravagna lab, Nat Commun 2024), from bioconda r-tinc 0.1.0: Strelka 2 PASS SNVs through
  `TINC::load_VCF_Strelka`, with Battenberg clone-A copy number in CNAqc columns when given, then `TINC::autofit`.
