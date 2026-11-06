# Changelog

All notable changes to the `capture_kit` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-11-05

This release was authored by Brett Collinge.

- Calls the exome capture kit of a BAM/CRAM from coverage. For every ordered pair of
  configured kits (A, B) the panel holds regions of A clear of B by `panel.gap` bp
  (segdups removed), plus a core shared by all kits and an off-target floor. Depth per
  set is divided by core depth.
- For kit A, `own` is the median over B of ratio(A not B) and `against` the max over B of
  ratio(B not A); margin = own - against. A kit passes when own >= `hi`, against <= `lo`
  and margin >= `min_margin`. A passing kit covered by another passing kit is dropped,
  which resolves nested designs (e.g. V5 inside V5+UTR inside V6+UTR).
- Labels: `kit`, `mixed` (two kits pass, or no kit passes but the best kit's own >= hi
  and the kit setting its against has its own exclusive territory >= lo), `none_known`,
  `too_shallow` (core depth < `min_core_depth`), `genome_like` (floor/core >=
  `genome_floor`). The component names of a mixture can name a nested design's sibling.
- Default thresholds were set on a known-answer set of Agilent V5, V5+UTR, V6+UTR, V7,
  IDT xGen v1 and v2 exomes, kit-mixed merges and genomes.
- Kits whose target files are identical cannot be separated by coverage.
- Cost is one single-threaded `samtools bedcov` per sample, dominated by reading the file.
