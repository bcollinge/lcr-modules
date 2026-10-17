# gatk_tangent

## Purpose

`gatk_tangent` is a Level 2 module that denoises whole-genome read counts against panels of normals with the GATK4
somatic CNV tools: `CreateReadCountPanelOfNormals` (SVD of the standardized panel counts) and `DenoiseReadCounts`
(projection of a sample onto the panel's eigensamples, the Tangent approach of Tabak et al. 2019 as implemented in
GATK). It bins whole contigs, so it suits shallow genomes where `gatk_cnv` (capture targets, allelic counts, PureCN)
does not apply. `ModelSegments` on the denoised copy ratios alone is optional.

Version 1.0 runs GATK 4.1.8.1 (bioconda `gatk4=4.1.8.1`), the env of `gatk_cnv` 1.0.

## Rules

| Rule | Command | Output |
|---|---|---|
| `_gatk_tangent_input_bam` | symlink BAM/CRAM and index | `00-inputs/bam/{seq_type}--{genome_build}/{sample_id}.{bam,cram}` |
| `_gatk_tangent_preprocess_intervals` | `PreprocessIntervals -L <contig>... --bin-length --padding` | `01-intervals/{genome_build}/preprocessed.interval_list` |
| `_gatk_tangent_annotate_intervals` | `AnnotateIntervals` (GC content) | `01-intervals/{genome_build}/annotated.tsv` |
| `_gatk_tangent_collect_read_counts` | `CollectReadCounts -R <genome.fa> --format HDF5` | `02-counts/{seq_type}--{genome_build}/{sample_id}.counts.hdf5` |
| `_gatk_tangent_create_pon` | `CreateReadCountPanelOfNormals --annotated-intervals` | `03-pon/{genome_build}/{panel}.pon.hdf5` |
| `_gatk_tangent_denoise` | `DenoiseReadCounts --count-panel-of-normals`, once per `options.denoise` label | `04-denoised/{variant}/.../{sample_id}--{panel}.{standardizedCR,denoisedCR}.tsv` |
| `_gatk_tangent_model_segments` | `ModelSegments --denoised-copy-ratios`, for the labels in `options.model_segments` | `05-model/{variant}/.../{sample_id}--{panel}.cr.seg` |
| `_gatk_tangent_output_*` | symlinks | `99-outputs/{denoisedCR,cr_seg}/{variant}/{seq_type}--{genome_build}/`, `99-outputs/pon/{genome_build}/` |

`PreprocessIntervals` drops bins made only of N. `AnnotateIntervals` adds GC content only; no mappability or
segmental-duplication track is passed. CRAM input is decoded with the `genome_build` fasta from `reference_files`;
GATK 4.1.8.1 (htsjdk 2.23) reads CRAM 3.0 but not 3.1.

## Configuration

Panels and denoising runs are declared in the module config, not derived from pairing:

```yaml
lcr-modules:
    gatk_tangent:
        inputs:
            sample_bam: "data/{sample_id}.cram"
            sample_bai: "data/{sample_id}.cram.crai"
        panels:
            ffpe_normals: [N01, N02, N03, N04, N05, N06]
            ffpe_normals_minus_N01: [N02, N03, N04, N05, N06]
        denoise:
            T01: [ffpe_normals]
            N01: [ffpe_normals_minus_N01]
```

Every panel member and denoised sample must be in the samples table. The module stops before running if a panel
spans genome builds or holds a sample of the denoised sample's patient.

`options.create_pon` holds the `CreateReadCountPanelOfNormals` arguments (the GATK defaults). The requested
`--number-of-eigensamples` is capped at the number of panel samples left after filtering. `options.denoise` maps a
label to `DenoiseReadCounts` arguments, one run per label (default `tangent: ""`, every eigensample in the panel); a
second label such as `tangent_k2: "--number-of-eigensamples 2"` denoises against the same panel with fewer
eigensamples. With a small panel the
sample filter `--extreme-sample-median-percentile` can be set to 0 to keep every member; the panel log reports the
samples and intervals kept and the number of eigensamples.

The Spark-backed tools run with pure-Java BLAS (`F2jBLAS`), as in `gatk_cnv` 1.0.
