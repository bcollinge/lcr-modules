# wisecondorx

## Purpose

`wisecondorx` is a Level 2 module that calls copy number from shallow whole-genome BAM/CRAM files with
[WisecondorX](https://github.com/CenterForMedicalGeneticsGhent/WisecondorX) (Raman et al. 2019,
[doi:10.1093/nar/gky1263](https://doi.org/10.1093/nar/gky1263)). Each sample is normalised against a
reference set of negative-control samples built by the module, so no matched normal is needed.

Version 1.0 runs WisecondorX 1.2.9 (bioconda `wisecondorx=1.2.9`). Settings follow the
[WisecondorX 1.2.9 README](https://github.com/CenterForMedicalGeneticsGhent/WisecondorX/blob/v1.2.9/README.md).

## Rules

| Rule | Command | Output |
|---|---|---|
| `_wisecondorx_input_bam` | symlink BAM/CRAM and index, keeping the extension | `00-inputs/bam/{seq_type}--{genome_build}/{sample_id}.{bam,cram}` |
| `_wisecondorx_convert` | `WisecondorX convert --reference <genome.fa> --binsize 5000` | `01-convert/{seq_type}--{genome_build}/{sample_id}.npz` |
| `_wisecondorx_newref` | `WisecondorX newref --binsize 100000 --refsize 300 --cpus N [--yfrac x]` | `02-newref/{genome_build}/{reference}.npz` |
| `_wisecondorx_predict` | `WisecondorX predict --bed --zscore 5 --alpha 0.0001 --maskrepeats 5 --minrefbins 150 --seed 1` | `03-predict/.../{sample_id}--{reference}_{bins,segments,aberrations}.bed`, `_statistics.txt` |
| `_wisecondorx_output` | symlinks | `99-outputs/bed/{seq_type}--{genome_build}/{sample_id}--{reference}_{bins,segments,aberrations}.bed`, `99-outputs/txt/.../{sample_id}--{reference}_statistics.txt` |

`convert` reads CRAM with the `genome_build` fasta from `reference_files`; WisecondorX picks the reader from the file
extension, so the input symlink keeps `.bam` or `.cram`. `convert` keeps reads with mapping quality >= 1 and drops
improper pairs and its own position duplicates (no `--normdup`); the README asks that no other read filtering be applied
beforehand. The `convert` bin size (5 kb) does not set the resolution; `newref --binsize` (100 kb) does, and must be a
multiple of it.

## Configuration

Reference sets and predictions are declared in the module config, not derived from pairing:

```yaml
lcr-modules:
    wisecondorx:
        inputs:
            sample_bam: "data/{sample_id}.cram"
            sample_bai: "data/{sample_id}.cram.crai"
        references:
            ffpe_normals: [N01, N02, N03, N04, N05, N06, N07, N08, N09, N10, N11]
        predict:
            T01: [ffpe_normals]
            T02: [ffpe_normals]
```

- `references`: name -> sample_ids. Members must be in the samples table at one `genome_build`, under one `seq_type`.
  A sample may sit in several sets, and a set may be built for one sample only (for example a set without a tumour's
  own matched normal, or a leave-one-out set for scoring a normal).
- `predict`: sample_id -> reference names. The module stops if a reference set holds a sample with the predicted
  sample's `patient_id`, so a tumour is never normalised against its own patient's normal.
- `options.yfrac`: `null` lets `newref` fit its Gaussian mixture to the chrY read fractions to assign sex; a number
  sets the cut-off by hand (`1` treats every sample as female).

`--refsize` (newref) and `--minrefbins` (predict) count within-sample reference *bins* per target bin, not reference
samples; the README says neither should generally be tweaked.

`pairing_config` only has to list each `seq_type` so `oncopipe` keeps the samples; runs are not used.

## Reference set size

The WisecondorX README asks for at least 50 reference samples per set ("the more the better", with little gain beyond
500), all negative controls from the same sequencer, mapper, reference genome and material as the test samples.
`newref` refuses fewer than 10 samples, and the module checks this before building. Sets between 10 and 50 run, but
fall below the authors' recommendation; with small sets `predict` may log "reference too small" and fall back to
unweighted segmentation and z-scores.

Sex: the Gaussian mixture is described as unreliable with fewer than 20 samples or without both sexes. `newref`
builds a female gonosomal reference only from >= 5 female samples and a male one only from >= 5 male samples;
otherwise samples of the missing sex are normalised on the other sex's chrX / chrY reference (`predict` warns). The
mixture fit can also fail outright (`IndexError` in `train_gender_model`) when it finds no local minimum on 0-0.02;
set `options.yfrac` then. Autosomes are unaffected.

## Outputs

- `_bins.bed`: per bin, log2 ratio and z-score against the within-sample reference bins; `nan` = masked bin.
- `_segments.bed`: weighted CBS segments, log2 ratio and between-sample z-score.
- `_aberrations.bed`: segments beyond `--zscore` (or `--beta` when given).
- `_statistics.txt`: per-chromosome ratio and z-score, sex call, read count, median segment variance and CPA score.

Chromosomes are written without a `chr` prefix (`1`-`22`, `X`, `Y`), starts 1-based.
