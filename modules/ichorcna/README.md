# ichorcna

# Purpose

The `ichorcna` is a Level 2 module that operates on `BAM/CRAM` files to use probabilistic model, implemented as a hidden Markov model (ichorCNA), and perform CNV calling for `genome` data. It generates `SEG` files as outputs.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    ichorcna:
        inputs:
            sample_bam: "data/{sample_id}.bam"
            sample_bai: "data/{sample_id}.bam.bai"
        options:
            bin_sizes: [500000]
            panel_of_normals:
                members:
                    my_panel: ["normal_1", "normal_2", "normal_3"]
                assign:
                    _default: "my_panel"
                exclude:
                    tumour_1: ["normal_2"]
```

Version 1.2 builds its panels of normals from samples in the sample table: every sample with a panel in
`assign` (or `_default`) is run against that panel, and every panel member is counted with readCounter.
A sample is left out of its own panel, and `exclude` removes further members for one sample (e.g. its
matched normal); either gives that sample its own panel file.

The example snakefile:

```python
#!/usr/bin/env snakemake

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")

subworkflow reference_files:
    workdir:
        "reference/"
    snakefile:
        "../workflows/reference_files/2.4/reference_files.smk"
    configfile:
        "../workflows/reference_files/2.4/config/default.yaml"

configfile: "../modules/ichorcna/1.2/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/ichorcna/1.2/ichorcna.smk"

rule all:
    input:
        rules._ichorcna_all.input
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/ichorcna/CHANGELOG.md)
