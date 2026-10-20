#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Brett Collinge
# Module Author:    Brett Collinge
# Contributors:     N/A


##### SETUP #####


import sys

import oncopipe as op
import pandas as pd

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
from importlib.metadata import version as pkg_version
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

current_version = pkg_version("oncopipe")
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section


# Setup module and store module-specific configuration in `CFG`
CFG = op.setup_module(
    name = "mrna_ashm_tally",
    version = "1.0",
    subdirectories = ["inputs", "tally", "outputs"]
)

# `runs` is popped from CFG by op.cleanup_module, so anything needed after that
# point is captured here
_SAMPLES = CFG["samples"]
_REGIONS = CFG["regions"]
_REGIONS_BED = CFG["regions_bed"]
_MERGE = CFG["options"]["merge"]

# One region file per genome build, so a build whose coordinates differ needs
# only its own entry in the config rather than a separate run
_BUILDS = sorted(set(_SAMPLES["genome_build"]))
_missing = [b for b in _BUILDS if b not in _REGIONS and b not in _REGIONS_BED]
assert not _missing, (
    "mrna_ashm_tally: no regions configured for genome build(s): "
    + ", ".join(_missing)
    + ". Set `regions` or `regions_bed` in the project config, keyed by genome build.")


# Define rules to be run locally when using a compute cluster
localrules:
    _mrna_ashm_tally_input_bam,
    _mrna_ashm_tally_regions,
    _mrna_ashm_tally_output_tsv,
    _mrna_ashm_tally_all,


wildcard_constraints:
    seq_type = r"[^/]+",
    genome_build = r"[^/]+",
    sample_id = r"[^/]+",


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mrna_ashm_tally_input_bam:
    input:
        bam = ancient(CFG["inputs"]["sample_bam"]),
        bai = ancient(CFG["inputs"]["sample_bai"])
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)


# Writes the regions for one genome build as a table the tally script reads,
# from the bed when one is configured for the build, else from `regions`.
# Overlapping and abutting intervals are merged and the table is sorted, so no
# base is reported twice and the per-sample output can be tabix-indexed.
def _mrna_ashm_tally_regions_bed(wildcards):
    bed = _REGIONS_BED.get(wildcards.genome_build)
    return [bed] if bed else []


rule _mrna_ashm_tally_regions:
    input:
        bed = _mrna_ashm_tally_regions_bed
    output:
        regions = CFG["dirs"]["inputs"] + "regions/{genome_build}.tsv"
    run:
        # 0-based half-open intervals
        ivs = []
        if input.bed:
            with open(input.bed[0]) as fh:
                for line in fh:
                    if not line.strip() or line.startswith(("#", "track", "browser")):
                        continue
                    f = line.rstrip("\n").split("\t")
                    name = f[3] if len(f) > 3 and f[3] else f"{f[0]}:{int(f[1]) + 1}-{f[2]}"
                    ivs.append((f[0], int(f[1]), int(f[2]), name))
        else:
            for name, coord in _REGIONS[wildcards.genome_build].items():
                chrom, span = coord.split(":")
                start, end = span.split("-")
                ivs.append((chrom, int(start) - 1, int(end), name))
        ivs.sort(key = lambda x: (x[0], x[1], x[2]))
        merged = []
        for chrom, start, end, name in ivs:
            if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
                m = merged[-1]
                m[2] = max(m[2], end)
                if name not in m[3]:
                    m[3].append(name)
            else:
                merged.append([chrom, start, end, [name]])
        pd.DataFrame(
            [{"region": ",".join(n), "chrom": c, "start": s + 1, "end": e}
             for c, s, e, n in merged]
        ).to_csv(output.regions, sep="\t", index=False)


rule _mrna_ashm_tally_run:
    input:
        bam = str(rules._mrna_ashm_tally_input_bam.output.bam),
        bai = str(rules._mrna_ashm_tally_input_bam.output.bai),
        regions = str(rules._mrna_ashm_tally_regions.output.regions),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        tally = CFG["dirs"]["tally"] + "{seq_type}--{genome_build}/{sample_id}.tally.tsv.gz",
        tbi = CFG["dirs"]["tally"] + "{seq_type}--{genome_build}/{sample_id}.tally.tsv.gz.tbi"
    log:
        stderr = CFG["logs"]["tally"] + "{seq_type}--{genome_build}/{sample_id}.tally.stderr.log"
    benchmark:
        CFG["dirs"]["tally"] + "{seq_type}--{genome_build}/{sample_id}.tally.benchmark.tsv"
    params:
        script = CFG["inputs"]["tally_script"],
        min_mapq = CFG["options"]["min_mapq"],
        min_baseq = CFG["options"]["min_baseq"],
        max_depth = CFG["options"]["max_depth"]
    conda:
        CFG["conda_envs"]["pysam"]
    threads:
        CFG["threads"]["tally"]
    resources:
        **CFG["resources"]["tally"]
    shell:
        op.as_one_line("""
        python {params.script}
        --bam {input.bam}
        --index {input.bai}
        --reference {input.fasta}
        --regions {input.regions}
        --sample {wildcards.sample_id}
        --out {output.tally}
        --min-mapq {params.min_mapq}
        --min-baseq {params.min_baseq}
        --max-depth {params.max_depth}
        2> {log.stderr}
        """)


# One table per seq_type and build, so a downstream consumer reads a single file
# rather than one per sample
def _mrna_ashm_tally_for_build(wildcards):
    keep = _SAMPLES[(_SAMPLES["seq_type"] == wildcards.seq_type) &
                    (_SAMPLES["genome_build"] == wildcards.genome_build)]
    return expand(str(rules._mrna_ashm_tally_run.output.tally),
                  seq_type = wildcards.seq_type,
                  genome_build = wildcards.genome_build,
                  sample_id = keep["sample_id"])


rule _mrna_ashm_tally_aggregate:
    input:
        tally = _mrna_ashm_tally_for_build
    output:
        tally = CFG["dirs"]["outputs"] + "merged/{seq_type}--{genome_build}.tally.tsv.gz"
    resources:
        **CFG["resources"]["aggregate"]
    shell:
        op.as_one_line("""
        gzip -dc {input.tally[0]} | head -n 1 | gzip > {output.tally}
        &&
        for f in {input.tally}; do gzip -dc $f | tail -n +2; done | gzip >> {output.tally}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _mrna_ashm_tally_output_tsv:
    input:
        tally = str(rules._mrna_ashm_tally_run.output.tally),
        tbi = str(rules._mrna_ashm_tally_run.output.tbi)
    output:
        tally = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.tally.tsv.gz",
        tbi = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.tally.tsv.gz.tbi"
    run:
        op.relative_symlink(input.tally, output.tally, in_module = True)
        op.relative_symlink(input.tbi, output.tbi, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _mrna_ashm_tally_all:
    input:
        expand(
            str(rules._mrna_ashm_tally_output_tsv.output.tally),
            zip,  # Run expand() with zip(), not product()
            seq_type = _SAMPLES["seq_type"],
            genome_build = _SAMPLES["genome_build"],
            sample_id = _SAMPLES["sample_id"]),
        expand(
            str(rules._mrna_ashm_tally_aggregate.output.tally),
            zip,
            seq_type = _SAMPLES.drop_duplicates(["seq_type", "genome_build"])["seq_type"],
            genome_build = _SAMPLES.drop_duplicates(["seq_type", "genome_build"])["genome_build"])
        if _MERGE else []


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
