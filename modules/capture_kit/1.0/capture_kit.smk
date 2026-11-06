#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Brett Collinge
# Module Author:    Brett Collinge
# Contributors:     N/A


##### SETUP #####

import oncopipe as op

# `CFG` is a shortcut to `config["lcr-modules"]["capture_kit"]`
CFG = op.setup_module(
    name = "capture_kit",
    version = "1.0",
    subdirectories = ["inputs", "kits", "panel", "bedcov", "call", "outputs"],
)

localrules:
    _capture_kit_input_bam,
    _capture_kit_input_references,
    _capture_kit_kits_tsv,
    _capture_kit_merge,
    _capture_kit_output,
    _capture_kit_output_merged,
    _capture_kit_all


def _capture_kit_build_key(genome_build):
    return "hg38" if "38" in str(genome_build) else "grch37"


##### RULES #####


wildcard_constraints:
    seq_type = "[^/]+",
    genome_build = "[^/]+",
    sample_id = "[^/]+"


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _capture_kit_input_bam:
    input:
        bam = ancient(CFG["inputs"]["sample_bam"]),
        bai = ancient(CFG["inputs"]["sample_bai"])
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# The reference is needed to decode CRAM input
rule _capture_kit_input_references:
    input:
        genome_fa = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        genome_fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        genome_fa = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.fa",
        genome_fai = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.fa.fai"
    run:
        op.absolute_symlink(input.genome_fa, output.genome_fa)
        op.absolute_symlink(input.genome_fai, output.genome_fai)


# Kit table (kit, build, path) for one build, from the `kits` config
rule _capture_kit_kits_tsv:
    output:
        tsv = CFG["dirs"]["kits"] + "{genome_build}/kits.tsv"
    params:
        build = lambda w: _capture_kit_build_key(w.genome_build)
    run:
        kits = config["lcr-modules"]["capture_kit"]["kits"]
        rows = [(k, params.build, b[params.build]) for k, b in kits.items() if params.build in b]
        if len(rows) < 2:
            raise ValueError(f"capture_kit: fewer than two kits configured for {params.build}")
        with open(output.tsv, "w") as f:
            f.write("kit\tbuild\tpath\n")
            for r in rows:
                f.write("\t".join(r) + "\n")


# Pairwise kit-exclusive sets, the shared core and the off-target floor, one BED4
rule _capture_kit_panel:
    input:
        tsv = str(rules._capture_kit_kits_tsv.output.tsv)
    output:
        panel = CFG["dirs"]["panel"] + "{genome_build}/panel.bed",
        sizes = CFG["dirs"]["panel"] + "{genome_build}/panel.sizes.tsv"
    log:
        stderr = CFG["logs"]["panel"] + "{genome_build}/build_panel.stderr.log"
    params:
        script = CFG["inputs"]["build_panel"],
        build = lambda w: _capture_kit_build_key(w.genome_build),
        annotation = CFG["panel"]["annotation_dir"],
        outdir = CFG["dirs"]["panel"] + "{genome_build}",
        gap = CFG["panel"]["gap"],
        min_region = CFG["panel"]["min_region"],
        panel_size = CFG["panel"]["panel_size"],
        floor_slop = CFG["panel"]["floor_slop"],
        floor_window = CFG["panel"]["floor_window"]
    conda:
        CFG["conda_envs"]["capture_kit"]
    threads:
        CFG["threads"]["panel"]
    resources:
        **CFG["resources"]["panel"]
    shell:
        op.as_one_line("""
        GAP={params.gap} MIN_REGION={params.min_region} PANEL_SIZE={params.panel_size}
        FLOOR_SLOP={params.floor_slop} FLOOR_WIN={params.floor_window}
        bash {params.script} {input.tsv} {params.build} {params.annotation} {params.outdir}
        2> {log.stderr}
        """)


# One samtools bedcov pass over the panel
rule _capture_kit_bedcov:
    input:
        bam = ancient(str(rules._capture_kit_input_bam.output.bam)),
        panel = str(rules._capture_kit_panel.output.panel),
        fasta = str(rules._capture_kit_input_references.output.genome_fa)
    output:
        bedcov = CFG["dirs"]["bedcov"] + "{seq_type}--{genome_build}/{sample_id}.tsv"
    log:
        stderr = CFG["logs"]["bedcov"] + "{seq_type}--{genome_build}/{sample_id}.stderr.log"
    params:
        script = CFG["inputs"]["fingerprint"]
    conda:
        CFG["conda_envs"]["capture_kit"]
    threads:
        CFG["threads"]["bedcov"]
    resources:
        **CFG["resources"]["bedcov"]
    shell:
        "bash {params.script} {input.bam} {input.panel} {input.fasta} {output.bedcov} 2> {log.stderr}"


# Per-set depth over the core, per-kit own/against/margin, and the label
rule _capture_kit_call:
    input:
        bedcov = str(rules._capture_kit_bedcov.output.bedcov)
    output:
        sets = CFG["dirs"]["call"] + "{seq_type}--{genome_build}/{sample_id}.capture_kit.tsv",
        call = CFG["dirs"]["call"] + "{seq_type}--{genome_build}/{sample_id}.capture_kit.call.tsv"
    log:
        stderr = CFG["logs"]["call"] + "{seq_type}--{genome_build}/{sample_id}.stderr.log"
    params:
        script = CFG["inputs"]["call_kit"],
        prefix = CFG["dirs"]["call"] + "{seq_type}--{genome_build}/{sample_id}",
        opts = op.as_one_line("""
            --hi {hi} --lo {lo} --min-margin {min_margin} --min-core-depth {min_core_depth}
            --min-set-bp {min_set_bp} --genome-floor {genome_floor}
            """.format(**CFG["call"]))
    conda:
        CFG["conda_envs"]["capture_kit"]
    threads:
        CFG["threads"]["call"]
    resources:
        **CFG["resources"]["call"]
    shell:
        "python3 {params.script} {input.bedcov} {wildcards.sample_id} {params.prefix} {params.opts} 2> {log.stderr}"


def _capture_kit_merge_inputs(wildcards):
    CFG = config["lcr-modules"]["capture_kit"]
    s = CFG["samples"]
    s = s[(s["seq_type"] == wildcards.seq_type) & (s["genome_build"] == wildcards.genome_build)]
    return expand(
        str(rules._capture_kit_call.output.call),
        zip,
        seq_type=s["seq_type"],
        genome_build=s["genome_build"],
        sample_id=s["sample_id"])


# One call table per seq_type and build
rule _capture_kit_merge:
    input:
        _capture_kit_merge_inputs
    output:
        CFG["dirs"]["call"] + "{seq_type}--{genome_build}.capture_kit.tsv"
    shell:
        op.as_one_line("""
        head -n 1 {input[0]} > {output} &&
        tail -n +2 -q {input} >> {output}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _capture_kit_output:
    input:
        sets = str(rules._capture_kit_call.output.sets),
        call = str(rules._capture_kit_call.output.call)
    output:
        sets = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.capture_kit_sets.tsv",
        call = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.capture_kit.tsv"
    run:
        op.relative_symlink(input.sets, output.sets, in_module = True)
        op.relative_symlink(input.call, output.call, in_module = True)


rule _capture_kit_output_merged:
    input:
        str(rules._capture_kit_merge.output)
    output:
        CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}.capture_kit.tsv"
    run:
        op.relative_symlink(input, output, in_module = True)


rule _capture_kit_all:
    input:
        expand(
            [
                str(rules._capture_kit_output.output.sets),
                str(rules._capture_kit_output.output.call)
            ],
            zip,
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"]),
        expand(
            str(rules._capture_kit_output_merged.output),
            zip,
            **CFG["samples"][["seq_type", "genome_build"]].drop_duplicates().to_dict("list"))


##### CLEANUP #####


op.cleanup_module(CFG)
