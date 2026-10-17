#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Broad Institute (GATK4 somatic CNV: CreateReadCountPanelOfNormals, DenoiseReadCounts)
# Module Author:    Brett Collinge
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import re
import sys

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(f"oncopipe version installed: {current_version}")
    logger.warning(f"This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment")
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["gatk_tangent"]`
CFG = op.setup_module(
    name = "gatk_tangent",
    version = "1.0",
    subdirectories = ["inputs", "intervals", "counts", "pon", "denoised", "model", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _gatk_tangent_input_bam,
    _gatk_tangent_output_denoised,
    _gatk_tangent_output_seg,
    _gatk_tangent_output_pon,
    _gatk_tangent_all,


##### PANELS AND DENOISING RUNS #####

# Kept as module globals: op.cleanup_module removes `CFG` and its sample table.
_GT_SAMPLES = CFG["samples"][["sample_id", "patient_id", "seq_type", "genome_build"]].drop_duplicates()
_GT_BAM = CFG["inputs"]["sample_bam"]
_GT_INPUTS_DIR = CFG["dirs"]["inputs"]
_GT_VARIANTS = dict(CFG["options"]["denoise"])
_GT_MODEL = list(CFG["options"]["model_segments"] or [])
for _v in list(_GT_VARIANTS) + _GT_MODEL:
    if not re.fullmatch(r"[A-Za-z0-9_]+", _v):
        sys.exit(f"gatk_tangent: denoise variant {_v!r} may use letters, digits and '_' only.")
    if _v not in _GT_VARIANTS:
        sys.exit(f"gatk_tangent: model_segments names {_v!r}, which is not an options.denoise variant.")

# Pure-Java BLAS for the Spark-backed SVD: the cluster nodes' system BLAS lacks cblas_dspr, so netlib-native fails.
_GT_F2J_BLAS = (
    "-Dcom.github.fommil.netlib.BLAS=com.github.fommil.netlib.F2jBLAS "
    "-Dcom.github.fommil.netlib.LAPACK=com.github.fommil.netlib.F2jLAPACK "
    "-Dcom.github.fommil.netlib.ARPACK=com.github.fommil.netlib.F2jARPACK"
)


def _gatk_tangent_sample_row(sample_id, context):
    rows = _GT_SAMPLES[_GT_SAMPLES["sample_id"] == sample_id]
    if len(rows) != 1:
        sys.exit(
            f"gatk_tangent: {context} names {sample_id!r}, which appears {len(rows)} times in the "
            "samples table (after the pairing_config seq_type filter); it must appear exactly once."
        )
    return rows.iloc[0]


# name -> {"genome_build": str, "members": [(seq_type, genome_build, sample_id)], "patients": set}
_GT_PANELS = {}
for _name, _members in (CFG["panels"] or {}).items():
    if not re.fullmatch(r"[A-Za-z0-9_.]+", _name):
        sys.exit(f"gatk_tangent: panel name {_name!r} may use letters, digits, '_' and '.' only.")
    _members = list(dict.fromkeys(_members))
    _rows = [_gatk_tangent_sample_row(m, f"panel {_name!r}") for m in _members]
    _builds = sorted({r["genome_build"] for r in _rows})
    if len(_builds) != 1:
        sys.exit(f"gatk_tangent: panel {_name!r} spans genome builds {_builds}.")
    _GT_PANELS[_name] = {
        "genome_build": _builds[0],
        "members": [(r["seq_type"], r["genome_build"], r["sample_id"]) for r in _rows],
        "patients": {r["patient_id"] for r in _rows},
    }

# One row per (sample, panel) denoising run
_GT_RUNS = {"seq_type": [], "genome_build": [], "sample_id": [], "panel": []}
for _sid, _panels in (CFG["denoise"] or {}).items():
    _row = _gatk_tangent_sample_row(_sid, "denoise")
    for _panel in _panels:
        if _panel not in _GT_PANELS:
            sys.exit(f"gatk_tangent: denoise {_sid!r} names panel {_panel!r}, which is not defined.")
        if _GT_PANELS[_panel]["genome_build"] != _row["genome_build"]:
            sys.exit(f"gatk_tangent: {_sid!r} is {_row['genome_build']} but panel {_panel!r} is {_GT_PANELS[_panel]['genome_build']}.")
        if _row["patient_id"] in _GT_PANELS[_panel]["patients"]:
            sys.exit(f"gatk_tangent: panel {_panel!r} holds a sample of {_sid!r}'s patient; exclude it from the panel.")
        _GT_RUNS["seq_type"].append(_row["seq_type"])
        _GT_RUNS["genome_build"].append(_row["genome_build"])
        _GT_RUNS["sample_id"].append(_sid)
        _GT_RUNS["panel"].append(_panel)
_GT_RUN_PANELS = sorted({(b, p) for b, p in zip(_GT_RUNS["genome_build"], _GT_RUNS["panel"])})


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/'),
# keeping the .bam / .cram extension and exposing the index under both names
rule _gatk_tangent_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.{ext}",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.{ext}.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.{ext}.crai"
    wildcard_constraints:
        ext = "bam|cram"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


def _gatk_tangent_get_reads(wildcards):
    w = dict(wildcards.items())
    src = str(_GT_BAM(wildcards) if callable(_GT_BAM) else _GT_BAM.format(**w))
    if src.endswith(".cram"):
        ext = "cram"
    elif src.endswith(".bam"):
        ext = "bam"
    else:
        raise ValueError(f"gatk_tangent: {src} is neither .bam nor .cram")
    return f"{_GT_INPUTS_DIR}bam/{w['seq_type']}--{w['genome_build']}/{w['sample_id']}.{ext}"


def _gatk_tangent_get_index(wildcards):
    reads = _gatk_tangent_get_reads(wildcards)
    return reads + (".crai" if reads.endswith(".cram") else ".bai")


# Fixed-length bins over the configured contigs; bins made only of N are dropped by the tool
rule _gatk_tangent_preprocess_intervals:
    input:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
        dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output:
        interval_list = CFG["dirs"]["intervals"] + "{genome_build}/preprocessed.interval_list"
    log:
        stderr = CFG["logs"]["intervals"] + "{genome_build}/preprocess_intervals.stderr.log"
    params:
        contigs = " ".join(f"-L {c}" for c in CFG["options"]["contigs"]),
        bin_length = CFG["options"]["bin_length"],
        padding = CFG["options"]["padding"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["preprocess_intervals"]
    resources:
        **CFG["resources"]["preprocess_intervals"]
    shell:
        op.as_one_line("""
        gatk PreprocessIntervals
            -R {input.fasta}
            {params.contigs}
            --bin-length {params.bin_length}
            --padding {params.padding}
            --interval-merging-rule OVERLAPPING_ONLY
            -O {output.interval_list}
            2> {log.stderr}
        """)


# GC content per bin, used by CreateReadCountPanelOfNormals for explicit GC-bias correction
rule _gatk_tangent_annotate_intervals:
    input:
        interval_list = str(rules._gatk_tangent_preprocess_intervals.output.interval_list),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output:
        annotated = CFG["dirs"]["intervals"] + "{genome_build}/annotated.tsv"
    log:
        stderr = CFG["logs"]["intervals"] + "{genome_build}/annotate_intervals.stderr.log"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["annotate_intervals"]
    resources:
        **CFG["resources"]["annotate_intervals"]
    shell:
        op.as_one_line("""
        gatk AnnotateIntervals
            -R {input.fasta}
            -L {input.interval_list}
            --interval-merging-rule OVERLAPPING_ONLY
            -O {output.annotated}
            2> {log.stderr}
        """)


# Read counts per bin; CRAM is decoded with the genome_build fasta
rule _gatk_tangent_collect_read_counts:
    input:
        bam = _gatk_tangent_get_reads,
        bai = _gatk_tangent_get_index,
        interval_list = str(rules._gatk_tangent_preprocess_intervals.output.interval_list),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        counts = CFG["dirs"]["counts"] + "{seq_type}--{genome_build}/{sample_id}.counts.hdf5"
    log:
        stderr = CFG["logs"]["counts"] + "{seq_type}--{genome_build}/{sample_id}.collect_read_counts.stderr.log"
    params:
        java_mem = lambda w, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["collect_read_counts"]
    resources:
        **CFG["resources"]["collect_read_counts"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}m" CollectReadCounts
            -I {input.bam}
            -R {input.fasta}
            -L {input.interval_list}
            --interval-merging-rule OVERLAPPING_ONLY
            --format HDF5
            -O {output.counts}
            2> {log.stderr}
        """)


def _gatk_tangent_get_panel_counts(wildcards):
    return [
        str(rules._gatk_tangent_collect_read_counts.output.counts).format(seq_type=s, genome_build=b, sample_id=i)
        for (s, b, i) in _GT_PANELS[wildcards.panel]["members"]
    ]


# Builds one panel of normals (SVD of the standardized panel counts)
rule _gatk_tangent_create_pon:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        counts = _gatk_tangent_get_panel_counts,
        annotated = str(rules._gatk_tangent_annotate_intervals.output.annotated)
    output:
        pon = CFG["dirs"]["pon"] + "{genome_build}/{panel}.pon.hdf5"
    log:
        stderr = CFG["logs"]["pon"] + "{genome_build}/{panel}.create_pon.stderr.log"
    params:
        counts = lambda w, input: " ".join(f"-I {p}" for p in input.counts),
        opts = CFG["options"]["create_pon"],
        java_mem = lambda w, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["create_pon"]
    resources:
        **CFG["resources"]["create_pon"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}m %s" CreateReadCountPanelOfNormals
            {params.counts}
            --annotated-intervals {input.annotated}
            {params.opts}
            -O {output.pon}
            2> {log.stderr}
        """ % _GT_F2J_BLAS)


# Denoises one sample against one panel, once per options.denoise variant
rule _gatk_tangent_denoise:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+",
        variant = r"[A-Za-z0-9_]+"
    input:
        counts = str(rules._gatk_tangent_collect_read_counts.output.counts),
        pon = str(rules._gatk_tangent_create_pon.output.pon)
    output:
        std_cr = CFG["dirs"]["denoised"] + "{variant}/{seq_type}--{genome_build}/{sample_id}--{panel}.standardizedCR.tsv",
        denoised_cr = CFG["dirs"]["denoised"] + "{variant}/{seq_type}--{genome_build}/{sample_id}--{panel}.denoisedCR.tsv"
    log:
        stderr = CFG["logs"]["denoised"] + "{variant}/{seq_type}--{genome_build}/{sample_id}--{panel}.denoise.stderr.log"
    params:
        opts = lambda w: _GT_VARIANTS[w.variant],
        java_mem = lambda w, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["denoise"]
    resources:
        **CFG["resources"]["denoise"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}m %s" DenoiseReadCounts
            -I {input.counts}
            --count-panel-of-normals {input.pon}
            {params.opts}
            --standardized-copy-ratios {output.std_cr}
            --denoised-copy-ratios {output.denoised_cr}
            2> {log.stderr}
        """ % _GT_F2J_BLAS)


# Segments the denoised copy ratios alone (no allelic counts)
rule _gatk_tangent_model_segments:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+",
        variant = r"[A-Za-z0-9_]+"
    input:
        denoised_cr = str(rules._gatk_tangent_denoise.output.denoised_cr)
    output:
        cr_seg = CFG["dirs"]["model"] + "{variant}/{seq_type}--{genome_build}/{sample_id}--{panel}/{sample_id}--{panel}.cr.seg"
    log:
        stderr = CFG["logs"]["model"] + "{variant}/{seq_type}--{genome_build}/{sample_id}--{panel}.model_segments.stderr.log"
    params:
        out_dir = CFG["dirs"]["model"] + "{variant}/{seq_type}--{genome_build}/{sample_id}--{panel}",
        opts = CFG["options"]["model_segments_opts"],
        java_mem = lambda w, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["model_segments"]
    resources:
        **CFG["resources"]["model_segments"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}m %s" ModelSegments
            --denoised-copy-ratios {input.denoised_cr}
            {params.opts}
            --output {params.out_dir}
            --output-prefix {wildcards.sample_id}--{wildcards.panel}
            2> {log.stderr}
        """ % _GT_F2J_BLAS)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _gatk_tangent_output_denoised:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+",
        variant = r"[A-Za-z0-9_]+"
    input:
        denoised_cr = str(rules._gatk_tangent_denoise.output.denoised_cr)
    output:
        denoised_cr = CFG["dirs"]["outputs"] + "denoisedCR/{variant}/{seq_type}--{genome_build}/{sample_id}--{panel}.denoisedCR.tsv"
    run:
        op.relative_symlink(input.denoised_cr, output.denoised_cr, in_module=True)


rule _gatk_tangent_output_seg:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+",
        variant = r"[A-Za-z0-9_]+"
    input:
        cr_seg = str(rules._gatk_tangent_model_segments.output.cr_seg)
    output:
        cr_seg = CFG["dirs"]["outputs"] + "cr_seg/{variant}/{seq_type}--{genome_build}/{sample_id}--{panel}.cr.seg"
    run:
        op.relative_symlink(input.cr_seg, output.cr_seg, in_module=True)


rule _gatk_tangent_output_pon:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        pon = str(rules._gatk_tangent_create_pon.output.pon)
    output:
        pon = CFG["dirs"]["outputs"] + "pon/{genome_build}/{panel}.pon.hdf5"
    run:
        op.relative_symlink(input.pon, output.pon, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _gatk_tangent_all:
    input:
        [
            f
            for variant in _GT_VARIANTS
            for f in expand(
                [str(rules._gatk_tangent_output_denoised.output.denoised_cr)]
                + ([str(rules._gatk_tangent_output_seg.output.cr_seg)] if variant in _GT_MODEL else []),
                zip,
                seq_type = _GT_RUNS["seq_type"],
                genome_build = _GT_RUNS["genome_build"],
                sample_id = _GT_RUNS["sample_id"],
                panel = _GT_RUNS["panel"],
                variant = [variant] * len(_GT_RUNS["sample_id"])
            )
        ],
        expand(
            str(rules._gatk_tangent_output_pon.output.pon),
            zip,
            genome_build = [b for b, p in _GT_RUN_PANELS],
            panel = [p for b, p in _GT_RUN_PANELS]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
