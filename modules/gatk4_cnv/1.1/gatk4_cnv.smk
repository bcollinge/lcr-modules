#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Broad Institute (GATK4 somatic CNV: CreateReadCountPanelOfNormals, DenoiseReadCounts);
#                   Olshen et al. 2004 (CBS, DNAcopy); CRUK CI (Rascal)
# Module Author:    Brett Collinge
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import pandas as pd
import re
import shlex
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
# `CFG` is a shortcut to `config["lcr-modules"]["gatk4_cnv"]`
CFG = op.setup_module(
    name = "gatk4_cnv",
    version = "1.1",
    subdirectories = ["inputs", "intervals", "counts", "pon", "denoised", "cbs", "fit", "corrected", "projected", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _gatk4_cnv_input_bam,
    _gatk4_cnv_capture_bed,
    _gatk4_cnv_output_seg,
    _gatk4_cnv_output_ratio,
    _gatk4_cnv_output_purity,
    _gatk4_cnv_output_plot,
    _gatk4_cnv_output_pon,
    _gatk4_cnv_output_projected,
    _gatk4_cnv_all,


##### PANELS AND RUNS #####

# Kept as module globals: op.cleanup_module removes `CFG` and its sample table.
_GC_COLUMNS = ["sample_id", "patient_id", "seq_type", "genome_build"]
if "capture_space" in CFG["samples"].columns:
    _GC_COLUMNS.append("capture_space")
_GC_SAMPLES = CFG["samples"][_GC_COLUMNS].drop_duplicates()
_GC_BAM = CFG["inputs"]["sample_bam"]
_GC_INPUTS_DIR = CFG["dirs"]["inputs"]
_GC_CONTIGS_CFG = CFG["options"]["contigs"]
_GC_CBS = CFG["options"]["cbs"]
_GC_RASCAL = CFG["options"]["rascal"]
_GC_CAPTURE = CFG["options"]["capture"]
_GC_MODES = CFG["switches"]["mode"]
_GC_BIN_LENGTH = CFG["options"]["bin_length"]
_GC_FIT_SCRIPTS = {"genome": CFG["inputs"]["fit_script"], "capture": CFG["inputs"]["fit_capture_script"]}

# Pure-Java BLAS for the Spark-backed SVD: the cluster nodes' system BLAS lacks cblas_dspr, so netlib-native fails.
_GC_F2J_BLAS = (
    "-Dcom.github.fommil.netlib.BLAS=com.github.fommil.netlib.F2jBLAS "
    "-Dcom.github.fommil.netlib.LAPACK=com.github.fommil.netlib.F2jLAPACK "
    "-Dcom.github.fommil.netlib.ARPACK=com.github.fommil.netlib.F2jARPACK"
)


# contigs: one list for every build, or a dict keyed by genome_build
def _gatk4_cnv_contigs(genome_build):
    c = _GC_CONTIGS_CFG.get(genome_build) if isinstance(_GC_CONTIGS_CFG, dict) else _GC_CONTIGS_CFG
    if not c:
        sys.exit(f"gatk4_cnv: options.contigs has no entry for genome_build {genome_build!r}.")
    return [str(x) for x in c]


# "genome" (fixed-length bins) or "capture" (capture-space targets), from switches.mode keyed by seq_type
def _gatk4_cnv_mode(seq_type):
    m = _GC_MODES.get(seq_type, _GC_MODES["_default"])
    if m not in ("genome", "capture"):
        sys.exit(f"gatk4_cnv: switches.mode gives {m!r} for seq_type {seq_type!r}; use 'genome' or 'capture'.")
    return m


# A capture-mode sample must name its capture space; there is no default
def _gatk4_cnv_capture_space(row, context):
    cs = row["capture_space"] if "capture_space" in row.index else None
    if cs is None or pd.isna(cs) or str(cs).strip().lower() in ("", "none", "na", "nan"):
        sys.exit(f"gatk4_cnv: {context} {row['sample_id']!r} is in capture mode but has no capture_space.")
    return str(cs)


def _gatk4_cnv_sample_row(sample_id, context):
    rows = _GC_SAMPLES[_GC_SAMPLES["sample_id"] == sample_id]
    if len(rows) != 1:
        sys.exit(
            f"gatk4_cnv: {context} names {sample_id!r}, which appears {len(rows)} times in the "
            "samples table (after the pairing_config seq_type filter); it must appear exactly once."
        )
    return rows.iloc[0]


# Capture space of each capture-mode sample named in a panel or a denoise: (seq_type, genome_build, sample_id) -> name
_GC_SPACE = {}

# Declared panels: name -> {"genome_build", "mode", "capture_space", "members": [(seq_type, genome_build, sample_id)]}
_GC_DECLARED = {}
for _name, _members in (CFG["panels"] or {}).items():
    if not re.fullmatch(r"[A-Za-z0-9_.]+", _name):
        sys.exit(f"gatk4_cnv: panel name {_name!r} may use letters, digits, '_' and '.' only.")
    _rows = [_gatk4_cnv_sample_row(m, f"panel {_name!r}") for m in dict.fromkeys(_members)]
    _builds = sorted({r["genome_build"] for r in _rows})
    if len(_builds) != 1:
        sys.exit(f"gatk4_cnv: panel {_name!r} spans genome builds {_builds}.")
    _modes = sorted({_gatk4_cnv_mode(r["seq_type"]) for r in _rows})
    if len(_modes) != 1:
        sys.exit(f"gatk4_cnv: panel {_name!r} mixes modes {_modes}.")
    _space = None
    if _modes[0] == "capture":
        for r in _rows:
            _GC_SPACE[(r["seq_type"], r["genome_build"], r["sample_id"])] = _gatk4_cnv_capture_space(r, f"panel {_name!r} member")
        _spaces = sorted({_GC_SPACE[(r["seq_type"], r["genome_build"], r["sample_id"])] for r in _rows})
        if len(_spaces) != 1:
            sys.exit(f"gatk4_cnv: panel {_name!r} spans capture spaces {_spaces}.")
        _space = _spaces[0]
    _GC_DECLARED[_name] = {
        "genome_build": _builds[0],
        "mode": _modes[0],
        "capture_space": _space,
        "members": [(r["seq_type"], r["genome_build"], r["sample_id"]) for r in _rows],
    }

# Panels built: every declared panel a sample is denoised against, with all its members.
_GC_PANELS = {}
_GC_RUNS = {"seq_type": [], "genome_build": [], "sample_id": [], "panel": []}
for _sid, _panel in (CFG["denoise"] or {}).items():
    _row = _gatk4_cnv_sample_row(_sid, "denoise")
    if _panel not in _GC_DECLARED:
        sys.exit(f"gatk4_cnv: denoise {_sid!r} names panel {_panel!r}, which is not defined.")
    _decl = _GC_DECLARED[_panel]
    if _decl["genome_build"] != _row["genome_build"]:
        sys.exit(f"gatk4_cnv: {_sid!r} is {_row['genome_build']} but panel {_panel!r} is {_decl['genome_build']}.")
    _mode = _gatk4_cnv_mode(_row["seq_type"])
    if _mode != _decl["mode"]:
        sys.exit(f"gatk4_cnv: {_sid!r} is in {_mode} mode but panel {_panel!r} is {_decl['mode']}.")
    if _mode == "capture":
        _space = _gatk4_cnv_capture_space(_row, "denoise")
        if _space != _decl["capture_space"]:
            sys.exit(f"gatk4_cnv: {_sid!r} is capture space {_space!r} but panel {_panel!r} is {_decl['capture_space']!r}.")
        _GC_SPACE[(_row["seq_type"], _row["genome_build"], _sid)] = _space
    _GC_PANELS[_panel] = _decl
    _GC_RUNS["seq_type"].append(_row["seq_type"])
    _GC_RUNS["genome_build"].append(_row["genome_build"])
    _GC_RUNS["sample_id"].append(_sid)
    _GC_RUNS["panel"].append(_panel)
_GC_RUN_PANELS = sorted({(b, p) for b, p in zip(_GC_RUNS["genome_build"], _GC_RUNS["panel"])})
_GC_CAPTURE_RUNS = {k: [v for v, s in zip(_GC_RUNS[k], _GC_RUNS["seq_type"]) if _gatk4_cnv_mode(s) == "capture"] for k in _GC_RUNS}


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/'),
# keeping the .bam / .cram extension and exposing the index under both names
rule _gatk4_cnv_input_bam:
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


def _gatk4_cnv_get_reads(wildcards):
    w = dict(wildcards.items())
    src = str(_GC_BAM(wildcards) if callable(_GC_BAM) else _GC_BAM.format(**w))
    if src.endswith(".cram"):
        ext = "cram"
    elif src.endswith(".bam"):
        ext = "bam"
    else:
        raise ValueError(f"gatk4_cnv: {src} is neither .bam nor .cram")
    return f"{_GC_INPUTS_DIR}bam/{w['seq_type']}--{w['genome_build']}/{w['sample_id']}.{ext}"


def _gatk4_cnv_get_index(wildcards):
    reads = _gatk4_cnv_get_reads(wildcards)
    return reads + (".crai" if reads.endswith(".cram") else ".bai")


# Fixed-length bins over the configured contigs; bins made only of N are dropped by the tool
rule _gatk4_cnv_preprocess_intervals:
    input:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
        dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output:
        interval_list = CFG["dirs"]["intervals"] + "{genome_build}/preprocessed.interval_list"
    log:
        stderr = CFG["logs"]["intervals"] + "{genome_build}/preprocess_intervals.stderr.log"
    params:
        contigs = lambda w: " ".join(f"-L {c}" for c in _gatk4_cnv_contigs(w.genome_build)),
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
rule _gatk4_cnv_annotate_intervals:
    input:
        interval_list = str(rules._gatk4_cnv_preprocess_intervals.output.interval_list),
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


# Capture space BED (unpadded, cols 1-3) restricted to the configured contigs and sorted in their order
rule _gatk4_cnv_capture_bed:
    input:
        bed = reference_files("genomes/{genome_build}/capture_space/{capture_space}.bed")
    output:
        bed = CFG["dirs"]["intervals"] + "{genome_build}/capture/{capture_space}/targets.bed"
    params:
        contigs = lambda w: ",".join(_gatk4_cnv_contigs(w.genome_build))
    shell:
        op.as_one_line("""
        awk -F '\\t' -v OFS='\\t' -v c={params.contigs:q}
            'BEGIN {{ n = split(c, a, ","); for (i = 1; i <= n; i++) o[a[i]] = i }} ($1 in o) {{ print o[$1], $1, $2, $3 }}'
            {input.bed} | sort -k1,1n -k3,3n -k4,4n | cut -f 2-4 > {output.bed}
            && test -s {output.bed}
        """)


# One interval per target (bin length 0); overlapping targets are merged
rule _gatk4_cnv_preprocess_targets:
    input:
        bed = str(rules._gatk4_cnv_capture_bed.output.bed),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
        dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output:
        interval_list = CFG["dirs"]["intervals"] + "{genome_build}/capture/{capture_space}/targets.interval_list"
    log:
        stderr = CFG["logs"]["intervals"] + "{genome_build}/capture/{capture_space}/preprocess_targets.stderr.log"
    params:
        bin_length = _GC_CAPTURE["bin_length"],
        padding = _GC_CAPTURE["padding"]
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
            -L {input.bed}
            --bin-length {params.bin_length}
            --padding {params.padding}
            --interval-merging-rule OVERLAPPING_ONLY
            -O {output.interval_list}
            2> {log.stderr}
        """)


rule _gatk4_cnv_annotate_targets:
    input:
        interval_list = str(rules._gatk4_cnv_preprocess_targets.output.interval_list),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output:
        annotated = CFG["dirs"]["intervals"] + "{genome_build}/capture/{capture_space}/targets.annotated.tsv"
    log:
        stderr = CFG["logs"]["intervals"] + "{genome_build}/capture/{capture_space}/annotate_targets.stderr.log"
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


# Fixed grid over each contig (bins start at 0) with the number of targets whose midpoint falls in each bin
rule _gatk4_cnv_target_grid:
    input:
        interval_list = str(rules._gatk4_cnv_preprocess_targets.output.interval_list),
        script = CFG["inputs"]["grid_script"]
    output:
        grid = CFG["dirs"]["intervals"] + "{genome_build}/capture/{capture_space}/targets.grid.tsv"
    log:
        stderr = CFG["logs"]["intervals"] + "{genome_build}/capture/{capture_space}/target_grid.stderr.log"
    params:
        contigs = lambda w: ",".join(_gatk4_cnv_contigs(w.genome_build)),
        grid_bin = _GC_CAPTURE["grid_bin"]
    conda:
        CFG["conda_envs"]["rascal"]
    threads:
        CFG["threads"]["target_grid"]
    resources:
        **CFG["resources"]["target_grid"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {input.script}
            {input.interval_list} {output.grid} {params.contigs:q} {params.grid_bin}
            2> {log.stderr}
        """)


def _gatk4_cnv_space(wildcards):
    return _GC_SPACE[(wildcards.seq_type, wildcards.genome_build, wildcards.sample_id)]


# Genome mode: the build's bins; capture mode: the sample's capture space targets
def _gatk4_cnv_intervals(wildcards):
    if _gatk4_cnv_mode(wildcards.seq_type) == "capture":
        return str(rules._gatk4_cnv_preprocess_targets.output.interval_list).format(
            genome_build=wildcards.genome_build, capture_space=_gatk4_cnv_space(wildcards))
    return str(rules._gatk4_cnv_preprocess_intervals.output.interval_list).format(genome_build=wildcards.genome_build)


def _gatk4_cnv_annotated(wildcards):
    decl = _GC_DECLARED[wildcards.panel]
    if decl["mode"] == "capture":
        return str(rules._gatk4_cnv_annotate_targets.output.annotated).format(
            genome_build=wildcards.genome_build, capture_space=decl["capture_space"])
    return str(rules._gatk4_cnv_annotate_intervals.output.annotated).format(genome_build=wildcards.genome_build)


def _gatk4_cnv_grid(wildcards):
    return str(rules._gatk4_cnv_target_grid.output.grid).format(
        genome_build=wildcards.genome_build, capture_space=_gatk4_cnv_space(wildcards))


# Read counts per bin with the tool's default read filters (duplicates and MAPQ < 30 removed);
# CRAM is decoded with the genome_build fasta
rule _gatk4_cnv_collect_read_counts:
    input:
        bam = _gatk4_cnv_get_reads,
        bai = _gatk4_cnv_get_index,
        interval_list = _gatk4_cnv_intervals,
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


def _gatk4_cnv_get_panel_counts(wildcards):
    return [
        str(rules._gatk4_cnv_collect_read_counts.output.counts).format(seq_type=s, genome_build=b, sample_id=i)
        for (s, b, i) in _GC_PANELS[wildcards.panel]["members"]
    ]


# Builds one panel of normals (SVD of the standardized panel counts)
rule _gatk4_cnv_create_pon:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        counts = _gatk4_cnv_get_panel_counts,
        annotated = _gatk4_cnv_annotated
    output:
        pon = CFG["dirs"]["pon"] + "{genome_build}/{panel}.pon.hdf5"
    log:
        stderr = CFG["logs"]["pon"] + "{genome_build}/{panel}.create_pon.stderr.log"
    params:
        counts = lambda w, input: " ".join(f"-I {p}" for p in input.counts),
        percentile = CFG["options"]["pon_percentile"],
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
            --minimum-interval-median-percentile {params.percentile}
            {params.opts}
            -O {output.pon}
            2> {log.stderr}
        """ % _GC_F2J_BLAS)


# Denoises one sample against its panel
rule _gatk4_cnv_denoise:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        counts = str(rules._gatk4_cnv_collect_read_counts.output.counts),
        pon = str(rules._gatk4_cnv_create_pon.output.pon)
    output:
        std_cr = CFG["dirs"]["denoised"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.standardizedCR.tsv",
        denoised_cr = CFG["dirs"]["denoised"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.denoisedCR.tsv"
    log:
        stderr = CFG["logs"]["denoised"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.denoise.stderr.log"
    params:
        opts = CFG["options"]["denoise"],
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
        """ % _GC_F2J_BLAS)


# CBS of the floored log2 denoised copy ratio (seeded per run); per bin: Ratio (denoised copy ratio),
# MedianRatio (median copy ratio of the bin's segment), CopyNumber (round half up of 2 x MedianRatio)
rule _gatk4_cnv_cbs:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        denoised_cr = str(rules._gatk4_cnv_denoise.output.denoised_cr),
        script = CFG["inputs"]["cbs_script"]
    output:
        ratio = CFG["dirs"]["cbs"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.ratio.tsv"
    log:
        stderr = CFG["logs"]["cbs"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.cbs.stderr.log"
    params:
        contigs = lambda w: ",".join(_gatk4_cnv_contigs(w.genome_build)),
        seed_key = lambda w: _GC_CBS["seed_key"].format(**dict(w.items())),
        seed_base = _GC_CBS["seed_base"],
        floor = _GC_CBS["floor"],
        alpha = _GC_CBS["alpha"],
        undo_sd = _GC_CBS["undo_sd"]
    conda:
        CFG["conda_envs"]["rascal"]
    threads:
        CFG["threads"]["cbs"]
    resources:
        **CFG["resources"]["cbs"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {input.script}
            {input.denoised_cr} {output.ratio} {params.contigs:q} {params.seed_key:q} {params.seed_base}
            {params.floor} {params.alpha} {params.undo_sd}
            2> {log.stderr}
        """)


# Rascal fit on the profile scaled to its median MedianRatio; the solution of record is the lowest
# distance x (1 + |ploidy - 2|)^penploidy / cellularity^penalty (ties to the higher cellularity)
rule _gatk4_cnv_fit:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        ratio = str(rules._gatk4_cnv_cbs.output.ratio),
        script = lambda w: _GC_FIT_SCRIPTS[_gatk4_cnv_mode(w.seq_type)]
    output:
        purity = CFG["dirs"]["fit"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.purity.tsv",
        solutions = CFG["dirs"]["fit"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.solutions.tsv"
    log:
        stderr = CFG["logs"]["fit"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.fit.stderr.log"
    params:
        run = "{sample_id}--{panel}",
        contigs = lambda w: ",".join(_gatk4_cnv_contigs(w.genome_build)),
        # genome: bin length (each bin weighs bin length bp); capture: options.capture.segment_weight
        bin_length = lambda w: _GC_CAPTURE["segment_weight"] if _gatk4_cnv_mode(w.seq_type) == "capture" else _GC_BIN_LENGTH,
        distance = _GC_RASCAL["distance_function"],
        penploidy = _GC_RASCAL["penploidy"],
        penalty = _GC_RASCAL["penalty"],
        amb_distance = _GC_RASCAL["ambiguous_distance"],
        amb_purity = _GC_RASCAL["ambiguous_purity"]
    conda:
        CFG["conda_envs"]["rascal"]
    threads:
        CFG["threads"]["fit"]
    resources:
        **CFG["resources"]["fit"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {input.script}
            {input.ratio} {output.purity} {output.solutions} {params.run:q} {params.contigs:q} {params.bin_length}
            {params.distance} {params.penploidy} {params.penalty} {params.amb_distance} {params.amb_purity}
            2> {log.stderr}
        """)


# Purity correction r' = max((r / m - (1 - p)) / p, 0), m the median MedianRatio, p the fitted purity
# (no fit: the profile scaled to m, not purity corrected); the lcr seg (log2 of the corrected segment
# value, floored as for CBS) and the profile plot
rule _gatk4_cnv_correct:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        ratio = str(rules._gatk4_cnv_cbs.output.ratio),
        purity = str(rules._gatk4_cnv_fit.output.purity),
        script = CFG["inputs"]["correct_script"]
    output:
        ratio = CFG["dirs"]["corrected"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.ratio.tsv",
        seg = CFG["dirs"]["corrected"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.seg",
        plot = CFG["dirs"]["corrected"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.png"
    log:
        stderr = CFG["logs"]["corrected"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.correct.stderr.log"
    params:
        seg_id = "{sample_id}",
        contigs = lambda w: ",".join(_gatk4_cnv_contigs(w.genome_build)),
        bin_length = lambda w: _GC_CAPTURE["bin_length"] if _gatk4_cnv_mode(w.seq_type) == "capture" else _GC_BIN_LENGTH,
        floor = _GC_CBS["floor"]
    conda:
        CFG["conda_envs"]["rascal"]
    threads:
        CFG["threads"]["correct"]
    resources:
        **CFG["resources"]["correct"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {input.script}
            {input.ratio} {input.purity} {output.ratio} {output.seg} {output.plot} {params.seg_id:q}
            {params.contigs:q} {params.bin_length} {params.floor}
            2> {log.stderr}
        """)


# Capture mode: the corrected profile on the target grid. Per bin: n_called (targets with a denoised ratio, by midpoint),
# Ratio (their median), MedianRatio (the segment holding most of them; ties to the first), CopyNumber; NA below
# min_targets_per_bin
rule _gatk4_cnv_project:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        ratio = str(rules._gatk4_cnv_correct.output.ratio),
        grid = _gatk4_cnv_grid,
        script = CFG["inputs"]["project_script"]
    output:
        grid = CFG["dirs"]["projected"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.grid.tsv"
    log:
        stderr = CFG["logs"]["projected"] + "{seq_type}--{genome_build}/{sample_id}--{panel}.project.stderr.log"
    params:
        min_targets = _GC_CAPTURE["min_targets_per_bin"]
    conda:
        CFG["conda_envs"]["rascal"]
    threads:
        CFG["threads"]["project"]
    resources:
        **CFG["resources"]["project"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {input.script}
            {input.ratio} {input.grid} {output.grid} {params.min_targets}
            2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _gatk4_cnv_output_seg:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        seg = str(rules._gatk4_cnv_correct.output.seg)
    output:
        seg = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{sample_id}--{panel}.seg"
    run:
        op.relative_symlink(input.seg, output.seg, in_module=True)


rule _gatk4_cnv_output_ratio:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        ratio = str(rules._gatk4_cnv_correct.output.ratio)
    output:
        ratio = CFG["dirs"]["outputs"] + "ratio/{seq_type}--{genome_build}/{sample_id}--{panel}.ratio.tsv"
    run:
        op.relative_symlink(input.ratio, output.ratio, in_module=True)


rule _gatk4_cnv_output_purity:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        purity = str(rules._gatk4_cnv_fit.output.purity)
    output:
        purity = CFG["dirs"]["outputs"] + "purity/{seq_type}--{genome_build}/{sample_id}--{panel}.purity.tsv"
    run:
        op.relative_symlink(input.purity, output.purity, in_module=True)


rule _gatk4_cnv_output_plot:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        plot = str(rules._gatk4_cnv_correct.output.plot)
    output:
        plot = CFG["dirs"]["outputs"] + "plot/{seq_type}--{genome_build}/{sample_id}--{panel}.png"
    run:
        op.relative_symlink(input.plot, output.plot, in_module=True)


rule _gatk4_cnv_output_pon:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        pon = str(rules._gatk4_cnv_create_pon.output.pon)
    output:
        pon = CFG["dirs"]["outputs"] + "pon/{genome_build}/{panel}.pon.hdf5"
    run:
        op.relative_symlink(input.pon, output.pon, in_module=True)


rule _gatk4_cnv_output_projected:
    wildcard_constraints:
        panel = r"[A-Za-z0-9_.]+"
    input:
        grid = str(rules._gatk4_cnv_project.output.grid)
    output:
        grid = CFG["dirs"]["outputs"] + "projected/{seq_type}--{genome_build}/{sample_id}--{panel}.grid.tsv"
    run:
        op.relative_symlink(input.grid, output.grid, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _gatk4_cnv_all:
    input:
        expand(
            [
                str(rules._gatk4_cnv_output_seg.output.seg),
                str(rules._gatk4_cnv_output_ratio.output.ratio),
                str(rules._gatk4_cnv_output_purity.output.purity),
                str(rules._gatk4_cnv_output_plot.output.plot),
            ],
            zip,
            seq_type = _GC_RUNS["seq_type"],
            genome_build = _GC_RUNS["genome_build"],
            sample_id = _GC_RUNS["sample_id"],
            panel = _GC_RUNS["panel"]
        ),
        expand(
            str(rules._gatk4_cnv_output_pon.output.pon),
            zip,
            genome_build = [b for b, p in _GC_RUN_PANELS],
            panel = [p for b, p in _GC_RUN_PANELS]
        ),
        expand(
            str(rules._gatk4_cnv_output_projected.output.grid),
            zip,
            seq_type = _GC_CAPTURE_RUNS["seq_type"],
            genome_build = _GC_CAPTURE_RUNS["genome_build"],
            sample_id = _GC_CAPTURE_RUNS["sample_id"],
            panel = _GC_CAPTURE_RUNS["panel"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
