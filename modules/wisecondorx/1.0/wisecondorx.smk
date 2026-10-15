#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Lennart Raman (WisecondorX)
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
# `CFG` is a shortcut to `config["lcr-modules"]["wisecondorx"]`
CFG = op.setup_module(
    name = "wisecondorx",
    version = "1.0",
    subdirectories = ["inputs", "convert", "newref", "predict", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _wisecondorx_input_bam,
    _wisecondorx_output,
    _wisecondorx_all,


##### REFERENCE SETS AND PREDICTIONS #####

# Kept as module globals: op.cleanup_module removes `CFG` and its sample table.
_WCX_SAMPLES = CFG["samples"][["sample_id", "patient_id", "seq_type", "genome_build"]].drop_duplicates()
_WCX_BAM = CFG["inputs"]["sample_bam"]
_WCX_BLACKLIST = CFG["inputs"].get("blacklist") or []
_WCX_INPUTS_DIR = CFG["dirs"]["inputs"]


def _wisecondorx_sample_row(sample_id, context):
    rows = _WCX_SAMPLES[_WCX_SAMPLES["sample_id"] == sample_id]
    if len(rows) != 1:
        sys.exit(
            f"wisecondorx: {context} names {sample_id!r}, which appears {len(rows)} times in the "
            "samples table (after the pairing_config seq_type filter); it must appear exactly once."
        )
    return rows.iloc[0]


# name -> {"genome_build": str, "members": [(seq_type, genome_build, sample_id)], "patients": set}
_WCX_REFS = {}
for _name, _members in (CFG["references"] or {}).items():
    if not re.fullmatch(r"[A-Za-z0-9_.]+", _name):
        sys.exit(f"wisecondorx: reference name {_name!r} may use letters, digits, '_' and '.' only.")
    _members = list(dict.fromkeys(_members))
    if len(_members) < 10:
        sys.exit(f"wisecondorx: reference {_name!r} has {len(_members)} members; newref needs at least 10.")
    _rows = [_wisecondorx_sample_row(m, f"reference {_name!r}") for m in _members]
    _builds = sorted({r["genome_build"] for r in _rows})
    if len(_builds) != 1:
        sys.exit(f"wisecondorx: reference {_name!r} spans genome builds {_builds}.")
    _WCX_REFS[_name] = {
        "genome_build": _builds[0],
        "members": [(r["seq_type"], r["genome_build"], r["sample_id"]) for r in _rows],
        "patients": {r["patient_id"] for r in _rows},
    }

# One row per (sample, reference) prediction
_WCX_PRED = {"seq_type": [], "genome_build": [], "sample_id": [], "reference": []}
for _sid, _refs in (CFG["predict"] or {}).items():
    _row = _wisecondorx_sample_row(_sid, "predict")
    for _ref in _refs:
        if _ref not in _WCX_REFS:
            sys.exit(f"wisecondorx: predict {_sid!r} names reference {_ref!r}, which is not defined.")
        if _WCX_REFS[_ref]["genome_build"] != _row["genome_build"]:
            sys.exit(f"wisecondorx: {_sid!r} is {_row['genome_build']} but reference {_ref!r} is {_WCX_REFS[_ref]['genome_build']}.")
        if _row["patient_id"] in _WCX_REFS[_ref]["patients"]:
            sys.exit(f"wisecondorx: reference {_ref!r} holds a sample of {_sid!r}'s patient; exclude it from the set.")
        _WCX_PRED["seq_type"].append(_row["seq_type"])
        _WCX_PRED["genome_build"].append(_row["genome_build"])
        _WCX_PRED["sample_id"].append(_sid)
        _WCX_PRED["reference"].append(_ref)


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/'),
# keeping the .bam / .cram extension WisecondorX uses to pick the reader
rule _wisecondorx_input_bam:
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


def _wisecondorx_get_reads(wildcards):
    w = dict(wildcards.items())
    src = str(_WCX_BAM(wildcards) if callable(_WCX_BAM) else _WCX_BAM.format(**w))
    if src.endswith(".cram"):
        ext = "cram"
    elif src.endswith(".bam"):
        ext = "bam"
    else:
        raise ValueError(f"wisecondorx: {src} is neither .bam nor .cram")
    return f"{_WCX_INPUTS_DIR}bam/{w['seq_type']}--{w['genome_build']}/{w['sample_id']}.{ext}"


# Counts reads per bin (WisecondorX convert)
rule _wisecondorx_convert:
    input:
        bam = _wisecondorx_get_reads,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        npz = CFG["dirs"]["convert"] + "{seq_type}--{genome_build}/{sample_id}.npz"
    log:
        stdout = CFG["logs"]["convert"] + "{seq_type}--{genome_build}/{sample_id}/convert.stdout.log",
        stderr = CFG["logs"]["convert"] + "{seq_type}--{genome_build}/{sample_id}/convert.stderr.log"
    params:
        binsize = CFG["options"]["convert_binsize"],
        opts = CFG["options"]["convert"]
    conda:
        CFG["conda_envs"]["wisecondorx"]
    threads:
        CFG["threads"]["convert"]
    resources:
        **CFG["resources"]["convert"]
    shell:
        op.as_one_line("""
        WisecondorX convert {input.bam} {output.npz}
            --reference {input.fasta}
            --binsize {params.binsize}
            {params.opts}
            > {log.stdout}
            2> {log.stderr}
        """)


def _wisecondorx_get_reference_npz(wildcards):
    return [
        str(rules._wisecondorx_convert.output.npz).format(seq_type=s, genome_build=b, sample_id=i)
        for (s, b, i) in _WCX_REFS[wildcards.reference]["members"]
    ]


# Builds one reference set (WisecondorX newref)
rule _wisecondorx_newref:
    wildcard_constraints:
        reference = r"[A-Za-z0-9_.]+"
    input:
        npz = _wisecondorx_get_reference_npz
    output:
        npz = CFG["dirs"]["newref"] + "{genome_build}/{reference}.npz"
    log:
        stdout = CFG["logs"]["newref"] + "{genome_build}/{reference}/newref.stdout.log",
        stderr = CFG["logs"]["newref"] + "{genome_build}/{reference}/newref.stderr.log"
    params:
        binsize = CFG["options"]["newref_binsize"],
        yfrac = "" if CFG["options"]["yfrac"] is None else "--yfrac " + str(CFG["options"]["yfrac"]),
        opts = CFG["options"]["newref"]
    conda:
        CFG["conda_envs"]["wisecondorx"]
    threads:
        CFG["threads"]["newref"]
    resources:
        **CFG["resources"]["newref"]
    shell:
        op.as_one_line("""
        WisecondorX newref {input.npz} {output.npz}
            --binsize {params.binsize}
            --cpus {threads}
            {params.yfrac}
            {params.opts}
            > {log.stdout}
            2> {log.stderr}
        """)


# Calls copy number against one reference set (WisecondorX predict)
rule _wisecondorx_predict:
    wildcard_constraints:
        reference = r"[A-Za-z0-9_.]+"
    input:
        npz = str(rules._wisecondorx_convert.output.npz),
        reference = str(rules._wisecondorx_newref.output.npz),
        blacklist = _WCX_BLACKLIST
    output:
        bins = CFG["dirs"]["predict"] + "{seq_type}--{genome_build}/{sample_id}--{reference}/{sample_id}--{reference}_bins.bed",
        segments = CFG["dirs"]["predict"] + "{seq_type}--{genome_build}/{sample_id}--{reference}/{sample_id}--{reference}_segments.bed",
        aberrations = CFG["dirs"]["predict"] + "{seq_type}--{genome_build}/{sample_id}--{reference}/{sample_id}--{reference}_aberrations.bed",
        statistics = CFG["dirs"]["predict"] + "{seq_type}--{genome_build}/{sample_id}--{reference}/{sample_id}--{reference}_statistics.txt"
    log:
        stdout = CFG["logs"]["predict"] + "{seq_type}--{genome_build}/{sample_id}--{reference}/predict.stdout.log",
        stderr = CFG["logs"]["predict"] + "{seq_type}--{genome_build}/{sample_id}--{reference}/predict.stderr.log"
    params:
        outid = CFG["dirs"]["predict"] + "{seq_type}--{genome_build}/{sample_id}--{reference}/{sample_id}--{reference}",
        blacklist = lambda w, input: f"--blacklist {input.blacklist}" if input.blacklist else "",
        opts = CFG["options"]["predict"]
    conda:
        CFG["conda_envs"]["wisecondorx"]
    threads:
        CFG["threads"]["predict"]
    resources:
        **CFG["resources"]["predict"]
    shell:
        op.as_one_line("""
        WisecondorX predict {input.npz} {input.reference} {params.outid}
            --bed
            {params.blacklist}
            {params.opts}
            > {log.stdout}
            2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _wisecondorx_output:
    wildcard_constraints:
        reference = r"[A-Za-z0-9_.]+"
    input:
        bins = str(rules._wisecondorx_predict.output.bins),
        segments = str(rules._wisecondorx_predict.output.segments),
        aberrations = str(rules._wisecondorx_predict.output.aberrations),
        statistics = str(rules._wisecondorx_predict.output.statistics)
    output:
        bins = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}/{sample_id}--{reference}_bins.bed",
        segments = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}/{sample_id}--{reference}_segments.bed",
        aberrations = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}/{sample_id}--{reference}_aberrations.bed",
        statistics = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{sample_id}--{reference}_statistics.txt"
    run:
        op.relative_symlink(input.bins, output.bins, in_module=True)
        op.relative_symlink(input.segments, output.segments, in_module=True)
        op.relative_symlink(input.aberrations, output.aberrations, in_module=True)
        op.relative_symlink(input.statistics, output.statistics, in_module=True)


# Generates the target sentinels for each prediction, which generate the symlinks
rule _wisecondorx_all:
    input:
        expand(
            [
                str(rules._wisecondorx_output.output.bins),
                str(rules._wisecondorx_output.output.segments),
                str(rules._wisecondorx_output.output.aberrations),
                str(rules._wisecondorx_output.output.statistics),
            ],
            zip,
            seq_type = _WCX_PRED["seq_type"],
            genome_build = _WCX_PRED["genome_build"],
            sample_id = _WCX_PRED["sample_id"],
            reference = _WCX_PRED["reference"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
