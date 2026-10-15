#!/usr/bin/env snakemake


# ---------------------------------------------------------------------------- #
##### ATTRIBUTION #####
# ---------------------------------------------------------------------------- #

# Original snakemake author: Jasper Wong
# Module author: Jasper Wong
# Additional contributors: Brett Collinge


# ---------------------------------------------------------------------------- #
##### SETUP #####
# ---------------------------------------------------------------------------- #

### Modules ###

import os
import sys
import oncopipe as op

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

### Directories ###
CFG = op.setup_module(
    name = "ichorcna",
    version = "1.2",
    subdirectories = ["inputs", "bam", "wig", "pon", "run", "outputs"]
)

wildcard_constraints:
    bin_size = r"\d+"

localrules:
    _ichorcna_input_bam,
    _ichorcna_pon_filelist,
    _ichorcna_output,
    _ichorcna_all


# ---------------------------------------------------------------------------- #
##### PANELS OF NORMALS #####
# ---------------------------------------------------------------------------- #

# Resolved at parse time from options.panel_of_normals, per sample_id: the panel
# name used in output paths and the panel file the run reads (the shared panel,
# or one with members removed).
_ICHOR_OPTS = CFG["options"]
_ICHOR_PON = _ICHOR_OPTS["panel_of_normals"]
_ICHOR_MEMBERS = {k: list(map(str, v)) for k, v in (_ICHOR_PON["members"] or {}).items()}
_ICHOR_ASSIGN = {str(k): v for k, v in (_ICHOR_PON["assign"] or {}).items()}
_ICHOR_EXCLUDE = {str(k): set(map(str, v)) for k, v in (_ICHOR_PON["exclude"] or {}).items()}

for _p in set(_ICHOR_ASSIGN.values()) - set(_ICHOR_MEMBERS):
    raise ValueError(f"ichorcna: panel {_p!r} is assigned but has no members.")

_ICHOR_PANEL_FILES = {}  # panel file id -> member sample_ids
_ICHOR_SAMPLE_PANEL = {}  # sample_id -> (panel, panel file id)
for _s in CFG["samples"]["sample_id"].astype(str):
    _p = _ICHOR_ASSIGN.get(_s, _ICHOR_ASSIGN.get("_default"))
    if _p is None:
        continue
    _drop = ({_s} | _ICHOR_EXCLUDE.get(_s, set())) & set(_ICHOR_MEMBERS[_p])
    _fid = f"{_p}--minus--{_s}" if _drop else _p
    _ICHOR_PANEL_FILES[_fid] = [m for m in _ICHOR_MEMBERS[_p] if m not in _drop]
    assert _ICHOR_PANEL_FILES[_fid], f"ichorcna: panel {_fid!r} has no members left."
    _ICHOR_SAMPLE_PANEL[_s] = (_p, _fid)

_ICHOR_BINS = [int(b) for b in _ICHOR_OPTS["bin_sizes"]]


def _ichorcna_bin_label(bin_size):
    return f"{int(bin_size) // 1000}kb"


def _ichorcna_build(wildcards):
    return _ICHOR_OPTS["ichorcna_build"][wildcards.genome_build]


def _ichorcna_style(wildcards):
    return _ICHOR_OPTS["genome_style"][wildcards.genome_build]


def _ichorcna_extdata(key):
    def _f(wildcards):
        return _ICHOR_OPTS[key].format(
            ichorcna_build = _ichorcna_build(wildcards),
            bin_label = _ichorcna_bin_label(wildcards.bin_size))
    return _f


def _ichorcna_centromere(wildcards):
    return _ICHOR_OPTS["centromere"][_ichorcna_build(wildcards)]


def _ichorcna_chromosomes(wildcards):
    chroms = [str(c) for c in range(1, 23)] + ["X", "Y"]
    if _ichorcna_style(wildcards) == "UCSC":
        chroms = ["chr" + c for c in chroms]
    return ",".join(chroms)


def _ichorcna_r_args(options):
    return " ".join(f"--{k} '{v}'" for k, v in options.items())


# ---------------------------------------------------------------------------- #
##### RULES #####
# ---------------------------------------------------------------------------- #

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _ichorcna_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)


# readCounter reads BAM only; CRAM (or BAM) input is written to an indexed BAM
rule _ichorcna_bam:
    input:
        aln = str(rules._ichorcna_input_bam.output.bam),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        bam = temp(CFG["dirs"]["bam"] + "{seq_type}--{genome_build}/{sample_id}.bam"),
        bai = temp(CFG["dirs"]["bam"] + "{seq_type}--{genome_build}/{sample_id}.bam.bai")
    log:
        stderr = CFG["logs"]["bam"] + "{seq_type}--{genome_build}/{sample_id}.stderr.log"
    conda:
        CFG["conda_envs"]["ichorcna"]
    threads:
        CFG["threads"]["bam"]
    resources:
        **CFG["resources"]["bam"]
    shell:
        op.as_one_line("""
        samtools view -b -@ {threads} -T {input.fasta} -o {output.bam} {input.aln} 2> {log.stderr} &&
        samtools index -@ {threads} {output.bam} 2>> {log.stderr}
        """)


rule _ichorcna_read_counter:
    input:
        bam = str(rules._ichorcna_bam.output.bam),
        bai = str(rules._ichorcna_bam.output.bai)
    output:
        wig = CFG["dirs"]["wig"] + "{seq_type}--{genome_build}/{bin_size}/{sample_id}.wig"
    log:
        stderr = CFG["logs"]["wig"] + "{seq_type}--{genome_build}/{bin_size}/{sample_id}.stderr.log"
    params:
        quality = CFG["options"]["readcounter"]["quality"],
        chroms = _ichorcna_chromosomes
    conda:
        CFG["conda_envs"]["ichorcna"]
    threads:
        CFG["threads"]["readcounter"]
    resources:
        **CFG["resources"]["readcounter"]
    shell:
        op.as_one_line("""
        readCounter --window {wildcards.bin_size} --quality {params.quality}
        --chromosome "{params.chroms}" {input.bam} > {output.wig} 2> {log.stderr}
        """)


def _ichorcna_pon_wigs(wildcards):
    return expand(
        str(rules._ichorcna_read_counter.output.wig),
        seq_type = wildcards.seq_type,
        genome_build = wildcards.genome_build,
        bin_size = wildcards.bin_size,
        sample_id = _ICHOR_PANEL_FILES[wildcards.panel_file])


rule _ichorcna_pon_filelist:
    input:
        wigs = _ichorcna_pon_wigs
    output:
        filelist = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{bin_size}/{panel_file}.files.txt"
    run:
        with open(output.filelist, "w") as fh:
            fh.write("".join(os.path.abspath(w) + "\n" for w in input.wigs))


rule _ichorcna_pon:
    input:
        filelist = str(rules._ichorcna_pon_filelist.output.filelist)
    output:
        rds = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{bin_size}/{panel_file}_" + CFG["options"]["panel_of_normals"]["method"] + ".rds"
    log:
        stdout = CFG["logs"]["pon"] + "{seq_type}--{genome_build}/{bin_size}/{panel_file}.stdout.log",
        stderr = CFG["logs"]["pon"] + "{seq_type}--{genome_build}/{bin_size}/{panel_file}.stderr.log"
    params:
        script = CFG["options"]["scripts"]["pon"],
        prefix = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{bin_size}/{panel_file}",
        method = CFG["options"]["panel_of_normals"]["method"],
        gc = _ichorcna_extdata("gc_wig"),
        map = _ichorcna_extdata("map_wig"),
        centromere = _ichorcna_centromere,
        build = _ichorcna_build,
        style = _ichorcna_style,
        opts = _ichorcna_r_args(CFG["options"]["pon"])
    conda:
        CFG["conda_envs"]["ichorcna"]
    threads:
        CFG["threads"]["pon"]
    resources:
        **CFG["resources"]["pon"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {params.script} --filelist {input.filelist} --outfile {params.prefix}
        --method {params.method} --gcWig {params.gc} --mapWig {params.map} --centromere {params.centromere}
        --genomeBuild {params.build} --genomeStyle {params.style} {params.opts}
        > {log.stdout} 2> {log.stderr}
        """)


def _ichorcna_panel_rds(wildcards):
    panel, panel_file = _ICHOR_SAMPLE_PANEL[wildcards.tumour_id]
    assert panel == wildcards.panel, f"{wildcards.tumour_id} is assigned panel {panel!r}, not {wildcards.panel!r}"
    return str(rules._ichorcna_pon.output.rds).format(
        seq_type = wildcards.seq_type, genome_build = wildcards.genome_build,
        bin_size = wildcards.bin_size, panel_file = panel_file)


rule _ichorcna_run:
    input:
        wig = CFG["dirs"]["wig"] + "{seq_type}--{genome_build}/{bin_size}/{tumour_id}.wig",
        pon = _ichorcna_panel_rds
    output:
        cna = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{bin_size}/{panel}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cna.seg",
        seg = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{bin_size}/{panel}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.seg",
        seg_txt = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{bin_size}/{panel}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.seg.txt",
        params = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{bin_size}/{panel}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.params.txt",
        depth = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{bin_size}/{panel}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.correctedDepth.txt"
    log:
        stdout = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{bin_size}/{panel}/{tumour_id}--{normal_id}--{pair_status}.stdout.log",
        stderr = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{bin_size}/{panel}/{tumour_id}--{normal_id}--{pair_status}.stderr.log"
    params:
        script = CFG["options"]["scripts"]["run"],
        out_dir = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{bin_size}/{panel}/{tumour_id}--{normal_id}--{pair_status}",
        gc = _ichorcna_extdata("gc_wig"),
        map = _ichorcna_extdata("map_wig"),
        centromere = _ichorcna_centromere,
        build = _ichorcna_build,
        style = _ichorcna_style,
        opts = _ichorcna_r_args(CFG["options"]["run"])
    conda:
        CFG["conda_envs"]["ichorcna"]
    threads:
        CFG["threads"]["run"]
    resources:
        **CFG["resources"]["run"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {params.script} --id {wildcards.tumour_id} --WIG {input.wig}
        --normalPanel {input.pon} --gcWig {params.gc} --mapWig {params.map} --centromere {params.centromere}
        --genomeBuild {params.build} --genomeStyle {params.style} --cores {threads} {params.opts}
        --outDir {params.out_dir} > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _ichorcna_output:
    input:
        cna = str(rules._ichorcna_run.output.cna),
        seg = str(rules._ichorcna_run.output.seg),
        params = str(rules._ichorcna_run.output.params)
    output:
        cna = CFG["dirs"]["outputs"] + "cna_seg/{seq_type}--{genome_build}/{bin_size}/{panel}/{tumour_id}--{normal_id}--{pair_status}.cna.seg",
        seg = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{bin_size}/{panel}/{tumour_id}--{normal_id}--{pair_status}.seg",
        params = CFG["dirs"]["outputs"] + "params/{seq_type}--{genome_build}/{bin_size}/{panel}/{tumour_id}--{normal_id}--{pair_status}.params.txt"
    run:
        op.relative_symlink(input.cna, output.cna, in_module = True)
        op.relative_symlink(input.seg, output.seg, in_module = True)
        op.relative_symlink(input.params, output.params, in_module = True)


# Runs for tumours with a panel, at every bin size
_ICHOR_RUNS = CFG["runs"][CFG["runs"]["tumour_sample_id"].astype(str).isin(_ICHOR_SAMPLE_PANEL)]

rule _ichorcna_all:
    input:
        [
            str(out).format(
                seq_type = r.tumour_seq_type,
                genome_build = r.tumour_genome_build,
                bin_size = b,
                panel = _ICHOR_SAMPLE_PANEL[str(r.tumour_sample_id)][0],
                tumour_id = r.tumour_sample_id,
                normal_id = r.normal_sample_id,
                pair_status = r.pair_status)
            for r in _ICHOR_RUNS.itertuples()
            for b in _ICHOR_BINS
            for out in rules._ichorcna_output.output
        ]


##### CLEANUP #####

op.cleanup_module(CFG)
