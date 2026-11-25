#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Giulio Caravagna (TINC)
# Module Author:    Brett Collinge
# Contributors:     N/A


##### SETUP #####

import oncopipe as op

# `CFG` is a shortcut to `config["lcr-modules"]["tinc"]`
CFG = op.setup_module(
    name = "tinc",
    version = "1.0",
    subdirectories = ["inputs", "vcf", "tinc", "outputs"],
)

localrules:
    _tinc_input_maf,
    _tinc_input_subclones,
    _tinc_output_tsv,
    _tinc_all


##### RULES #####


wildcard_constraints:
    seq_type = "[^/]+",
    genome_build = "[^/]+",
    tumour_id = "[^/]+",
    normal_id = "[^/]+",
    pair_status = "[^/]+"


_TINC_PAIR = "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}"
_TINC_RUNS = CFG["runs"]

assert CFG["options"]["cna_mode"] in ("clonal", "majority", "state1"), "options.cna_mode"
assert CFG["options"]["cna_filler"] in ("drop", "as_written"), "options.cna_filler"


rule _tinc_input_maf:
    input:
        maf = CFG["inputs"]["maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/" + _TINC_PAIR + ".maf"
    run:
        op.absolute_symlink(input.maf, output.maf)


rule _tinc_input_subclones:
    input:
        subclones = CFG["inputs"]["subclones"]
    output:
        subclones = CFG["dirs"]["inputs"] + "subclones/" + _TINC_PAIR + ".subclones.txt"
    run:
        op.absolute_symlink(input.subclones, output.subclones)


# SNVs of the MAF as a two-sample VCF: ##tumor_sample / ##normal_sample header lines, FORMAT AD from the MAF's
# ref and alt counts, FILTER PASS.
rule _tinc_maf_to_vcf:
    input:
        maf = str(rules._tinc_input_maf.output.maf)
    output:
        vcf = CFG["dirs"]["vcf"] + _TINC_PAIR + ".vcf"
    log:
        stderr = CFG["logs"]["vcf"] + _TINC_PAIR + ".maf_to_vcf.stderr.log"
    conda:
        CFG["conda_envs"]["tinc"]
    threads:
        CFG["threads"]["maf_to_vcf"]
    resources:
        **CFG["resources"]["maf_to_vcf"]
    script:
        "src/maf_to_vcf.R"


rule _tinc_run:
    input:
        vcf = str(rules._tinc_maf_to_vcf.output.vcf),
        subclones = str(rules._tinc_input_subclones.output.subclones)
    output:
        tsv = CFG["dirs"]["tinc"] + _TINC_PAIR + "/tinc.tsv",
        rds = CFG["dirs"]["tinc"] + _TINC_PAIR + "/tinc.rds",
        pdf = CFG["dirs"]["tinc"] + _TINC_PAIR + "/tinc_report.pdf"
    log:
        stderr = CFG["logs"]["tinc"] + _TINC_PAIR + ".tinc.stderr.log"
    params:
        cna_mode = CFG["options"]["cna_mode"],
        cna_filler = CFG["options"]["cna_filler"],
        seed = CFG["options"]["seed"],
        fast = bool(CFG["options"]["fast"]),
        call_estimate = CFG["options"]["call"]["estimate"],
        call_above = CFG["options"]["call"]["above"],
        call_purity = CFG["options"]["call"]["purity_at_least"]
    conda:
        CFG["conda_envs"]["tinc"]
    threads:
        CFG["threads"]["tinc"]
    resources:
        mem_mb = lambda wildcards, attempt: CFG["resources"]["tinc"]["mem_mb"] * attempt
    script:
        "src/tinc_run.R"


rule _tinc_output_tsv:
    input:
        tsv = str(rules._tinc_run.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/" + _TINC_PAIR + ".tinc.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module = True)


def _tinc_build_tsvs(wildcards):
    runs = _TINC_RUNS[(_TINC_RUNS["tumour_seq_type"] == wildcards.seq_type) &
                      (_TINC_RUNS["tumour_genome_build"] == wildcards.genome_build)]
    return expand(str(rules._tinc_output_tsv.output.tsv), zip,
                  seq_type = runs["tumour_seq_type"],
                  genome_build = runs["tumour_genome_build"],
                  tumour_id = runs["tumour_sample_id"],
                  normal_id = runs["normal_sample_id"],
                  pair_status = runs["pair_status"])


# One row per pair, per seq_type and genome build.
rule _tinc_merge:
    input:
        tsvs = _tinc_build_tsvs
    output:
        tsv = CFG["dirs"]["outputs"] + "merged/{seq_type}--{genome_build}.tinc.tsv"
    resources:
        **CFG["resources"]["merge"]
    run:
        with open(output.tsv, "w") as out:
            for k, path in enumerate(input.tsvs):
                with open(path) as f:
                    header = f.readline()
                    if k == 0:
                        out.write(header)
                    out.writelines(f)


rule _tinc_all:
    input:
        expand(str(rules._tinc_merge.output.tsv), zip,
               seq_type = _TINC_RUNS[["tumour_seq_type", "tumour_genome_build"]].drop_duplicates()["tumour_seq_type"],
               genome_build = _TINC_RUNS[["tumour_seq_type", "tumour_genome_build"]].drop_duplicates()["tumour_genome_build"])


##### CLEANUP #####


op.cleanup_module(CFG)
