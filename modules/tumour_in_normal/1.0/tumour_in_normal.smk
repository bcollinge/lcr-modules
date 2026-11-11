#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Amaro Taylor-Weiner (deTiN), Giulio Caravagna (TINC)
# Module Author:    Brett Collinge
# Contributors:     N/A


##### SETUP #####

import os
import oncopipe as op

# `CFG` is a shortcut to `config["lcr-modules"]["tumour_in_normal"]`
CFG = op.setup_module(
    name = "tumour_in_normal",
    version = "1.0",
    subdirectories = ["inputs", "callstats", "population", "detin", "tinc", "outputs"],
)

localrules:
    _tumour_in_normal_detin_source,
    _tumour_in_normal_strip_headers,
    _tumour_in_normal_callstats_positions,
    _tumour_in_normal_detin_summary,
    _tumour_in_normal_output,
    _tumour_in_normal_report,
    _tumour_in_normal_all


##### RULES #####


wildcard_constraints:
    seq_type = "[^/]+",
    genome_build = "[^/]+",
    tumour_id = "[^/]+",
    normal_id = "[^/]+",
    pair_status = "[^/]+",
    arm = "ssnv|joint|ascna",
    pset = "[^/]+",
    cna_mode = "none|clonal|majority|state1"


_TIN_PAIR = "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}"
# Scripts run in the deTiN environment take arguments: its pandas cannot read the snakemake object.
_TIN_SRC = os.path.join(workflow.current_basedir, "src") + "/"
_TIN_OPTS = CFG["options"]
_TIN_INPUTS = CFG["inputs"]
_TIN_RUNS = CFG["runs"]


def _tin_given(key):
    return _TIN_INPUTS.get(key) not in ("", None)


def _tin_path(key, wildcards):
    value = _TIN_INPUTS[key]
    return value(wildcards) if callable(value) else value.format(**wildcards)


# deTiN runs: SNVs only (ssnv), or with copy number both the joint estimate (joint) and its copy-number model
# alone (ascna, deTiN's --only_ascnas). deTiN needs the SNV table in every run.
_TIN_SNV = _tin_given("detin_mutect2_vcf")
_TIN_CN = all(_tin_given(k) for k in ("detin_seg", "detin_tumour_hets", "detin_normal_hets"))
_TIN_ARMS = (["joint", "ascna"] if _TIN_CN else ["ssnv"]) if _TIN_SNV else []

_TIN_PSETS = list(_TIN_OPTS["detin_params"] or {})
if _TIN_ARMS:
    assert _TIN_PSETS, "options.detin_params must name at least one set of deTiN settings"
    for _pset, _vals in _TIN_OPTS["detin_params"].items():
        for _key in ("TiN_prior", "mutation_prior", "SSNV_af_threshold"):
            assert (_vals or {}).get(_key) is not None, f"options.detin_params.{_pset}.{_key} must be set"
    assert _TIN_OPTS["detin_call"].get("ci_low_above") is not None, "options.detin_call.ci_low_above must be set"
    assert _tin_given("detin_population_vcf"), "inputs.detin_population_vcf must be set with detin_mutect2_vcf"

_TIN_CNA_MODES = list(_TIN_OPTS["tinc_cna_modes"] or []) if _tin_given("tinc_cna") else ["none"]
if _tin_given("tinc_vcf"):
    assert _TIN_CNA_MODES, "options.tinc_cna_modes must be set when inputs.tinc_cna is given"
    assert "none" not in _TIN_CNA_MODES or not _tin_given("tinc_cna"), "tinc_cna_modes 'none' runs without tinc_cna"
    if set(_TIN_CNA_MODES) - {"clonal", "none"}:
        assert _TIN_OPTS["tinc_cna_filler"] in ("drop", "as_written"), "options.tinc_cna_filler must be set"


# deTiN source at the configured commit; deTiN.py imports its siblings as top-level modules.
rule _tumour_in_normal_detin_source:
    output:
        script = CFG["dirs"]["detin"] + "source/" + _TIN_OPTS["detin_dir"] + "/deTiN/deTiN.py"
    params:
        url = _TIN_OPTS["detin_url"],
        dest = CFG["dirs"]["detin"] + "source/"
    conda:
        CFG["conda_envs"]["detin"]
    shell:
        "curl -fsSL {params.url} | tar -xz -C {params.dest}"


# deTiN's call-stats table from the matched Mutect2 call (src/mutect2_to_callstats.py).
rule _tumour_in_normal_callstats:
    input:
        vcf = _TIN_INPUTS["detin_mutect2_vcf"] or []
    output:
        callstats = CFG["dirs"]["callstats"] + _TIN_PAIR + "/callstats.tsv"
    params:
        candidate_filters = _TIN_OPTS["detin_candidate_filters"]
    threads:
        CFG["threads"]["population"]
    resources:
        **CFG["resources"]["population"]
    script:
        "src/mutect2_to_callstats.py"


# Copy-number inputs without the "@" header lines GATK writes; deTiN skips only "#" lines.
rule _tumour_in_normal_strip_headers:
    input:
        seg = _TIN_INPUTS["detin_seg"] or [],
        tumour_hets = _TIN_INPUTS["detin_tumour_hets"] or [],
        normal_hets = _TIN_INPUTS["detin_normal_hets"] or []
    output:
        seg = CFG["dirs"]["inputs"] + "detin/" + _TIN_PAIR + "/seg.tsv",
        tumour_hets = CFG["dirs"]["inputs"] + "detin/" + _TIN_PAIR + "/tumour_hets.tsv",
        normal_hets = CFG["dirs"]["inputs"] + "detin/" + _TIN_PAIR + "/normal_hets.tsv"
    shell:
        op.as_one_line("""
        grep -v '^@' {input.seg} > {output.seg} &&
        grep -v '^@' {input.tumour_hets} > {output.tumour_hets} &&
        grep -v '^@' {input.normal_hets} > {output.normal_hets}
        """)


rule _tumour_in_normal_callstats_positions:
    input:
        callstats = str(rules._tumour_in_normal_callstats.output.callstats)
    output:
        positions = temp(CFG["dirs"]["population"] + _TIN_PAIR + "/positions.tsv")
    params:
        script = _TIN_SRC + "callstats_positions.py"
    conda:
        CFG["conda_envs"]["detin"]
    shell:
        "python {params.script} {input.callstats} {output.positions}"


# Population records at the call-stats positions.
rule _tumour_in_normal_population_sites:
    input:
        positions = str(rules._tumour_in_normal_callstats_positions.output.positions),
        vcf = _TIN_INPUTS["detin_population_vcf"] or []
    output:
        sites = CFG["dirs"]["population"] + _TIN_PAIR + "/sites.tsv"
    log:
        stderr = CFG["logs"]["population"] + _TIN_PAIR + "/population_sites.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["population"]
    resources:
        **CFG["resources"]["population"]
    shell:
        op.as_one_line("""
        bcftools query -T {input.positions} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\n' {input.vcf}
            > {output.sites} 2> {log.stderr}
        """)


# deTiN's population pickle (build_exac_pickle layout), restricted to the pair's positions.
rule _tumour_in_normal_population_pickle:
    input:
        sites = str(rules._tumour_in_normal_population_sites.output.sites)
    output:
        pickle = CFG["dirs"]["population"] + _TIN_PAIR + "/population.pickle"
    params:
        min_af = _TIN_OPTS["population_min_af"],
        script = _TIN_SRC + "population_pickle.py"
    conda:
        CFG["conda_envs"]["detin"]
    threads:
        CFG["threads"]["population"]
    resources:
        **CFG["resources"]["population"]
    shell:
        "python {params.script} {input.sites} {params.min_af} {output.pickle}"


def _tin_detin_inputs(wildcards):
    inputs = {
        "script": str(rules._tumour_in_normal_detin_source.output.script),
        "callstats": str(rules._tumour_in_normal_callstats.output.callstats),
        "pickle": str(rules._tumour_in_normal_population_pickle.output.pickle)
    }
    if wildcards.arm in ("joint", "ascna"):
        inputs["seg"] = str(rules._tumour_in_normal_strip_headers.output.seg)
        inputs["tumour_hets"] = str(rules._tumour_in_normal_strip_headers.output.tumour_hets)
        inputs["normal_hets"] = str(rules._tumour_in_normal_strip_headers.output.normal_hets)
    return inputs


def _tin_detin_args(wildcards, input):
    args = ["--mutation_data_path", input.callstats, "--exac_data_path", input.pickle]
    if "seg" in input.keys():
        args += ["--cn_data_path", input.seg,
                 "--tumor_het_data_path", input.tumour_hets,
                 "--normal_het_data_path", input.normal_hets]
    if wildcards.arm == "ascna":
        args.append("--only_ascnas")
    args += [f"--{key} {value}" for key, value in _TIN_OPTS["detin_params"][wildcards.pset].items()]
    return " ".join(args)


# deTiN as published. Its stdout is kept as an output: it holds each model's own estimate and site count.
rule _tumour_in_normal_detin:
    input:
        unpack(_tin_detin_inputs)
    output:
        tin = CFG["dirs"]["detin"] + _TIN_PAIR + "/{pset}/{arm}/{tumour_id}--{normal_id}.TiN_estimate.txt",
        ci = CFG["dirs"]["detin"] + _TIN_PAIR + "/{pset}/{arm}/{tumour_id}--{normal_id}.TiN_estimate_CI.txt",
        added = CFG["dirs"]["detin"] + _TIN_PAIR + "/{pset}/{arm}/{tumour_id}--{normal_id}.number_of_SSNVs_added.txt",
        stdout = CFG["dirs"]["detin"] + _TIN_PAIR + "/{pset}/{arm}/deTiN.stdout.txt"
    log:
        stderr = CFG["logs"]["detin"] + _TIN_PAIR + "/{pset}/{arm}.stderr.log"
    params:
        args = _tin_detin_args,
        opts_seq_type = lambda w: _TIN_OPTS["detin_seq_type"].get(w.seq_type, ""),
        outdir = CFG["dirs"]["detin"] + _TIN_PAIR + "/{pset}/{arm}"
    conda:
        CFG["conda_envs"]["detin"]
    threads:
        CFG["threads"]["detin"]
    resources:
        **CFG["resources"]["detin"]
    shell:
        op.as_one_line("""
        MPLBACKEND=Agg python {input.script} {params.args}
            --output_name {wildcards.tumour_id}--{wildcards.normal_id} --output_dir {params.outdir}
            {params.opts_seq_type}
            > {output.stdout} 2> {log.stderr}
        """)


rule _tumour_in_normal_detin_summary:
    input:
        tin = expand(str(rules._tumour_in_normal_detin.output.tin), pset = _TIN_PSETS, arm = _TIN_ARMS, allow_missing = True),
        ci = expand(str(rules._tumour_in_normal_detin.output.ci), pset = _TIN_PSETS, arm = _TIN_ARMS, allow_missing = True),
        added = expand(str(rules._tumour_in_normal_detin.output.added), pset = _TIN_PSETS, arm = _TIN_ARMS, allow_missing = True),
        stdout = expand(str(rules._tumour_in_normal_detin.output.stdout), pset = _TIN_PSETS, arm = _TIN_ARMS, allow_missing = True)
    output:
        tsv = CFG["dirs"]["detin"] + _TIN_PAIR + "/detin.tsv"
    params:
        runs = [(p, a) for p in _TIN_PSETS for a in _TIN_ARMS]
    script:
        "src/detin_summary.py"


def _tin_existing(key):
    # An input given as a path pattern is used when the file exists; a callable is taken as produced upstream.
    def resolve(wildcards):
        if not _tin_given(key):
            return ""
        path = _tin_path(key, wildcards)
        return path if callable(_TIN_INPUTS[key]) or os.path.exists(path) else ""
    return resolve


# TINC as published: Strelka PASS SNVs through TINC's own parser, Battenberg copy number per cna_mode.
rule _tumour_in_normal_tinc:
    input:
        vcf = lambda w: _tin_existing("tinc_vcf")(w) or [],
        cna = lambda w: (_tin_existing("tinc_cna")(w) or []) if w.cna_mode != "none" else []
    output:
        tsv = CFG["dirs"]["tinc"] + _TIN_PAIR + "/{cna_mode}/tinc.tsv",
        rds = CFG["dirs"]["tinc"] + _TIN_PAIR + "/{cna_mode}/tinc.rds",
        pdf = CFG["dirs"]["tinc"] + _TIN_PAIR + "/{cna_mode}/tinc_report.pdf"
    log:
        stderr = CFG["logs"]["tinc"] + _TIN_PAIR + "/{cna_mode}.tinc.stderr.log"
    params:
        vcf = _tin_existing("tinc_vcf"),
        cna = _tin_existing("tinc_cna"),
        filler = _TIN_OPTS["tinc_cna_filler"] or "",
        seed = _TIN_OPTS["tinc_seed"],
        fast = bool(_TIN_OPTS["tinc_fast"])
    conda:
        CFG["conda_envs"]["tinc"]
    threads:
        CFG["threads"]["tinc"]
    resources:
        **CFG["resources"]["tinc"]
    script:
        "src/tinc_run.R"


rule _tumour_in_normal_output:
    input:
        tsv = str(rules._tumour_in_normal_detin_summary.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "detin/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.detin.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module = True)


def _tin_pairs(tool):
    runs = _TIN_RUNS[_TIN_RUNS["tumour_seq_type"].isin(CFG["tools"].get(tool, []))]
    return [dict(seq_type = r.tumour_seq_type, genome_build = r.tumour_genome_build, tumour_id = r.tumour_sample_id,
                 normal_id = r.normal_sample_id, pair_status = r.pair_status) for r in runs.itertuples()]


# deTiN's coordinates are hg19 linear positions; it runs only on the builds listed in options.detin_builds.
_TIN_DETIN_RUN = [p for p in _tin_pairs("detin") if _TIN_ARMS and p["genome_build"] in _TIN_OPTS["detin_builds"]]
_TIN_DETIN_NOT_RUN = [dict(p, tool = "detin", reason = "deTiN not run: " + (
    "no SNV input" if not _TIN_ARMS else f"build {p['genome_build']} not in options.detin_builds"))
    for p in _tin_pairs("detin") if p not in _TIN_DETIN_RUN]
_TIN_TINC_RUN = _tin_pairs("tinc") if _tin_given("tinc_vcf") else []


rule _tumour_in_normal_report:
    input:
        detin = [str(rules._tumour_in_normal_detin_summary.output.tsv).format(**p) for p in _TIN_DETIN_RUN],
        tinc = [str(rules._tumour_in_normal_tinc.output.tsv).format(cna_mode = m, **p)
                for p in _TIN_TINC_RUN for m in _TIN_CNA_MODES]
    output:
        pairs = CFG["dirs"]["outputs"] + "report/pairs.tsv",
        normals = CFG["dirs"]["outputs"] + "report/normals.tsv"
    params:
        opts = dict(_TIN_OPTS, tinc_cna_modes = _TIN_CNA_MODES, detin_psets = _TIN_PSETS),
        not_run = _TIN_DETIN_NOT_RUN
    script:
        "src/report.py"


rule _tumour_in_normal_all:
    input:
        [str(rules._tumour_in_normal_output.output.tsv).format(**p) for p in _TIN_DETIN_RUN],
        str(rules._tumour_in_normal_report.output.pairs),
        str(rules._tumour_in_normal_report.output.normals)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
