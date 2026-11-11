#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Amaro Taylor-Weiner (deTiN), Giulio Caravagna (TINC)
# Module Author:    Brett Collinge
# Contributors:     N/A


##### SETUP #####

import oncopipe as op

# `CFG` is a shortcut to `config["lcr-modules"]["tumour_in_normal"]`
CFG = op.setup_module(
    name = "tumour_in_normal",
    version = "1.0",
    subdirectories = ["inputs", "population", "detin", "tinc", "outputs"],
)

localrules:
    _tumour_in_normal_detin_source,
    _tumour_in_normal_strip_headers,
    _tumour_in_normal_detin_summary,
    _tumour_in_normal_output,
    _tumour_in_normal_all


##### RULES #####


wildcard_constraints:
    seq_type = "[^/]+",
    genome_build = "[^/]+",
    tumour_id = "[^/]+",
    normal_id = "[^/]+",
    pair_status = "[^/]+",
    arm = "ssnv|ascna|joint"


_TIN_PAIR = "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}"
_TIN_OPTS = CFG["options"]
_TIN_INPUTS = CFG["inputs"]


def _tin_given(key):
    return CFG["inputs"].get(key, "") not in ("", None)


def _tin_arms():
    arms = []
    if _tin_given("detin_callstats"):
        arms.append("ssnv")
    if all(_tin_given(k) for k in ("detin_seg", "detin_tumour_hets", "detin_normal_hets")):
        arms.append("ascna")
    if len(arms) == 2:
        arms.append("joint")
    return arms


_TIN_ARMS = _tin_arms()


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


# Copy-number inputs without the "@" header lines GATK writes; deTiN skips only "#" lines.
rule _tumour_in_normal_strip_headers:
    input:
        seg = CFG["inputs"]["detin_seg"] or [],
        tumour_hets = CFG["inputs"]["detin_tumour_hets"] or [],
        normal_hets = CFG["inputs"]["detin_normal_hets"] or []
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


# Population records at the call-stats positions.
rule _tumour_in_normal_population_sites:
    input:
        callstats = CFG["inputs"]["detin_callstats"] or [],
        vcf = CFG["inputs"]["detin_population_vcf"] or []
    output:
        positions = temp(CFG["dirs"]["population"] + _TIN_PAIR + "/positions.tsv"),
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
        awk -F '\\t' 'NR == 1 {{ for (i = 1; i <= NF; i++) col[$i] = i; next }}
            {{ print $col["contig"] "\\t" $col["position"] }}' {input.callstats}
            | sort -k1,1 -k2,2n -u > {output.positions}
        &&
        bcftools query -T {output.positions} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\n' {input.vcf}
            > {output.sites} 2> {log.stderr}
        """)


# deTiN's population pickle (build_exac_pickle layout), restricted to the pair's positions.
rule _tumour_in_normal_population_pickle:
    input:
        sites = str(rules._tumour_in_normal_population_sites.output.sites)
    output:
        pickle = CFG["dirs"]["population"] + _TIN_PAIR + "/population.pickle"
    params:
        min_af = _TIN_OPTS["population_min_af"]
    conda:
        CFG["conda_envs"]["detin"]
    threads:
        CFG["threads"]["population"]
    resources:
        **CFG["resources"]["population"]
    script:
        "src/population_pickle.py"


def _tin_path(key, wildcards):
    value = _TIN_INPUTS[key]
    return value(wildcards) if callable(value) else value.format(**wildcards)


def _tin_detin_inputs(wildcards):
    inputs = {"script": str(rules._tumour_in_normal_detin_source.output.script)}
    if wildcards.arm in ("ssnv", "joint"):
        inputs["callstats"] = _tin_path("detin_callstats", wildcards)
        inputs["pickle"] = str(rules._tumour_in_normal_population_pickle.output.pickle)
    if wildcards.arm in ("ascna", "joint"):
        inputs["seg"] = str(rules._tumour_in_normal_strip_headers.output.seg)
        inputs["tumour_hets"] = str(rules._tumour_in_normal_strip_headers.output.tumour_hets)
        inputs["normal_hets"] = str(rules._tumour_in_normal_strip_headers.output.normal_hets)
    return inputs


def _tin_detin_args(wildcards, input):
    args = []
    if "callstats" in input.keys():
        args += ["--mutation_data_path", input.callstats, "--exac_data_path", input.pickle]
    if "seg" in input.keys():
        args += ["--cn_data_path", input.seg,
                 "--tumor_het_data_path", input.tumour_hets,
                 "--normal_het_data_path", input.normal_hets]
    return " ".join(args)


# deTiN as published, once per arm: SNVs only, copy number only, and both (its joint estimate).
rule _tumour_in_normal_detin:
    input:
        unpack(_tin_detin_inputs)
    output:
        tin = CFG["dirs"]["detin"] + _TIN_PAIR + "/{arm}/{tumour_id}--{normal_id}.TiN_estimate.txt",
        ci = CFG["dirs"]["detin"] + _TIN_PAIR + "/{arm}/{tumour_id}--{normal_id}.TiN_estimate_CI.txt",
        added = CFG["dirs"]["detin"] + _TIN_PAIR + "/{arm}/{tumour_id}--{normal_id}.number_of_SSNVs_added.txt"
    log:
        stdout = CFG["logs"]["detin"] + _TIN_PAIR + "/{arm}.stdout.log",
        stderr = CFG["logs"]["detin"] + _TIN_PAIR + "/{arm}.stderr.log"
    params:
        args = _tin_detin_args,
        opts = _TIN_OPTS["detin"],
        opts_seq_type = lambda w: _TIN_OPTS["detin_seq_type"].get(w.seq_type, ""),
        outdir = CFG["dirs"]["detin"] + _TIN_PAIR + "/{arm}"
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
            {params.opts} {params.opts_seq_type}
            > {log.stdout} 2> {log.stderr}
        """)


# One row per arm: the estimate, its interval and the SNVs deTiN would recover.
rule _tumour_in_normal_detin_summary:
    input:
        tin = expand(str(rules._tumour_in_normal_detin.output.tin), arm = _TIN_ARMS, allow_missing = True),
        ci = expand(str(rules._tumour_in_normal_detin.output.ci), arm = _TIN_ARMS, allow_missing = True),
        added = expand(str(rules._tumour_in_normal_detin.output.added), arm = _TIN_ARMS, allow_missing = True)
    output:
        tsv = CFG["dirs"]["detin"] + _TIN_PAIR + "/detin.tsv"
    run:
        with open(output.tsv, "w") as out:
            out.write("tumour_id\tnormal_id\tseq_type\tgenome_build\tarm\tTiN\tTiN_CI_low\tTiN_CI_high\tn_SSNVs_added\n")
            for arm, tin, ci, added in zip(_TIN_ARMS, input.tin, input.ci, input.added):
                low, high = [x.strip() for x in open(ci).read().split(" - ")]
                out.write("\t".join([wildcards.tumour_id, wildcards.normal_id, wildcards.seq_type,
                                     wildcards.genome_build, arm, open(tin).read().strip(), low, high,
                                     open(added).read().strip()]) + "\n")


# TINC as published: Strelka PASS SNVs through TINC's own parser, Battenberg copy number when given.
rule _tumour_in_normal_tinc:
    input:
        vcf = CFG["inputs"]["tinc_vcf"] or [],
        cna = CFG["inputs"]["tinc_cna"] or []
    output:
        tsv = CFG["dirs"]["tinc"] + _TIN_PAIR + "/tinc.tsv",
        rds = CFG["dirs"]["tinc"] + _TIN_PAIR + "/tinc.rds",
        pdf = CFG["dirs"]["tinc"] + _TIN_PAIR + "/tinc_report.pdf"
    log:
        stderr = CFG["logs"]["tinc"] + _TIN_PAIR + "/tinc.stderr.log"
    params:
        fast = bool(_TIN_OPTS["tinc_fast"])
    conda:
        CFG["conda_envs"]["tinc"]
    threads:
        CFG["threads"]["tinc"]
    resources:
        **CFG["resources"]["tinc"]
    script:
        "src/tinc_run.R"


def _tin_output_source(wildcards):
    if wildcards.tool == "detin":
        return str(rules._tumour_in_normal_detin_summary.output.tsv)
    return str(rules._tumour_in_normal_tinc.output.tsv)


rule _tumour_in_normal_output:
    input:
        tsv = _tin_output_source
    output:
        tsv = CFG["dirs"]["outputs"] + "{tool}/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{tool}.tsv"
    wildcard_constraints:
        tool = "detin|tinc"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module = True)


def _tin_targets():
    runs = CFG["runs"]
    targets = []
    for tool in ("detin", "tinc"):
        if tool == "detin" and not _TIN_ARMS:
            continue
        if tool == "tinc" and not _tin_given("tinc_vcf"):
            continue
        keep = runs["tumour_seq_type"].isin(CFG["tools"].get(tool, []))
        sub = runs[keep]
        targets += expand(str(rules._tumour_in_normal_output.output.tsv), zip,
                          tool = [tool] * len(sub),
                          seq_type = sub["tumour_seq_type"],
                          genome_build = sub["tumour_genome_build"],
                          tumour_id = sub["tumour_sample_id"],
                          normal_id = sub["normal_sample_id"],
                          pair_status = sub["pair_status"])
    return targets


rule _tumour_in_normal_all:
    input:
        _tin_targets()


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
