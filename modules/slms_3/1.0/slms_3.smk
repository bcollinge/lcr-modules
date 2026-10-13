#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package oncopipe and inspect packages. 
import oncopipe as op
import inspect

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
from importlib.metadata import version as pkg_version
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

current_version = pkg_version("oncopipe")
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 


# Setup module and store module-specific configuration in `CFG_SLMS3`
# `CFG_SLMS3` is a shortcut to `config["lcr-modules"]["slms_3"]`
# This is not assigned to `CFG` as in other modules because the sub-modules here will use `CFG`. 
CFG_SLMS3 = op.setup_module(
    name = "slms_3",
    version = "1.0",
    subdirectories = ["inputs", "strelka_gnomad", "lofreq_gnomad", "strelka_lofreq_union", "sage_gnomad", "mutect2_depth_filt", "union"]
        + (["gnomad_lookup"] if config["lcr-modules"]["slms_3"].get("options", {}).get("gnomad_callside", False) else [])
        + ["outputs"],
)

# Ensure each of the submodule versions meets requirements. 

assert float(CFG_SLMS3["module_versions"]["manta"]) >= 2.0, (
    f"The current manta module version is {CFG_SLMS3['module_versions']['manta']}. "
    "SLMS-3 requires manta module version 2.0 or higher. "
)

assert float(CFG_SLMS3["module_versions"]["starfish"]) >= 2.0, (
    f"The current starfish module version is {CFG_SLMS3['module_versions']['starfish']}. "
    "SLMS-3 requires starfish module version 2.0 or higher. "
)

assert float(CFG_SLMS3["module_versions"]["mutect2"]) >= 2.0, (
    f"The current mutect2 module version is {CFG_SLMS3['module_versions']['mutect2']}. "
    "SLMS-3 requires mutect2 module version 2.0 or higher. "
)

assert float(CFG_SLMS3["module_versions"]["strelka"]) >= 1.1, (
    f"The current strelka module version is {CFG_SLMS3['module_versions']['strelka']}. "
    "SLMS-3 requires strelka module version 1.1 or higher. "
)

# Define rules to be run locally when using a compute cluster
localrules:
    _slms_3_input_sage_vcf,
    _slms_3_input_strelka_vcf,
    _slms_3_input_mutect_vcf, 
    _slms_3_input_lofreq_vcf, 
    _slms_3_mutect2_samples_table, 
    _slms_3_output_vcf, 
    _slms_3_all,

##### MODULE CONFIGFILES #####

# Load default module configs and update with values from slms_3 config

for tool in CFG_SLMS3["module_versions"].keys():
    lcr_modules = config["lcr-modules"]["_shared"]["lcr-modules"] + "modules/"
    configfile: lcr_modules + "/" +  tool + "/" + CFG_SLMS3["module_versions"][tool] + "/config/default.yaml"
    snakemake.utils.update_config(
        config["lcr-modules"][tool], 
        config["lcr-modules"]["slms_3"][tool]
    )
    config["lcr-modules"][tool]["inputs"]["sample_bam"] = CFG_SLMS3["inputs"]["sample_bam"]
    config["lcr-modules"][tool]["inputs"]["sample_bai"] = CFG_SLMS3["inputs"]["sample_bai"]    


##### FIRST PASS VARIANT CALLING MODULE SNAKEFILES #####

# Load Manta first
include: "../../manta/" + CFG_SLMS3["module_versions"]["manta"] + "/manta.smk"

# Update Strelka config to use Manta Candidate Small Indels    
config["lcr-modules"]["strelka"]["inputs"]["candidate_small_indels"] = expand(str(rules._manta_output_vcf.output.vcf), vcf_name = "candidateSmallIndels", allow_missing=True)

# Load Strelka, SAGE, and Lofreq
include: "../../strelka/" + CFG_SLMS3["module_versions"]["strelka"] + "/strelka.smk"
include: "../../sage/" + CFG_SLMS3["module_versions"]["sage"] + "/sage.smk"
include: "../../lofreq/" + CFG_SLMS3["module_versions"]["lofreq"] + "/lofreq.smk"
 
##### FIRST PASS VARIANT CALLING RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')

rule _slms_3_input_strelka_vcf: 
    input:
        vcf = str(rules._strelka_output_filtered_vcf.output.vcf), 
        tbi = str(rules._strelka_output_filtered_vcf.output.vcf_tbi)
    output:
        vcf = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka.combined.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module = True)
        op.relative_symlink(input.tbi, output.tbi, in_module = True)

rule _slms_3_input_sage_vcf: 
    input:
        vcf = str(rules._sage_output_vcf.output.combined), 
        tbi = str(rules._sage_output_vcf.output.combined_tbi)
    output:
        vcf = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.combined.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module = True)
        op.relative_symlink(input.tbi, output.tbi, in_module = True)

rule _slms_3_input_lofreq_vcf: 
    input:
        vcf = str(rules._lofreq_output_vcf.output.vcf_all), 
    output:
        vcf = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lofreq.snvs.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lofreq.snvs.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module = True)
        op.relative_symlink(input.vcf + ".tbi", output.tbi, in_module = True)

# gnomAD sites VCF stem per genome build (options.gnomad_resource); unlisted builds use af-only-gnomad.
def _slms_3_gnomad(wildcards, stems=dict(CFG_SLMS3["options"].get("gnomad_resource", {}))):
    stem = stems.get(wildcards.genome_build, "af-only-gnomad")
    return reference_files(f"genomes/{wildcards.genome_build}/variation/{stem}.{wildcards.genome_build}.vcf.gz")

_SLMS3_CALLSIDE = CFG_SLMS3["options"].get("gnomad_callside", False)
_SLMS3_DROPPED = CFG_SLMS3["options"].get("gnomad_dropped", False)
_SLMS3_GATE_AF = str(CFG_SLMS3["options"].get("gnomad_gate_af", 0.0001))

# Call-side lookup (options.gnomad_callside): a split, left-aligned copy of the caller VCF (MNVs also split with
# options.gnomad_atomize) is looked up by exact allele in the normalized resource; INFO/AF is written back onto the
# original records, 0 when absent.
def _slms_3_lookup_vcf(wildcards):
    return {"strelka": str(rules._slms_3_input_strelka_vcf.output.vcf),
            "lofreq": str(rules._slms_3_input_lofreq_vcf.output.vcf),
            "sage": str(rules._slms_3_input_sage_vcf.output.vcf)}[wildcards.caller]

_SLMS3_LOOKUP_VCF = CFG_SLMS3["dirs"]["gnomad_lookup"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}.gnomad_af.vcf"

if _SLMS3_CALLSIDE:
    rule _slms_3_gnomad_lookup:
        input:
            vcf = _slms_3_lookup_vcf,
            gnomad = _slms_3_gnomad,
            fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
        output:
            vcf = _SLMS3_LOOKUP_VCF + ".gz" if _SLMS3_DROPPED else temp(_SLMS3_LOOKUP_VCF)
        log:
            stderr = CFG_SLMS3["logs"]["gnomad_lookup"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}.gnomad_lookup.stderr.log"
        wildcard_constraints:
            caller = "strelka|lofreq|sage"
        params:
            script = CFG_SLMS3["inputs"]["gnomad_lookup"],
            atomize = int(bool(CFG_SLMS3["options"].get("gnomad_atomize", False))),
            bgzip = "| bgzip" if _SLMS3_DROPPED else ""
        conda:
            CFG_SLMS3["conda_envs"]["bcftools_lookup"]
        container:
            CFG_SLMS3["container_envs"]["bcftools_lookup"]
        threads:
            CFG_SLMS3["threads"]["gnomad_lookup"]
        resources:
            **CFG_SLMS3["resources"]["gnomad_lookup"]
        shell:
            op.as_one_line("""
            bash {params.script} {input.vcf} {input.fasta} {input.gnomad} AF {params.atomize} {output.vcf}.tmp
            2> {log.stderr} {params.bgzip} > {output.vcf}
            """)

def _slms_3_af_input(caller):
    def _f(wildcards):
        if not _SLMS3_CALLSIDE:
            return []
        return expand(str(rules._slms_3_gnomad_lookup.output.vcf), caller = caller, **dict(wildcards.items()))
    return _f

# Writes the caller VCF with INFO/AF from the resource on stdout (AF=0 when absent).
def _slms_3_af_cmd(caller):
    def _f(wildcards, input, threads):
        if _SLMS3_CALLSIDE:
            return f"bcftools view {input.af}"
        return (f"bcftools annotate --threads {threads} -a {input.gnomad} -c INFO/AF {input.vcf} | "
                "awk 'BEGIN {FS=OFS=\"\\t\"} { if ($1 !~ /^#/ && $8 !~ \";AF=\") $8=$8\";AF=0\"; print $0; }'")
    return _f

# options.gnomad_dropped: each gnomAD filter also writes the records it removes, with the looked-up AF and the failed
# clause, to {caller}.gnomad_dropped.tsv.gz. Off, the outputs and commands are unchanged.
def _slms_3_dropped_out(subdir, caller):
    if not _SLMS3_DROPPED:
        return {}
    return {"dropped": CFG_SLMS3["dirs"][subdir] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/" + caller + ".gnomad_dropped.tsv.gz"}

# Saves the AF-annotated stream so the kept and the dropped records are both read from it.
def _slms_3_dropped_split(wildcards, output):
    if not _SLMS3_DROPPED:
        return ""
    return f"bcftools view -Ob -o {output.vcf}.af.bcf && bcftools view {output.vcf}.af.bcf |"

# Records failing `keep` from the saved stream: CHROM POS REF ALT AF failed (awk_body prints the row from the query
# columns, extra ones from `fields`).
def _slms_3_dropped_cmd(keep, awk_body, fields = ""):
    def _f(wildcards, input, output):
        if not _SLMS3_DROPPED:
            return ""
        return (f"&& {{ echo '# gnomAD resource: {input.gnomad}'; printf 'CHROM\\tPOS\\tREF\\tALT\\tAF\\tfailed\\n'; "
                f"bcftools view -e '{keep}' {output.vcf}.af.bcf | "
                f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF{fields}\\n' | "
                f"awk 'BEGIN {{FS=OFS=\"\\t\"}} {{{awk_body}}}'; }} | "
                f"gzip > {output.dropped} && rm {output.vcf}.af.bcf")
    return _f

_SLMS3_KEEP_AF = "INFO/AF < " + _SLMS3_GATE_AF
_SLMS3_KEEP_STRELKA = "INFO/AF < " + _SLMS3_GATE_AF + " && INFO/SomaticEVS >= 10 && FMT/DP[1] >= 10"
_SLMS3_PRINT_AF = 'print $1, $2, $3, $4, $5, "af"'
# Strelka: af when no ALT is under the cut; evs and dp from SomaticEVS and the tumour DP
_SLMS3_PRINT_STRELKA = ('split($7, d, ","); m = 1; n = split($5, a, ","); for (i = 1; i <= n; i++) if (a[i] + 0 < m) m = a[i] + 0; '
                        'f = (m >= ' + _SLMS3_GATE_AF + ') ? "af" : ""; if ($6 == "." || $6 + 0 < 10) f = f (f == "" ? "" : "+") "evs"; '
                        'if (d[2] == "." || d[2] + 0 < 10) f = f (f == "" ? "" : "+") "dp"; print $1, $2, $3, $4, $5, f')

# Annotate Strelka VCF and remove common GnomAD variants

rule _slms_3_annotate_strelka_gnomad:
    input:
        vcf = str(rules._slms_3_input_strelka_vcf.output.vcf),
        tbi = str(rules._slms_3_input_strelka_vcf.output.tbi),
        gnomad = _slms_3_gnomad,
        af = _slms_3_af_input("strelka")
    output:
        vcf = CFG_SLMS3["dirs"]["strelka_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka.combined.gnomad.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["strelka_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka.combined.gnomad.vcf.gz.tbi",
        **_slms_3_dropped_out("strelka_gnomad", "strelka")
    log:
        stderr = CFG_SLMS3["logs"]["strelka_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_gnomad.stderr.log"
    params:
        af_cmd = _slms_3_af_cmd("strelka"),
        split = _slms_3_dropped_split,
        dropped = _slms_3_dropped_cmd(_SLMS3_KEEP_STRELKA, _SLMS3_PRINT_STRELKA, "\\t%INFO/SomaticEVS\\t[%DP,]"),
        keep = _SLMS3_KEEP_STRELKA
    conda:
        CFG_SLMS3["conda_envs"]["bcftools"]
    container:
        CFG_SLMS3["container_envs"]["bcftools"]
    threads:
        CFG_SLMS3["threads"]["strelka_gnomad"]
    resources:
        **CFG_SLMS3["resources"]["strelka_gnomad"]
    shell:
        op.as_one_line("""
        {params.af_cmd} | {params.split}
        bcftools view -i '{params.keep}' -Oz -o {output.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        {params.dropped}
        """)


# Annotate LoFreq VCF and remove common GnomAD variants

rule _slms_3_annotate_lofreq_gnomad:
    input: 
        vcf = str(rules._slms_3_input_lofreq_vcf.output.vcf), 
        tbi = str(rules._slms_3_input_lofreq_vcf.output.tbi),
        gnomad = _slms_3_gnomad,
        af = _slms_3_af_input("lofreq")
    output: 
        vcf = CFG_SLMS3["dirs"]["lofreq_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq.snvs.gnomad.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["lofreq_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq.snvs.gnomad.vcf.gz.tbi",
        **_slms_3_dropped_out("lofreq_gnomad", "lofreq")
    log:
        stderr = CFG_SLMS3["logs"]["lofreq_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq_gnomad.stderr.log"
    params:
        af_cmd = _slms_3_af_cmd("lofreq"),
        split = _slms_3_dropped_split,
        dropped = _slms_3_dropped_cmd(_SLMS3_KEEP_AF, _SLMS3_PRINT_AF),
        keep = _SLMS3_KEEP_AF
    conda:
        CFG_SLMS3["conda_envs"]["bcftools"]
    container:
        CFG_SLMS3["container_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["lofreq_gnomad"]
    threads: 
        CFG_SLMS3["threads"]["lofreq_gnomad"]
    shell: 
        op.as_one_line("""
        {params.af_cmd} | {params.split}
        bcftools view -i '{params.keep}' -Oz -o {output.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        {params.dropped}
        """)

# Create a union of Strelka and Lofreq as candidates for Mutect2

rule _slms_3_strelka_lofreq_union: 
    input: 
        strelka = str(rules._slms_3_annotate_strelka_gnomad.output.vcf),
        strelka_tbi = str(rules._slms_3_annotate_strelka_gnomad.output.tbi), 
        lofreq = str(rules._slms_3_annotate_lofreq_gnomad.output.vcf), 
        lofreq_tbi = str(rules._slms_3_annotate_lofreq_gnomad.output.tbi)
    output: 
        vcf = CFG_SLMS3["dirs"]["strelka_lofreq_union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka_lofreq_union_gnomad.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["strelka_lofreq_union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka_lofreq_union_gnomad.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["strelka_lofreq_union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_lofreq_union_gnomad.stderr.log"
    conda:
        CFG_SLMS3["conda_envs"]["bcftools"]
    container:
        CFG_SLMS3["container_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["strelka_lofreq_union"]
    threads: 
        CFG_SLMS3["threads"]["strelka_lofreq_union"]
    shell: 
        op.as_one_line("""
        bcftools merge --force-samples -m all --threads {threads} 
        -Oz -o {output.vcf} {input.strelka} {input.lofreq} 2> {log.stderr} && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

# Modify the SAGE VCF file to use `NORMAL` and `TUMOR` as sample IDs
rule _slms_3_annotate_sage_gnomad: 
    input: 
        vcf = str(rules._slms_3_input_sage_vcf.output.vcf), 
        gnomad = _slms_3_gnomad,
        af = _slms_3_af_input("sage")
    output: 
        vcf = CFG_SLMS3["dirs"]["sage_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sage.renamed.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["sage_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sage.renamed.vcf.gz.tbi",
        **_slms_3_dropped_out("sage_gnomad", "sage")
    log:
        stderr = CFG_SLMS3["logs"]["sage_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_lofreq_union_gnomad.stderr.log"
    params:
        af_cmd = _slms_3_af_cmd("sage"),
        split = _slms_3_dropped_split,
        dropped = _slms_3_dropped_cmd(_SLMS3_KEEP_AF, _SLMS3_PRINT_AF),
        keep = _SLMS3_KEEP_AF
    conda:
        CFG_SLMS3["conda_envs"]["bcftools"]
    container:
        CFG_SLMS3["container_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["sage_gnomad"]
    threads: 
        CFG_SLMS3["threads"]["sage_gnomad"]
    shell: 
        op.as_one_line("""
        {params.af_cmd} | {params.split}
        perl -ne '$norm="{wildcards.normal_id}";
                  $tum="{wildcards.tumour_id}";
                  s/(\s)$tum(\s)/$1TUMOR$2/;
                  s/(\s)$norm(\s)/$1NORMAL$2/;
                  print;' |
        bcftools view -s "NORMAL,TUMOR" -i '{params.keep}' -Oz -o {output.vcf} 2> {log.stderr} 
        &&
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        {params.dropped}
        """)


##### SECOND PASS VARIANT CALLING MODULE SNAKEFILES #####

# Upate the Mutect2 config with the _strelka_lofreq_vcf

config["lcr-modules"]["mutect2"]["inputs"]["candidate_positions"] = str(rules._slms_3_strelka_lofreq_union.output.vcf)

include: "../../mutect2/" + CFG_SLMS3["module_versions"]["mutect2"] + "/mutect2.smk"

##### SECOND PASS VARIANT CALLING RULES #####

rule _slms_3_input_mutect_vcf: 
    input:
        vcf = str(rules._mutect2_output_vcf.output.vcf), 
        tbi = str(rules._mutect2_output_vcf.output.tbi)
    output:
        vcf = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.mutect2.combined.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.mutect2.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module = True)
        op.relative_symlink(input.tbi, output.tbi, in_module = True)


# Depth/VAF filter the Mutect2 output. This must be completed here and not in the mutect2 module so that the TUMOR column can be properly indexed. 
rule _slms_3_mutect2_samples_table: 
    output: 
        table = CFG_SLMS3["dirs"]["mutect2_depth_filt"] + "samples.txt"
    shell: 
        """echo "TUMOR" > {output.table}"""

rule _slms_3_mutect2_depth_filt: 
    input: 
        vcf = str(rules._slms_3_input_mutect_vcf.output.vcf), 
        table = str(rules._slms_3_mutect2_samples_table.output.table)
    output: 
        vcf = CFG_SLMS3["dirs"]["mutect2_depth_filt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.depthfilt.mutect2.combined.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["mutect2_depth_filt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.depthfilt.mutect2.combined.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["mutect2_depth_filt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_depth_filt.stderr.log"
    conda:
        CFG_SLMS3["conda_envs"]["bcftools"]
    container:
        CFG_SLMS3["container_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["mutect2_depth_filt"]
    threads: 
        CFG_SLMS3["threads"]["mutect2_depth_filt"]
    shell: 
        op.as_one_line("""
        bcftools view {input.vcf} | 
        perl -ne 'if(/^\#\#normal_sample=(.+)$/){{$norm=$1;}}if(/tumor_sample=(.+)$/){{$tum = $1;}}s/(\s)$tum(\s)/$1TUMOR$2/;s/(\s)$norm(\s)/$1NORMAL$2/;print;'|
        sed 's/##INFO=<ID=AS_FilterStatus,Number=A/##INFO=<ID=AS_FilterStatus,Number=1/' |   
        bcftools view  -s "NORMAL,TUMOR" -i 'FMT/DP[@{input.table}] >= 10 && FMT/AD[@{input.table}:1] >= 4 && FMT/AF[@{input.table}:0] >= 0.1' 
        -Oz -o {output.vcf} 2> {log.stderr} && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

# Update Starfish config to use outputs generated here

snakemake.utils.update_config(config["lcr-modules"]["starfish"], {
    "dirs": {
        "_parent": CFG_SLMS3["dirs"]["_parent"] + "starfish-" + CFG_SLMS3["module_versions"]["starfish"]
    },
    "inputs": {
        "names": 
            ["mutect", "sage", "strelka", "lofreq"], 
        "paths": 
            [
                str(rules._slms_3_mutect2_depth_filt.output.vcf),
                str(rules._slms_3_annotate_sage_gnomad.output.vcf),
                str(rules._slms_3_annotate_strelka_gnomad.output.vcf), 
                str(rules._slms_3_annotate_lofreq_gnomad.output.vcf)
            ]
    }
})

include: "../../starfish/" + CFG_SLMS3["module_versions"]["starfish"] + "/starfish.smk"

# Create a final union VCF file summarizing which variants were called by which variant callers
# First rename the sample columns in each VCF as TUMOR_{caller} or NORMAL_{caller}

rule _slms_3_rename_samples_all: 
    input: 
        vcf = lambda w: config["lcr-modules"]["starfish"]["inputs"]["vcf"][w.caller]
    output:
        vcf = CFG_SLMS3["dirs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}.tmp.vcf.gz",
        tbi = CFG_SLMS3["dirs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}.tmp.vcf.gz.tbi",
        samples = CFG_SLMS3["dirs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}.samples.txt"
    log:
        stderr = CFG_SLMS3["logs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/rename_samples_{caller}.stderr.log"
    conda:
        CFG_SLMS3["conda_envs"]["bcftools"]
    container:
        CFG_SLMS3["container_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["rename_all"]
    threads: 
        CFG_SLMS3["threads"]["rename_all"]
    shell: 
        op.as_one_line("""
        printf "TUMOR TUMOR_{wildcards.caller}\\nNORMAL NORMAL_{wildcards.caller}\\n" > {output.samples} 
        && 
        bcftools reheader -s {output.samples} -o {output.vcf} {input.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

rule _slms_3_union_vcf: 
    input: 
        vcf = ancient(expand(
            rules._slms_3_rename_samples_all.output.vcf,
            caller = config["lcr-modules"]["starfish"]["inputs"]["names"],
            allow_missing = True
        )),
        tbi = ancient(expand(
            rules._slms_3_rename_samples_all.output.tbi,
            caller = config["lcr-modules"]["starfish"]["inputs"]["names"],
            allow_missing = True
        ))
    output: 
        vcf = CFG_SLMS3["dirs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/union.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/union.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/union.stderr.log"
    conda:
        CFG_SLMS3["conda_envs"]["bcftools"]
    container:
        CFG_SLMS3["container_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["union_vcf"]
    threads: 
        CFG_SLMS3["threads"]["union_vcf"]
    shell: 
        op.as_one_line("""
        bcftools merge --threads {threads} -m both -Oz -o {output.vcf} {input.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)
    

# Symlinks the final output files into the module results directory (under '99-outputs/')


rule _slms_3_output_vcf:
    input:
        union_vcf = str(rules._slms_3_union_vcf.output.vcf), 
        union_tbi = str(rules._slms_3_union_vcf.output.tbi),
        dispatched = str(rules._starfish_dispatch.output)
    output:
        isec_vcf = CFG_SLMS3["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3.final.vcf.gz", 
        isec_tbi = CFG_SLMS3["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3.final.vcf.gz.tbi", 
        union_vcf = CFG_SLMS3["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3.union.vcf.gz",
        union_tbi = CFG_SLMS3["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3.union.vcf.gz.tbi"
    params: 
        vcf_in = config["lcr-modules"]["starfish"]["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.3+.vcf.gz", 
        rm_files = CFG_SLMS3["options"]["cleanup_vcfs"], 
        files_to_rm = [
            str(rules._slms_3_annotate_strelka_gnomad.output.vcf), 
            str(rules._slms_3_annotate_lofreq_gnomad.output.vcf), 
            str(rules._slms_3_strelka_lofreq_union.output.vcf), 
            str(rules._slms_3_mutect2_depth_filt.output.vcf), 
            str(rules._slms_3_annotate_sage_gnomad.output.vcf)
        ] 
    run:
        op.relative_symlink(params.vcf_in, output.isec_vcf, in_module = True)
        op.relative_symlink(params.vcf_in + ".tbi", output.isec_tbi, in_module = True)
        op.relative_symlink(input.union_vcf, output.union_vcf, in_module = True)
        op.relative_symlink(input.union_tbi, output.union_tbi, in_module = True)
        if params.rm_files: 
            for file in params.files_to_rm: 
                if os.path.exists(file): 
                    os.remove(file)




# Generates the target sentinels for each run, which generate the symlinks
rule _slms_3_all:
    input:
        rules._manta_all.input,
        rules._strelka_all.input, 
        rules._lofreq_all.input, 
        rules._sage_all.input,
        rules._mutect2_all.input,
        rules._starfish_all.input, 
        expand(
            [
                str(rules._slms_3_output_vcf.output.isec_vcf),
                str(rules._slms_3_output_vcf.output.union_vcf)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG_SLMS3["runs"]["tumour_seq_type"],
            genome_build=CFG_SLMS3["runs"]["tumour_genome_build"],
            tumour_id=CFG_SLMS3["runs"]["tumour_sample_id"],
            normal_id=CFG_SLMS3["runs"]["normal_sample_id"],
            pair_status=CFG_SLMS3["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
CFG=CFG_SLMS3
op.cleanup_module(CFG)
