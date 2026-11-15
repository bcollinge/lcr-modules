#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     Brett Collinge


##### SETUP #####


import oncopipe as op

# `CFG` is a shortcut to `config["lcr-modules"]["somalier"]`
CFG = op.setup_module(
    name = "somalier",
    version = "1.0",
    subdirectories = ["inputs", "extract", "relate", "outputs"],
)

localrules:
    _somalier_input_bam,
    _somalier_output_relate,
    _somalier_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _somalier_input_bam:
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


def _somalier_get_sites(wildcards):
    CFG = config["lcr-modules"]["somalier"]
    return CFG["inputs"]["sites"][wildcards.genome_build]


# Allele counts at the somalier sites. The digest's sample name is "{seq_type}--{genome_build}--{sample_id}:<read-group SM>";
# a file carrying more than one read-group SM fails rather than yielding several digests.
rule _somalier_extract:
    input:
        bam = ancient(str(rules._somalier_input_bam.output.bam)),
        bai = ancient(str(rules._somalier_input_bam.output.bai)),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        sites = ancient(_somalier_get_sites)
    output:
        digest = CFG["dirs"]["extract"] + "{seq_type}--{genome_build}/{sample_id}.somalier"
    log:
        stderr = CFG["logs"]["extract"] + "{seq_type}--{genome_build}/{sample_id}.extract.stderr.log"
    params:
        opts = CFG["options"]["extract"],
        tmp = CFG["dirs"]["extract"] + "{seq_type}--{genome_build}/{sample_id}.tmp"
    conda:
        CFG["conda_envs"]["somalier"]
    container:
        CFG["container_envs"]["somalier"]
    threads:
        CFG["threads"]["extract"]
    resources:
        **CFG["resources"]["extract"]
    shell:
        op.as_one_line("""
        rm -rf {params.tmp} && mkdir -p {params.tmp} &&
        somalier extract {params.opts} -d {params.tmp} --sites {input.sites} -f {input.fasta}
        --sample-prefix "{wildcards.seq_type}--{wildcards.genome_build}--{wildcards.sample_id}:" {input.bam} 2> {log.stderr} &&
        n=$(find {params.tmp} -name "*.somalier" | wc -l) &&
        if [ "$n" -ne 1 ]; then echo "expected 1 digest, found $n" >> {log.stderr}; exit 1; fi &&
        mv "$(find {params.tmp} -name '*.somalier')" {output.digest} && rm -rf {params.tmp}
        """)


def _somalier_relate_inputs(wildcards):
    CFG = config["lcr-modules"]["somalier"]
    s = CFG["samples"]
    return expand(str(rules._somalier_extract.output.digest), zip,
                  seq_type=s["seq_type"], genome_build=s["genome_build"], sample_id=s["sample_id"])


# Relatedness, IBS0/IBS2 and per-sample stats (including X heterozygous sites and Y depth) over every digest in the run.
rule _somalier_relate:
    input:
        digests = _somalier_relate_inputs
    output:
        pairs = CFG["dirs"]["relate"] + "relate.pairs.tsv",
        samples = CFG["dirs"]["relate"] + "relate.samples.tsv",
        groups = CFG["dirs"]["relate"] + "relate.groups.tsv"
    log:
        stderr = CFG["logs"]["relate"] + "relate.stderr.log"
    params:
        opts = CFG["options"]["relate"],
        prefix = CFG["dirs"]["relate"] + "relate"
    conda:
        CFG["conda_envs"]["somalier"]
    container:
        CFG["container_envs"]["somalier"]
    threads:
        CFG["threads"]["relate"]
    resources:
        **CFG["resources"]["relate"]
    shell:
        "somalier relate {params.opts} -o {params.prefix} {input.digests} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _somalier_output_relate:
    input:
        pairs = str(rules._somalier_relate.output.pairs),
        samples = str(rules._somalier_relate.output.samples)
    output:
        pairs = CFG["dirs"]["outputs"] + "relate.pairs.tsv",
        samples = CFG["dirs"]["outputs"] + "relate.samples.tsv"
    run:
        op.relative_symlink(input.pairs, output.pairs, in_module= True)
        op.relative_symlink(input.samples, output.samples, in_module= True)


# Every digest is a target, so a sample added to the run is extracted and its newer digest re-runs relate
rule _somalier_all:
    input:
        expand(str(rules._somalier_extract.output.digest), zip,
               seq_type=CFG["samples"]["seq_type"], genome_build=CFG["samples"]["genome_build"],
               sample_id=CFG["samples"]["sample_id"]),
        str(rules._somalier_output_relate.output.pairs),
        str(rules._somalier_output_relate.output.samples)


##### CLEANUP #####


op.cleanup_module(CFG)
