#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostia Dreval
# Module Author:    Kostia Dreval
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op

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
    print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
    print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["qc"]`
CFG = op.setup_module(
    name = "qc",
    version = "1.1",
    subdirectories = ["inputs", "samtools", "gatk", "aggregated_metrics", "contamination", "complexity", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _qc_input_bam,
    _qc_input_references,
    _qc_download_baits,
    _qc_output_tsv,
    _qc_output_contamination,
    _qc_output_complexity,
    _qc_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _qc_input_bam:
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

# symlink the reference files to ensure all index/dictionaries are available for GATK tools
rule _qc_input_references:
    input:
        genome_fa = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        genome_fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
        genome_dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output:
        genome_fa = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.fa",
        genome_fai = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.fa.fai",
        genome_dict = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.dict"
    shell:
        op.as_one_line("""
        ln -s {input.genome_fa} {output.genome_fa} &&
        ln -s {input.genome_fai} {output.genome_fai} &&
        ln -s {input.genome_dict} {output.genome_dict}
        """)


# Collect samtools stats
rule _qc_samtools_stat:
    input:
        bam = ancient(str(rules._qc_input_bam.output.bam))
    output:
        samtools_stat = CFG["dirs"]["samtools"] + "{seq_type}--{genome_build}/{sample_id}.{genome_build}.stat"
    log:
        stderr = CFG["logs"]["samtools"] + "{seq_type}--{genome_build}/{sample_id}.run_samtools.stderr.log"
    params:
        opts = CFG["options"]["samtools_stat"]
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        CFG["container_envs"]["samtools"]
    threads:
        CFG["threads"]["samtools_stat"]
    resources:
        **CFG["resources"]["samtools_stat"]
    shell:
        op.as_one_line("""
        samtools stats
        {params.opts}
        --threads {threads}
        {input.bam}
        >> {output.samtools_stat}
        2>> {log.stderr}
        """)

# Collect samtools per-chromosome coverage
rule _qc_samtools_coverage:
    input:
        bam = ancient(str(rules._qc_input_bam.output.bam))
    output:
        samtools_coverage = CFG["dirs"]["samtools"] + "{seq_type}--{genome_build}/{sample_id}.{genome_build}.coverage"
    log:
        stderr = CFG["logs"]["samtools"] + "{seq_type}--{genome_build}/{sample_id}.run_samtools_coverage.stderr.log"
    params:
        opts = CFG["options"]["samtools_coverage"]
    conda:
        CFG["conda_envs"]["samtools_cov"]
    container:
        CFG["container_envs"]["samtools_cov"]
    threads:
        CFG["threads"]["samtools_coverage"]
    resources:
        **CFG["resources"]["samtools_coverage"]
    shell:
        op.as_one_line("""
        samtools coverage
        {input.bam}
        -o {output.samtools_coverage}
        {params.opts}
        2>> {log.stderr}
        """)


# Collecting GATK base quality metrics
rule _qc_gatk_basequality:
    input:
        bam = ancient(str(rules._qc_input_bam.output.bam)),
        fasta = str(rules._qc_input_references.output.genome_fa)
    output:
        gatk_basequal = CFG["dirs"]["gatk"] + "QualityScoreDistribution/{seq_type}--{genome_build}/{sample_id}.{genome_build}.QualityScoreDistribution.txt",
        gatk_basequal_chart = CFG["dirs"]["gatk"] + "QualityScoreDistribution/{seq_type}--{genome_build}/{sample_id}.{genome_build}.QualityScoreDistribution.pdf"
    log:
        stdout = CFG["logs"]["gatk"] + "QualityScoreDistribution/{seq_type}--{genome_build}/{sample_id}.QualityScoreDistribution.stdout.log",
        stderr = CFG["logs"]["gatk"] + "QualityScoreDistribution/{seq_type}--{genome_build}/{sample_id}.QualityScoreDistribution.stderr.log"
    params:
        opts = CFG["options"]["QualityScoreDistribution"],
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatkR"]
    container:
        CFG["container_envs"]["gatkR"]
    threads:
        CFG["threads"]["QualityScoreDistribution"]
    resources:
        **CFG["resources"]["QualityScoreDistribution"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        gatk
        --java-options "-Xmx{params.jvmheap}m -XX:ConcGCThreads=1"
        QualityScoreDistribution
        {params.opts}
        -I {input.bam}
        -O {output.gatk_basequal}
        -CHART {output.gatk_basequal_chart}
        -R {input.fasta}
        >> {log.stdout}
        2>> {log.stderr} &&
        echo "DONE {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        """)


# Collect GATK QC metrics (separately for genomes/exomes)
rule _qc_gatk_wgs:
    input:
        bam = ancient(str(rules._qc_input_bam.output.bam)),
        fasta = str(rules._qc_input_references.output.genome_fa),
        samtools_stats = str(rules._qc_samtools_stat.output.samtools_stat)
    output:
        gatk_wgs = CFG["dirs"]["gatk"] + "CollectMetrics/{seq_type}--{genome_build}/{sample_id}.{genome_build}.CollectWgsMetrics.txt"
    log:
        stdout = CFG["logs"]["gatk"] + "CollectMetrics/{seq_type}--{genome_build}/{sample_id}.CollectWgsMetrics.stdout.log",
        stderr = CFG["logs"]["gatk"] + "CollectMetrics/{seq_type}--{genome_build}/{sample_id}.CollectWgsMetrics.stderr.log"
    wildcard_constraints:
        seq_type = "genome|lpwgs"
    params:
        opts = CFG["options"]["CollectWgsMetrics"],
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatkR"]
    container:
        CFG["container_envs"]["gatkR"]
    threads:
        CFG["threads"]["CollectMetrics"]
    resources:
        **CFG["resources"]["CollectMetrics"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        max_read_length=$(grep ^SN {input.samtools_stats} | cut -f 2- | grep "maximum length:" | cut -f 2); 
        gatk
        --java-options "-Xmx{params.jvmheap}m -XX:ConcGCThreads=1"
        CollectWgsMetrics
        {params.opts}
        -I {input.bam}
        -O {output.gatk_wgs}
        -R {input.fasta}
        --READ_LENGTH $((max_read_length + 5)) 
        >> {log.stdout}
        2>> {log.stderr} &&
        echo "DONE {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        """)


def _qc_get_baits(wildcards):
    CFG = config["lcr-modules"]["qc"]
    # simplify config options by re-using the capture space bertween e.g grch37 and hs37d5
    if "38" in str(wildcards.genome_build):
        this_genome_build = "hg38"
    else:
        this_genome_build = "grch37"

    if "baits_regions" in  CFG:
        # If user specified path to bed file in the sample table, use it from dictionary in config
        if wildcards.baits_regions in list(CFG["baits_regions"][this_genome_build].keys()):
            these_baits = str(wildcards.baits_regions)
        else:
            logger.warning(
                f"WARNING: the baits regions were specified in the sample table, but were not found in config "
                f"for {wildcards.baits_regions} ({this_genome_build}) and {wildcards.seq_type} combination. Using the default space ..."
            )
            these_baits = "_default"
    else:
        # use default
        these_baits = "_default"

    return(str(CFG["baits_regions"][this_genome_build][these_baits]))


# Downloads the baits BED file from a URL or copies it from a local path.
# Uses a run: block so it executes in the main Python process — no container
# needed, avoiding issues with curl/wget absence in bioinformatics containers.
rule _qc_download_baits:
    output:
        raw_baits = temp(CFG["dirs"]["inputs"] + "references/{genome_build}/raw/{baits_regions}.bed")
    wildcard_constraints:
        genome_build = "[^/]+",
        baits_regions = "[^/]+"
    params:
        baits = _qc_get_baits
    run:
        import urllib.request, shutil, os
        src = params.baits
        if os.path.isfile(src):
            shutil.copy(src, output.raw_baits)
        else:
            urllib.request.urlretrieve(src, output.raw_baits)


# Subset contigs in bed to those present in the reference and sort the output
# This is needed to pass GATK validation
rule _qc_sort_baits:
    input:
        fai = str(rules._qc_input_references.output.genome_fai),
        raw_baits = str(rules._qc_download_baits.output.raw_baits)
    output:
        baits = CFG["dirs"]["inputs"] + "references/{genome_build}/{baits_regions}.bed",
        intermediate_baits = temp(CFG["dirs"]["inputs"] + "references/{genome_build}/{baits_regions}.INTERMEDIATE.bed")
    wildcard_constraints:
        genome_build = "[^/]+",
        baits_regions = "[^/]+"
    conda:
        CFG["conda_envs"]["bedtools"]
    container:
        CFG["container_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        cut -f 1-3 {input.raw_baits} > {output.intermediate_baits}
            &&
        QC_REF_PREFIXED=$(head -1 {input.fai} | cut -f 1)
            &&
        BED_PREFIXED=$(head -1 {output.intermediate_baits} | cut -f 1)
            &&
        if [[ $QC_REF_PREFIXED == *"chr"* && ! $BED_PREFIXED == *"chr"* ]]; then
            awk '{{print "chr" $0}}' {output.intermediate_baits} > {output.intermediate_baits}.fix_prefix;
            rm {output.intermediate_baits};
            mv {output.intermediate_baits}.fix_prefix {output.intermediate_baits};
        fi
            &&
        cut -f 1-2  {input.fai} |
        awk '{{print $1 "\\t0\\t" $2}}' |
        bedtools intersect -wa -a {output.intermediate_baits} -b stdin |
        bedtools sort |
        bedtools merge >
        {output.baits}
        """)

# Create interval list
rule _qc_baits_to_interval_list:
    input:
        baits = str(rules._qc_sort_baits.output.baits),
        genome_dict = str(rules._qc_input_references.output.genome_dict)
    output:
        interval_list = CFG["dirs"]["inputs"] + "references/{genome_build}/{baits_regions}.interval_list"
    conda:
        CFG["conda_envs"]["gatkR"]
    container:
        CFG["container_envs"]["gatkR"]
    shell:
        op.as_one_line("""
        gatk BedToIntervalList
        --INPUT {input.baits}
        -SD {input.genome_dict}
        -O {output.interval_list}
        """)

# get the proper interval list corresponding to each sample, if provided
def _qc_get_intervals(wildcards):
    CFG = config["lcr-modules"]["qc"]
    samples = CFG["samples"]
    this_sample = samples.loc[(samples['sample_id'] == wildcards.sample_id) &
            (samples['genome_build'] == wildcards.genome_build) &
            (samples['seq_type'] == wildcards.seq_type)]

    # ensure the sample is uniquely mapped
    if len(this_sample) != 1:
        raise AssertionError("Found %s matches when examining the sample table for \'%s\' \'%s\' \'%s\'" % (len(sample), sample_id, genome_build, seq_type))

    # simplify config options by re-using the capture space bertween e.g grch37 and hs37d5
    if "38" in str(wildcards.genome_build):
        this_genome_build = "hg38"
    else:
        this_genome_build = "grch37"

    if "baits_regions" in  this_sample:
        # If user specified path to bed file in the sample table, use it from dictionary in config
        if str(this_sample.iloc[0]['baits_regions']) in list(CFG["baits_regions"][this_genome_build].keys()):
            these_baits = str(this_sample.iloc[0]['baits_regions'])
        elif str(this_sample.iloc[0]['baits_regions']) is None:
            these_baits = "_default"
        elif str(this_sample.iloc[0]['baits_regions']) in ["none"]:
            these_baits = "_default"
        else:
            logger.warning(
                f"WARNING: the baits regions were specified in the sample table, but were not found in config "
                f"for {str(this_sample.iloc[0]['sample_id'])} ({str(this_sample.iloc[0]['baits_regions'])}, {this_genome_build} and {wildcards.seq_type} combination). Using the default space ..."
            )
            these_baits = "_default"
    else:
        # use default
        these_baits = "_default"
    return(expand(str(rules._qc_baits_to_interval_list.output.interval_list), baits_regions=these_baits,allow_missing=True))

# Collect metrics on WES samples
rule _qc_gatk_wes:
    input:
        bam = ancient(str(rules._qc_input_bam.output.bam)),
        fasta = str(rules._qc_input_references.output.genome_fa),
        intervals = _qc_get_intervals
    output:
        gatk_wes = CFG["dirs"]["gatk"] + "CollectMetrics/{seq_type}--{genome_build}/{sample_id}.{genome_build}.CollectHsMetrics.txt"
    log:
        stdout = CFG["logs"]["gatk"] + "CollectMetrics/{seq_type}--{genome_build}/{sample_id}.CollectHsMetrics.stdout.log",
        stderr = CFG["logs"]["gatk"] + "CollectMetrics/{seq_type}--{genome_build}/{sample_id}.CollectHsMetrics.stderr.log"
    wildcard_constraints:
        seq_type = "capture"
    params:
        opts = CFG["options"]["CollectHsMetrics"],
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatkR"]
    container:
        CFG["container_envs"]["gatkR"]
    threads:
        CFG["threads"]["CollectMetrics"]
    resources:
        **CFG["resources"]["CollectMetrics"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        gatk
        --java-options "-Xmx{params.jvmheap}m -XX:ConcGCThreads=1"
        CollectHsMetrics
        {params.opts}
        -I {input.bam}
        -O {output.gatk_wes}
        -R {input.fasta}
        -BI {input.intervals}
        -TI {input.intervals}
        >> {log.stdout}
        2>> {log.stderr} &&
        echo "DONE {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        """)


def _qc_get_stats(wildcards):
    if wildcards.seq_type in ["capture"]:
        output = str(rules._qc_gatk_wes.output.gatk_wes)
    else:
        output = str(rules._qc_gatk_wgs.output.gatk_wgs)
    return output


# Collect required metrix into a tidy table
rule _qc_collect_metrics:
    input:
        stat = str(rules._qc_samtools_stat.output.samtools_stat),
        base_qual = str(rules._qc_gatk_basequality.output.gatk_basequal),
        metrics = _qc_get_stats
    output:
        stat = CFG["dirs"]["aggregated_metrics"] + "{seq_type}--{genome_build}/{sample_id}.metrix.tsv"
    conda:
        CFG["conda_envs"]["gatkR"]
    container:
        CFG["container_envs"]["Rscript"]
    threads:
        CFG["threads"]["collect"]
    resources:
        **CFG["resources"]["collect"]
    script:
        "src/R/aggregate_metrics.R"


# Combine all individual metrics into one single file
def _qc_merge_metrics_inputs(wildcards):
    # Restrict to samples of THIS seq_type before zipping. CFG["samples"] can be a mixed
    # genome/capture/lpwgs table; leaving seq_type as a free wildcard would request e.g.
    # capture--<build>/<sample>.metrix.tsv for genome-only samples (which don't exist ->
    # _qc_input_bam KeyError). Zipping seq_type too keeps each input fully resolved.
    # Re-fetch CFG from config: cleanup_module() deletes the module-level `CFG` name, so
    # at DAG-eval time (input functions) it must be read back off config (as at L216/297).
    CFG = config["lcr-modules"]["qc"]
    s = CFG["samples"][CFG["samples"]["seq_type"] == wildcards.seq_type]
    return expand(
        str(rules._qc_collect_metrics.output.stat),
        zip,
        seq_type=s["seq_type"],
        genome_build=s["genome_build"],
        sample_id=s["sample_id"])

rule _qc_merge_metrics:
    input:
        _qc_merge_metrics_inputs
    output:
        CFG["dirs"]["aggregated_metrics"] + "{seq_type}/master_qc_metrics_output.tsv"
    shell:
        op.as_one_line("""
        head -n1 {input[0]} > {output};
        tail -n+2 -q {input} >> {output};
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _qc_symlink_output:
    input:
        str(rules._qc_merge_metrics.output)
    output:
        CFG["dirs"]["outputs"] + "{seq_type}.qc_metrics.tsv",
    run:
        op.relative_symlink(input, output, in_module= True)

rule _qc_output_tsv:
    input:
        stat = str(rules._qc_samtools_stat.output.samtools_stat),
        coverage = str(rules._qc_samtools_coverage.output.samtools_coverage),
        base_qual = str(rules._qc_gatk_basequality.output.gatk_basequal),
        metrics = _qc_get_stats
    output:
        stat = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.stat.tsv",
        coverage = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.coverage.tsv",
        base_qual = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.QualityScoreDistribution.tsv",
        metrics = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.CollectMetrics.tsv",
    run:
        op.relative_symlink(input.stat, output.stat, in_module= True)
        op.relative_symlink(input.coverage, output.coverage, in_module= True)
        op.relative_symlink(input.base_qual, output.base_qual, in_module= True)
        op.relative_symlink(input.metrics, output.metrics, in_module= True)


# Common-SNP sites for GetPileupSummaries: CFG["contamination_sites"][genome_build][key], where key is the sample's
# baits_regions for capture samples and "_genome" otherwise. The VCF contigs must match the genome_build reference.
def _qc_get_contamination_sites(wildcards):
    CFG = config["lcr-modules"]["qc"]
    s = CFG["samples"]
    row = s.loc[(s["sample_id"] == wildcards.sample_id) & (s["genome_build"] == wildcards.genome_build) &
                (s["seq_type"] == wildcards.seq_type)]
    assert len(row) == 1, f"{len(row)} sample table rows for {wildcards.sample_id} {wildcards.seq_type} {wildcards.genome_build}"
    key = str(row.iloc[0]["baits_regions"]) if wildcards.seq_type == "capture" and "baits_regions" in row else "_genome"
    return CFG["contamination_sites"][wildcards.genome_build][key]



_QC_CONTAM = CFG["samples"][CFG["samples"]["seq_type"].isin(CFG["contamination_seq_types"])]

# Allele counts at common SNPs, for CalculateContamination (alone here, or against a matched normal downstream)
rule _qc_gatk_pileup_summaries:
    input:
        bam = ancient(str(rules._qc_input_bam.output.bam)),
        fasta = str(rules._qc_input_references.output.genome_fa),
        sites = ancient(_qc_get_contamination_sites)
    output:
        pileup = CFG["dirs"]["contamination"] + "{seq_type}--{genome_build}/{sample_id}.pileups.table"
    log:
        stdout = CFG["logs"]["contamination"] + "{seq_type}--{genome_build}/{sample_id}.pileups.stdout.log",
        stderr = CFG["logs"]["contamination"] + "{seq_type}--{genome_build}/{sample_id}.pileups.stderr.log"
    params:
        opts = CFG["options"]["GetPileupSummaries"],
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    container:
        CFG["container_envs"]["gatk"]
    threads:
        CFG["threads"]["contamination"]
    resources:
        **CFG["resources"]["contamination"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.jvmheap}m" GetPileupSummaries {params.opts}
        -I {input.bam} -R {input.fasta} -V {input.sites} -L {input.sites} -O {output.pileup}
        > {log.stdout} 2> {log.stderr}
        """)


# Contamination of each sample on its own pileups
rule _qc_gatk_calc_contamination:
    input:
        pileup = str(rules._qc_gatk_pileup_summaries.output.pileup)
    output:
        contamination = CFG["dirs"]["contamination"] + "{seq_type}--{genome_build}/{sample_id}.contamination.table",
        segments = CFG["dirs"]["contamination"] + "{seq_type}--{genome_build}/{sample_id}.segments.table"
    log:
        stdout = CFG["logs"]["contamination"] + "{seq_type}--{genome_build}/{sample_id}.contamination.stdout.log",
        stderr = CFG["logs"]["contamination"] + "{seq_type}--{genome_build}/{sample_id}.contamination.stderr.log"
    params:
        opts = CFG["options"]["CalculateContamination"],
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    container:
        CFG["container_envs"]["gatk"]
    threads:
        CFG["threads"]["contamination"]
    resources:
        **CFG["resources"]["contamination"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.jvmheap}m" CalculateContamination {params.opts}
        -I {input.pileup} -O {output.contamination} --tumor-segmentation {output.segments}
        > {log.stdout} 2> {log.stderr}
        """)


rule _qc_output_contamination:
    input:
        pileup = str(rules._qc_gatk_pileup_summaries.output.pileup),
        contamination = str(rules._qc_gatk_calc_contamination.output.contamination)
    output:
        pileup = CFG["dirs"]["outputs"] + "contamination/{seq_type}--{genome_build}/{sample_id}.pileups.table",
        contamination = CFG["dirs"]["outputs"] + "contamination/{seq_type}--{genome_build}/{sample_id}.contamination.table"
    run:
        op.relative_symlink(input.pileup, output.pileup, in_module= True)
        op.relative_symlink(input.contamination, output.contamination, in_module= True)


_QC_COMPLEX = CFG["samples"][CFG["samples"]["seq_type"].isin(CFG["complexity_seq_types"])]


# Read length and paired share from the first 200,000 primary reads
rule _qc_read_profile:
    input:
        bam = ancient(str(rules._qc_input_bam.output.bam)),
        fasta = str(rules._qc_input_references.output.genome_fa)
    output:
        profile = CFG["dirs"]["complexity"] + "{seq_type}--{genome_build}/{sample_id}.read_profile.tsv"
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        CFG["container_envs"]["samtools"]
    shell:
        op.as_one_line("""
        (samtools view -T {input.fasta} -F 0x904 {input.bam} || true) | head -n 200000 |
        awk 'BEGIN {{ OFS = "\\t"; print "reads", "mean_read_length", "paired_share" }}
             {{ n++; l += length($10); if (and($2, 1)) p++ }} END {{ print n, (n ? l / n : "NA"), (n ? p / n : "NA") }}'
        > {output.profile}
        """)


# Coverage regions: the capture targets for capture samples, fixed-size bins otherwise
def _qc_mosdepth_by(wildcards):
    if wildcards.seq_type == "capture":
        return _qc_get_intervals(wildcards)[0].replace(".interval_list", ".bed")
    return []


# Coverage distributions and per-region means without ("dedup", default flags) and with ("withdup") duplicates
rule _qc_mosdepth:
    input:
        bam = ancient(str(rules._qc_input_bam.output.bam)),
        bai = ancient(str(rules._qc_input_bam.output.bai)),
        fasta = str(rules._qc_input_references.output.genome_fa),
        regions = _qc_mosdepth_by
    output:
        dist = CFG["dirs"]["complexity"] + "{seq_type}--{genome_build}/{sample_id}.{dups}.mosdepth.global.dist.txt",
        regions = CFG["dirs"]["complexity"] + "{seq_type}--{genome_build}/{sample_id}.{dups}.regions.bed.gz"
    log:
        CFG["logs"]["complexity"] + "{seq_type}--{genome_build}/{sample_id}.{dups}.mosdepth.log"
    wildcard_constraints:
        dups = "dedup|withdup"
    params:
        prefix = CFG["dirs"]["complexity"] + "{seq_type}--{genome_build}/{sample_id}.{dups}",
        flag = lambda w: "1796" if w.dups == "dedup" else "772",
        by = lambda w, input: input.regions if w.seq_type == "capture" else str(CFG["complexity_bin_size"]),
        opts = CFG["options"]["mosdepth"]
    conda:
        CFG["conda_envs"]["mosdepth"]
    container:
        CFG["container_envs"]["mosdepth"]
    threads:
        CFG["threads"]["mosdepth"]
    resources:
        **CFG["resources"]["mosdepth"]
    shell:
        op.as_one_line("""
        mosdepth {params.opts} -n -t {threads} --flag {params.flag} --by {params.by} --fasta {input.fasta}
        {params.prefix} {input.bam} > {log} 2>&1
        """)


# GC content and N count of each fixed-size bin of the genome
rule _qc_bin_gc:
    input:
        fasta = str(rules._qc_input_references.output.genome_fa),
        fai = str(rules._qc_input_references.output.genome_fai)
    output:
        gc = CFG["dirs"]["complexity"] + "bin_gc/{genome_build}.bins.gc.tsv"
    params:
        size = CFG["complexity_bin_size"]
    conda:
        CFG["conda_envs"]["bedtools"]
    container:
        CFG["container_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        bedtools makewindows -g {input.fai} -w {params.size} |
        bedtools nuc -fi {input.fasta} -bed stdin |
        awk 'BEGIN {{ OFS = "\\t"; print "chrom", "start", "end", "gc", "n_share" }} NR > 1 {{ print $1, $2, $3, $5, $10 / $12 }}'
        > {output.gc}
        """)


# GC bias over the whole genome, for every seq_type
rule _qc_gc_bias:
    input:
        bam = ancient(str(rules._qc_input_bam.output.bam)),
        fasta = str(rules._qc_input_references.output.genome_fa)
    output:
        detail = CFG["dirs"]["complexity"] + "{seq_type}--{genome_build}/{sample_id}.gc_bias.detail.txt",
        summary = CFG["dirs"]["complexity"] + "{seq_type}--{genome_build}/{sample_id}.gc_bias.summary.txt",
        chart = CFG["dirs"]["complexity"] + "{seq_type}--{genome_build}/{sample_id}.gc_bias.pdf"
    log:
        CFG["logs"]["complexity"] + "{seq_type}--{genome_build}/{sample_id}.gc_bias.log"
    params:
        opts = CFG["options"]["CollectGcBiasMetrics"],
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatkR"]
    container:
        CFG["container_envs"]["gatkR"]
    threads: 1
    resources:
        **CFG["resources"]["gc_bias"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.jvmheap}m" CollectGcBiasMetrics {params.opts}
        -I {input.bam} -R {input.fasta} -O {output.detail} -S {output.summary} -CHART {output.chart} > {log} 2>&1
        """)


# Complexity curve: expected distinct fragments (paired) or reads against reads sequenced. A library that the default
# fit rejects is refitted in defects mode (-D); a library neither fits gets a header-only curve.
rule _qc_preseq:
    input:
        bam = ancient(str(rules._qc_input_bam.output.bam)),
        fasta = str(rules._qc_input_references.output.genome_fa),
        profile = str(rules._qc_read_profile.output.profile)
    output:
        curve = CFG["dirs"]["complexity"] + "{seq_type}--{genome_build}/{sample_id}.preseq.lc_extrap.txt"
    log:
        CFG["logs"]["complexity"] + "{seq_type}--{genome_build}/{sample_id}.preseq.log"
    params:
        opts = CFG["options"]["preseq"]
    conda:
        CFG["conda_envs"]["preseq"]
    container:
        CFG["container_envs"]["preseq"]
    threads: 2
    resources:
        **CFG["resources"]["preseq"]
    shell:
        op.as_one_line("""
        P=$(awk 'NR == 2 && $3 != "NA" && $3 > 0.5 {{ print "-P" }}' {input.profile});
        run() {{ samtools view -u -T {input.fasta} -F 0x904 {input.bam} |
                 preseq lc_extrap {params.opts} -B $P "$@" -o {output.curve} /dev/stdin; }};
        run > {log} 2>&1 || {{ echo "refit in defects mode" >> {log}; run -D >> {log} 2>&1; }} ||
        {{ echo "no fit" >> {log}; printf "TOTAL_READS\\tEXPECTED_DISTINCT\\tLOWER_0.95CI\\tUPPER_0.95CI\\n" > {output.curve}; }}
        """)


def _qc_complexity_inputs(wildcards):
    i = {
        "profile": str(rules._qc_read_profile.output.profile),
        "dedup_dist": expand(str(rules._qc_mosdepth.output.dist), dups="dedup", allow_missing=True)[0],
        "dedup_regions": expand(str(rules._qc_mosdepth.output.regions), dups="dedup", allow_missing=True)[0],
        "withdup_dist": expand(str(rules._qc_mosdepth.output.dist), dups="withdup", allow_missing=True)[0],
        "withdup_regions": expand(str(rules._qc_mosdepth.output.regions), dups="withdup", allow_missing=True)[0],
        "gc_summary": str(rules._qc_gc_bias.output.summary),
        "preseq": str(rules._qc_preseq.output.curve)
    }
    if wildcards.seq_type != "capture":
        i["bin_gc"] = str(rules._qc_bin_gc.output.gc)
    return i


# One row of evenness, GC and complexity metrics per sample
rule _qc_complexity_metrics:
    input:
        unpack(_qc_complexity_inputs)
    output:
        tsv = CFG["dirs"]["complexity"] + "{seq_type}--{genome_build}/{sample_id}.complexity.tsv"
    params:
        fixed_reads = lambda w: CFG["preseq_fixed_reads"][w.seq_type]
    conda:
        CFG["conda_envs"]["gatkR"]
    container:
        CFG["container_envs"]["Rscript"]
    resources:
        **CFG["resources"]["collect"]
    script:
        "src/R/complexity_metrics.R"


rule _qc_output_complexity:
    input:
        tsv = str(rules._qc_complexity_metrics.output.tsv),
        gc_detail = str(rules._qc_gc_bias.output.detail),
        preseq = str(rules._qc_preseq.output.curve)
    output:
        tsv = CFG["dirs"]["outputs"] + "complexity/{seq_type}--{genome_build}/{sample_id}.complexity.tsv",
        gc_detail = CFG["dirs"]["outputs"] + "complexity/{seq_type}--{genome_build}/{sample_id}.gc_bias.detail.txt",
        preseq = CFG["dirs"]["outputs"] + "complexity/{seq_type}--{genome_build}/{sample_id}.preseq.lc_extrap.txt"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)
        op.relative_symlink(input.gc_detail, output.gc_detail, in_module= True)
        op.relative_symlink(input.preseq, output.preseq, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _qc_all:
    input:
        expand(
            [
                str(rules._qc_output_tsv.output.stat),
                str(rules._qc_output_tsv.output.coverage),
                str(rules._qc_output_tsv.output.base_qual),
                str(rules._qc_output_tsv.output.metrics)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"]),
        expand(
            [
                str(rules._qc_symlink_output.output)
            ],
            seq_type=list(CFG["samples"]["seq_type"].unique())),
        expand(
            [
                str(rules._qc_output_contamination.output.contamination)
            ],
            zip,
            seq_type=_QC_CONTAM["seq_type"],
            genome_build=_QC_CONTAM["genome_build"],
            sample_id=_QC_CONTAM["sample_id"]),
        expand(
            [
                str(rules._qc_output_complexity.output.tsv)
            ],
            zip,
            seq_type=_QC_COMPLEX["seq_type"],
            genome_build=_QC_COMPLEX["genome_build"],
            sample_id=_QC_COMPLEX["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
