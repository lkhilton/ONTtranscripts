# ---------------------------------------------------------------------------------------------- #
# GLOBAL
# ---------------------------------------------------------------------------------------------- #

### PACKAGES ###

import pandas as pd
import numpy as np

### CONFIGURATION ###

configfile: "config/config.yaml"

### VARIABLES ###

FASTQ_DIR = "data/fastq"
FAST5_DIR = "data/fast5"
METADATA="data/metadata/ONT_Metadata.txt"

wildcard_constraints:
    IX = "IX731[12]"


# Lists here will be used in Snakemake rules -----------------------------

MD = pd.read_csv(METADATA, sep="\t")
SAMPLE = MD['External_Identifier'].tolist()
SAMPLE_IX = (MD['External_Identifier'] + "_" + MD['IX']).tolist()
GENES = ["FCGR2C", "FCGR2B"]


rule all:
    input:
      expand("results/minimap2/{sample}.hg19a.bam", sample=SAMPLE),
      expand("results/minimap2/03-subsampled-fastq/{sample_IX}.index.complete", sample_IX=SAMPLE_IX),
      expand("results/nanopolish/05-phased-reads/{sample}.{gene}.phased.bam", sample=SAMPLE, gene=GENES),
      "results/nanopolish/06-phased-fasta/ONTtranscripts.fasta",
      expand("results/nanopolish/07-clustal/{gene}.clustal", gene=GENES),
      expand("results/nanopolish/08-translated/all.{gene}.protein.clustal", gene=GENES)


# ---------------------------------------------------------------------------------------------- #
# MINIMAP2
# ---------------------------------------------------------------------------------------------- #

### VARIABLES ###

MINIMAP_DIR = "results/minimap2"
MINIMAP_HG19 = f"{MINIMAP_DIR}/01-hg19"
MINIMAP_TX = f"{MINIMAP_DIR}/02-transcript-aligned"
MINIMAP_FQ = f"{MINIMAP_DIR}/03-subsampled-fastq"

rule minimap_to_hg19:
    input:
        fastq = f"{FASTQ_DIR}""/{sample}.fastq"
    output:
        bam = f"{MINIMAP_DIR}""/{sample}.hg19a.bam",
        idx = f"{MINIMAP_DIR}""/{sample}.hg19a.bam.idxstats"
    params:
        ref = config["reference"]["grch37"]["minimap"]
    threads: 4
    conda: "config/envs/minimap2.yaml"
    shell:
        'minimap2 -ax splice -t {threads} {params.ref} {input.fastq} | '
        'samtools view -b -@ {threads} | samtools sort -@ {threads} > {output.bam} && '
        'samtools index {output.bam} && '
        'samtools idxstats {output.bam} > {output.idx}'

rule minimap_to_transcripts:
    input:
        fastq = f"{FASTQ_DIR}""/{sample}.fastq"
    output:
        bam = f"{MINIMAP_TX}""/{sample}.transcripts.bam",
        idx = f"{MINIMAP_TX}""/{sample}.transcripts.bam.idxstats"
    params:
        ref = config["reference"]["grch37"]["transcripts"]
    threads: 4
    conda: "config/envs/minimap2.yaml"
    shell:
        'minimap2 -ax map-ont -t {threads} {params.ref} {input.fastq} | '
        'samtools view -b -m 1050 -s 0.1 -@ {threads} | samtools sort -@ {threads} > {output.bam} && '
        'samtools index {output.bam} && '
        'samtools idxstats {output.bam} > {output.idx}'

rule bam_to_fastq:
    input:
        bam = f"{MINIMAP_TX}""/{sample}.transcripts.bam"
    output:
        fastq = f"{MINIMAP_FQ}""/{sample}.subsampled.fastq"
    threads: 4
    conda: "config/envs/minimap2.yaml"
    shell:
        'samtools fastq {input.bam} | '
        'awk \'BEGIN {{OFS = \"\\n\"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 1050) {{print header, seq, qheader, qseq}} }}\' '
        '> {output.fastq}'

# ---------------------------------------------------------------------------------------------- #
# NANOPOLISH
# ---------------------------------------------------------------------------------------------- #

### VARIABLES ###

NANOPOLISH_DIR = "results/nanopolish"
NP_CONSENSUS = f"{NANOPOLISH_DIR}/01-consensus"
NP_VARIANTS = f"{NANOPOLISH_DIR}/02-variants"
NP_MERGED = f"{NANOPOLISH_DIR}/03-merged-variants"
NP_PHASED_VCF = f"{NANOPOLISH_DIR}/04-phased-vcf"
NP_PHASED = f"{NANOPOLISH_DIR}/05-phased-reads"
NP_PHASED_FASTA = f"{NANOPOLISH_DIR}/06-phased-fasta"
NP_CLUSTAL = f"{NANOPOLISH_DIR}/07-clustal"
NP_TRANSLATE = f"{NANOPOLISH_DIR}/08-translated"

rule nanopolish_index:
    input:
        fast5 = f"{FAST5_DIR}""/{IX}",
        fastq = f"{MINIMAP_FQ}""/{sample}.subsampled.fastq",
        seq_summary = "data/sequencing_reports/{IX}/sequencing_summary.txt"
    output:
        complete = f"{MINIMAP_FQ}""/{sample}_{IX}.index.complete"
    conda: "config/envs/nanopolish.yaml"
    shell:
        'nanopolish index -s {input.seq_summary} -d {input.fast5} {input.fastq} && '
        'touch {output.complete}'

rule nanopolish_consensus:
    input:
        bam = f"{MINIMAP_TX}""/{sample}.transcripts.bam",
        fastq = f"{MINIMAP_FQ}""/{sample}.subsampled.fastq"
    output:
        vcf = f"{NP_CONSENSUS}""/{sample}.{gene}.nanopolish.vcf"
    params:
        ref = config["reference"]["grch37"]["transcripts"]
    conda: "config/envs/nanopolish.yaml"
    threads: 8
    shell:
        'nanopolish variants -t {threads} --consensus --ploidy 2 -w {wildcards.gene}:1-2000 '
        '--reads {input.fastq} --bam {input.bam} --genome {params.ref} -o {output.vcf}'

rule nanopolish_variants:
    input:
        bam = f"{MINIMAP_TX}""/{sample}.transcripts.bam",
        fastq = f"{MINIMAP_FQ}""/{sample}.subsampled.fastq"
    output:
        vcf = f"{NP_VARIANTS}""/{sample}.{gene}.nanopolish.vcf"
    params:
        ref = config["reference"]["grch37"]["transcripts"]
    conda: "config/envs/nanopolish.yaml"
    threads: 8
    shell:
        'nanopolish variants -t {threads} --ploidy 2 -w {wildcards.gene}:1-2000 '
        '--reads {input.fastq} --bam {input.bam} --genome {params.ref} -o {output.vcf}'

rule split_bam:
    input:
        bam = f"{MINIMAP_TX}""/{sample}.transcripts.bam"
    output:
        bam = f"{MINIMAP_TX}""/{sample}.{gene}.transcripts.bam"
    conda: "config/envs/nanopolish.yaml"
    shell:
        'samtools view -b {input.bam} {wildcards.gene} > {output.bam} && '
        'samtools index {output.bam}'

rule merge_variants:
    input:
        vars = f"{NP_VARIANTS}""/{sample}.{gene}.nanopolish.vcf",
        cons = f"{NP_CONSENSUS}""/{sample}.{gene}.nanopolish.vcf"
    output:
        vcf = f"{NP_MERGED}""/{sample}.{gene}.merged.nanopolish.vcf"
    conda: "config/envs/vcftools.yaml"
    shell:
        'vcfcombine {input.vars} {input.cons} > {output.vcf}'

rule whatshap_gt:
    input:
        vcf = f"{NP_VARIANTS}""/{sample}.{gene}.nanopolish.vcf",
        bam = f"{MINIMAP_TX}""/{sample}.{gene}.transcripts.bam"
    output:
        vcf = f"{NP_PHASED_VCF}""/{sample}.{gene}.whatshap_gt.vcf"
    params:
        ref = config["reference"]["grch37"]["transcripts"]
    conda: "config/envs/whatshap.yaml"
    shell:
        'whatshap genotype --indels --ignore-read-groups -o {output.vcf} --reference {params.ref} {input.vcf} {input.bam}'

rule phase_vcf:
    input:
        vcf = f"{NP_PHASED_VCF}""/{sample}.{gene}.whatshap_gt.vcf",
        bam = f"{MINIMAP_TX}""/{sample}.{gene}.transcripts.bam"
    output:
        vcf = f"{NP_PHASED_VCF}""/{sample}.{gene}.phased.vcf"
    params:
        ref = config["reference"]["grch37"]["transcripts"]
    conda: "config/envs/whatshap.yaml"
    shell:
        "whatshap phase --indels --ignore-read-groups --reference {params.ref} -o {output.vcf} {input.vcf} {input.bam}"

rule phase_bam:
    input:
        vcf = f"{NP_PHASED_VCF}""/{sample}.{gene}.phased.vcf",
        bam = f"{MINIMAP_TX}""/{sample}.{gene}.transcripts.bam"
    output:
        vcf = f"{NP_PHASED_VCF}""/{sample}.{gene}.phased.vcf.gz",
        bam = f"{NP_PHASED}""/{sample}.{gene}.phased.bam"
    params:
        ref = config["reference"]["grch37"]["transcripts"]
    conda: "config/envs/whatshap.yaml"
    shell:
        'bgzip -c {input.vcf} > {output.vcf} && tabix {output.vcf} && '
        'whatshap haplotag --ignore-read-groups -o {output.bam} --reference {params.ref} {output.vcf} {input.bam} && '
        'samtools index {output.bam}'

rule consensus_fasta:
    input:
        vcf = f"{NP_PHASED_VCF}""/{sample}.{gene}.phased.vcf.gz"
    output:
        fasta = f"{NP_PHASED_FASTA}""/{sample}.{gene}.hap{hap}.fasta"
    params:
        ref = config["reference"]["grch37"]["transcripts"]
    conda: "config/envs/whatshap.yaml"
    shell:
        'samtools faidx {params.ref} {wildcards.gene} | '
        'bcftools consensus -H {wildcards.hap} {input.vcf}  | '
        'sed \'s|>|>{wildcards.sample}.hap{wildcards.hap}.|g\' > {output.fasta}'

rule concatenate_fasta:
    input:
        fasta = expand(f"{NP_PHASED_FASTA}""/{sample}.{gene}.hap{hap}.fasta", sample=SAMPLE, gene=GENES, hap=[1,2])
    output:
        fasta = f"{NP_PHASED_FASTA}""/ONTtranscripts.fasta"
    shell:
        'cat {input.fasta} >> {output.fasta}'

rule clustalo:
    input:
        fasta = expand(f"{NP_PHASED_FASTA}""/{sample}.{{gene}}.hap{hap}.fasta", sample=SAMPLE, hap=[1,2])
    output:
        clustal = f"{NP_CLUSTAL}""/{gene}.clustal"
    params:
        ref = config["reference"]["grch37"]["transcripts"]
    conda: "config/envs/clustalo.yaml"
    shell:
        'samtools faidx {params.ref} {wildcards.gene} | '
        'cat - {input.fasta} | '
        'clustalo -i - -o {output.clustal} --outfmt=clu '
        '--residuenumber --wrap=100 --output-order=input-order'

rule concatenate_by_gene:
    input:
        fasta = expand(f"{NP_PHASED_FASTA}""/{sample}.{{gene}}.hap{hap}.fasta", sample=SAMPLE, hap=[1,2])
    output:
        fasta = f"{NP_PHASED_FASTA}""/all.{gene}.fasta"
    shell:
        'cat {input.fasta} >> {output.fasta}'

rule translate:
    input:
        FCGR2B = f"{NP_PHASED_FASTA}""/all.FCGR2B.fasta",
        FCGR2C = f"{NP_PHASED_FASTA}""/all.FCGR2C.fasta"
    output:
        FCGR2B = f"{NP_TRANSLATE}""/all.FCGR2B.translated.fasta",
        FCGR2C = f"{NP_TRANSLATE}""/all.FCGR2C.translated.fasta"
    params:
        FCGR2B = "128-1200",
        FCGR2C = "100-1200"
    conda: "config/envs/emboss.yaml"
    shell:
        "transeq -region {params.FCGR2B} -sequence {input.FCGR2B} -outseq {output.FCGR2B} && "
        "transeq -region {params.FCGR2C} -sequence {input.FCGR2C} -outseq {output.FCGR2C}"

rule clustalo_protein:
    input:
        fasta = f"{NP_TRANSLATE}""/all.{gene}.translated.fasta"
    output:
        clustal = f"{NP_TRANSLATE}""/all.{gene}.protein.clustal"
    conda: "config/envs/clustalo.yaml"
    shell:
        'clustalo -i {input.fasta} -o {output.clustal} --outfmt=clu '
        '--residuenumber --wrap=50 --output-order=input-order'
