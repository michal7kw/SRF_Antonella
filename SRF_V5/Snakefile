import os
from glob import glob

# Configuration
configfile: "config.yaml"

# NOTE: SES-V5ChIP-Seq1_S5 sample is excluded from analysis due to corrupted R1 fastq file
# Only SES-V5ChIP-Seq2_S6 is being analyzed as a representative sample
SAMPLES = ["SES-V5ChIP-Seq2_S6"]  # Manually specified due to data availability
INPUT_SAMPLE = "InputSES-V5ChIP-Seq_S2"  # The input control sample

rule all:
    input:
        # Quality control outputs
        expand("qc/fastqc/{sample}_R{read}_001_fastqc.html", sample=SAMPLES + [INPUT_SAMPLE], read=[1,2]),
        "qc/multiqc/multiqc_report.html",
        # Alignment outputs
        expand("alignment/{sample}.sorted.bam", sample=SAMPLES + [INPUT_SAMPLE]),
        expand("alignment/{sample}.sorted.bam.bai", sample=SAMPLES + [INPUT_SAMPLE]),
        # Peak calling outputs
        expand("peaks/{sample}_peaks.narrowPeak", sample=SAMPLES),
        # BigWig files for visualization
        expand("bigwig/{sample}.bw", sample=SAMPLES + [INPUT_SAMPLE])

rule fastqc:
    input:
        r1="u251_epigenomics/{sample}_R1_001.fastq.gz",
        r2="u251_epigenomics/{sample}_R2_001.fastq.gz"
    output:
        html1="qc/fastqc/{sample}_R1_001_fastqc.html",
        html2="qc/fastqc/{sample}_R2_001_fastqc.html"
    threads: 4
    resources:
        mem_mb=8000
    shell:
        """
        mkdir -p qc/fastqc
        fastqc -t {threads} -o qc/fastqc {input.r1} {input.r2}
        """

rule multiqc:
    input:
        expand("qc/fastqc/{sample}_R{read}_001_fastqc.html", sample=SAMPLES + [INPUT_SAMPLE], read=[1,2])
    output:
        "qc/multiqc/multiqc_report.html"
    resources:
        mem_mb=4000
    shell:
        """
        mkdir -p qc/multiqc
        multiqc qc/fastqc -o qc/multiqc
        """

rule align:
    input:
        r1="u251_epigenomics/{sample}_R1_001.fastq.gz",
        r2="u251_epigenomics/{sample}_R2_001.fastq.gz"
    output:
        bam="alignment/{sample}.sorted.bam",
        bai="alignment/{sample}.sorted.bam.bai"
    params:
        index=config["genome_index"],
        tmp_dir="alignment/tmp_{sample}"
    log:
        "logs/align/{sample}.log"
    threads: 64
    resources:
        mem_mb=64000,
        runtime=1440
    shell:
        """
        mkdir -p alignment
        mkdir -p {params.tmp_dir}
        
        # Run alignment with error logging
        (bowtie2 -p {threads} -x {params.index} \
            -1 {input.r1} -2 {input.r2} | \
            samtools sort -@ {threads} \
            -T {params.tmp_dir} \
            -m 3G \
            -o {output.bam}) 2> {log}
            
        # Index BAM file
        samtools index {output.bam}
        
        # Clean up temporary directory
        rm -rf {params.tmp_dir}
        """

rule call_peaks:
    input:
        treatment="alignment/{sample}.sorted.bam",
        control=f"alignment/{INPUT_SAMPLE}.sorted.bam"
    output:
        "peaks/{sample}_peaks.narrowPeak"
    params:
        genome_size=config["genome_size"]
    threads: 8
    resources:
        mem_mb=32000
    shell:
        """
        mkdir -p peaks
        macs2 callpeak -t {input.treatment} -c {input.control} \
            -f BAMPE -g {params.genome_size} -n {wildcards.sample} \
            --outdir peaks
        """

rule create_bigwig:
    input:
        bam="alignment/{sample}.sorted.bam",
        bai="alignment/{sample}.sorted.bam.bai"
    output:
        "bigwig/{sample}.bw"
    threads: 8
    resources:
        mem_mb=32000
    shell:
        """
        mkdir -p bigwig
        bamCoverage -b {input.bam} -o {output} \
            --binSize 10 --normalizeUsing RPKM \
            --numberOfProcessors {threads}
        """
