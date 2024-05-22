configfile: "config.yaml"

import os
def create_dirs(*dirs):
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
create_dirs(config["output_dir"], f"{config['output_dir']}/mapped_reads",
               f"{config['output_dir']}/epic2", f"{config['output_dir']}/macs3"
    )


rule all:
    input:
        f"{config['output_dir']}/epic2_default.bedgraph",
        f"{config['output_dir']}/epic2_bs2000.bedgraph",
        f"{config['output_dir']}/chip_vs_input_bs200.bw",
        f"{config['output_dir']}/chip_vs_input_bs2000.bw",
        f"{config['output_dir']}/macs3/macs3_peaks.broadPeak",
        f"{config['output_dir']}/macs3_peaks.bedgraph"

rule bowtie2_build:
    input:
        config["genome_fasta"]
    output:
        expand("{genome}.{ext}", genome=config["genome_fasta"], ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    conda: "envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build {input} {input}
        """

rule bowtie2_align:
    input:
        fasta=expand("{genome}.{ext}",genome=config["genome_fasta"],ext=["1.bt2"]),
        fastq=lambda wildcards: config["samples"][wildcards.sample],
        genome=config["genome_fasta"]
    output:
        bam="{output_dir}/mapped_reads/{sample}.bam"
    conda: "envs/bowtie2.yaml"
    threads: workflow.cores * 0.8
    shell:
        """
        bowtie2 -p {threads}  -x {input.genome} -U {input.fastq} | samtools view -Sb - > {output.bam}
        """

rule bowtie2_sort:
    input:
        "{output_dir}/mapped_reads/{sample}.bam"
    output:
        "{output_dir}/mapped_reads/{sample}.sorted.bam"
    conda: "envs/bowtie2.yaml"
    shell:
        """
        samtools sort -o {output} {input}
        """

rule bowtie2_index:
    input:
        "{output_dir}/mapped_reads/{sample}.sorted.bam"
    output:
        "{output_dir}/mapped_reads/{sample}.sorted.bam.csi"

    conda: "envs/bowtie2.yaml"
    shell:
        """
        samtools index {input} -c 2> log.txt
        """
rule make_genome_index:
    input:
        config["genome_fasta"]
    output:
        config["genome_fasta"] + ".fai"
    conda: "envs/bowtie2.yaml"
    shell:
        """
        samtools faidx {input}
        """
rule make_chromsizes:
    input:
        config["genome_fasta"] + ".fai"
    output:
        config["genome_fasta"] + ".chromsizes"
    shell:
        """
        cut -f1,2 {input} > {output}
        """


rule get_genome_size:
    input:
        config["genome_fasta"] + ".chromsizes"
    output:
        config["genome_fasta"] + ".size"
    run:
        with open(input[0]) as f, open(output[0], "w") as out:
            out.write(str(sum(int(line.strip().split()[1]) for line in f)))


rule run_epic2:
    input:
        input_bam="{output_dir}/mapped_reads/input.sorted.bam",
        input_csi="{output_dir}/mapped_reads/input.sorted.bam.csi",
        chip_bam="{output_dir}/mapped_reads/chip.sorted.bam",
        chip_csi="{output_dir}/mapped_reads/chip.sorted.bam.csi",
        chromsizes=config["genome_fasta"] + ".chromsizes"

    output:
        epic2_default="{output_dir}/epic2/epic2_default.csv",
        epic2_bs2000="{output_dir}/epic2/epic2_bs2000.csv"
    log:
        "{output_dir}/epic2.log"
    conda: "envs/epic2.yaml"
    shell:
        """
        epic2 -t {input.chip_bam} -c {input.input_bam} -cs {input.chromsizes} > {output.epic2_default} 2>>{log}
        epic2 -t {input.chip_bam} -c {input.input_bam} -cs {input.chromsizes} -bin 2000 > {output.epic2_bs2000} 2>> {log}
        """
rule epic2csv_to_bedgraph:
    input:
        csv1="{output_dir}/epic2/epic2_default.csv",
        csv2="{output_dir}/epic2/epic2_bs2000.csv"
    output:
        bedgraph1="{output_dir}/epic2_default.bedgraph",
        bedgraph2="{output_dir}/epic2_bs2000.bedgraph"
    shell:
        """
        cut -f 1-3,10 {input.csv1} > {output.bedgraph1}
        cut -f 1-3,10 {input.csv2} > {output.bedgraph2}
        """


rule bamCompare:
    input:
        input=config["output_dir"] + "/mapped_reads/input.sorted.bam",
        input_csi=config["output_dir"] + "/mapped_reads/input.sorted.bam.csi",
        chip=config["output_dir"] + "/mapped_reads/chip.sorted.bam",
        chip_csi=config["output_dir"] + "/mapped_reads/chip.sorted.bam.csi"
    output:
        bw200=config["output_dir"] + "/chip_vs_input_bs200.bw",
        bw2000=config["output_dir"] + "/chip_vs_input_bs2000.bw"
    log:
        config["output_dir"] + "/bamCompare.log"
    conda: "envs/deeptools.yaml"
    shell:
        """
        bamCompare -b1 {input.chip} -b2 {input.input} -o {output.bw200} --binSize 200 2>>{log}
        bamCompare -b1 {input.chip} -b2 {input.input} -o {output.bw2000} --binSize 2000 2>>{log}
        """
rule macs3:
    input:
        chip=config["output_dir"] + "/mapped_reads/chip.sorted.bam",
        chip_csi=config["output_dir"] + "/mapped_reads/chip.sorted.bam.csi",
        input=config["output_dir"] + "/mapped_reads/input.sorted.bam",
        input_csi=config["output_dir"] + "/mapped_reads/input.sorted.bam.csi",
        genome_size=config["genome_fasta"] + ".size"
    output:
        dir=directory(config["output_dir"] + "/macs3"),
        macs3_peaks=config["output_dir"] + "/macs3/macs3_peaks.broadPeak",
        macs3_bedgraph=config["output_dir"] + "/macs3_peaks.bedgraph"
    conda: "envs/macs3.yaml"
    log:
        config["output_dir"] + "/macs3.log"
    shell:
        """
        GS=$(cat {input.genome_size})
        macs3 callpeak -t {input.chip} -c {input.input} -f BAM -g $GS -n macs3 --outdir {output.dir} --broad --nomodel --extsize 200 2>> {log}
        # make bedgraph
        cut -f 1-3,7 {output.dir}/macs3_peaks.broadPeak > {output.macs3_bedgraph}
        """
