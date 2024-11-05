configfile: "config.yaml"

import os
def create_dirs(*dirs):
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
create_dirs(config["output_dir"], f"{config['output_dir']}/mapped_reads",
               f"{config['output_dir']}/epic2", f"{config['output_dir']}/macs3",
                f"{config['output_dir']}/peakBeast"
    )

# list output dir to see if subdirectories are created
print(os.listdir(config["output_dir"]))

# add symbolic link of genome fasta to output dir (if not already there)
# use new name genome.fasta
if not os.path.exists(f"{config['output_dir']}/genome.fasta"):
    # use absolute path
    genome_abs = os.path.abspath(config["genome_fasta"])
    os.symlink(genome_abs, f"{config['output_dir']}/genome.fasta")






rule all:
    input:
        f"{config['output_dir']}/epic2_all.default.bedgraph",
        f"{config['output_dir']}/epic2_all.bs2000.bedgraph",
        f"{config['output_dir']}/epic2_unique.default.bedgraph",
        f"{config['output_dir']}/epic2_unique.bs2000.bedgraph",
        f"{config['output_dir']}/chip_vs_input_all.bs200.bw",
        f"{config['output_dir']}/chip_vs_input_all.bs2000.bw",
        f"{config['output_dir']}/chip_vs_input_unique.bs200.bw",
        f"{config['output_dir']}/chip_vs_input_unique.bs2000.bw",
        f"{config['output_dir']}/macs3/macs3_all_peaks.broadPeak",
        f"{config['output_dir']}/macs3_all_peaks.bedgraph",
        f"{config['output_dir']}/macs3/macs3_unique_peaks.broadPeak",
        f"{config['output_dir']}/macs3_unique_peaks.bedgraph",
        f"{config['output_dir']}/input_coverage_bs2000.bw",
        f"{config['output_dir']}/chip_coverage_bs2000.bw",
        f"{config['output_dir']}/peakBeast/peakBeast_10k_norm.bw",
        f"{config['output_dir']}/summary_plot.png"



rule bowtie2_build:
    input:
        config["genome_fasta"]
    output:
        config["genome_fasta"] + ".bowtie2_index_check"
    conda: "envs/bowtie2.yaml"
    shell:
        """
        # check if {input}.bt2 OR {input}.bt2l files exists (one is enough)
        # if not, build index
        if [ -f {input}.1.bt2 ] || [ -f {input}.1.bt2l ]; then
            echo "Index already exists"
            touch {output}
            exit 0
        else
            echo "Building index"
            bowtie2-build {input} {input}
            if [ -f {input}.1.bt2 ] || [ -f {input}.1.bt2l ]; then
                touch {output}
            else
                echo "Index building failed"
                exit 1
            fi
        fi
        """

rule bowtie2_align:
    input:
        bowtie2_index_check=config["genome_fasta"] + ".bowtie2_index_check",
        fastq=lambda wildcards: config["samples"][wildcards.sample],
        genome=config["genome_fasta"]
    output:
        bam="{output_dir}/mapped_reads/{sample}.all.bam"
    conda: "envs/bowtie2.yaml"
    threads: workflow.cores * 1
    shell:
        """
        bowtie2 -p {threads} -x {input.genome} -U {input.fastq} | samtools view -Sb - > {output.bam}
        """

rule bowtie2_sort:
    input:
        "{output_dir}/mapped_reads/{sample}.all.bam"
    output:
        "{output_dir}/mapped_reads/{sample}.all.sorted.bam"
    conda: "envs/bowtie2.yaml"
    shell:
        """
        samtools sort -o {output} {input}
        """

rule bowtie2_index:
    # make index of all sorted bam files
    input:
        "{output_dir}/mapped_reads/{sample}.all.sorted.bam",
        "{output_dir}/mapped_reads/{sample}.unique.sorted.bam"
    output:
        "{output_dir}/mapped_reads/{sample}.all.sorted.bam.csi",
        "{output_dir}/mapped_reads/{sample}.unique.sorted.bam.csi"

    conda: "envs/bowtie2.yaml"
    shell:
        """
        samtools index {input[0]} -c 
        samtools index {input[1]} -c 
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


rule get_unique_mapped_bam:
    input:
        input=f"{config['output_dir']}/mapped_reads/input.all.sorted.bam",
        chip=f"{config['output_dir']}/mapped_reads/chip.all.sorted.bam"
    output:
        input=f"{config['output_dir']}/mapped_reads/input.unique.sorted.bam",
        chip=f"{config['output_dir']}/mapped_reads/chip.unique.sorted.bam"
    conda: "envs/sambamba.yaml"
    threads: workflow.cores * 0.8

    shell:
        """
        sambamba view -t {threads} -f bam -h -F "[XS] == null and not unmapped  and not duplicate" {input.input} > {output.input}
        sambamba view -t {threads} -f bam -h -F "[XS] == null and not unmapped  and not duplicate" {input.chip} > {output.chip}
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
        input_bam="{output_dir}/mapped_reads/input.all.sorted.bam",
        input_csi="{output_dir}/mapped_reads/input.all.sorted.bam.csi",
        chip_bam="{output_dir}/mapped_reads/chip.all.sorted.bam",
        chip_csi="{output_dir}/mapped_reads/chip.all.sorted.bam.csi",
        chromsizes=config["genome_fasta"] + ".chromsizes"

    output:
        epic2_default="{output_dir}/epic2/epic2_all.default.csv",
        epic2_bs2000="{output_dir}/epic2/epic2_all.bs2000.csv"
    conda: "envs/epic2.yaml"
    shell:
        """
        epic2 -t {input.chip_bam} -c {input.input_bam} -cs {input.chromsizes} > {output.epic2_default} 
        epic2 -t {input.chip_bam} -c {input.input_bam} -cs {input.chromsizes} -bin 2000 > {output.epic2_bs2000}
        """

rule run_epic2_on_unique:
    input:
        input_bam="{output_dir}/mapped_reads/input.unique.sorted.bam",
        input_csi="{output_dir}/mapped_reads/input.unique.sorted.bam.csi",
        chip_bam="{output_dir}/mapped_reads/chip.unique.sorted.bam",
        chip_csi="{output_dir}/mapped_reads/chip.unique.sorted.bam.csi",
        chromsizes=config["genome_fasta"] + ".chromsizes"

    output:
        epic2_default="{output_dir}/epic2/epic2_unique.default.csv",
        epic2_bs2000="{output_dir}/epic2/epic2_unique.bs2000.csv"
    conda: "envs/epic2.yaml"
    shell:
        """
        epic2 -t {input.chip_bam} -c {input.input_bam} -cs {input.chromsizes} > {output.epic2_default}
        epic2 -t {input.chip_bam} -c {input.input_bam} -cs {input.chromsizes} -bin 2000 > {output.epic2_bs2000}
        """


rule epic2csv_to_bedgraph:
    input:
        csv1="{output_dir}/epic2/epic2_all.default.csv",
        csv2="{output_dir}/epic2/epic2_all.bs2000.csv",
        csv3="{output_dir}/epic2/epic2_unique.default.csv",
        csv4="{output_dir}/epic2/epic2_unique.bs2000.csv"
    output:
        bedgraph1="{output_dir}/epic2_all.default.bedgraph",
        bedgraph2="{output_dir}/epic2_all.bs2000.bedgraph",
        bedgraph3="{output_dir}/epic2_unique.default.bedgraph",
        bedgraph4="{output_dir}/epic2_unique.bs2000.bedgraph"
    shell:
        """
        cut -f 1-3,10 {input.csv1} > {output.bedgraph1}
        cut -f 1-3,10 {input.csv2} > {output.bedgraph2}
        cut -f 1-3,10 {input.csv3} > {output.bedgraph3}
        cut -f 1-3,10 {input.csv4} > {output.bedgraph4}
        """


rule bamCompare:
    input:
        input=config["output_dir"] + "/mapped_reads/input.all.sorted.bam",
        input_csi=config["output_dir"] + "/mapped_reads/input.all.sorted.bam.csi",
        chip=config["output_dir"] + "/mapped_reads/chip.all.sorted.bam",
        chip_csi=config["output_dir"] + "/mapped_reads/chip.all.sorted.bam.csi"
    output:
        bw200=config["output_dir"] + "/chip_vs_input_all.bs200.bw",
        bw2000=config["output_dir"] + "/chip_vs_input_all.bs2000.bw"

    conda: "envs/deeptools.yaml"
    shell:
        """
        bamCompare -b1 {input.chip} -b2 {input.input} -o {output.bw200} --binSize 200 
        bamCompare -b1 {input.chip} -b2 {input.input} -o {output.bw2000} --binSize 2000 
        """

rule bamCompare_on_unique:
    input:
        input=config["output_dir"] + "/mapped_reads/input.unique.sorted.bam",
        input_csi=config["output_dir"] + "/mapped_reads/input.unique.sorted.bam.csi",
        chip=config["output_dir"] + "/mapped_reads/chip.unique.sorted.bam",
        chip_csi=config["output_dir"] + "/mapped_reads/chip.unique.sorted.bam.csi"
    output:
        bw200=config["output_dir"] + "/chip_vs_input_unique.bs200.bw",
        bw2000=config["output_dir"] + "/chip_vs_input_unique.bs2000.bw"
    conda: "envs/deeptools.yaml"
    shell:
        """
        bamCompare -b1 {input.chip} -b2 {input.input} -o {output.bw200} --binSize 200 
        bamCompare -b1 {input.chip} -b2 {input.input} -o {output.bw2000} --binSize 2000
        """

rule bamCoverage:
    input:
        input=config["output_dir"] + "/mapped_reads/input.all.sorted.bam",
        input_csi=config["output_dir"] + "/mapped_reads/input.all.sorted.bam.csi",
        chip=config["output_dir"] + "/mapped_reads/chip.all.sorted.bam",
        chip_csi=config["output_dir"] + "/mapped_reads/chip.all.sorted.bam.csi"
    output:
        input_bw2000=config["output_dir"] + "/input_coverage_bs2000.bw",
        chip_bw2000=config["output_dir"] + "/chip_coverage_bs2000.bw"
    conda: "envs/deeptools.yaml"
    threads: workflow.cores
    shell:
        """
        bamCoverage -b {input.input} -o {output.input_bw2000} --binSize 2000 -p {threads} 
        bamCoverage -b {input.chip} -o {output.chip_bw2000} --binSize 2000 -p {threads}
        """


rule macs3:
    input:
        chip=config["output_dir"] + "/mapped_reads/chip.all.sorted.bam",
        chip_csi=config["output_dir"] + "/mapped_reads/chip.all.sorted.bam.csi",
        input=config["output_dir"] + "/mapped_reads/input.all.sorted.bam",
        input_csi=config["output_dir"] + "/mapped_reads/input.all.sorted.bam.csi",
        genome_size=config["genome_fasta"] + ".size"
    output:
        macs3_peaks=config["output_dir"] + "/macs3/macs3_all_peaks.broadPeak",
        macs3_bedgraph=config["output_dir"] + "/macs3_all_peaks.bedgraph"
    conda: "envs/macs3.yaml"
    shell:
        """
        macs3_dir=$(dirname {output.macs3_peaks})
        GS=$(cat {input.genome_size})
        macs3 callpeak -t {input.chip} -c {input.input} -f BAM -g $GS -n macs3_all --outdir $macs3_dir --broad --nomodel --extsize 200
        # make bedgraph
        cut -f 1-3,7 {output.macs3_peaks} > {output.macs3_bedgraph}
        """

rule macs3_on_unique:
    input:
        chip=config["output_dir"] + "/mapped_reads/chip.unique.sorted.bam",
        chip_csi=config["output_dir"] + "/mapped_reads/chip.unique.sorted.bam.csi",
        input=config["output_dir"] + "/mapped_reads/input.unique.sorted.bam",
        input_csi=config["output_dir"] + "/mapped_reads/input.unique.sorted.bam.csi",
        genome_size=config["genome_fasta"] + ".size",
    output:
        macs3_peaks=config["output_dir"] + "/macs3/macs3_unique_peaks.broadPeak",
        macs3_bedgraph=config["output_dir"] + "/macs3_unique_peaks.bedgraph"
    conda: "envs/macs3.yaml"
    shell:
        """
        macs3_dir=$(dirname {output.macs3_peaks})
        GS=$(cat {input.genome_size})
        macs3 callpeak -t {input.chip} -c {input.input} -f BAM -g $GS -n macs3_unique --outdir $macs3_dir --broad --nomodel --extsize 200
        # make bedgraph
        cut -f 1-3,7 {output.macs3_peaks} > {output.macs3_bedgraph}
        """


rule peakBeast:
    input:
        chip=config["output_dir"] + "/mapped_reads/chip.all.sorted.bam",
        chip_csi=config["output_dir"] + "/mapped_reads/chip.all.sorted.bam.csi",
        input=config["output_dir"] + "/mapped_reads/input.all.sorted.bam",
        input_csi=config["output_dir"] + "/mapped_reads/input.all.sorted.bam.csi"
    output:
        bwn = config["output_dir"] + "/peakBeast/peakBeast_10k_norm.bw"
    params:
        prefix = config["output_dir"] + "/peakBeast/peakBeast",
        basedir = workflow.current_basedir
    conda: "envs/peakBeast.yaml"
    threads: workflow.cores
    shell:
        """
        # get absolute path of scripts directory
        scripts_dir={params.basedir}/scripts
        echo "scripts dir: $scripts_dir"
        echo "---------------------------------"
        export PATH=$scripts_dir:$PATH
        peakBeast.R --input {input.input} --chip {input.chip} --prefix {params.prefix} --threads {threads} --normalized_only -S 
        """


rule plot_summary:
    input:
        config["output_dir"]+"/chip_coverage_bs2000.bw",
        config["output_dir"]+"/input_coverage_bs2000.bw",
        config["output_dir"]+"/chip_vs_input_all.bs2000.bw",
        config["output_dir"]+"/chip_vs_input_unique.bs2000.bw",
        config["output_dir"]+"/epic2_all.bs2000.bedgraph",
        config["output_dir"]+"/epic2_unique.bs2000.bedgraph",
        config["output_dir"]+"/macs3_all_peaks.bedgraph",
        config["output_dir"]+"/macs3_unique_peaks.bedgraph",
        config["output_dir"]+"/peakBeast/peakBeast_10k_norm.bw",
        chrom_sizes = config["genome_fasta"] + ".chromsizes"
    output:
        config["output_dir"]+"/summary_plot.png"
    params:
        basedir = workflow.current_basedir,
        outdir = config["output_dir"]
    conda: "envs/peakBeast.yaml"
    threads: 1
    shell:
        """
        # get absolute path of scripts directory
        scripts_dir={params.basedir}/scripts
        echo "scripts dir: $scripts_dir"
        echo "---------------------------------"
        export PATH=$scripts_dir:$PATH
        plot_summary.R --dir {params.outdir} --chrom_sizes {input.chrom_sizes} --output {output}
        """
