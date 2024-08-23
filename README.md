# CHIP seq analysis pipeline


This pipeline is designed to analyze ChIP-seq data using a series of tools like epic2, Macs3 and deeptools. The pipeline is intended to be used for detection of broad peaks. The pipeline is encapsulated in a Singularity container.

## Requirements 
Singularity is required to use the container. Singularity can be installed using conda environment. 

```bash
conda create -n singularity3 -c conda-forge "singularity>=3.6"
conda activate singularity3
```

## Quick Start
Singularity image (.sif file) can be downloaded from https://github.com/kavonrtep/cenh3_chip_seq_pipeline/releases 


Input files for the pipeline are provided in `config.yaml` file. The file contains the following fields:

```yaml
genome_fasta: /path/to/genome.fasta
samples:
   input: /path/to/input.fastq.gz
   chip: /path/to/chip.fastq.gz
output_dir: /path/to/output
```


Input and ChIP could be `fastq` or `fastq.gz` files. To run the pipeline, execute the following command:

```bash
singularity run -B /path/to/ -B $PWD chipseq_pipeline.sif -c config.yaml -t 20
``` 
Parameter `-t` specifies the number of threads to use. Singularity parameter `-B` is used to bind the input and output directories to the container. Without this parameter, the container will not be able to access the input and output files. File `config.yaml` must be also in directory which is accessible to the container. In the example above this is the current directory `$PWD`. 

If you have files in different directories, you can specify multiple `-B` parameters. For example if your `config.yaml` is:
```yaml
genome_fasta: /mnt/data/genomes/genome.fasta
samples:
   input: /mnt/data/fastq_reads/input.fastq.gz
   chip: /mnt/data/fastq_reads/chip.fastq.gz
output_dir: /mnt/data/outputs/output
```

Then you can use:
```bash
singularity run -B /mnt/data/genomes \
-B /mnt/data/fastq_reads -B /mnt/data/outputs \
-B $PWD chipseq_pipeline.sif -c config.yaml -t 20
```
Or you can use the following command to bind the whole directory `/mnt/data` to the container
```bash
singularity run -B /mnt/data -B $PWD chipseq_pipeline.sif -c config.yaml -t 20
```


## Output structure
The pipeline will create the following files and folders in the output directory:

```bash
output/
├── epic2/
│   ├── epic2_all.bs2000.csv
│   ├── epic2_all.default.csv
│   ├── epic2_unique.bs2000.csv
│   └── epic2_unique.default.csv
├── macs3/
│   ├── macs3_all_peaks.broadPeak
│   ├── macs3_all_peaks.gappedPeak
│   ├── macs3_all_peaks.xls
│   ├── macs3_unique_peaks.broadPeak
│   ├── macs3_unique_peaks.gappedPeak
│   └── macs3_unique_peaks.xls
├── mapped_reads/
│   ├── chip.all.bam
│   ├── chip.all.sorted.bam
│   ├── chip.all.sorted.bam.csi
│   ├── chip.unique.sorted.bam
│   ├── chip.unique.sorted.bam.csi
│   ├── input.all.bam
│   ├── input.all.sorted.bam
│   ├── input.all.sorted.bam.csi
│   ├── input.unique.sorted.bam
│   └── input.unique.sorted.bam.csi
├── peakBeast/
│   ├── peakBeast_10k_norm.bw       # output from peakBest program - enrichment score in 10k windows
│   ├── peakBeast_2k_norm.bw        # output from peakBest program - enrichment score in 2k windows
│   └── ...
├── summary_plot.png                # graphical summary of ChIP-seq analysis
├── chip_vs_input_all.bs2000.bw     # output from bamCompare, bin witdh 2000bp
├── chip_vs_input_all.bs200.bw      # output from bamCompare, bin witdh 200bp
├── chip_vs_input_unique.bs2000.bw  # output from bamCompare on uniquely mapped reads, bin with 2000bp
├── chip_vs_input_unique.bs200.bw   # output from bamCompare on uniquely mapped reads, bin with 200bp
├── epic2_all.bs2000.bedgraph       # output from epic2, bin width 2000bp
├── epic2_all.default.bedgraph      # output from epic2, default bin width (200bp)
├── epic2_unique.bs2000.bedgraph    # output from epic2 on uniquely mapped reads, bin width 2000bp
├── epic2_unique.default.bedgraph   # output from epic2 on uniquely mapped reads, default bin width (200bp)
├── macs3_all_peaks.bedgraph        # output from macs3, broadPeaks as bedgraph
└── macs3_unique_peaks.bedgraph     # output from macs3 on uniquely mapped reads, broadPeaks as bedgraph

```


## Build the container

To build the container, run the following command:

```bash
SINGULARITY=`which singularity`
sudo $SINGULARITY build chipseq_pipeline_0.2.0.sif Singularity
```