# CHIP seq analysis pipeline


This pipeline is designed to analyze ChIP-seq data using a series of tools like epic2, Macs3 and deeptools. The pipeline is intended to be used for detection of broad peaks. The pipeline is encapsulated in a Singularity container.

## Requirements 
Singularity is required to use the contained. Singularity can be installed using conda environment. 

```bash
conda install -c conda-forge singularity
```

## Quick Start

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


## Build the container

To build the container, run the following command:

```bash
SINGULARITY=`which singularity`
sudo $SINGULARITY build chipseq_pipeline.sif Singularity
```