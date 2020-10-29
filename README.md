# To use this pipeline you need to:

## 1 - Create the expected folder structure

Go to the folder containing the Snakefile (via command line) and run:

```bash
mkdir -p ./Protocol/{data/{assembled,collapsed,removal/{index,reference},raw,trimmed},alignment/{db,index},taxonomic,functional}
```

## 2 - Install Snakemake

The recommended way to install Snakemake is via the conda package manager. Read the [installation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to get a conda environment with the latest Snakemake version.

## 3 - Move the input files to the expected location

Move your raw fastq files to the inputDIR specified in the Snakefile. By default, the inputDIR is "Protocol/data/raw" and paired-end filenames are expected to present the suffixes "_1.fastq" and "_2.fastq". Alternatively, you may change the inputDIR editing the Snakefile. It is worth to mention that all paths in the Snakefile are interpreted relative to the directory Snakemake is executed in.

## 4 - Run snakemake

To start this pipeline go to the folder containing the Snakefile (via command line) and run:

```bash
snakemake --cores --use-conda
```

This will use all available cores whenever is possible. Alternatively, you may define the number of available cores.
