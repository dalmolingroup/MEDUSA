# To execute the workflows from this repository you need to:

## 1 - Download this repository

Download this repository to get the files describing the pipeline rules:

```bash
wget https://github.com/arthurvinx/Medusa/archive/refs/heads/main.zip
unzip main.zip
cd Medusa-main/
```

## 2 - Create the expected folder structure

Go to the folder containing the Snakefile (via command line) and create the expected folder structure with:

```bash
mkdir -p ./Protocol/{data/{assembled,collapsed,removal/{index,reference},raw,trimmed},alignment/{db,index},taxonomic/db,functional/db}
```

## 3 - Install the conda package manager

An Anaconda environment was created to ease the installation of the software required by this pipeline. The environment recipe is available at the Anaconda cloud: https://anaconda.org/arthurvinx/medusapipeline.

A conda package manager, such as miniconda, must be installed to get this environment and to install the Snakemake workflow management system.

Download and install the [latest miniconda release](https://docs.conda.io/en/latest/miniconda.html) for Linux (adapt the commands if needed):

```bash
cd ~/Downloads
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh
chmod u+x Miniconda3-py39_4.9.2-Linux-x86_64.sh
./Miniconda3-py39_4.9.2-Linux-x86_64.sh
```

Close the terminal after the installation, open it again, and then check if the installation was successful with:

```bash
conda -V
```

## 4 - Get the pipeline environment

Get the pipeline environment from the Anaconda cloud:

```bash
conda activate base
conda install anaconda-client -y
conda env create arthurvinx/medusaPipeline
conda activate medusaPipeline
pip3 install -U plyvel --no-cache-dir --no-deps --force-reinstall
```

## 5 - Install Snakemake

The recommended way to install Snakemake is via the conda package manager. The following commands will create a conda environment with the full version of Snakemake. More details can be found at the [Snakemake installation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

```bash
conda activate base
conda install -c conda-forge mamba -y
mamba create -c bioconda -c conda-forge -n snakemake snakemake
```

Check whether your installation succeeded by typing:

```bash
conda activate snakemake
snakemake --help
```

## 6 - Move the input files to the expected location

Move your raw fastq files to the inputDIR specified in the Snakefile. By default, the inputDIR is "Protocol/data/raw" and paired-end filenames are expected to present the suffixes "_1.fastq" and "_2.fastq". Alternatively, you may change the inputDIR editing the Snakefile. It is worth to mention that all paths in the Snakefile are interpreted relative to the directory Snakemake is executed in.

## 7 - Run snakemake

To start this pipeline, go to the folder containing the Snakefile (via command line) and run:

```bash
snakemake --cores --use-conda
```

This will use all available cores whenever is possible. Alternatively, you may define the number of available cores.
