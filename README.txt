To use this pipeline you need to:

1 - Create the expected folder structure

Go to the folder containing the
Snakefile (via command line) and run:

mkdir -p ./Protocol/{data/{assembled,collapsed,removal/{index,reference},raw,trimmed},alignment/{db,index},taxonomic,functional}

2 - Install Snakemake

The recommended way to install Snakemake is via the conda
package manager. Read the installation guide presented at

https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

to get a conda environment with the latest Snakemake version.

3 - Move the input files to the expected location

Download/Move your raw fastq files to the inputDIR specified
in the Snakefile. By default, the inputDIR is setted as
"Protocol/data/raw". Alternatively, you may change the inputDIR
variable. It is worth to mention that all paths in the Snakefile
are interpreted relative to the directory Snakemake is executed
in.

4 - Start the pipeline

To start this pipeline go to the folder containing the
Snakefile (via command line) and run:

snakemake --cores --use-conda

This will use all available cores whenever is possible.
Alternatively, you may define the number of available cores.
