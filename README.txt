To use this pipeline you need to:

1 - Install Snakemake

The recommended way to install Snakemake is via the conda
package manager. Read the installation guide presented at

https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

to get a conda environment with the latest Snakemake version.

2 - Move the input files to the expected location

Download/Move your raw fastq files to the inputDIR specified
in the Snakefile. By default, the inputDIR is setted as
"Protocol/data/raw". Alternatively, you may change the inputDIR
variable. It is worth to mention that all paths in the Snakefile
are interpreted relative to the directory Snakemake is executed
in.

3 - Start the pipeline

To start this pipeline (i) go to the folder containing the
Snakefile (via command line) and (ii) run:

snakemake --cores --use-conda

This will use all available cores whenever is possible.
Alternatively, you may define the number of available cores.
