from os.path import join

inputDIR = "Protocol/data/raw"

IDs = glob_wildcards(join(inputDIR, "{id, [A-Za-z0-9]+}{suffix, (_1)?|(_2)?}.fastq"))

# env created at .snakemake/conda
env = "envs/protocolMeta.yaml"

rule all:
    input:
        #expand("Protocol/data/collapsed/{id}.fasta", id = IDs.id)
        expand("{id}.fasta", id = IDs.id)

rule fastpSingle:
    input: join(inputDIR, "{id}.fastq")
    output: "{id}.fasta"
    conda: env
    #threads: workflow.cores
    shell: "touch {output}"

rule fastpPaired:
    input:
        f = join(inputDIR, "{id}_1.fastq"),
        r = join(inputDIR, "{id}_2.fastq"),
    output: "{id}.fasta"
    conda: env
    shell: "touch {output}"