inputDIR = "fastq"
outputDIR = "results"

IDS, = glob_wildcards(os.path.join(inputDIR, "{id}.txt"))
# env created at .snakemake/conda
env = "envs/protocolMeta.yaml"

rule all:
    input: expand(os.path.join(outputDIR, "{id}_resultado.txt"), id = IDS)

ruleorder: teste_paired > teste_single

rule teste_single:
    input: os.path.join(inputDIR, "{id}.txt")
    output: os.path.join(outputDIR, "{id}_resultado.txt")
    conda: env
    threads: workflow.cores
    shell: "echo {threads} > {output}"

rule teste_paired:
    input:
        r1 = os.path.join(inputDIR, "{id}_1.txt"),
        r2 = os.path.join(inputDIR, "{id}_2.txt"),
    output: os.path.join(outputDIR, "{id}_resultado.txt")
    conda: env
    shell: "fastp -v > {output}"