from os.path import join

# default args
inputDIR = "Protocol/data/raw"
phredQuality = "20"
trimmedDIR = "Protocol/data/trimmed" 

IDs = glob_wildcards(join(inputDIR, "{id, [A-Za-z0-9]+}{suffix, (_[12])?}.fastq"))

# env created at .snakemake/conda
env = "envs/protocolMeta.yaml"

rule all:
    input:
        expand("Protocol/data/trimmed/{id}_trim.fastq", id = IDs.id)

rule singleQC:
    input: join(inputDIR, "{id}.fastq")
    output:
        trimmed = "{trimmedDIR}/{id}_trim.fastq",
        html = "{trimmedDIR}/{id}_report.html",
        json = "{trimmedDIR}/{id}_report.json"
    conda: env
    threads: workflow.cores
    shell: "fastp -i {input} -o {output.trimmed} -q {phredQuality} -w {threads} -h {output.html} -j {output.json}"

rule pairedQC:
    input:
        f = join(inputDIR, "{id}_1.fastq"),
        r = join(inputDIR, "{id}_2.fastq"),
        outputDIR = "Protocol/data/trimmed"
    output:
        f = "{trimmedDIR}/{id}_1_trim.fastq",
        r = "{trimmedDIR}/{id}_2_trim.fastq",
        html = "{trimmedDIR}/{id}_report.html",
        json = "{trimmedDIR}/{id}_report.json"
    conda: env
    threads: workflow.cores
    shell: "fastp -i {input.f} -I {input.r} -o {output.f} -O {output.r} -q {phredQuality} -w {threads} --detect_adapter_for_pe -h {output.html} -j {output.json}"