from os.path import join

# default args
inputDIR = "Protocol/data/raw"
phredQuality = "20"
trimmedDIR = "Protocol/data/trimmed"
referenceDIR = "Protocol/data/removal/reference"
ensemblRelease = "101"
bowtie2IndexDIR = "Protocol/data/removal/index"

# env created at .snakemake/conda
env = "envs/protocolMeta.yaml"

IDs = glob_wildcards(join(inputDIR,
    "{id, [A-Za-z0-9]+}{suffix, (_[12])?}.fastq"))

rule all:
    input:
        expand("Protocol/data/trimmed/{id}_trim.fastq", id = IDs.id)

rule qualityControlSingle:
    input: join(inputDIR, "{id}.fastq")
    output:
        trimmed = "{trimmedDIR}/{id}_trim.fastq",
        html = "{trimmedDIR}/{id}_report.html",
        json = "{trimmedDIR}/{id}_report.json"
    conda: env
    threads: workflow.cores
    shell: "fastp -i {input} -o {output.trimmed} -q {phredQuality} \
        -w {threads} -h {output.html} -j {output.json}"

rule qualityControlPaired:
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
    shell: "fastp -i {input.f} -I {input.r} -o {output.f} \
        -O {output.r} -q {phredQuality} -w {threads} \
        --detect_adapter_for_pe -h {output.html} -j {output.json}"

rule downloadHumanPrimaryAssembly:
    output: "{referenceDIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    conda: env
    threads: workflow.cores
    shell: "wget -P {referenceDIR} \
    ftp://ftp.ensembl.org/pub/release-{ensemblRelease}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    && pigz -d {output} -p {threads}"

rule bowtie2BuildHumanIndex:
    input: "{referenceDIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output: "{bowtie2IndexDIR}/hostHS"
    conda: env
    threads: workflow.cores
    shell: "bowtie2-build {input} {output} --threads {threads}"
