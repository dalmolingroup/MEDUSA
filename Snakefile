from os.path import join

# input args
preprocessingDIR = "Protocol/data"
alignmentDIR = "Protocol/alignment"
taxonomicDIR = "Protocol/taxonomic"
functionalDIR = "Protocol/functional"
inputDIR = join(preprocessingDIR, "raw")
phredQuality = "20"
trimmedDIR = join(preprocessingDIR, "trimmed")
removalDIR = join(preprocessingDIR, "removal")
referenceDIR = join(removalDIR, "reference")
ensemblRelease = "101"
bowtie2IndexDIR = join(removalDIR, "index")
assembledDIR = join(preprocessingDIR, "assembled")
collapsedDIR = join(preprocessingDIR, "collapsed")
NRDIR = join(alignmentDIR, "db")
diamondIndexDIR = join(alignmentDIR, "index")

# env created at .snakemake/conda
env = "envs/protocolMeta.yaml"

IDs = glob_wildcards(join(inputDIR,
    "{id, [A-Za-z0-9]+}{suffix, (_[1-2])?}.fastq"))

ruleorder: qualityControlPaired > qualityControlSingle
ruleorder: removeHumanContaminantsPaired > removeHumanContaminantsSingle
ruleorder: deduplicatePaired > deduplicateSingle

rule all:
    input:
        expand(join(taxonomicDIR, "{id}_tax.txt"), id = IDs.id),
        expand(join(functionalDIR, "{id}_functional.txt"), id = IDs.id)

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
    params: indexPrefix = join(bowtie2IndexDIR, "hostHS")
    input: join(referenceDIR, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    output:
        expand(join(bowtie2IndexDIR, "hostHS.{index}.bt2l"), index = range(1, 5)),
        expand(join(bowtie2IndexDIR, "hostHS.rev.{index}.bt2l"), index = range(1, 3))
    conda: env
    threads: workflow.cores
    shell: "bowtie2-build --large-index {input} {params.indexPrefix} --threads {threads}"

rule removeHumanContaminantsSingle:
    params: indexPrefix = join(bowtie2IndexDIR, "hostHS"),
    input:
        expand(join(bowtie2IndexDIR, "hostHS.{index}.bt2l"), index = range(1, 5)),
        expand(join(bowtie2IndexDIR, "hostHS.rev.{index}.bt2l"), index = range(1, 3)),
        trimmed = join(trimmedDIR, "{id}_trim.fastq")
    output: "{removalDIR}/{id}_unaligned.fastq"
    conda: env
    threads: workflow.cores
    shell: "bowtie2 -x {params.indexPrefix} -U {input.trimmed} -S {removalDIR}/{wildcards.id}.sam -p {threads} \
        && samtools view -bS {removalDIR}/{wildcards.id}.sam > {removalDIR}/{wildcards.id}.bam \
        && samtools view -b -f 4 -F 256 {removalDIR}/{wildcards.id}.bam > {removalDIR}/{wildcards.id}_unaligned.bam \
        && samtools sort -n {removalDIR}/{wildcards.id}_unaligned.bam -o {removalDIR}/{wildcards.id}_unaligned_sorted.bam \
        && samtools bam2fq {removalDIR}/{wildcards.id}_unaligned_sorted.bam > {removalDIR}/{wildcards.id}_unaligned.fastq"

rule removeHumanContaminantsPaired:
    params: indexPrefix = join(bowtie2IndexDIR, "hostHS"),
    input:
        expand(join(bowtie2IndexDIR, "hostHS.{index}.bt2l"), index = range(1, 5)),
        expand(join(bowtie2IndexDIR, "hostHS.rev.{index}.bt2l"), index = range(1, 3)),
        f = join(trimmedDIR, "{id}_1_trim.fastq"),
        r = join(trimmedDIR, "{id}_2_trim.fastq")
    output:
        f = "{removalDIR}/{id}_unaligned_1.fastq",
        r = "{removalDIR}/{id}_unaligned_2.fastq"
    conda: env
    threads: workflow.cores
    shell: "bowtie2 -x {params.indexPrefix} -1 {input.f} -2 {input.r} -S {removalDIR}/{wildcards.id}.sam -p {threads} \
        && samtools view -bS {removalDIR}/{wildcards.id}.sam > {removalDIR}/{wildcards.id}.bam \
        && samtools view -b -f 12 -F 256 {removalDIR}/{wildcards.id}.bam > {removalDIR}/{wildcards.id}_unaligned.bam \
        && samtools sort -n {removalDIR}/{wildcards.id}_unaligned.bam -o {removalDIR}/{wildcards.id}_unaligned_sorted.bam \
        && samtools bam2fq {removalDIR}/{wildcards.id}_unaligned_sorted.bam > {removalDIR}/{wildcards.id}_unaligned.fastq \
        && cat {removalDIR}/{wildcards.id}_unaligned.fastq | grep '^@.*/1$' -A 3 --no-group-separator > {removalDIR}/{wildcards.id}_unaligned_1.fastq \
        && cat {removalDIR}/{wildcards.id}_unaligned.fastq | grep '^@.*/2$' -A 3 --no-group-separator > {removalDIR}/{wildcards.id}_unaligned_2.fastq"

rule mergePaired:
    input:
        f = join(removalDIR, "{id}_unaligned_1.fastq"),
        r = join(removalDIR, "{id}_unaligned_2.fastq")
    output:
        merged = "{assembledDIR}/{id}_assembled.fastq",
        html = "{assembledDIR}/{id}_report2.html",
        json = "{assembledDIR}/{id}_report2.json"
    conda: env
    threads: workflow.cores
    shell: "fastp -i {input.f} -I {input.r} -o {assembledDIR}/{wildcards.id}_unassembled_1.fastq \
        -O {assembledDIR}/{wildcards.id}_unassembled_2.fastq -q {phredQuality} -w {threads} \
        --detect_adapter_for_pe -h {output.html} -j {output.json} \
        -m --merged_out {output.merged}"

rule deduplicateSingle:
    input: join(removalDIR, "{id}_unaligned.fastq")
    output: "{collapsedDIR}/{id}_collapsed.fasta"
    conda: env
    shell: "fastx_collapser -i {input} -o {output}"

rule deduplicatePaired:
    input: join(assembledDIR, "{id}_assembled.fastq")
    output: "{collapsedDIR}/{id}_collapsed.fasta"
    conda: env
    shell: "fastx_collapser -i {input} -o {output}"

rule downloadNR:
    output: "{NRDIR}/nr"
    conda: env
    threads: workflow.cores
    shell: "wget -P {NRDIR} \
        ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz \
        && pigz -d {output} -p {threads}"

rule diamondMakeDB:
    input: join(NRDIR, "nr")
    output: "{diamondIndexDIR}/nr.dmnd"
    conda: env
    threads: workflow.cores
    shell: "diamond makedb --in {input} -d {output} --threads {threads}"

rule alignment:
    input:
        index = join(diamondIndexDIR, "nr.dmnd"),
        reads = join(collapsedDIR, "{id}_collapsed.fasta")
    output:
        matches = "{alignmentDIR}/{id}.m8",
        unaligned = "{alignmentDIR}/{id}_unaligned.fasta"
    conda: env
    threads: workflow.cores
    shell: "touch {output.unaligned} \
        && diamond blastx -d {input.index} -q {input.reads} -o {output.matches} --top 3 --un {output.unaligned} --threads {threads}"

rule downloadTaxonomy:
    output:
        "{taxonomicDIR}/db/complete_taxa.db",
        "{taxonomicDIR}/db/prot_mapping.db"
    conda: env
    shell: "basta taxonomy -d {taxonomicDIR}/db \
        && basta download prot -d {taxonomicDIR}/db"

rule taxonomicClassification:
    input:
        matches = join(alignmentDIR, "{id}.m8"),
        levedb = "{taxonomicDIR}/db/prot_mapping.db"
    output:
        out = "{taxonomicDIR}/{id}_tax.txt",
        lca = "{taxonomicDIR}/{id}_lca.html",
        best = "{taxonomicDIR}/{id}_best.html"
    conda: env
    threads: workflow.cores
    shell: "basta sequence {input.matches} {output.out} prot -d {taxonomicDIR}/db -l 1 -m 1 -b True \
        && awk -F \"\t\" '{{print $1\"\t\"$2}}' {output.out} > {taxonomicDIR}/{wildcards.id}_lca.txt \
        && awk -F \"\t\" '{{print $1\"\t\"$3}}' {output.out} > {taxonomicDIR}/{wildcards.id}_best.txt \
        && basta2krona.py {taxonomicDIR}/{wildcards.id}_lca.txt {output.lca} \
        && basta2krona.py {taxonomicDIR}/{wildcards.id}_best.txt {output.best}"

rule downloadUniprotMapping:
    output: "{functionalDIR}/idmapping_selected.tab"
    conda: env
    threads: workflow.cores
    shell: "wget -P {functionalDIR} \
        ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz \
        && pigz -d {output} -p {threads}"

rule createDictionaryNR2GO:
    input: join(functionalDIR, "idmapping_selected.tab") 
    output: "{functionalDIR}/db/NR2GO_mapping.db"
    conda: env
    threads: workflow.cores
    shell: "awk -F \"\t\" '{{if(($7!=\"\") && ($18!=\"\")){{print $18\"\t\"$7}}}}' {input} > {functionalDIR}/genbank2GO.txt \
        && awk -F \"\t\" '{{if(($4!=\"\") && ($7!=\"\")){{print $4\"\t\"$7}}}}' {input} > {functionalDIR}/refseq2GO.txt \
        && Rscript createDictionary.R {functionalDIR}/dictionary.txt {functionalDIR}/genbank2GO.txt {functionalDIR}/refseq2GO.txt {threads} \
        && basta create_db {functionalDIR}/dictionary.txt NR2GO_mapping.db 0 1 -d {functionalDIR}/db"

rule annotate:
    input:
        matches = join(alignmentDIR, "{id}.m8"),
        leveldb = "{functionalDIR}/db/NR2GO_mapping.db"
    output: "{functionalDIR}/{id}_functional.txt"
    conda: env
    threads: workflow.cores
    shell: "annotate {input.matches} {output} NR2GO -l -1 -d {functionalDIR}/db"
