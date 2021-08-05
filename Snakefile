from os.path import join
import os

# input args
preprocessingDIR = "Pipeline/data"
alignmentDIR = "Pipeline/alignment"
taxonomicDIR = "Pipeline/taxonomic"
functionalDIR = "Pipeline/functional"
resultDIR = "Pipeline/result"
inputDIR = join(preprocessingDIR, "raw")
phredQuality = "20"
trimmedDIR = join(preprocessingDIR, "trimmed")
removalDIR = join(preprocessingDIR, "removal")
referenceDIR = join(removalDIR, "reference")
ensemblRelease = "104"
bowtie2IndexDIR = join(removalDIR, "index")
# max memory in byte to be used in SdBG construction (if set between 0-1, fraction of the machine's total memory)
megahitMemory = "0.9"
# max number of threads for kaiju-mkbwt
kaijumkbwtThreads = 10
assembledDIR = join(preprocessingDIR, "assembled")
mergedDIR = join(preprocessingDIR, "merged")
collapsedDIR = join(preprocessingDIR, "collapsed")
NRDIR = join(alignmentDIR, "db")
diamondIndexDIR = join(alignmentDIR, "index")

IDs = glob_wildcards(join(inputDIR,
    "{id, [A-Za-z0-9]+}{suffix, (_[1-2])?}.fastq"))

ruleorder: qualityControlPaired > qualityControlSingle
ruleorder: removeHumanContaminantsPaired > removeHumanContaminantsSingle
ruleorder: deduplicatePaired > deduplicateSingle
ruleorder: diamondMakeDB > kaijuPrepare
ruleorder: assemblyPaired > assemblySingle
ruleorder: taxonomicClassificationPaired > taxonomicClassificationSingle

rule all:
    input:
        expand(join(resultDIR, "{id}_kaiju.names"), id = IDs.id),
        expand(join(resultDIR, "{id}_contigs_kaiju.names"), id = IDs.id),
        expand(join(resultDIR, "{id}_functional_GO.txt"), id = IDs.id),
        expand(join(resultDIR, "{id}_functional_contigs_GO.txt"), id = IDs.id)

rule qualityControlSingle:
    input: join(inputDIR, "{id}.fastq")
    output:
        trimmed = "{trimmedDIR}/{id}_trim.fastq",
        html = "{trimmedDIR}/{id}_report.html",
        json = "{trimmedDIR}/{id}_report.json"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && fastp -i {input} -o {output.trimmed} -q {phredQuality} -w {threads} -h {output.html} -j {output.json}"

rule qualityControlPaired:
    input:
        f = join(inputDIR, "{id}_1.fastq"),
        r = join(inputDIR, "{id}_2.fastq"),
    output:
        f = "{trimmedDIR}/{id}_1_trim.fastq",
        r = "{trimmedDIR}/{id}_2_trim.fastq",
        html = "{trimmedDIR}/{id}_report.html",
        json = "{trimmedDIR}/{id}_report.json"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && fastp -i {input.f} -I {input.r} -o {output.f} \
        -O {output.r} -q {phredQuality} -w {threads} \
        --detect_adapter_for_pe -h {output.html} -j {output.json}"

rule downloadHumanPrimaryAssembly:
    output: "{referenceDIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && wget -P {referenceDIR} \
        ftp://ftp.ensembl.org/pub/release-{ensemblRelease}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
        && pigz -d {output} -p {threads}"

rule bowtie2BuildHumanIndex:
    params: indexPrefix = join(bowtie2IndexDIR, "hostHS")
    input: join(referenceDIR, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    output:
        expand(join(bowtie2IndexDIR, "hostHS.{index}.bt2l"), index = range(1, 5)),
        expand(join(bowtie2IndexDIR, "hostHS.rev.{index}.bt2l"), index = range(1, 3))
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && bowtie2-build --large-index {input} {params.indexPrefix} --threads {threads} \
        && rm {input}"

rule removeHumanContaminantsSingle:
    params: indexPrefix = join(bowtie2IndexDIR, "hostHS"),
    input:
        expand(join(bowtie2IndexDIR, "hostHS.{index}.bt2l"), index = range(1, 5)),
        expand(join(bowtie2IndexDIR, "hostHS.rev.{index}.bt2l"), index = range(1, 3)),
        trimmed = join(trimmedDIR, "{id}_trim.fastq")
    output: "{removalDIR}/{id}_unaligned.fastq"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && bowtie2 -x {params.indexPrefix} -U {input.trimmed} -S {removalDIR}/{wildcards.id}.sam -p {threads} \
        && samtools view -bS {removalDIR}/{wildcards.id}.sam > {removalDIR}/{wildcards.id}.bam \
        && samtools view -b -f 4 -F 256 {removalDIR}/{wildcards.id}.bam > {removalDIR}/{wildcards.id}_unaligned.bam \
        && samtools sort -n {removalDIR}/{wildcards.id}_unaligned.bam -o {removalDIR}/{wildcards.id}_unaligned_sorted.bam \
        && samtools bam2fq {removalDIR}/{wildcards.id}_unaligned_sorted.bam > {removalDIR}/{wildcards.id}_unaligned.fastq \
        && rm {removalDIR}/{wildcards.id}.sam {removalDIR}/{wildcards.id}_unaligned.bam {removalDIR}/{wildcards.id}_unaligned_sorted.bam {input.trimmed}"

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
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && bowtie2 -x {params.indexPrefix} -1 {input.f} -2 {input.r} -S {removalDIR}/{wildcards.id}.sam -p {threads} \
        && samtools view -bS {removalDIR}/{wildcards.id}.sam > {removalDIR}/{wildcards.id}.bam \
        && samtools view -b -f 12 -F 256 {removalDIR}/{wildcards.id}.bam > {removalDIR}/{wildcards.id}_unaligned.bam \
        && samtools sort -n {removalDIR}/{wildcards.id}_unaligned.bam -o {removalDIR}/{wildcards.id}_unaligned_sorted.bam \
        && samtools bam2fq {removalDIR}/{wildcards.id}_unaligned_sorted.bam > {removalDIR}/{wildcards.id}_unaligned.fastq \
        && cat {removalDIR}/{wildcards.id}_unaligned.fastq | grep '^@.*/1$' -A 3 --no-group-separator > {removalDIR}/{wildcards.id}_unaligned_1.fastq \
        && cat {removalDIR}/{wildcards.id}_unaligned.fastq | grep '^@.*/2$' -A 3 --no-group-separator > {removalDIR}/{wildcards.id}_unaligned_2.fastq \
        && rm {removalDIR}/{wildcards.id}.sam {removalDIR}/{wildcards.id}_unaligned.bam {removalDIR}/{wildcards.id}_unaligned_sorted.bam {removalDIR}/{wildcards.id}_unaligned.fastq {input.f} {input.r}"

rule mergePaired:
    input:
        f = join(removalDIR, "{id}_unaligned_1.fastq"),
        r = join(removalDIR, "{id}_unaligned_2.fastq")
    output:
        merged = "{mergedDIR}/{id}_merged.fastq",
        html = "{mergedDIR}/{id}_report.html",
        json = "{mergedDIR}/{id}_report.json"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && fastp -i {input.f} -I {input.r} -o {mergedDIR}/{wildcards.id}_unmerged_1.fastq \
        -O {mergedDIR}/{wildcards.id}_unmerged_2.fastq -q {phredQuality} -w {threads} \
        --detect_adapter_for_pe -h {output.html} -j {output.json} \
        -m --merged_out {output.merged}"

rule assemblySingle:
    input:
        reads = join(removalDIR, "{id}_unaligned.fastq")
    output:
        contigs = "{assembledDIR}/{id}/final.contigs.fa"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && megahit -r {input.reads} --force -o {assembledDIR}/{wildcards.id} -t {threads} -m {megahitMemory}"

rule assemblyPaired:
    input:
        f = join(removalDIR, "{id}_unaligned_1.fastq"),
        r = join(removalDIR, "{id}_unaligned_2.fastq")
    output:
        contigs = "{assembledDIR}/{id}/final.contigs.fa"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && megahit -1 {input.f} -2 {input.r} --force -o {assembledDIR}/{wildcards.id} -t {threads} -m {megahitMemory}"

rule taxonomicClassificationSingle:
    input:
        reads = join(collapsedDIR, "{id}_collapsed.fasta"),
        fmi = join(taxonomicDIR, "db/kaijuNR.fmi"),
        names = join(taxonomicDIR, "db/names.dmp"),
        nodes = join(taxonomicDIR, "db/nodes.dmp")
    output:
        ranks = "{resultDIR}/{id}_kaiju.names",
        html = "{resultDIR}/{id}_krona.html"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && kaiju -t {input.nodes} -f {input.fmi} -i {input.reads} -o {resultDIR}/kaiju.out -z {threads} \
        && kaiju-addTaxonNames -t {input.nodes} -n {input.names} -r superkingdom,phylum,class,order,family,genus,species -i {resultDIR}/kaiju.out -o {output.ranks} \
        && kaiju2krona -t {input.nodes} -n {input.names} -i {resultDIR}/kaiju.out -o {resultDIR}/kaiju_krona \
        && ktImportText -o {output.html} {resultDIR}/kaiju_krona \
        && rm {resultDIR}/kaiju.out {resultDIR}/kaiju_krona"

rule taxonomicClassificationPaired:
    input:
        f = join(removalDIR, "{id}_unaligned_1.fastq"),
        r = join(removalDIR, "{id}_unaligned_2.fastq"),
        fmi = join(taxonomicDIR, "db/kaijuNR.fmi"),
        names = join(taxonomicDIR, "db/names.dmp"),
        nodes = join(taxonomicDIR, "db/nodes.dmp")
    output:
        ranks = "{resultDIR}/{id}_kaiju.names",
        html = "{resultDIR}/{id}_krona.html"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && kaiju -t {input.nodes} -f {input.fmi} -i {input.f} -j {input.r} -o {resultDIR}/kaiju.out -z {threads} \
        && kaiju-addTaxonNames -t {input.nodes} -n {input.names} -r superkingdom,phylum,class,order,family,genus,species -i {resultDIR}/kaiju.out -o {output.ranks} \
        && kaiju2krona -t {input.nodes} -n {input.names} -i {resultDIR}/kaiju.out -o {resultDIR}/kaiju_krona \
        && ktImportText -o {output.html} {resultDIR}/kaiju_krona \
        && rm {resultDIR}/kaiju.out {resultDIR}/kaiju_krona"

rule taxonomicClassificationContigs:
    input:
        reads = join(assembledDIR, "{id}/final.contigs.fa"),
        fmi = join(taxonomicDIR, "db/kaijuNR.fmi"),
        names = join(taxonomicDIR, "db/names.dmp"),
        nodes = join(taxonomicDIR, "db/nodes.dmp")
    output:
        ranks = "{resultDIR}/{id}_contigs_kaiju.names",
        html = "{resultDIR}/{id}_contigs_krona.html"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && kaiju -t {input.nodes} -f {input.fmi} -i {input.reads} -o {resultDIR}/kaiju.out -z {threads} \
        && kaiju-addTaxonNames -t {input.nodes} -n {input.names} -r superkingdom,phylum,class,order,family,genus,species -i {resultDIR}/kaiju.out -o {output.ranks} \
        && kaiju2krona -t {input.nodes} -n {input.names} -i {resultDIR}/kaiju.out -o {resultDIR}/kaiju_krona \
        && ktImportText -o {output.html} {resultDIR}/kaiju_krona \
        && rm {resultDIR}/kaiju.out {resultDIR}/kaiju_krona"

rule deduplicateSingle:
    input: join(removalDIR, "{id}_unaligned.fastq")
    output: "{collapsedDIR}/{id}_collapsed.fasta"
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && fastx_collapser -i {input} -o {output}"

rule deduplicatePaired:
    input: join(mergedDIR, "{id}_merged.fastq")
    output: "{collapsedDIR}/{id}_collapsed.fasta"
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && fastx_collapser -i {input} -o {output}"

rule downloadNR:
    output: "{NRDIR}/nr"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && wget -P {NRDIR} \
        ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz \
        && pigz -d {output} -p {threads}"

rule diamondMakeDB:
    input: join(NRDIR, "nr")
    output: "{diamondIndexDIR}/nr.dmnd"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && diamond makedb --in {input} -d {output} --threads {threads}"

rule kaijuPrepare:
    input: join(NRDIR, "nr")
    output:
        convertedNR = "{taxonomicDIR}/db/kaijuNR.fasta",
        names = "{taxonomicDIR}/db/names.dmp",
        nodes = "{taxonomicDIR}/db/nodes.dmp"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && wget -P {taxonomicDIR}/db ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz \
        && tar -xf {taxonomicDIR}/db/taxdump.tar.gz nodes.dmp names.dmp \
        && mv nodes.dmp {taxonomicDIR}/db && mv names.dmp {taxonomicDIR}/db \
        && rm {taxonomicDIR}/db/taxdump.tar.gz \
        && wget -P {taxonomicDIR}/db ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz \
        && pigz -d {taxonomicDIR}/db/prot.accession2taxid.gz -p {threads} \
        && kaiju-convertNR -t {output.nodes} -g {taxonomicDIR}/db/prot.accession2taxid -e $(conda info --base)/envs/medusaPipeline/bin/kaiju-excluded-accessions.txt -a -o {output.convertedNR} -i {input} \
        && rm {taxonomicDIR}/db/prot.accession2taxid {input}"

rule kaijuMakeBWT:
    input: join(taxonomicDIR, "db/kaijuNR.fasta")
    output:
        bwt = "{taxonomicDIR}/db/kaijuNR.bwt",
        sa = "{taxonomicDIR}/db/kaijuNR.sa"
    threads: kaijumkbwtThreads
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && kaiju-mkbwt -n {threads} -a ACDEFGHIKLMNPQRSTVWY -o {taxonomicDIR}/db/kaijuNR {input} \
        && rm {input}"

rule kaijuMakeFMI:
    input:
        bwt = join(taxonomicDIR, "db/kaijuNR.bwt"),
        sa = join(taxonomicDIR, "db/kaijuNR.sa")
    output:
        fmi = "{taxonomicDIR}/db/kaijuNR.fmi"
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && kaiju-mkfmi {taxonomicDIR}/db/kaijuNR \
        && rm {input.bwt} {input.sa}"

rule alignment:
    input:
        index = join(diamondIndexDIR, "nr.dmnd"),
        reads = join(collapsedDIR, "{id}_collapsed.fasta")
    output:
        matches = "{alignmentDIR}/{id}.m8",
        unaligned = "{alignmentDIR}/{id}_unaligned.fasta"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && touch {output.unaligned} \
        && diamond blastx -d {input.index} -q {input.reads} -o {output.matches} --top 3 --un {output.unaligned} --threads {threads}"

rule alignmentContigs:
    input:
        index = join(diamondIndexDIR, "nr.dmnd"),
        reads = join(assembledDIR, "{id}/final.contigs.fa")
    output:
        matches = "{alignmentDIR}/{id}_contigs.m8",
        unaligned = "{alignmentDIR}/{id}_contigs_unaligned.fasta"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && touch {output.unaligned} \
        && diamond blastx -d {input.index} -q {input.reads} -o {output.matches} --top 3 -F 15 --range-culling --un {output.unaligned} --threads {threads}"

rule downloadUniprotMapping:
    output: "{functionalDIR}/idmapping_selected.tab"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && wget -P {functionalDIR} \
        ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz \
        && pigz -d {output} -p {threads}"

rule createDictionaries:
    input: join(functionalDIR, "idmapping_selected.tab") 
    output:
        directory("{functionalDIR}/db/NR2GO.ldb")
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && awk -F \"\t\" '{{if(($7!=\"\") && ($18!=\"\")){{print $18\"\t\"$7}}}}' {input} > {functionalDIR}/genbank2GO.txt \
        && awk -F \"\t\" '{{if(($4!=\"\") && ($7!=\"\")){{print $4\"\t\"$7}}}}' {input} > {functionalDIR}/refseq2GO.txt \
        && Rscript createDictionary.R {functionalDIR}/NR2GO.txt {functionalDIR}/genbank2GO.txt {functionalDIR}/refseq2GO.txt {threads} \
        && annotate createdb {functionalDIR}/NR2GO.txt NR2GO 0 1 -d {functionalDIR}/db \
        && rm {functionalDIR}/genbank2GO.txt {functionalDIR}/refseq2GO.txt"

rule annotateGO:
    input:
        matches = join(alignmentDIR, "{id}.m8"),
        NR2GO = join(functionalDIR, "db/NR2GO.ldb")
    output:
        GO = "{resultDIR}/{id}_functional_GO.txt",
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && annotate idmapping {input.matches} {output.GO} NR2GO -l 1 -d {functionalDIR}/db"

rule annotateGOContigs:
    input:
        matchesContigs = join(alignmentDIR, "{id}_contigs.m8"),
        NR2GO = join(functionalDIR, "db/NR2GO.ldb")
    output:
        contigsGO = "{resultDIR}/{id}_functional_contigs_GO.txt"
    threads: workflow.cores
    shell: "set +eu \
        && . $(conda info --base)/etc/profile.d/conda.sh && conda activate medusaPipeline \
        && annotate idmapping {input.matchesContigs} {output.contigsGO} NR2GO -l 1 -d {functionalDIR}/db"