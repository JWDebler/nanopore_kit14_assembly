
// Files need to be duplex called and split into "sampleID.duplex.fastq.gz" and "sampleID.simplex.fastq.gz"


def helpMessage() {
    log.info"""
    # Nanopore Q20+ genome assembly and polishing
    A pipeline for assembly and polishing of fungal genomes from Oxford Nanopore Q20+ reads from Kit14 using Flye and Medaka.

    ## Examples
    nextflow run jwdebler/nanopore_kit14_assembly -resume -latest -profile docker,nimbus --reads "reads/"

    ## Parameters
    --reads <glob>
        Required
        A folder containing two files per sample. 
        The basename of the file is used as the sample ID and must contain 
        `duplex` and `simplex`. 
        Example of file names: `Sample1.duplex.fastq.gz`, `Sample1.simplex.fastq.gz`.
        (Default: a folder called `reads/`)

    --genomeSize <glob>
        not required
        Size of genome, for example "42m" (Default: 42m)

    --medakaModel <glob>
        not required
        Which basecaller model was used?
        r1041_e82_400bps_sup_v5.0.0 (kit114, sup, 5 kHz)
        r1041_e82_400bps_sup_v4.1.0 (kit114, sup, 4 kHz)
        (Default: r1041_e82_400bps_sup_v4.3.0)

    --outdir <path>
        The directory to store the results in.
        (Default: `assembly`)

    --minlen
        Min read length to keep for assembly
        (Default: 1000)

    --quality
        Min read q-score to keep for read filtering
        (Default: 8)

    ## Exit codes
    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """
}

//by default script looks for reads in the current directory

params.reads=""
params.size="42m"
params.medakaModel="r1041_e82_400bps_sup_v5.0.0"
params.minlen="1000"
params.quality="8"

if (params.help) {
    helpMessage()
    exit 0
}

if ( params.reads ) {
    nanoporeReads = Channel
    .fromFilePairs(params.reads + "*{sim,du}plex.fastq.gz")
    .map { sampleID, reads -> [sampleID.tokenize('.')[0], reads] }
    .map { sampleID -> [sampleID[0]] + sampleID[1] }
    .tap { ReadsDuplexForChopper }
    .tap { ReadsSimplexForChopper }
    .tap { ReadsDuplexForQC }
    .tap { ReadsSimplexForQC }
    .view()

} else {
    nanoporeReads = Channel
    .fromFilePairs("*{sim,du}plex.fastq.gz")
    .map { sampleID, reads -> [sampleID.tokenize('.')[0], reads] }
    .map { sampleID -> [sampleID[0]] + sampleID[1] }
    .tap { ReadsDuplexForChopper }
    .tap { ReadsSimplexForChopper }
    .tap { ReadsDuplexForQC }
    .tap { ReadsSimplexForQC }
    .view()
}


process version_medaka {

    label "medaka"

    output:
    path 'versions.txt' into medaka_version

    """
    echo medaka: >> versions.txt
    medaka --version >> versions.txt
    echo --------------- >> versions.txt
    """
}


process version_seqkit {

    label "seqkit"

    output:
    path 'versions.txt' into seqkit_version

    """
    echo seqkit: >> versions.txt
    seqkit version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process version_flye{

    label "flye"

    output:
    path 'versions.txt' into flye_version

    """
    echo flye: >> versions.txt
    flye --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process version_nextdenovo {

    label "nextdenovo"

    output:
    path 'versions.txt' into nextdenovo_version

    """
    echo nextdenovo: >> versions.txt
    nextDenovo --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process version_chopper {

    label "chopper"

    output:
    path 'versions.txt' into chopper_version

    """
    echo chopper: >> versions.txt
    chopper --version >> versions.txt
    echo --------------- >> versions.txt
    """
}


process versions {

    input:
    path "medaka.txt" from medaka_version
    path "seqkit.txt" from seqkit_version
    path "flye.txt" from flye_version
    path "nextdenovo.txt" from nextdenovo_version
    path "chopper.txt" from chopper_version

    publishDir "${params.outdir}/", mode: 'copy', pattern: 'versions.txt'

    output:
    path "versions.txt"

    script:
    """
    cat medaka.txt seqkit.txt flye.txt nextdenovo.txt chopper.txt> versions.txt
    """
}

// filtering reads
// checks for duplicate read IDs and removes DCS reads

process QC_chopper_Duplex {

    label "chopper"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.fastq.gz'

    input:
    tuple sampleID, 'duplex.fastq.gz', 'simplex.fastq.gz' from ReadsDuplexForChopper

    output:
    tuple sampleID, "${sampleID}.duplex.chopper.200bp.q${params.quality}.fastq.gz" into FilteredDuplex200
    tuple sampleID, "${sampleID}.duplex.chopper.${params.minlen}bp.q${params.quality}.fastq.gz" into FilteredDuplex1000
    

    """
    wget https://raw.githubusercontent.com/JWDebler/nanopore_kit14_assembly/main/data/DCS.fasta
    minimap2 -d dcs.mmi DCS.fasta
    seqkit rmdup -n duplex.fastq.gz | minimap2 -t "${task.cpus}" -ax lr:hq dcs.mmi - | samtools view -O fastq -@ "${task.cpus}" - | chopper -t ${task.cpus}  -q ${params.quality} -l 200 | pigz -9 > ${sampleID}.duplex.chopper.200bp.q${params.quality}.fastq.gz 
    chopper -i ${sampleID}.duplex.chopper.200bp.q${params.quality}.fastq.gz -t ${task.cpus} -l ${params.minlen} | pigz -9 > ${sampleID}.duplex.chopper.${params.minlen}bp.q${params.quality}.fastq.gz
    """
}

process QC_chopper_Simplex {

    label "chopper"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.fastq.gz'

    input:
    tuple sampleID, 'duplex.fastq.gz', 'simplex.fastq.gz' from ReadsSimplexForChopper

    output:
    tuple sampleID, "${sampleID}.simplex.chopper.200bp.q${params.quality}.fastq.gz" into FilteredSimplex200
    tuple sampleID, "${sampleID}.simplex.chopper.${params.minlen}bp.q${params.quality}.fastq.gz" into FilteredSimplex1000
    
    """
    wget https://raw.githubusercontent.com/JWDebler/nanopore_kit14_assembly/main/data/DCS.fasta
    minimap2 -d dcs.mmi DCS.fasta
    seqkit rmdup -n simplex.fastq.gz | minimap2 -t "${task.cpus}" -ax lr:hq dcs.mmi - | samtools view -O fastq -@ "${task.cpus}" - | chopper -t ${task.cpus}  -q ${params.quality} -l 200 | pigz -9 > ${sampleID}.simplex.chopper.200bp.q${params.quality}.fastq.gz 
    chopper -i ${sampleID}.simplex.chopper.200bp.q${params.quality}.fastq.gz -t ${task.cpus} -l ${params.minlen} | pigz -9 > ${sampleID}.simplex.chopper.${params.minlen}bp.q${params.quality}.fastq.gz
    """
}

// read QC
process QC_nanoplot_Raw_Duplex {

    label "nanoplot"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/01-QC", pattern: '*.html'

    input:
    tuple sampleID, "${sampleID}.duplex.fastq.gz", "${sampleID}.simplex.fastq.gz" from ReadsDuplexForQC

    output:
    path "*.html"
    
    """
    NanoPlot \
    --fastq ${sampleID}.duplex.fastq.gz \
    -o output && \
    cp output/NanoPlot-report.html ${sampleID}.nanoplot.duplex.html
    """
}
process QC_nanoplot_Raw_Simplex {

    label "nanoplot"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/01-QC", pattern: '*.html'

    input:
    tuple sampleID, "${sampleID}.duplex.fastq.gz", "${sampleID}.simplex.fastq.gz" from ReadsSimplexForQC

    output:
    path "*.html"
    
    
    """
    NanoPlot \
    --fastq ${sampleID}.simplex.fastq.gz \
    -o output && \
    cp output/NanoPlot-report.html ${sampleID}.nanoplot.simplex.html
    """
}

process QC_nanoplot_Chopper_Duplex {

    label "nanoplot"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/01-QC", pattern: '*.html'

    input:
    tuple sampleID, "${sampleID}.duplex.fastq.gz" from FilteredDuplex1000

    output:
    path "*.html"
    tuple sampleID, "${sampleID}.duplex.fastq.gz" into FilterdForAssemblyDuplex
    
    
    """
    NanoPlot \
    --fastq ${sampleID}.duplex.fastq.gz \
    -o output && \
    cp output/NanoPlot-report.html ${sampleID}.nanoplot.chopper.duplex.html
    """
}

process QC_nanoplot_Chopper_Simplex {

    label "nanoplot"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/01-QC", pattern: '*.html'

    input:
    tuple sampleID, "${sampleID}.simplex.fastq.gz" from FilteredSimplex1000

    output:
    path "*.html"
    tuple sampleID, "${sampleID}.simplex.fastq.gz" into FilterdForAssemblySimplex
    
    """
    NanoPlot \
    --fastq ${sampleID}.simplex.fastq.gz \
    -o output && \
    cp output/NanoPlot-report.html ${sampleID}.nanoplot.chopper.simplex.html
    """
}

FilterdForAssemblyDuplex.join(FilterdForAssemblySimplex)
.tap { FilteredForFlye }
//.tap { FilteredForNextdenovo }

process mergeFilteredReads {

    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*mergedSimplexDuplex.fastq.gz'

    input:
    tuple sampleID, "duplex.fastq.gz", "simplex.fastq.gz" from FilteredDuplex200.join(FilteredSimplex200)

    output:
    tuple sampleID, "${sampleID}.mergedSimplexDuplex.fastq.gz" into MergedFilteredForMedakaFlye
    tuple sampleID, "${sampleID}.mergedSimplexDuplex.fastq.gz" into MergedFilteredForMedakaNextdenovo

    """
    cat duplex.fastq.gz simplex.fastq.gz > ${sampleID}.mergedSimplexDuplex.fastq.gz
    """
}

// flye assembly
process Assembly_flye {

    label "flye"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/03-assembly", pattern: '*_flye.*'

    input:
    tuple sampleID, "${sampleID}.duplex.chopper.fastq.gz", "${sampleID}.simplex.chopper.fastq.gz" from FilteredForFlye

    output:
    tuple sampleID, "${sampleID}_flye.fasta" into MedakaFlye
    tuple sampleID, "${sampleID}.duplex.chopper.fastq.gz", "${sampleID}.simplex.chopper.fastq.gz" into FilteredForNextdenovo
    file "${sampleID}_flye.assembly_info.txt"
    """
     flye \
    --nano-hq ${sampleID}.duplex.chopper.fastq.gz ${sampleID}.simplex.chopper.fastq.gz \
    --read-error 0.03 \
    --genome-size ${params.size} \
    --asm-coverage 50 \
    --threads ${task.cpus} \
    --out-dir ${sampleID}.flye 

    cp ${sampleID}.flye/assembly.fasta ${sampleID}_flye.fasta
    cp ${sampleID}.flye/assembly_info.txt ${sampleID}_flye.assembly_info.txt
    """
}

process Assembly_nextdenovo {

    label "nextdenovo"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/03-assembly", pattern: '*_nextdenovo.fasta'
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*corredted.fasta'

    input:
    tuple sampleID, "duplex.fastq.gz", "simplex.fastq.gz" from FilteredForNextdenovo

    output:
    tuple sampleID, "${sampleID}_nextdenovo.fasta", "simplex.fastq.gz" into MedakaNextdenovo
    tuple sampleID, "${sampleID}.nextdenovo.corredted.fasta"
    

    """
    ls duplex.fastq.gz simplex.fastq.gz > ${sampleID}.fofn

    echo '''
    [General]
    job_type = local
    job_prefix = ${sampleID}.nextdenovo
    task = all
    rewrite = yes
    deltmp = yes
    parallel_jobs = 5
    input_type = raw
    read_type = ont # clr, ont, hifi
    input_fofn = ${sampleID}.fofn
    workdir = ${sampleID}.nextdenovo

    [correct_option]
    read_cutoff = 1k
    genome_size = ${params.size} # estimated genome size
    sort_options = -m 4g -t 12
    minimap2_options_raw = -t 12
    pa_correction = 5
    correction_options = -p 10

    [assemble_option]
    minimap2_options_cns = -t 12
    nextgraph_options = -a 1
    ''' > ${sampleID}.config

    nextDenovo ${sampleID}.config

    cp ${sampleID}.nextdenovo/03.ctg_graph/nd.asm.fasta ${sampleID}_nextdenovo.fasta

    cat ${sampleID}.nextdenovo/02.cns_align/01.seed_cns.sh.work/seed_cns*/cns.fasta > ${sampleID}.nextdenovo.corredted.fasta
    """
}


process Polishing_medaka_flye {

    label "medaka"
    tag {sampleID}

    input:
    tuple sampleID, "flye.fasta", "mergedReads.fastq.gz" from MedakaFlye.join(MergedFilteredForMedakaFlye)

    output:
    tuple sampleID, "${sampleID}_flye_medaka.fasta" into SeqkitFlye

    """
    
    medaka_consensus \
    -i mergedReads.fastq.gz \
    -d flye.fasta \
    -o ${sampleID}_medaka_output \
    -m ${params.medakaModel} \
    -t ${task.cpus}

    cp ${sampleID}_medaka_output/consensus.fasta ${sampleID}_flye_medaka.fasta
    """
}

process Polishing_medaka_nextdenovo {

    label "medaka"
    tag {sampleID}

    input:
    tuple sampleID, "nextdenovo.fasta", "simplex.fastq.gz", "mergedReads.fastq.gz"  from MedakaNextdenovo.join(MergedFilteredForMedakaNextdenovo)

    output:
    tuple sampleID, "${sampleID}_nextenovo_medaka.fasta" into SeqkitNextdenovo
    tuple sampleID, "simplex.fastq.gz" into Correction_dechat

    """
    medaka_consensus \
    -i mergedReads.fastq.gz \
    -d nextdenovo.fasta \
    -o ${sampleID}_medaka_output \
    -m ${params.medakaModel} \
    -t ${task.cpus}

    cp ${sampleID}_medaka_output/consensus.fasta ${sampleID}_nextenovo_medaka.fasta
    """
}

process Cleanup_seqkitFlye {

    label "seqkit"
    tag {sampleID}

    publishDir "${params.outdir}/${sampleID}/04-medaka-polished", pattern: '*.fasta'

    input:
    tuple sampleID, "${sampleID}_flye_medaka.unsorted.fasta" from SeqkitFlye

    output:
    tuple sampleID, "${sampleID}_flye_medaka.fasta" into FlyeForRagtag

    """
    seqkit sort -lr ${sampleID}_flye_medaka.unsorted.fasta > ${sampleID}_flye_medaka.sorted.fasta
    seqkit replace -p '.+' -r '${sampleID}_ctg_{nr}' --nr-width 2 ${sampleID}_flye_medaka.sorted.fasta > ${sampleID}_flye_medaka.fasta
    """
}

process Cleanup_seqkitNextdenovo {

    label "seqkit"
    tag {sampleID}

    publishDir "${params.outdir}/${sampleID}/04-medaka-polished", pattern: '*.fasta'

    input:
    tuple sampleID, "${sampleID}_nextdenovo_medaka.unsorted.fasta" from SeqkitNextdenovo

    output:
    tuple sampleID, "${sampleID}_nextdenovo_medaka.fasta" into NextDenovoForRagtag

    """
    seqkit sort -lr ${sampleID}_nextdenovo_medaka.unsorted.fasta > ${sampleID}_nextdenovo_medaka.sorted.fasta
    seqkit replace -p '.+' -r '${sampleID}_ctg_{nr}' --nr-width 2 ${sampleID}_nextdenovo_medaka.sorted.fasta > ${sampleID}_nextdenovo_medaka.fasta
    """
}

// compare assemblies with RagTag and order according to Nextdenovo

process Cleanup_ragtag {

    label "ragtag"
    tag {sampleID}

    publishDir "${params.outdir}/${sampleID}/05-ragtag"

    input:
    tuple sampleID, "nextdenovo.fasta", "flye.fasta" from NextDenovoForRagtag.join(FlyeForRagtag)

    output:
    path "ragtag.*"

    """
    ragtag.py scaffold nextdenovo.fasta flye.fasta
    cp ragtag_output/* .
    """
}

process Correction_dechat {

    label "dechat"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.fasta.gz'

    input:
    tuple sampleID,  "${sampleID}.simplex.fastq.gz" from Correction_dechat

    output:
    tuple sampleID,  "${sampleID}.corrected.dechat.fasta.gz" 

    script:

    """
    dechat \
    -t ${task.cpus} \
    -o ${sampleID} \
    -i ${sampleID}.simplex.fastq.gz

    pigz -9 ${sampleID}.ec.fa

    mv ${sampleID}.ec.fa.gz ${sampleID}.corrected.dechat.fasta.gz
    """
}