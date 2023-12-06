
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
        r1041_e82_400bps_sup_v4.2.0 (kit114, sup, 5 kHz)
        r1041_e82_400bps_sup_v4.1.0 (kit114, sup, 4 kHz)
        (Default: r1041_e82_400bps_sup_v4.2.0)

    --canuSlow
        Disables canu fast mode.
        (Default: false)

    --outdir <path>
        The directory to store the results in.
        (Default: `assembly`)

    --minlen
        Min read length to keep for assembly
        (Default: 1000)

    --quality
        Min read q-score to keep for read filtering
        (Default: 10)

    ## Exit codes
    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """
}

//by default script looks for reads in the current directory

params.reads=""
params.size="42m"
params.medakaModel="r1041_e82_400bps_sup_v4.2.0"
params.minlen="1000"
params.quality="10"

if (params.help) {
    helpMessage()
    exit 0
}

if ( params.reads ) {
    nanoporeReads = Channel
    .fromFilePairs(params.reads + "*{sim,du}plex.fastq.gz")
    .map { sampleID, reads -> [sampleID.tokenize('.')[0], reads] }
    .map { sampleID -> [sampleID[0]] + sampleID[1] }
    .tap { ReadsForDCSQC }
    .view()

} else {
    nanoporeReads = Channel
    .fromFilePairs("*{sim,du}plex.fastq.gz")
    .map { sampleID, reads -> [sampleID.tokenize('.')[0], reads] }
    .map { sampleID -> [sampleID[0]] + sampleID[1] }
    .tap { ReadsForDCSQC }
    .view()
}

process version_canu {

    label "canu"

    output:
    path 'versions.txt' into canu_version

    """
    echo canu: >> versions.txt
    canu --version >> versions.txt
    echo --------------- >> versions.txt
    """
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

process version_minimap2 {

    label "minimap2"

    output:
    path 'versions.txt' into minimap2_version

    """
    echo minimap2: >> versions.txt
    minimap2 --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process version_samtools {

    label "samtools"

    output:
    path 'versions.txt' into samtools_version

    """
    echo samtools: >> versions.txt
    samtools version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process versions {

    input:
    path "canu.txt" from canu_version
    path "medaka.txt" from medaka_version
    path "seqkit.txt" from seqkit_version
    path "flye.txt" from flye_version
    path "nextdenovo.txt" from nextdenovo_version
    path "chopper.txt" from chopper_version
    path "minimap2.txt" from minimap2_version
    path "samtools.txt" from samtools_version

    publishDir "${params.outdir}/", mode: 'copy', pattern: 'versions.txt'

    output:
    path "versions.txt"

    script:
    """
    cat canu.txt medaka.txt seqkit.txt flye.txt nextdenovo.txt chopper.txt minimap2.txt samtools.txt > versions.txt
    """
}

process QC_DCS_minimap {

    label "minimap2"
    tag {sampleID}

    input:
    tuple sampleID, 'duplex.fastq.gz', 'simplex.fastq.gz' from ReadsForDCSQC

    output:
    tuple sampleID, 'duplex.sam', 'simplex.sam' into DCSalignments

    """
    wget https://raw.githubusercontent.com/JWDebler/nanopore_kit14_assembly/main/data/DCS.fasta
    minimap2 -d dcs.mmi DCS.fasta
    minimap2 -t "${task.cpus}" -ax map-ont dcs.mmi duplex.fastq.gz > duplex.sam
    minimap2 -t "${task.cpus}" -ax map-ont dcs.mmi simplex.fastq.gz > simplex.sam
    """

}

process QC_DCS_filtering_reads {

    label "samtools"
    tag {sampleID}

    input:
    tuple sampleID, 'duplex.sam', 'simplex.sam' from DCSalignments

    output:
    tuple sampleID, 'duplex.fastq', 'simplex.fastq' into DCSFilteredReads

    """
    samtools view -@ "${task.cpus}" -b -f 4 duplex.sam | samtools fastq -@ "${task.cpus}" - > duplex.fastq
    samtools view -@ "${task.cpus}" -b -f 4 simplex.sam | samtools fastq -@ "${task.cpus}" - > simplex.fastq
    """
}

DCSFilteredReads
.tap { ReadsDuplexForChopper }
.tap { ReadsSimplexForChopper }
.tap { ReadsDuplexForQC }
.tap { ReadsSimplexForQC }

// filtering reads

process QC_chopper_Simplex {

    label "chopper"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.fastq.gz'

    input:
    tuple sampleID, 'duplex.fastq', 'simplex.fastq' from ReadsSimplexForChopper

    output:
    //path "${sampleID}.simplex.chopper.fastq.gz"
    tuple sampleID, "${sampleID}.simplex.chopper.200bp.q${params.quality}.fastq.gz" into FilteredSimplex200
    tuple sampleID, "${sampleID}.simplex.chopper.${params.minlen}bp.q${params.quality}.fastq.gz" into FilteredSimplex1000
    //tuple sampleID, "${sampleID}.simplex.chopper.${params.minlen}bp.q${params.quality}.fastq.gz" into ReadsForCorrection

    """
    cat simplex.fastq | chopper -q ${params.quality} -l 200 > ${sampleID}.simplex.chopper.200bp.q${params.quality}.fastq
    cat ${sampleID}.simplex.chopper.200bp.q${params.quality}.fastq | chopper -l ${params.minlen} -q ${params.quality}| gzip -9 > ${sampleID}.simplex.chopper.${params.minlen}bp.q${params.quality}.fastq.gz
    gzip -9 ${sampleID}.simplex.chopper.200bp.q${params.quality}.fastq 

    """
}

process QC_chopper_Duplex {

    label "chopper"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.fastq.gz'

    input:
    tuple sampleID, 'duplex.fastq', 'simplex.fastq' from ReadsDuplexForChopper

    output:
    //path "${sampleID}.simplex.chopper.fastq.gz"
    tuple sampleID, "${sampleID}.duplex.chopper.200bp.q${params.quality}.fastq.gz" into FilteredDuplex200
    tuple sampleID, "${sampleID}.duplex.chopper.${params.minlen}bp.q${params.quality}.fastq.gz" into FilteredDuplex1000
    

    """
    cat duplex.fastq | chopper -q ${params.quality} -l 200 > ${sampleID}.duplex.chopper.200bp.q${params.quality}.fastq
    cat ${sampleID}.duplex.chopper.200bp.q${params.quality}.fastq | chopper -l ${params.minlen} | gzip -9 > ${sampleID}.duplex.chopper.${params.minlen}bp.q${params.quality}.fastq.gz
    gzip -9 ${sampleID}.duplex.chopper.200bp.q${params.quality}.fastq
    """
}

// read QC
process QC_nanoplot_Raw_Duplex {

    label "nanoplot"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/01-QC", pattern: '*.html'

    input:
    tuple sampleID, 'duplex.fastq', 'simplex.fastq' from ReadsDuplexForQC

    output:
    path "*.html"
    
    """
    NanoPlot \
    --fastq duplex.fastq \
    -o output && \
    cp output/NanoPlot-report.html ${sampleID}.nanoplot.duplex.html
    """
}

process QC_nanoplot_Raw_Simplex {

    label "nanoplot"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/01-QC", pattern: '*.html'

    input:
    tuple sampleID, 'duplex.fastq', 'simplex.fastq' from ReadsSimplexForQC

    output:
    path "*.html"
    
    
    """
    NanoPlot \
    --fastq simplex.fastq \
    -o output && \
    cp output/NanoPlot-report.html ${sampleID}.nanoplot.simplex.html
    """
}

process QC_nanoplot_Chopper_Duplex {

    label "nanoplot"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/01-QC", pattern: '*.html'

    input:
    tuple sampleID, 'reads.fastq.gz' from FilteredDuplex1000

    output:
    path "*.html"
    tuple sampleID, "reads.fastq.gz" into FilterdForAssemblyDuplex
    
    
    """
    NanoPlot \
    --fastq reads.fastq.gz \
    -o output && \
    cp output/NanoPlot-report.html ${sampleID}.nanoplot.chopper.duplex.html
    """
}

process QC_nanoplot_Chopper_Simplex {

    label "nanoplot"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/01-QC", pattern: '*.html'

    input:
    tuple sampleID, 'reads.fastq.gz' from FilteredSimplex1000

    output:
    path "*.html"
    tuple sampleID, "reads.fastq.gz" into FilterdForAssemblySimplex
    
    """
    NanoPlot \
    --fastq reads.fastq.gz \
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
    publishDir "${params.outdir}/${sampleID}/03-assembly", pattern: '*_flye.fasta'

    input:
    tuple sampleID, "${sampleID}.duplex.chopper.fastq.gz", "${sampleID}.simplex.chopper.fastq.gz" from FilteredForFlye

    output:
    tuple sampleID, "${sampleID}_flye.fasta" into MedakaFlye
    tuple sampleID, "${sampleID}.duplex.chopper.fastq.gz", "${sampleID}.simplex.chopper.fastq.gz" into FilteredForNextdenovo

    """
     flye \
    --nano-hq ${sampleID}.duplex.chopper.fastq.gz ${sampleID}.simplex.chopper.fastq.gz \
    --read-error 0.03 \
    --genome-size ${params.size} \
    --asm-coverage 50 \
    --threads "${task.cpus}" \
    --out-dir ${sampleID}.flye 

    cp ${sampleID}.flye/assembly.fasta ${sampleID}_flye.fasta
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
    tuple sampleID, "${sampleID}_nextdenovo.fasta" into MedakaNextdenovo
    tuple sampleID, "${sampleID}.nextdenovo.corredted.fasta"
    tuple sampleID, "duplex.fastq.gz", "simplex.fastq.gz" into ReadsForCorrection

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
    sort_options = -m 4g -t 3
    minimap2_options_raw = -t 3
    pa_correction = 5
    correction_options = -p 10

    [assemble_option]
    minimap2_options_cns = -t 8
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
    -d flye.fasta \
    -i mergedReads.fastq.gz \
    -o ${sampleID}_medaka_output \
    -t ${task.cpus} \
    -m ${params.medakaModel}

    cp ${sampleID}_medaka_output/consensus.fasta ${sampleID}_flye_medaka.fasta
    """
}

process Polishing_medaka_nextdenovo {

    label "medaka"
    tag {sampleID}

    input:
    tuple sampleID, "nextdenovo.fasta", "mergedReads.fastq.gz"  from MedakaNextdenovo.join(MergedFilteredForMedakaNextdenovo)

    output:
    tuple sampleID, "${sampleID}_nextenovo_medaka.fasta" into SeqkitNextdenovo

    """
    medaka_consensus \
    -d nextdenovo.fasta \
    -i mergedReads.fastq.gz \
    -o ${sampleID}_medaka_output \
    -t ${task.cpus} \
    -m ${params.medakaModel}

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
    seqkit replace -p '.+' -r 'flye_ctg_{nr}' --nr-width 2 ${sampleID}_flye_medaka.sorted.fasta > ${sampleID}_flye_medaka.fasta
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
    seqkit replace -p '.+' -r 'nd_ctg_{nr}' --nr-width 2 ${sampleID}_nextdenovo_medaka.sorted.fasta > ${sampleID}_nextdenovo_medaka.fasta
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

// read correction
process Correction_canu {

    label "canu"
    tag {sampleID}
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.fasta.gz'
    publishDir "${params.outdir}/${sampleID}/02-processed-reads", pattern: '*.report'

    input:
    tuple sampleID,  "duplex.fastq.gz", "${sampleID}.simplex.fastq.gz" from ReadsForCorrection

    output:
    path "${sampleID}.simplex.corrected.fasta.gz"
    path "${sampleID}.simplex.corrected.report"

    script:
    // See: https://groovy-lang.org/operators.html#_elvis_operator
    fast_option = params.canuSlow ? "" : "-fast "

    """
    canu \
    -correct \
    -p ${sampleID} \
    -d ${sampleID} \
    genomeSize=${params.size} \
    ${fast_option} \
    -nanopore ${sampleID}.simplex.fastq.gz

    cp ${sampleID}/*correctedReads.fasta.gz ${sampleID}.simplex.corrected.fasta.gz
    cp ${sampleID}/*.report ${sampleID}.simplex.corrected.report
    """
}
