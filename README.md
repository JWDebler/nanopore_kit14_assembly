# Nanopore assembly pipeline for kit 14 duplex called reads

- [Nanopore assembly pipeline for kit 14 duplex called reads](#nanopore-assembly-pipeline-for-kit-14-duplex-called-reads)
  - [Input](#input)
    - [Basecall with Dorado (non barcoded samples):](#basecall-with-dorado-non-barcoded-samples)
      - [Seperate simplex and duplex reads](#seperate-simplex-and-duplex-reads)
    - [Basecall with Dorado (barcoded samples as of Dorado 0.5.3)](#basecall-with-dorado-barcoded-samples-as-of-dorado-053)
  - [Running the pipeline](#running-the-pipeline)
  - [Profiles](#profiles)
  - [Parameters](#parameters)
  - [Known problems](#known-problems)


This is a Nanopore assembly pipeline aimed at Kit 14 generated reads which have been duplex called using the [dorado](https://github.com/nanoporetech/dorado/) basecaller. It requires that you have separated the reads into `sampleID.duplex.fastq.gz` and `sampleID.simplex.fastq.gz`.

It then does some basice QC by removing control DNA which is sometimes used during a run to debug potential problems, but which should not end up in the final assembly.

Assembly happens using two different assemblers, Flye and nextDenovo. Both are very fast and have different strenghts. I have found that the nextDenovo assembly overall is better with fewer contigs, but tends to trim telomeres and sometimes loses the mitogenome. Flye however is great at maintaining telomeres and the mitogenome tends to fall out as a single, non-concatenated contig (other assemblers create tandem copies as they don't expect a circular sequence).
Canu is used to generate corredted reads which I use to manually check and curate the assemblies in Geneious.

Currently this pipeline is optimised to run on a Nimbus instance with 16 cores and 64 GB of RAM.

## Input

The pipeline requires you to basecall your raw fast5 or pod5 files with dorado and then split the reads into simplex and duplex.

### Basecall with Dorado (non barcoded samples):

```
dorado duplex sup pod5s/ -r > sampleID.dorado.bam
```
#### Seperate simplex and duplex reads

Dorado generates a BAM file containing both the simplex and duplex reads. The respective reads can be extracted from that BAM file using samtools.

```
samtools view -O fastq -d dx:0 sampleID.dorado.bam |  gzip -9 >  sampleID.simplex.fastq.gz && \
samtools view -O fastq -d dx:1 sampleID.dorado.bam |  gzip -9 >  sampleID.duplex.fastq.gz
```

### Basecall with Dorado (barcoded samples as of Dorado 0.5.3)
As of version 0.5.3 Dorado can not yet do duplex calling and barcode demultiplexing in one go. The current workaround is to first call the reads in simplex and classification of barcodes, and then to extract a list of reads per barcode. Then you call the data again only using a subset of the already classified reads.
The commands below also demultiplex the pod5 files into individual files per barcode. For this [pod5 tools](https://github.com/nanoporetech/pod5-file-format) need to be installed.

```
location=/path/to/raw/data/containing/folder

dorado basecaller sup -r $location --kit-name SQK-NBD114-24 > all.bam && \
dorado demux --output-dir simplex --no-classify all.bam && \
rm simplex/*unclassified* && \
for file in simplex/*.bam; do id=$(echo "$file" | grep -oP 'barcode\d+'); samtools view -h $file | cut -f 1 > simplex/$id.readids.txt; done && \
for file in simplex/*readids.txt; do id=$(echo "$file" | grep -oP 'barcode\d+'); grep -v "@" $file > simplex/$id.clean.txt; done && \
for file in simplex/*clean.txt; do id=$(echo "$file" | grep -oP 'barcode\d+'); pod5 filter -r $location -i $file -t 10 --missing-ok --duplicate-ok --output $id.pod5; done && \
for file in *.pod5; do id=$(echo "$file" | grep -oP 'barcode\d+'); dorado duplex sup $file > $id.duplex.untrimmed.bam; done && \
for file in *untrimmed.bam; do id=$(echo "$file" | grep -oP 'barcode\d+'); dorado trim $file > $id.duplex.bam && \
rm *untrimmed.bam && \
for file in *duplex.bam; do id=$(echo "$file" | grep -oP 'barcode\d+'); samtools view  -@ 8 -O fastq -d dx:1 $file | pigz -9 > $id.duplex.fastq.gz; done && \
for file in *duplex.bam; do id=$(echo "$file" | grep -oP 'barcode\d+'); samtools view  -@ 8 -O fastq -d dx:0 $file | pigz -9 > $id.simplex.fastq.gz; done 
```

## Running the pipeline

After you have separated your reads into `sampleID.simplex.fastq.gz` and `sampleID.duplex.fastq.gz` you can run the pipeline by pointing `--reads` to the folder containing the fastq.gz files.

On a Nimbus instance with 64 Gb of ram and 16 Cores:
```
nextflow run jwdebler/nanopore_kit14_assembly -resume -latest -profile docker,nimbus --reads "reads/"
```
On the P2 Solo lab computer with 64 Gb of ram and 30 Cores:
```
nextflow run jwdebler/nanopore_kit14_assembly -resume -latest -profile docker,solo --reads "reads/"
```

## Profiles

We have a few profiles available to customise how the pipeline will run.

- `nimbus` sets the canu assembler to use 15 CPUs and 60GB RAM.
- `solo` sets paramaters to run on the lab computer running the P2 Solo sequencer.
- `docker` and `docker_sudo` sets it to use docker containers, `docker_sudo` is identical except that docker is run as root (required for some installations of docker).



## Parameters

```
--reads <glob>
    Required
    A folder containing two files per sample. The basename of the file is used as the sample ID and must contain `duplex` and `simplex`. Example of file name: `Sample1.duplex.fastq.gz`, `Sample1.simplex.fastq.gz`.

--genomeSize
    Default: 42m
    Used by the assemblers to calculate read coverage

--medakaModel <glob>
    Default: r1041_e82_400bps_sup_v4.3.0 (kit 14, sup, 5 kHz)
    The model that was used during basecalling.
    r1041_e82_400bps_sup_v4.1.0 (kit 14, sup, 4 kHz)
    r941_min_sup_g507 (LSK109, sup, 4kHz)

--minlen
    Min read length to keep for assembly
    (Default: 1000)

--quality
    Min read q-score to keep for read filtering
    (Default: 10)

--outdir <path>
    Default: `assembly`
    The directory to store the results in.
```

## Known problems
If the pipeline gets stopped or crashes during the medaka steps there is a chance it will get stuck on that forever if restarted due to an index file that is too large. Should that happen to you check the hash for the assembly processes "Assembly_flye" and "Assembly_nextdenovo". The hashes look like that at the start of the respective line:

[bc/a53374] process > Assembly_flye  
[a4/0cfae1] process > Assembly_nextdenovo

Go to the folders `work/bc/a53374...` and delete the `*.fai` and `*.mmi` files. After that resume the pipeline and it will run medaka properly.