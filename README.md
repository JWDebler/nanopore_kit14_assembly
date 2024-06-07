# Nanopore assembly pipeline for kit 14 duplex called reads

- [Nanopore assembly pipeline for kit 14 duplex called reads](#nanopore-assembly-pipeline-for-kit-14-duplex-called-reads)
  - [Input](#input)
    - [Basecall with Dorado - non barcoded samples:](#basecall-with-dorado---non-barcoded-samples)
      - [Trimming with dorado](#trimming-with-dorado)
      - [Seperate simplex and duplex reads](#seperate-simplex-and-duplex-reads)
    - [Basecall with Dorado - barcoded samples (Dorado 0.7.1)](#basecall-with-dorado---barcoded-samples-dorado-071)
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

### Basecall with Dorado - non barcoded samples:

```
pod5 view folder_with_pod5_files/*.pod5 -t $(nproc) --include "read_id, channel" --output channel.summary.tsv
pod5 subset folder_with_pod5_files/ -t $(nproc) --summary channel.summary.tsv --columns channel --output pod5_by_channel/ 
dorado duplex sup pod5_by_channel/  > sampleID.dorado.untrimmed.bam
```

#### Trimming with dorado

```
dorado trim -t $(nproc) sampleID.dorado.untrimmed.bam >  sampleID.dorado.trimmed.bam
```
#### Seperate simplex and duplex reads

Dorado generates a BAM file containing both the simplex and duplex reads. The respective reads can be extracted from that BAM file using samtools.

```
samtools view -O fastq -d dx:0 sampleID.dorado.trimmed.bam |  gzip -9 >  sampleID.simplex.fastq.gz && \
samtools view -O fastq -d dx:1 sampleID.dorado.trimmed.bam |  gzip -9 >  sampleID.duplex.fastq.gz
```

### Basecall with Dorado - barcoded samples (Dorado 0.7.1)
As of version 0.7.1 Dorado can not yet do duplex calling and barcode demultiplexing in one go. The current workaround is to first call the reads in simplex mode and do classification of barcodes, and then to extract a list of reads per barcode. Then you call the data again only using a subset of the already classified reads.
The script below also demultiplexes the pod5 files into individual files per barcode and then into individual channels. For this [pod5 tools](https://github.com/nanoporetech/pod5-file-format) need to be installed.

Use the `basecalling.sh` script to demultiplex and duplex call a folder of pod5 files. Your output will be a folder of demultiplexed and separated by channel pod5 files. Separating by channel significantly speeds up duplex calling.
The final output are three files per barcode called `barcodeX.duplex.fastq.gz`, `barcodeX.simplex.fastq.gz` and `barcodeX.simplex.corrected.fasta.gz`.
The two `fastq.gz` files are going to be the input files for this assembly pipeline. You can use the `corrected.fasta.gz` file for mapping back to your draft assembly and manually fixing it.
At this point I am using the 'raw' files as imput for the assemblers, this might change in the future. 

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

- `nimbus` sets parameters for a Nimbus cloud instance with 16 cores and 64 Gb of RAM.
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
    Default: r1041_e82_400bps_sup_v5.0.0 (kit 14, sup, 5 kHz)
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

