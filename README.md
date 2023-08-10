# Nanopore assembly pipeline for kit 14 duplex called reads

This is a Nanopore assembly pipeline aimed at Kit 14 generated reads which have been duplex called using the [dorado](https://github.com/nanoporetech/dorado/) basecaller. It requires that you have separated the reads into `sampleID.duplex.fastq.gz` and `sampleID.simplex.fastq.gz`.

It then does some basice QC by removing control DNA which is sometimes used during a run to debug potential problems, but which should not end up in the final assembly.

Assembly happens using two different assemblers, Flye and nextDenovo. Both are very fast and have different strenghts. I have found that the nextDenovo assembly overall is better with fewer contigs, but tends to trim telomeres and sometimes loses the mitogenome. Flye however is great at maintaining telomeres and the mitogenome tends to fall out as a single, non-concatenated contig (other assemblers create tandem copies as they don't expect a circular sequence).
Canu is used to generate corredted reads which I use to manually check and curate the assemblies in Geneious.

Currently this pipeline is optimised to run on a Nimbus instance with 16 cores and 64 GB of RAM.

## Input

The pipeline requires you to basecall your raw fast5 or pod5 files with dorado and then split the reads into simplex and duplex.

### Basecall with Dorado:

- First download the latest model, this one is for kit 14 using the new 5 kHz sampling rate.
```
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.2.0
```
- Then run basecalling and simplex/duplex separation:

```
dorado duplex dna_r10.4.1_e8.2_400bps_sup@v4.2.0 pod5s/ > sampleID.dorado.bam
```
### Seperate simplex and duplex reads

Dorado generates a BAM file containing both the simplex and duplex reads. The respective reads can be extracted from that BAM file using samtools.

```
samtools view -O fastq -d dx:0 sampleID.dorado.bam |  gzip -9 >  sampleID.simplex.fastq.gz &&
samtools view -O fastq -d dx:1 sampleID.dorado.bam |  gzip -9 >  sampleID.duplex.fastq.gz
```


## Running the pipeline

After you have separated your reads into `sampleID.simplex.fastq.gz` and `sampleID.duplex.fastq.gz` you can run the pipeline by pointing `--reads` to the folder containing the fastq.gz files.

```
nextflow run jwdebler/nanopore_kit14_assembly -resume -latest -profile docker,nimbus --reads "reads/"
```


## Profiles

We have a few profiles available to customise how the pipeline will run.

- `nimbus` sets the canu assembler to use 15 CPUs and 60GB RAM.
- `zeus` sets the canu assembler to use 14 CPUs and 64GB RAM, and sets some cluster specific options to use the slurm based scheduler at Pawsey.
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
    Default: r1041_e82_400bps_sup_v4.2.0 (kit 14, sup, 5 kHz)
    The model that was used during basecalling.
    r1041_e82_400bps_sup_v4.1.0 (kit 14, sup, 4 kHz)
    r941_min_sup_g507 (LSK109, sup, 4kHz)

--outdir <path>
    Default: `assembly`
    The directory to store the results in.
```
