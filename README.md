# Nanopore assembly pipeline

This is a draft nanopore assembly pipeline used by some groups at the CCDM.
More documentation to come.

## Examples

```bash
nextflow run ccdmb/nanopore_assembly -resume -latest --nanoporeReads "*.fastq.gz"
```


## Profiles

We have a few profiles available to customise how the pipeline will run.

- `nimbus` sets the canu assembler to use 8 CPUs and 30GB RAM.
- `zeus` sets the canu assembler to use 14 CPUs and 64GB RAM, and sets some cluster specific options to use the slurm based scheduler at Pawsey.
- `docker` and `docker_sudo` sets it to use docker containers, `docker_sudo` is identical except that docker is run as root (required for some installations of docker).
- `singularity` sets it to use singularity containers, these are automatically converted from docker containers.


```bash
# Using docker
nextflow run ccdmb/nanopore_assembly -profile docker -resume -latest --nanoporeReads "*.fastq.gz"

# Using docker, and also increasing CPUs for canu.
nextflow run ccdmb/nanopore_assembly -profile docker,nimbus -resume -latest --nanoporeReads "*.fastq.gz"
```


## Parameters

```
--nanoporeReads <glob>
    Required
    A glob of the fastq.gz files of the adapter and barcode trimmed reads.
    The basename of the file needs to match the basename of the respective genome.

--canuSlow
    Default: false
    Disables canu fast mode.

--outdir <path>
    Default: `assembly`
    The directory to store the results in.
```
