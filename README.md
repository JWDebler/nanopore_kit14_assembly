# Nanopore assembly pipeline for kit 14 duplex called reads

This is a draft nanopore assembly pipeline used by some groups at the CCDM.
More documentation to come.

## Examples

```bash
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
    A glob of the fastq.gz files of the adapter and barcode trimmed reads.
    The basename of the file needs to match the basename of the respective genome.

--canuSlow
    Default: false
    Disables canu fast mode.

--genomeSize
    Default: 42m
    Used by the assemblers to calculate read coverage

--medakaModel <glob>
        Default: r1041_e82_400bps_sup_v4.2.0 (kit114, sup, 5khz)
        The model that was used during basecalling.
        r1041_e82_400bps_sup_v4.1.0 (kit114, sup, 4khz)
        r941_min_sup_g507 (LSK109)

--outdir <path>
    Default: `assembly`
    The directory to store the results in.
```
