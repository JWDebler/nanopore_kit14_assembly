#!/bin/bash

##############################################
# This script is meant to duplex basecall a folder of multiplexed pod5 files and separate the raw data by barcode and channel.
# The output per barcode will be 4 files, 'barcodeX.duplex.fastq.gz', 'barcodeX.simplex.fastq.gz', 'barcodeX.simplex.corrected.fasta.gz' and 'barcodeX.pod5.tar.gz'.
# The script will also output a checksum file 'checksums.md5' for file integrity.
# Tools that need to be installed and in your path:
# - Dorado (min 0.7.1) (https://github.com/nanoporetech/dorado)
# - pod5 tools (https://github.com/nanoporetech/pod5-file-format)
# - samtools
# - pigz (sudo apt install pigz)
# - chopper (https://github.com/wdecoster/chopper)
#############################################
#
# Change these if needed
# barcoding kit used:
kit_name="SQK-NBD114-24"
sampling_rate="5kHz"
flowcell_type="R10.4.1"
duplex_model="dna_r10.4.1_e8.2_5khz_stereo@v1.3"   # Change to latest duplex model
simplex_model="dna_r10.4.1_e8.2_400bps_sup@v5.0.0" # Change to latest simplex model
#############################################
# quality filtering, currently disabled, you can always filter later using chopper
# quality="10"
#############################################

# Function to display help message
display_help() {
    echo "Usage: $0 [--rapid] <path_input> <path_output>"
    echo
    echo "Options:"
    echo "  --rapid       Use the rapid barcoding kit (SQK-RBK114-24)"
    echo
    echo "Arguments:"
    echo "  <path_input>  Path to the input folder containing pod5 files"
    echo "  <path_output> Path to the output folder where results will be stored"
    exit 1
}

# Parse command line arguments
if [ "$#" -eq 0 ]; then
    display_help
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --rapid) kit_name="SQK-RBK114-24"; shift ;;
        *) break ;;
    esac
done

# Check if the number of arguments is correct
if [ "$#" -ne 2 ]; then
    display_help
fi

# Function to check if a command is available
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check if required commands are installed
required_commands=("dorado" "samtools" "pod5" "pigz" "chopper")
missing_commands=()
for cmd in "${required_commands[@]}"; do
    if ! command_exists "$cmd"; then
        missing_commands+=("$cmd")
    fi
done

if [ ${#missing_commands[@]} -ne 0 ]; then
    echo "Error: The following required programs are not installed or not in your PATH: ${missing_commands[*]}"
    echo "Please install them before running this script."
    exit 1
fi

# Assign input arguments to variables
path_input="$1"
path_output="$2"

# Function to log execution
log_file() {
    echo "$(date) - $1" >> "$log_dir/$2.log"
}

# Print the paths
echo "============================================================================"
echo "Input files located:     $path_input"
echo "Using barcoding kit:     $kit_name"
echo "Output files written to: $path_output"
echo "============================================================================"

if [ ! -d $path_output ]; then
    mkdir -p $path_output
    mkdir -p "$path_output/logs"
fi

# Create output and log directories if they don't exist
log_dir="$path_output/logs"

echo "============================================================================"
echo "$(date) - Starting HAC basecalling"
echo "============================================================================"

dorado basecaller hac -r "$path_input" --kit-name "$kit_name" > "$path_output/all.bam"

echo "============================================================================"
echo "$(date) - Extracting split reads"
echo "============================================================================"

samtools view -@ $(nproc) -h "$path_output/all.bam" | grep -w "pi:Z" > "$path_output/splitreads.bam"

echo "============================================================================"
echo "$(date) - Starting read demultiplexing"
echo "============================================================================"

dorado demux -t $(nproc) --output-dir $path_output/simplex --no-classify $path_output/all.bam
dorado demux -t $(nproc) --output-dir $path_output/split_reads --no-classify $path_output/splitreads.bam
rm $path_output/simplex/*unclassified.bam
rm $path_output/split_reads/*unclassified.bam

find "$path_output/split_reads" -name "*.bam" | xargs -P $(nproc) -I {} sh -c '
  id=$(echo "$1" | grep -oP "barcode\d+");
  samtools view -@ $(nproc) -h -O fastq "$1" >> "'"$path_output"'/split_reads/$id.simplex.untrimmed.fastq" 2>&1 | tee "'"$log_dir"'/demultiplex_${id}.log"
' _ {}

echo "============================================================================"
echo "$(date) - Generating barcode ID lists"
echo "============================================================================"

mkdir -p $path_output/simplex/tmp

find "$path_output/simplex" -name "*.bam" | xargs -P $(nproc) -I {} sh -c '
  id=$(echo "$1" | grep -oP "barcode\d+");
  samtools view -@ $(nproc) -h "{}" | cut -f 1 >> "'"$path_output"'/simplex/tmp/$id.readids.txt" 2>&1 | tee "'"$log_dir"'/extract_readids_$id.log"
' _ {}

find "$path_output/simplex/tmp" -name "*readids.txt" | xargs -P $(nproc) -I {} sh -c '
  id=$(echo "$1" | grep -oP "barcode\d+");
  grep -v "@" "$1" > "'"$path_output"'/simplex/tmp/$id.clean.txt" 2>&1 | tee "'"$log_dir"'/clean_readids_$id.log"
' _ {}

echo "============================================================================"
echo "$(date) - Separating POD5 files by barcode"
echo "============================================================================"

mkdir -p $path_output/pod5_by_barcode

find "$path_output/simplex/tmp" -name "*clean.txt" | xargs -P $(nproc) -I {} sh -c '
  id=$(echo "$1" | grep -oP "barcode\d+");
  pod5 filter -r "'"$path_input"'" -i "$1" -t $(nproc) --missing-ok --duplicate-ok --output "'"$path_output"'/pod5_by_barcode/$id.pod5" 2>&1 | tee "'"$log_dir"'/pod5_separation_$id.log"
' _ {}

echo "============================================================================"
echo "$(date) - Extracting channel information from POD5 files"
echo "============================================================================"

find "$path_output/pod5_by_barcode" -name "*.pod5" | xargs -P $(nproc) -I {} sh -c '
  id=$(echo "$1" | grep -oP "barcode\d+");
  pod5 view "$1" -t $(nproc) --include "read_id, channel" --output "'"$path_output"'/simplex/tmp/$id.summary.tsv" 2>&1 | tee "'"$log_dir"'/channel_extraction_$id.log"
' _ {}

echo "============================================================================"
echo "$(date) - Extracting raw reads into channel specific POD5s"
echo "============================================================================"

find "$path_output/simplex/tmp" -name "*.summary.tsv" | xargs -P $(nproc) -I {} sh -c '
  id=$(echo "$1" | grep -oP "barcode\d+");
  pod5 subset "'"$path_output"'/pod5_by_barcode/$id.pod5" -t $(nproc) --summary "$1" --columns channel --output "'"$path_output"'/pod5_by_barcode/$id";
  rm "'"$path_output"'/pod5_by_barcode/$id.pod5" 2>&1 | tee "'"$log_dir"'/raw_read_extraction_$id.log"
' _ {}

echo "============================================================================"
echo "$(date) - Duplex calling of separated reads"
echo "============================================================================"

mkdir -p $path_output/duplex
mkdir -p $path_output/tmp_model

for folder in $path_output/pod5_by_barcode/*;
  do
  id=$(basename $folder);
  echo "============================================================================";
  echo "$(date) - Duplex calling $id";
  echo "============================================================================";
  dorado duplex sup -r $folder --models-directory $path_output/tmp_model > $path_output/duplex/$id.duplex.untrimmed.bam;
  done

echo "============================================================================"
echo "$(date) - Extracting Duplex Reads"
echo "============================================================================"

mkdir -p $path_output/duplex/tmp

find "$path_output/duplex" -name "*duplex.untrimmed.bam" | xargs -P $(nproc) -I {} sh -c '
  id=$(echo "$1" | grep -oP "barcode\d+");
  samtools view -@ $(nproc) -h -O fastq -d dx:1 "$1" > "'"$path_output"'/duplex/tmp/$id.duplex.untrimmed.fastq" 2>&1 | tee "'"$log_dir"'/duplex_read_extraction_$id.log"
' _ {}

echo "============================================================================"
echo "$(date) - Extracting Simplex Reads"
echo "============================================================================"

find "$path_output/duplex" -name "*duplex.untrimmed.bam" | xargs -P $(nproc) -I {} sh -c '
  id=$(echo "$1" | grep -oP "barcode\d+");
  samtools view -@ $(nproc) -h -O fastq -d dx:0 "$1" > "'"$path_output"'/simplex/tmp/$id.simplex.untrimmed.fastq" 2>&1 | tee "'"$log_dir"'/simplex_read_extraction_$id.log"
' _ {}

echo "============================================================================"
echo "$(date) - Merging Simplex reads with split Simplex reads"
echo "============================================================================"

find "$path_output/split_reads" -name "*.fastq" | xargs -P $(nproc) -I {} sh -c '
  id=$(echo "$1" | grep -oP "barcode\\d+");
  filename=$(basename "$1");
  log_file="'"$log_dir"'/merge_${id}.log"

  echo "Processing $id - $(date)" > "$log_file"

  if [ -f "'"$path_output"'/simplex/tmp/${filename}" ]; then
    echo "Merging ${filename} with tmp/${filename}" >> "$log_file"
    cat "$1" "'"$path_output"'/simplex/tmp/${filename}" > "'"$path_output"'/simplex/${id}.simplex.untrimmed.fastq"
  else
    echo "No matching file in tmp/, copying ${filename} to simplex/" >> "$log_file"
    cp "$1" "'"$path_output"'/simplex/"
  fi

  echo "Completed $id - $(date)" >> "$log_file"
' _ {}

echo "============================================================================"
echo "$(date) - Trimming All Reads"
echo "============================================================================"

# Had to revert to hard trimming with chopper because 'dorado trim' just lets too many reads with adapter
# and barcode sequences through and they mess with my assemblies. This might change in the future.

find "$path_output/duplex/tmp" -name "*duplex.untrimmed.fastq" | xargs -P $(nproc) -I {} sh -c '
  id=$(echo "$1" | grep -oP "barcode\d+");
  chopper --minlength 300 --headcrop 75 --tailcrop 75 -t $(nproc) -i "$1" | pigz -9 > "'"$path_output"'/${id}.duplex.fastq.gz" 2>&1 | tee "'"$log_dir"'/trimming_duplex_reads_$id.log"
' _ {}

find "$path_output/simplex" -name "*simplex.untrimmed.fastq" | xargs -P $(nproc) -I {} sh -c '
  id=$(echo "$1" | grep -oP "barcode\d+");
  chopper --minlength 300 --headcrop 75 --tailcrop 75 -t $(nproc) -i "$1" > "'"$path_output"'/${id}.simplex.fastq" 2>&1 | tee "'"$log_dir"'/trimming_simplex_reads_$id.log"
' _ {}

echo "============================================================================"
echo "$(date) - Simplex correction"
echo "============================================================================"

for file in $path_output/*.simplex.fastq;
  do
  id=$(echo "$file" | grep -oP 'barcode\d+');
  echo "============================================================================";
  echo "$(date) - Correcting $id";
  echo "============================================================================";
  dorado correct $file > $path_output/$id.simplex.corrected.fasta;
  pigz -9 $file;
  pigz -9 $path_output/$id.simplex.corrected.fasta
  done

echo "============================================================================"
echo "$(date) - Archiving Raw Reads"
echo "============================================================================"

find "$path_output/pod5_by_barcode" -mindepth 1 -maxdepth 1 -type d | xargs -P $(nproc) -I {} sh -c '
  id=$(basename "{}");
  tar -cf - "{}" | pigz -0 - > "'"$path_output"'/$id.'"$sampling_rate"'.'"$flowcell_type"'.pod5.tar.gz" 2>&1 | tee "'"$log_dir"'/compress_pod5_$id.log"
' _ {}


echo "============================================================================"
echo "$(date) - Cleanup"
echo "============================================================================"

rm -r $path_output/simplex/
rm -r $path_output/duplex/
rm -r $path_output/split_reads/
rm -r $path_output/pod5_by_barcode/
rm -r $path_output/tmp_model/
rm $path_output/all.bam $path_output/splitreads.bam $path_output/*.fai

echo "============================================================================"
echo "$(date) - Creating file integrity checksum files"
echo "============================================================================"

cd $path_output
find . -type f -print0 | xargs -0 -P $(nproc) md5sum >> checksums.md5
cd ..

echo "============================================================================"
echo "$(date) - FINISHED"
echo "                 .                 "
echo "                ':'                "
echo "              ___:____     |'\/'|  "
echo "            ,'        '.    \  /   "
echo "            |  O        \___/  |   "
echo "          ~^~^~^~^~^~^~^~^~^~^~^~^~"
echo "============================================================================"

# Completion
log_file "FINISHED" "process"
