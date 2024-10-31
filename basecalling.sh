#!/bin/bash

#############################################
# This script is meant to duplex basecall a folder of multiplexed pod5 files and separate the raw data by barcode and channel.
# The output per barcode will be 3 files, 'barcodeX.duplex.fastq.gz', 'barcodeX.simplex.fastq.gz', and 'barcodeX.simplex.corrected.fasta.gz'
# Tools that need to be installed and in your path:
# - Dorado (min 0.7.1) (https://github.com/nanoporetech/dorado)
# - pod5 tools (https://github.com/nanoporetech/pod5-file-format)
# - samtools (sudo apt install samtools)
# - pigz (sudo apt install pigz)
#############################################
#
# Change these if needed
# barcoding kit used:
kit_name="SQK-NBD114-24" 
# quality filtering, currently disabled, you can always filter later using chopper
# quality="10" 
#############################################



# Check if the number of arguments is correct
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path_input> <path_output>"
    exit 1
fi

# Function to check if a command is available
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check if required commands are installed
required_commands=("dorado" "samtools" "pod5" "pigz")
missing_commands=()
for cmd in "${required_commands[@]}"; do
    if ! command_exists "$cmd"; then
        missing_commands+=("$cmd")
    fi
done

if [ ${#missing_commands[@]} -ne 0 ]; then
    echo "Error: The following required programs are not installed: ${missing_commands[*]}"
    echo "Please install them before running this script."
    exit 1
fi

# Check if the number of arguments is correct
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path_input> <path_output>"
    exit 1
fi

# Assign input arguments to variables
path_input="$1"
path_output="$2"

# Print the paths
echo "============================================================================"
echo "Input files located:     $path_input"
echo "Output files written to: $path_output"
echo "============================================================================"

if [ ! -d $path_output ]; then
    mkdir -p $path_output
fi

echo "============================================================================"
echo "$(date) - Starting SUP basecalling"
echo "============================================================================"

dorado basecaller hac -r $path_input --kit-name $kit_name > $path_output/all.bam 

echo "============================================================================"
echo "$(date) - Extracting split reads"
echo "============================================================================"

samtools view -@ $(nproc) -h $path_output/all.bam | grep -w "pi:Z" > $path_output/splitreads.bam

echo "============================================================================"
echo "$(date) - Starting read demultiplexing"
echo "============================================================================"

dorado demux -t $(nproc) --output-dir $path_output/simplex --no-classify $path_output/all.bam
dorado demux -t $(nproc) --output-dir $path_output/split_reads --no-classify $path_output/splitreads.bam
rm $path_output/simplex/*unclassified.bam
rm $path_output/split_reads/*unclassified.bam

for file in $path_output/split_reads/*.bam;
  do
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  samtools view -@ $(nproc) -h -O fastq $file >  $path_output/split_reads/$id.simplex.untrimmed.fastq;
  done

echo "============================================================================"
echo "$(date) - Generating barcode ID lists"
echo "============================================================================"

mkdir $path_output/simplex/tmp

for file in $path_output/simplex/*.bam; 
  do 
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  samtools view -@ $(nproc) -h $file | cut -f 1 > $path_output/simplex/tmp/$id.readids.txt; 
  done

for file in $path_output/simplex/tmp/*readids.txt; 
  do 
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  grep -v "@" $file > $path_output/simplex/tmp/$id.clean.txt; 
  done

mkdir $path_output/pod5_by_barcode

echo "============================================================================"
echo "$(date) - Separating POD5 files by barcode"
echo "============================================================================"

for file in $path_output/simplex/tmp/*clean.txt; 
  do 
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  pod5 filter -r $path_input -i $file -t $(nproc) --missing-ok --duplicate-ok --output $path_output/pod5_by_barcode/$id.pod5; 
  done

echo "============================================================================"
echo "$(date) - Extracting channel information from POD5 files"
echo "============================================================================"

for file in $path_output/pod5_by_barcode/*.pod5; 
  do 
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  pod5 view $file -t $(nproc) --include "read_id, channel" --output $path_output/simplex/tmp/$id.summary.tsv; 
  done

echo "============================================================================"
echo "$(date) - Extracting raw reads into channel specific POD5s"
echo "============================================================================"

for file in $path_output/simplex/tmp/*.summary.tsv; 
  do 
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  pod5 subset $path_output/pod5_by_barcode/$id.pod5 -t $(nproc) --summary $file --columns channel --output $path_output/pod5_by_barcode/$id; 
  rm $path_output/pod5_by_barcode/$id.pod5;
  done

echo "============================================================================"
echo "$(date) - Duplex calling of separated reads"
echo "============================================================================"

mkdir $path_output/duplex

for folder in $path_output/pod5_by_barcode/*; 
  do 
  id=$(basename $folder); 
  echo "============================================================================";
  echo "$(date) - Duplex calling $id";
  echo "============================================================================";
  dorado duplex sup -r $folder > $path_output/duplex/$id.duplex.untrimmed.bam; 
  done

echo "============================================================================"
echo "$(date) - Extracting Duplex Reads"
echo "============================================================================"

mkdir $path_output/duplex/tmp

for file in $path_output/duplex/*duplex.untrimmed.bam; 
  do 
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  samtools view -@ $(nproc) -h -O fastq -d dx:1 $file >  $path_output/duplex/tmp/$id.duplex.untrimmed.fastq;
  done

echo "============================================================================"
echo "$(date) - Extracting Simplex Reads"
echo "============================================================================"

for file in $path_output/duplex/*duplex.untrimmed.bam; 
  do 
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  samtools view -@ $(nproc) -h -O fastq -d dx:0 $file > $path_output/simplex/tmp/$id.simplex.untrimmed.fastq; 
  done

echo "============================================================================"
echo "$(date) - Merging Simplex reads with split Simplex reads"
echo "============================================================================"

for file in $path_output/split_reads/*.fastq;
  do
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  filename=$(basename $file);
  if [ -f $path_output/simplex/tmp/$filename ]; 
    then cat $file $path_output/simplex/tmp/$filename > $path_output/simplex/$id.simplex.untrimmed.fastq;
  else
    cp $file $path_output/simplex/;
  fi;
  done


echo "============================================================================"
echo "$(date) - Trimming All Reads"
echo "============================================================================"

for file in $path_output/duplex/tmp/*duplex.untrimmed.fastq;
  do
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  dorado trim -t $(nproc) --emit-fastq $file | pigz -9 >  $path_output/$id.duplex.fastq.gz
  done

for file in $path_output/simplex/*simplex.untrimmed.fastq;
  do
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  dorado trim -t $(nproc) --emit-fastq $file >  $path_output/$id.simplex.fastq
  done

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
echo "$(date) - Cleanup"
echo "============================================================================"

rm -r $path_output/simplex/
rm -r $path_output/duplex/
rm -r $path_output/split_reads/
rm $path_output/all.bam $path_output/splitreads.bam $path_output/*.fai

echo "============================================================================"
echo "$(date) - FINISHED"
echo "                 .                 "
echo "                ':'                "
echo "              ___:____     |'\/'|  "
echo "            ,'        '.    \  /   "
echo "            |  O        \___/  |   "
echo "          ~^~^~^~^~^~^~^~^~^~^~^~^~"
echo "============================================================================"