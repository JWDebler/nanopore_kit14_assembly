#!/bin/bash

#############################################
# This script is meant to  basecall a folder of multiplexed pod5 files and separate the raw data by barcode and channel.
# The output per barcode will be 3 files, 'barcodeX.simplex.fastq.gz', 'barcodeX.simplex.corrected.fasta.gz' and 'barcodeX.pod5.tar.gz'.
# The script will also output a checksum file 'checksums.md5' for file integrity.
# Tools that need to be installed and in your path:
# - Dorado (min 1.3.0 (https://github.com/nanoporetech/dorado)
# - pod5 tools (https://github.com/nanoporetech/pod5-file-format)
# - samtools (sudo apt install samtools)
# - pigz (sudo apt install pigz)
# - chopper (https://github.com/wdecoster/chopper)
#############################################
#
# Change these if needed
# barcoding kit used:
kit_name="SQK-NBD114-96" 
sampling_rate="5kHz"
flowcell_type="R10.4.1"
simplex_model="dna_r10.4.1_e8.2_400bps_sup@v5.2.0" # Change to latest simplex model
dorado_path="/opt/dorado/current/bin/dorado"

export sampling_rate flowcell_type simplex_model

#############################################
# quality filtering, currently disabled, you can always filter later using chopper
# quality="10" 
#############################################

# Function to create checkpoint
create_checkpoint() {
    local step_name="$1"
    local checkpoint_file="$path_output/.checkpoint_${step_name}"
    touch "$checkpoint_file"
    echo "============================================================================"
    echo "$(date) - Checkpoint created: $step_name"
}

# Function to check if checkpoint exists
checkpoint_exists() {
    local step_name="$1"
    local checkpoint_file="$path_output/.checkpoint_${step_name}"
    [ -f "$checkpoint_file" ]
}

# Function to skip step if checkpoint exists
skip_if_done() {
    local step_name="$1"
    local description="$2"
    if checkpoint_exists "$step_name"; then
        echo "============================================================================"
        echo "$(date) - SKIPPING: $description (checkpoint exists)"
        echo "============================================================================"
        return 0
    else
        echo "============================================================================"
        echo "$(date) - Starting: $description"
        echo "============================================================================"
        return 1
    fi
}

# Function to display help message
display_help() {
    echo "Usage: $0 [--rapid] [--methylation] [--resume] [--clean-checkpoints] [--sample-map <file>] <path_input> <path_output>"
    echo
    echo "Options:"
    echo "  --rapid              Use the rapid barcoding kit (SQK-RBK114-24)"
    echo "  --methylation        Perform methylation calling using 5mC_5hmC,6mA models"
    echo "  --resume             Resume from existing checkpoints (default behavior)"
    echo "  --clean-checkpoints  Remove all checkpoint files and start fresh"
    echo "  --sample-map <file>  Text file mapping barcodeXX to sample names (two columns)"
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

clean_checkpoints=false
methylation=""
sample_map_file=""
declare -A barcode_to_sample

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --rapid)
            kit_name="SQK-RBK114-24"
            shift
            ;;
        --methylation)
            methylation=",5mC_5hmC,6mA"
            shift
            ;;
        --resume)
            shift ;; # Default behavior, just consume the flag
        --clean-checkpoints)
            clean_checkpoints=true
            shift
            ;;
        --sample-map)
            sample_map_file="$2"
            shift 2
            ;;
        *)
            break
            ;;
    esac
done

# Check if the number of arguments is correct
if [ "$#" -ne 2 ]; then
    display_help
fi

# Parse sample map if provided
if [ -n "$sample_map_file" ]; then
    if [ ! -f "$sample_map_file" ]; then
        echo "Error: Sample map file '$sample_map_file' not found!"
        exit 1
    fi
    
    echo "============================================================================"
    echo "$(date) - Parsing sample map file: $sample_map_file"
    echo "============================================================================"
    
    # Handle different line endings (DOS/Windows, Unix/Linux, old Mac)
    # Convert to Unix line endings and remove any trailing whitespace
    temp_map_file=$(mktemp)
    sed 's/\r$//' "$sample_map_file" | sed 's/[[:space:]]*$//' > "$temp_map_file"
    
    sample_count=0
    line_number=0
    
    while IFS=$'\t' read -r barcode sample rest || [ -n "$barcode" ]; do
        line_number=$((line_number + 1))
        
        # Skip empty lines
        [ -z "$barcode" ] && continue
        
        # Strip leading/trailing whitespace from both fields
        barcode=$(echo "$barcode" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        sample=$(echo "$sample" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        
        # Skip lines that start with # (comments)
        [[ "$barcode" =~ ^[[:space:]]*# ]] && continue
        
        # Validate barcode format and ensure sample is non-empty
        if [[ "$barcode" =~ ^barcode[0-9][0-9]$ ]] && [ -n "$sample" ]; then
            barcode_to_sample[$barcode]="$sample"
            echo "  Line $line_number: $barcode -> $sample"
            sample_count=$((sample_count + 1))
        else
            echo "  WARNING: Line $line_number skipped - invalid format or empty sample"
            echo "    Content: '$barcode' -> '$sample'"
            if [[ ! "$barcode" =~ ^barcode[0-9][0-9]$ ]]; then
                echo "    Issue: Barcode must be in format 'barcodeXX' (e.g., barcode01, barcode12)"
            fi
            if [ -z "$sample" ]; then
                echo "    Issue: Sample name cannot be empty"
            fi
        fi
    done < "$temp_map_file"
    
    # Clean up temporary file
    rm "$temp_map_file"
    
    if [ $sample_count -eq 0 ]; then
        echo "ERROR: No valid sample mappings found in $sample_map_file"
        echo "Expected format (tab-separated):"
        echo "barcode01<TAB>sample_name_1"
        echo "barcode02<TAB>sample_name_2"
        exit 1
    fi
    
    echo "============================================================================"
    echo "Successfully parsed $sample_count sample mappings from $sample_map_file"
    echo "============================================================================"
else
    echo "============================================================================"
    echo "No sample map file provided - barcodes will be used as identifiers"
    echo "============================================================================"
fi

# Function to check if a command is available
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check if required commands are installed
required_commands=("$dorado_path" "samtools" "pod5" "pigz")
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

export path_output

# Create output directory if it doesn't exist
if [ ! -d "$path_output" ]; then
    mkdir -p "$path_output"
fi

# Clean checkpoints if requested
if [ "$clean_checkpoints" = true ]; then
    echo "============================================================================"
    echo "$(date) - Cleaning all checkpoint files"
    echo "============================================================================"
    rm -f "$path_output"/.checkpoint_*
fi

# Input validation: Check input directory for .pod5 files
if [ ! -d "$path_input" ]; then
    echo "Error: Input directory '$path_input' does not exist!"
    exit 1
fi
if ! find "$path_input" -type f -name "*.pod5" | grep -q .; then
    echo "Error: No .pod5 files found in input directory '$path_input' (searched recursively)!"
    exit 1
fi

# Warn if output directory is not empty
if [ -d "$path_output" ] && [ "$(ls -A $path_output)" ]; then
    echo "Warning: Output directory '$path_output' is not empty. Files may be overwritten or mixed."
fi

# Print the paths
echo "============================================================================"
echo "Input files located:     $path_input"
echo "Using barcoding kit:     $kit_name"
echo "Calling methylation:     $methylation"
echo "Output files written to: $path_output"
echo "============================================================================"

# STEP 1: SUP basecalling
if ! skip_if_done "basecalling" "SUP basecalling"; then
    # Capture existing directories before basecalling
    existing_dirs=$(find "$path_output" -maxdepth 1 -type d ! -path "$path_output" -printf "%f\n" 2>/dev/null | sort)
    
    $dorado_path basecaller sup$methylation -r "$path_input" --kit-name "$kit_name" -o "$path_output"
    
    # Capture new directories after basecalling
    new_dirs=$(find "$path_output" -maxdepth 1 -type d ! -path "$path_output" -printf "%f\n" 2>/dev/null | sort)
    
    # Store the difference (new directories created by dorado)
    dorado_created_dirs=$(comm -13 <(echo "$existing_dirs") <(echo "$new_dirs"))
    echo "$dorado_created_dirs" > "$path_output/.dorado_created_dirs"
    
    create_checkpoint "basecalling"
fi

#STEP 2: Merging BAM files
if ! skip_if_done "merging_bam" "merging BAM files"; then

    # First, collect all unique barcode IDs across all FLOWCELLID directories
    unique_barcodes=$(find "$path_output" -type d -name "barcode[0-9][0-9]" | grep "/bam_pass/" | sed 's/.*\/\(barcode[0-9][0-9]\)$/\1/' | sort -u)
    
    # Process each unique barcode
    for barcode_id in $unique_barcodes; do
        echo "Processing barcode: $barcode_id"
        
        # Find ALL bam files for this barcode across ALL FLOWCELLID directories
        bam_files=$(find "$path_output" -path "*/bam_pass/${barcode_id}/*.bam" | tr '\n' ' ')
        
        if [ -n "$bam_files" ]; then
            echo "Merging bam files for $barcode_id from multiple flowcells:"
            echo "  Files: $bam_files"
            samtools merge -@ $(nproc) -o "$path_output/${barcode_id}.${simplex_model}${methylation}.bam" $bam_files
            echo "  Successfully merged $barcode_id"
        else
            echo "No bam files found for $barcode_id"
        fi
    done
    
    # Clean up the original bam files after merging
    find "$path_output" -name "*.bam" -path "*/bam_pass/*" -delete
    
    create_checkpoint "merging_bam"
fi

# STEP 3: Generating barcode ID lists
if ! skip_if_done "barcode_ids" "generating barcode ID lists"; then

    mkdir -p "$path_output/simplex/tmp"
    
    for file in "$path_output/"*.bam; do 
        id=$(echo "$file" | grep -oP 'barcode\d+'); 
        samtools view -@ $(nproc) -h "$file" | cut -f 1 > "$path_output/simplex/tmp/$id.readids.txt"; 
    done
    
    # Parallelise this bit to avoid bottlenecking
    ls "$path_output/simplex/tmp/"*readids.txt | xargs -P $(nproc) -I {} bash -c ' 
        file="$1"
        path_output="$2"
        id=$(basename "$file" | grep -oP "barcode\\d+")
        grep -v "@" "$file" > "$path_output/simplex/tmp/$id.clean.txt"
    ' _ {} "$path_output"
    
    create_checkpoint "barcode_ids"
fi

# STEP 4: Extracting Simplex Reads
if ! skip_if_done "simplex_extraction" "extracting Simplex Reads"; then

    for file in "$path_output/"*.bam; do 
        id=$(echo "$file" | grep -oP 'barcode\d+');  
        $dorado_path demux -t $(nproc) --output-dir "$path_output/simplex/" --no-classify "$file" --emit-fastq
    done

    # Find all fastq files with barcodes in any subdirectory and merge by barcode
    for barcode in $(find "$path_output/simplex" -name "*barcode[0-9][0-9]*" -name "*.fastq*" | grep -oP 'barcode\d+' | sort -u); do
        echo "Merging files for $barcode"
        find "$path_output" -name "*${barcode}*" -name "*.fastq*" -exec cat {} \; > "$path_output/${barcode}.${simplex_model}.fastq"
    done

    create_checkpoint "simplex_extraction"
fi

# STEP 5: Simplex correction
if ! skip_if_done "simplex_correction" "simplex correction"; then
    # Download model if not already present
    if ! checkpoint_exists "model_download"; then
        dorado download --model herro-v1
        create_checkpoint "model_download"
    fi
    
    # Create temporary directory
    mkdir -p "$path_output/tmp_correction"
    
    for file in "$path_output/"*.fastq; do
        id=$(echo "$file" | grep -oP 'barcode\d+');
        
        final_file="$path_output/${id}.${simplex_model}.simplex.corrected.fasta.gz"
        temp_file="$path_output/tmp_correction/${id}.simplex.corrected.fasta"
        temp_gz_file="$path_output/tmp_correction/${id}.simplex.corrected.fasta.gz"
        
        # Skip if final corrected file already exists and is not empty
        if [[ -f "$final_file" && -s "$final_file" ]]; then
            echo "$(date) - Skipping correction for $id (already exists)"
            continue
        fi
        
        # Clean up any existing temporary files for this barcode
        rm -f "$temp_file" "$temp_gz_file"
        
        echo "============================================================================";
        echo "$(date) - Correcting $id";
        echo "============================================================================";
        
        # dorado correct defaults have changed as of 0.9.5, so now we have to supply these as additional arguments in order not to lose almost all reads
        if $dorado_path correct -m herro-v1 "$file" --min-chain-score 4000 --mid-occ-frac 0.0002 > "$temp_file"; then
            # Only proceed if dorado correct succeeded
            if [[ -f "$temp_file" && -s "$temp_file" ]]; then
                # Compress the corrected file
                if pigz -9 "$temp_file"; then
                    # Move to final location only if compression succeeded
                    mv "$temp_gz_file" "$final_file"
                    echo "$(date) - Successfully completed correction for $id"
                    
                    # Compress original fastq file
                    pigz -9 "$file"
                else
                    echo "$(date) - Failed to compress corrected file for $id"
                    rm -f "$temp_file" "$temp_gz_file"
                fi
            else
                echo "$(date) - Correction produced empty file for $id"
                rm -f "$temp_file"
            fi
        else
            echo "$(date) - dorado correct failed for $id"
            rm -f "$temp_file"
        fi
    done
    
    # Clean up temporary directory
    rmdir "$path_output/tmp_correction" 2>/dev/null || true
    
    create_checkpoint "simplex_correction"
fi


# STEP 6: Separating POD5 files by barcode
if ! skip_if_done "pod5_separation" "separating POD5 files by barcode"; then

    mkdir -p "$path_output/pod5_by_barcode"
    
    for file in "$path_output/simplex/tmp/"*.clean.txt; do 
        id=$(echo "$file" | grep -oP 'barcode\d+')
        output_file="$path_output/pod5_by_barcode/$id.pod5"
        temp_file="$path_output/pod5_by_barcode/$id.pod5.tmp"
        
        # Check if the final POD5 file already exists and is not empty
        if [[ -f "$output_file" && -s "$output_file" ]]; then
            echo "POD5 file already exists for $id, skipping..."
            continue
        fi
        
        # Remove any existing temporary file (from interrupted run)
        if [[ -f "$temp_file" ]]; then
            echo "Removing incomplete temporary file for $id..."
            rm "$temp_file"
        fi
        
        echo "Processing $id..."
        if pod5 filter -r "$path_input" -i "$file" -t $(nproc) --missing-ok --duplicate-ok --output "$temp_file" > /dev/null 2>&1; then
            # Only rename to final name if pod5 filter succeeded
            mv "$temp_file" "$output_file"
            echo "Completed $id"
        else
            echo "Failed to process $id, removing temporary file..."
            rm -f "$temp_file"
        fi
    done
    
    create_checkpoint "pod5_separation"
fi

# STEP 7: Archiving Raw Reads
if ! skip_if_done "archiving" "archiving Raw Reads"; then

    for file in "$path_output/pod5_by_barcode/"*.pod5; do
        # Check if the glob matched any files
        if [ ! -f "$file" ]; then
            echo "No .pod5 files found in $path_output/pod5_by_barcode/"
            break
        fi
        
        echo "Processing: $file"
        id=$(echo "$file" | grep -oP 'barcode\d+'); 
        echo "Barcode ID: $id"
        
        # Skip if archive already exists
        if [ -f "$path_output/${id}.${sampling_rate}.${flowcell_type}.pod5.tar.gz" ]; then
            echo "$(date) - Skipping archive for $id (already exists)"
            continue
        fi
        
        tar -cf - "$file" | pigz -0 - > "$path_output/${id}.${sampling_rate}.${flowcell_type}.pod5.tar.gz"
        echo "Created archive: $path_output/${id}.${sampling_rate}.${flowcell_type}.pod5.tar.gz"
        echo "Removing: $file"
        rm "$file"
    done
    
    create_checkpoint "archiving"
fi

# STEP 8: Creating file integrity checksum files
if ! skip_if_done "checksums" "creating file integrity checksum files"; then

    # If sample mapping is provided, rename files before checksums
    if [ -n "$sample_map_file" ]; then
        for barcode in "${!barcode_to_sample[@]}"; do
            sample="${barcode_to_sample[$barcode]}"
            # Rename all relevant files in output directory
            for ext in ".bam" ".fastq.gz" ".simplex.corrected.fasta.gz" ".pod5.tar.gz"; do
                for f in "$path_output/"${barcode}*${ext}; do
                    [ -e "$f" ] || continue
                    newf="${f//$barcode/$sample}"
                    if [ "$f" != "$newf" ]; then
                        mv "$f" "$newf"
                    fi
                done
            done
        done
    fi

    cd "$path_output"
    # Remove existing checksums file to avoid duplicates
    rm -f checksums.md5
    find . -type f \( -name "*.bam" -o -name "*.gz" \) -print0 | xargs -0 md5sum >> checksums.md5
    cd - > /dev/null

    # Output file manifest: list all output files and their sample associations
    manifest_file="$path_output/output_manifest.tsv"
    echo -e "output_file\tsample_or_barcode" > "$manifest_file"
    for f in "$path_output"/*.{bam,fastq.gz,simplex.corrected.fasta.gz,pod5.tar.gz}; do
        [ -e "$f" ] || continue
        fname=$(basename "$f")
        # Try to extract sample or barcode from filename
        if [ -n "$sample_map_file" ]; then
            found=0
            for barcode in "${!barcode_to_sample[@]}"; do
                sample="${barcode_to_sample[$barcode]}"
                if [[ "$fname" == *"$sample"* ]]; then
                    echo -e "$fname\t$sample" >> "$manifest_file"
                    found=1
                    break
                fi
            done
            if [ $found -eq 0 ]; then
                # Fallback: try to extract barcode
                if [[ "$fname" =~ (barcode[0-9][0-9]) ]]; then
                    echo -e "$fname\t${BASH_REMATCH[1]}" >> "$manifest_file"
                else
                    echo -e "$fname\tUNKNOWN" >> "$manifest_file"
                fi
            fi
        else
            if [[ "$fname" =~ (barcode[0-9][0-9]) ]]; then
                echo -e "$fname\t${BASH_REMATCH[1]}" >> "$manifest_file"
            else
                echo -e "$fname\tUNKNOWN" >> "$manifest_file"
            fi
        fi
    done

    create_checkpoint "checksums"
fi

# STEP 9: Cleanup
if ! skip_if_done "cleanup" "cleanup"; then

    rm -r "$path_output/simplex/"
    rm -r "$path_output/pod5_by_barcode/"
    rm -r "herro-v1"
    rm -rf "$path_output/tmp_model/"
    rm -rf .temp_dorado_model*/
    rm "$path_output/"*.fai 

    create_checkpoint "cleanup"
fi

# Summary report
num_samples=0
if [ -n "$sample_map_file" ]; then
    num_samples=${#barcode_to_sample[@]}
    echo "============================================================================"
    echo "Summary: $num_samples samples processed (from sample map)"
    echo "Sample list:"
    for barcode in "${!barcode_to_sample[@]}"; do
        echo "  $barcode -> ${barcode_to_sample[$barcode]}"
    done
else
    # Count unique barcodes in output
    num_samples=$(ls "$path_output" | grep -oP 'barcode[0-9][0-9]' | sort -u | wc -l)
    echo "============================================================================"
    echo "Summary: $num_samples barcodes processed (no sample map)"
    echo "Barcodes:"
    ls "$path_output" | grep -oP 'barcode[0-9][0-9]' | sort -u
fi
if [ -f "$path_output/output_manifest.tsv" ]; then
    echo "Output file manifest written to: $path_output/output_manifest.tsv"
fi

echo "============================================================================"
echo "$(date) - FINISHED"
echo "                 .                 "
echo "                ':'                "
echo "              ___:____     |'\/'|  "
echo "            ,'        '.    \  /   "
echo "            |  O        \___/  |   "
echo "          ~^~^~^~^~^~^~^~^~^~^~^~^~"
echo "============================================================================"
echo "All checkpoint files are stored in: $path_output/.checkpoint_*"
echo "To restart from scratch, use: $0 --clean-checkpoints [other options]"
echo "============================================================================"