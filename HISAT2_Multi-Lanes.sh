#!/bin/bash

# Define the HISAT2 index path
REF_GENOME="/home/weilan/HISAT2/HISAT2_Indexes/GRCh38_reference_indexes"

# Directory containing your FASTQ files
FASTQ_DIR="/home/weilan/MECP2/0_Raw"

# Output directory for BAM files
OUTPUT_DIR="/home/weilan/MECP2/2_HISAT2/Untrimmed/"
mkdir -p $OUTPUT_DIR

# Number of cores to use
THREADS=16

# Loop through each unique sample ID in your FASTQ directory
for sample in $(ls $FASTQ_DIR | grep -oP 'IJ-\d+_S\d+' | sort -u); do
    echo "Processing sample: $sample"
    
    # Define input files for both lanes and both reads, taking into account all possible lanes
    R1_FILES=$(ls $FASTQ_DIR/${sample}_L00{1,2}_R1_*.fastq.gz 2>/dev/null | tr '\n' ',' | sed 's/,$//')
    R2_FILES=$(ls $FASTQ_DIR/${sample}_L00{1,2}_R2_*.fastq.gz 2>/dev/null | tr '\n' ',' | sed 's/,$//')

    # Check if files were found
    if [[ -z $R1_FILES ]] || [[ -z $R2_FILES ]]; then
        echo "No files found for $sample. Skipping..."
        continue
    fi

    # Define the output SAM and BAM file names
    SAM_FILE="$OUTPUT_DIR/${sample}.sam"
    BAM_FILE="$OUTPUT_DIR/${sample}.bam"
    
    # Run HISAT2 for alignment with multi-threading
    hisat2 -p $THREADS -x $REF_GENOME -1 $R1_FILES -2 $R2_FILES -S $SAM_FILE
    
    # Convert SAM to BAM, sort with multi-threading, and remove the original SAM file
    samtools view -@ $THREADS -bS $SAM_FILE | samtools sort -@ $THREADS -o $BAM_FILE
    rm $SAM_FILE
    
    echo "Finished processing sample: $sample"
done

echo "All samples processed."
