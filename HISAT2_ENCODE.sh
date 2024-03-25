#!/bin/bash

# Define the HISAT2 index path
REF_GENOME="/home/weilan/HISAT2/HISAT2_Indexes/GRCh38_reference_indexes"

# Directory containing your FASTQ files
FASTQ_DIR="/home/weilan/ENCODE/1_Fastp/HepG2_Trimed/Set1"

# Output directory for BAM files
OUTPUT_DIR="/home/weilan/ENCODE/2_HISAT2/"
mkdir -p "$OUTPUT_DIR"

# Number of cores to use
THREADS=16

# Extract unique base names without R1/R2 and the sequence part number
for filename in $(ls $FASTQ_DIR | grep -oP '.*(?=_R[12]_001.fastq.gz)' | sort -u); do
    echo "Processing sample: $filename"
    
    # Define input files for R1 and R2 reads
    R1_FILES="${FASTQ_DIR}/${filename}_R1_001.fastq.gz"
    R2_FILES="${FASTQ_DIR}/${filename}_R2_001.fastq.gz"
    
    # Define the output SAM and BAM file names
    SAM_FILE="${OUTPUT_DIR}/${filename}.sam"
    BAM_FILE="${OUTPUT_DIR}/${filename}.bam"
    
    # Run HISAT2 for alignment with multi-threading
    hisat2 -p $THREADS -x $REF_GENOME -1 $R1_FILES -2 $R2_FILES -S $SAM_FILE
    
    # Convert SAM to BAM, sort with multi-threading, and remove the original SAM file
    samtools view -@ $THREADS -b $SAM_FILE | samtools sort -@ $THREADS -o $BAM_FILE
    rm $SAM_FILE
    
    echo "Finished processing sample: $filename"
done

echo "All samples processed."
