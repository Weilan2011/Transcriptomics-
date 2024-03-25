#!/bin/bash

# Define the base directory for rMATS
RMATS_DIR="$HOME/rmats-turbo"

# Define the GTF file path
GTF_FILE="/home/weilan/HISAT2/HISAT2_Indexes/Homo_sapiens.GRCh38.111.gtf"

# Output directory for rMATS analysis
OUTPUT_DIR="/home/weilan/ENCODE/3_rMATS/Set2"

# Activate the rmats Conda environment
source activate rmats

# Iterate over all control BAM list files in the output directory
for control_bam_list in $OUTPUT_DIR/*_control_bams.txt; do
    # Extract the group name by removing the directory path and file suffix
    group_name=$(basename "$control_bam_list" _control_bams.txt)

    # Define the corresponding target BAM list
    target_bam_list="$OUTPUT_DIR/${group_name}_target_bams.txt"

    # Run rMATS for the group
    python ${RMATS_DIR}/rmats.py --b1 "$control_bam_list" \
                                 --b2 "$target_bam_list" \
                                 --gtf "${GTF_FILE}" \
                                 -t paired \
                                 --novelSS \
                                 --readLength 100 \
                                 --variable-read-length \
                                 --nthread 16 \
                                 --od "$OUTPUT_DIR/${group_name}_output" \
                                 --tmp "$OUTPUT_DIR/${group_name}_output/tmp"

    echo "Finished comparison: $group_name"
done

echo "All rMATS analyses are completed."