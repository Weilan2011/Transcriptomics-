#!/bin/bash

# Define the base directory for rMATS
RMATS_DIR="$HOME/rmats-turbo"

# Define the base directory for your group files
GROUP_DIR="/home/weilan/MECP2/3_rMATS"

# Define the GTF file path
GTF_FILE="/home/weilan/HISAT2/HISAT2_Indexes/Homo_sapiens.GRCh38.111.gtf"

# Activate the rmats Conda environment
source activate rmats

# Define your conditions and time points
conditions=("CTR" "MDS" "MRL")
time_points=("30" "60")

# Loop through each condition for comparison between time points
for condition in "${conditions[@]}"; do
    # Define input files for both time points
    GROUP1_FILE="${GROUP_DIR}/${condition}_30.txt"
    GROUP2_FILE="${GROUP_DIR}/${condition}_60.txt"
    
    # Define output directories
    OUTPUT_DIR="${GROUP_DIR}/${condition}_30vs60"
    TMP_DIR="${OUTPUT_DIR}/tmp"
    
    # Ensure output and temporary directories exist
    mkdir -p "${OUTPUT_DIR}"
    mkdir -p "${TMP_DIR}"
    
    # Execute rmats.py with current condition comparison between time points
    python ${RMATS_DIR}/rmats.py --b1 "${GROUP1_FILE}" \
                                 --b2 "${GROUP2_FILE}" \
                                 --gtf "${GTF_FILE}" \
                                 -t paired \
                                 --novelSS \
                                 --readLength 100 \
                                 --variable-read-length \
                                 --nthread 16 \
                                 --od "${OUTPUT_DIR}" \
                                 --tmp "${TMP_DIR}"
    
    echo "Finished comparison: ${condition} 30 mins vs. 60 mins"
done

# Deactivate Conda environment
conda deactivate

echo "All condition comparisons completed."
