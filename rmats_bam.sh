#!/bin/bash

# Define the base directory for rMATS
RMATS_DIR="$HOME/rmats-turbo"

# Activate the rmats Conda environment
source activate rmats

# Define the base directory for your group files
GROUP_DIR="/home/weilan/MECP2/3_rMATS"

# Define the GTF and BAM file path
GTF_FILE="/home/weilan/HISAT2/HISAT2_Indexes/Homo_sapiens.GRCh38.111.gtf"

GROUP1_FILE="/home/weilan/MECP2/3_rMATS/CTR_30.txt" 
GROUP2_FILE="/home/weilan/MECP2/3_rMATS/CTR_untrimmed_30.txt" 

# Directory where you want the output to be stored
OUTPUT_DIR="/home/weilan/MECP2/3_rMATS/comparison_CTR30_trimed_vs_umtrimmed"
TMP_DIR="/home/weilan/MECP2/3_rMATS/comparison_CTR30_trimed_vs_umtrimmed/tmp"


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
    
echo "Finished comparison: CTR_30 trimmed vs. untrimmed"


# Deactivate Conda environment
conda deactivate

echo "All condition comparisons completed."