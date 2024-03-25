#!/bin/bash

# Define paths
RMATS_DIR="$HOME/rmats-turbo"
GROUP_DIR="/home/weilan/MECP2/3_rMATS"
GTF_FILE="/home/weilan/HISAT2/HISAT2_Indexes/Homo_sapiens.GRCh38.111.gtf"
FILE_LIST="/home/weilan/ENCODE/2_HISAT2/ENCODE_rMATS_K562_pairs_rMATS.txt"
BAM_DIR="/home/weilan/ENCODE/2_HISAT2/K562"  # Path to the BAM files

# Activate rmats environment
source activate rmats

# Process each line in the file
while IFS=$'\t' read -r target_files control_files; do
    declare -A gene_combinations

    # Split target files, group by cell line + gene
    IFS=', ' read -r -a targets <<< "$target_files"
    for target in "${targets[@]}"; do
        cell_line_gene="${target%_rep*}"
        gene_combinations["$cell_line_gene"]+="$BAM_DIR/$target,"
    done

    for combo in "${!gene_combinations[@]}"; do
        targets="${gene_combinations[$combo]}"
        targets=${targets%,}  # Remove trailing comma

        # Create temporary files for target and control BAM lists
        targets_file=$(mktemp)
        controls_file=$(mktemp)
        echo "${targets}" | tr ',' '\n' > "$targets_file"
        echo "$control_files" | sed "s/,/\\n/g" | sed "s/^/$BAM_DIR\//" > "$controls_file"

        OUTPUT_DIR="$GROUP_DIR/output/${combo}_vs_controls"
        mkdir -p "$OUTPUT_DIR"

        echo "Running rMATS for $combo vs controls"
        python $RMATS_DIR/rmats.py --b1 "$targets_file" --b2 "$controls_file" -t paired --gtf "$GTF_FILE" --od "$OUTPUT_DIR" --tmp "$OUTPUT_DIR/tmp" --variable-read-lengt --readLength 100 --nthread 16 

        # Cleanup
        rm "$targets_file" "$controls_file"
    done
done < "$FILE_LIST"








