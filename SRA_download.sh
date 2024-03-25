#!/bin/bash

# Define a list of SRA accession numbers
SRA_ACCESSIONS=("SRR6342225" "SRR6342226" "SRR6342227" "SRR6342228" "SRR6342229" "SRR6342230" "SRR6342231" "SRR6342232" "SRR6342233" "SRR6342234" "SRR6342235" "SRR6342236" "SRR6342237" "SRR6342238" "SRR6342239" "SRR6342240")

# Loop through each accession number in the list
for ACCESSION in "${SRA_ACCESSIONS[@]}"
do
    echo "Downloading FASTQ for accession: $ACCESSION"
    # Call fastq-dump with each accession number
    fastq-dump --split-files $ACCESSION
    # Optionally, specify additional options to fastq-dump if needed
    # For example, to download split files (one per read in paired-end data), use: fastq-dump --split-files $ACCESSION
done

echo "Download completed."