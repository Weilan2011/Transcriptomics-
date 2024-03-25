#!/bin/bash

# Path to the directory where you want to check for files
DIRECTORY="/home/weilan/ENCODE/0_1_Raw"

# Path to the file containing the list of filenames to check
FILE_LIST="/home/weilan/ENCODE/0_1_Raw/ENCODE_set2_files.txt"

# Your output file paths
EXISTING_FILES="/home/weilan/ENCODE/0_1_Raw/existing_files.txt"
MISSING_FILES="/home/weilan/ENCODE/0_1_Raw/missing_files.txt"

# Clear previous outputs
> "$EXISTING_FILES"
> "$MISSING_FILES"

# Debugging: Echo directory being checked
echo "Checking in directory: $DIRECTORY"

# Read each line in FILE_LIST
while IFS= read -r filename; do
    full_path="$DIRECTORY/$filename"

    # Debugging: Echo full path being checked
    echo "Checking for file: $full_path"
    
    if [ -f "$full_path" ]; then
        echo "$filename" >> "$EXISTING_FILES"
    else
        echo "$filename" >> "$MISSING_FILES"
    fi
done < "$FILE_LIST"
