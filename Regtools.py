import os
import glob
import subprocess
import pandas as pd

# Define input directory containing SAM files
input_dir = "/home/weilan/STAR/subset/Regtools/HISAT2_regtools/SAM"

# Step 1: Convert SAM to BAM and index BAM files
sam_files = glob.glob(os.path.join(input_dir, "*.sam"))
for sam_file in sam_files:
    base = os.path.basename(sam_file).replace(".sam", "")
    bam_file = os.path.join(input_dir, base + ".bam")
    
    # Check if BAM file already exists
    if os.path.exists(bam_file):
        print(f"BAM file already exists for {sam_file}. Skipping conversion.")
    else:
        # Convert SAM to BAM
        os.system(f"samtools view -b {sam_file} > {bam_file}")
        
        # Index BAM file
        os.system(f"samtools index {bam_file}")
        print(f"Converted {sam_file} to BAM and indexed {bam_file}.")

print("SAM to BAM conversion and indexing completed.")

# Step 2: Extract junctions and annotate junctions from BAM files
bam_files = glob.glob(os.path.join(input_path, "*hisat2.dam"))
for bam_file in bam_files:
    base = os.path.basename(bam_file).replace("hisat2.bam", "")
    
    # Extract junctions
    os.system(f"/home/hari/Tools/regtools/build/regtools junctions extract {bam_file} -s XS > {base}_junctions.bed")
    
    # Annotate junctions
    os.system(f"/home/hari/Tools/regtools/build/regtools junctions annotate {base}_junctions.bed /home/weilan/HISAT2/Homo_sapiens.GRCh38.dna.primary_assembly.fasta /home/weilan/HISAT2/Homo_sapiens.GRCh38.111.gtf -o {base}hisat2_junctions_annotation.txt")

# Step 3: Summarize counts from text files
# Discover text files in a directory
dir_path = "/home/weilan/STAR/subset/Regtools/HISAT2_regtools/"  # Replace with the path to your directory containing text files
file_paths = glob.glob(os.path.join(dir_path, "*hisat2_junctions_annotation.txt"))

# Read and summarize counts from text files
def summarize_counts(file_path):
    # Read counts data from text file
    counts_data = pd.read_table(file_path)  # Adjust parameters based on file format
    
    # Calculate total counts for column named "Column4"
    total_counts = counts_data["name"].count()
    
    # Calculate counts for categories in two columns
    category_counts_col1 = counts_data["anchor"].value_counts()
    category_counts_col2 = counts_data["known_junction"].value_counts()
    
    # Return summarized counts data
    return pd.DataFrame({
        "File": [os.path.basename(file_path)],
        "Total_Counts": [total_counts],
        "Category_Counts_Column1": [category_counts_col1.to_dict()],
        "Category_Counts_Column2": [category_counts_col2.to_dict()]
    })

# Apply the summarize_counts function to each file path
summarized_counts = [summarize_counts(file_path) for file_path in file_paths]

# Step 4: Combine summarized counts data into one data frame
combined_counts = pd.concat(summarized_counts)

# Step 5: Write combined counts data to a new text file
output_file_path = os.path.join(dir_path, "hisat2_combined_counts_summary.txt")
combined_counts.to_csv(output_file_path, sep="\t", index=False)
