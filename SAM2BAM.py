import os
import glob

# Define input directory containing SAM files
input_dir = "/home/weilan/STAR/subset/Regtools/HISAT2_regtools/SAM"
# Define output directory for BAM files
bam_dir = "/home/weilan/STAR/subset/Regtools/HISAT2_regtools/BAM"

# Step 1: Convert SAM to BAM and index BAM files
sam_files = glob.glob(os.path.join(input_dir, "*.sam"))
for sam_file in sam_files:
    base = os.path.basename(sam_file).replace(".sam", "")
    bam_file = os.path.join(bam_dir, base + ".bam")
    
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

