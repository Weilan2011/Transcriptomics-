import os
import glob
import subprocess

# Define input directory containing BAM files
bam_dir = "/home/weilan/STAR/subset/Regtools/HISAT2_regtools/BAM/"

# Define annotation file path
annotation_file = "/home/weilan/HISAT2/Homo_sapiens.GRCh38.111.gtf"  # Replace with the path to your annotation file

# Define output directory for gene counts
output_dir = "/home/weilan/HISAT2/Subset_out/"

# Step 1: Count gene reads using FeatureCounts
bam_files = glob.glob(os.path.join(bam_dir, "*sorted.bam"))
for bam_file in bam_files:
    
    subprocess.run(["samtools", "index", bam_file])
    
    # Construct output file path
    output_file = os.path.join(output_dir, os.path.basename(bam_file).replace("sorted.bam", "_gene_counts.txt"))
    
    # Execute FeatureCounts command
    command = f"featureCounts -a {annotation_file} -o {output_file} {bam_file}"
    subprocess.run(command, shell=True)

print("Gene counting completed.")
