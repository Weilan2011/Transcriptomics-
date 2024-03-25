#HISTA2 Alignment Script for NGS Pipeline (Weilan)

import os, glob, shutil

# Specify input and output directories
inputDir = "/home/weilan/STAR/subset/fastp/Trimed"
outputDir = "/home/weilan/STAR/subset/HISAT2/"

# Specify HISAT2 index
hisat2_index = "/path/to/hisat2_index/index_prefix"

# Specify the number of CPU cores
cores = 10

# Create a list of FASTQ files
files = glob.glob(inputDir + "/*/*.fastq.gz")

# Processed samples
processed = []

# Run HISAT2 alignments
for file in files:
    # Extract sample name
    sample = "_".join(file.split("/")[-1].split("_")[:3])
    
    if sample in processed:
        continue
    
    processed.append(sample)
    print("\nProcessing: " + sample + "\n")

    # Create output directory for the sample
    sample_output_dir = os.path.join(outputDir, sample)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Find R1 and R2 FASTQ files
    r1 = glob.glob(inputDir + "*/" + sample + "*R1_001*.fastq.gz")
    r2 = glob.glob(inputDir + "*/" + sample + "*R2_001*.fastq.gz")
    
    r1.sort()
    r2.sort()

    # If multiple lanes for R1/R2, join them with commas
    if len(r1) > 1:
        r1 = ",".join(r1)
        r2 = ",".join(r2)
    else:
        r1 = r1[0]
        r2 = r2[0]

    # HISAT2 alignment command
    cmd = f'hisat2 -p {cores} -x {hisat2_index} -1 {r1} -2 {r2} -S {sample_output_dir}/{sample}.sam'
    print(cmd + "\n")
    os.system(cmd)

    # Convert SAM to BAM
    cmd = f'samtools view -bS {sample_output_dir}/{sample}.sam > {sample_output_dir}/{sample}.bam'
    os.system(cmd)

    # Sort and index BAM file
    cmd = f'samtools sort -@ {cores} -o {sample_output_dir}/{sample}.sorted.bam {sample_output_dir}/{sample}.bam'
    os.system(cmd)
    
    cmd = f'samtools index {sample_output_dir}/{sample}.sorted.bam'
    os.system(cmd)

    print("Finished processing: " + sample)

print("All samples processed.")