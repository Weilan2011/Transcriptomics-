import glob
import os
import subprocess
from argparse import ArgumentParser

def parse_arguments():
    parser = ArgumentParser(description="Automate FASTQ trimming and MultiQC report generation for paired-end reads")
    parser.add_argument("--inputDir", required=True, help="Input directory where FASTQ.gz files are located")
    parser.add_argument("--outputDir", required=True, help="Output directory for trimmed FASTQ.gz files and MultiQC report")
    parser.add_argument("--trim", type=int, default=1, help="Number indicating the trim iteration (1, 2, 3, etc.)")
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    print("Started trimming with Fastp")

    # Ensure the output directory exists
    os.makedirs(args.outputDir, exist_ok=True)

    # Process only R1 files and automatically infer R2 files
    r1_files = glob.glob(os.path.join(args.inputDir, "*_R1.fastq.gz"))

    for r1_path in r1_files:
        # Infer R2 path from R1 path
        r2_path = r1_path.replace("_R1.fastq.gz", "_R2.fastq.gz")
        # Ensure the R2 file exists
        if not os.path.exists(r2_path):
            print(f"Missing R2 for {r1_path}, skipping...")
            continue
        
        # Base name without R1/R2 for output naming
        base_name = os.path.basename(r1_path).replace("_R1.fastq.gz", "")
        
        # Output paths
        output_r1 = os.path.join(args.outputDir, f"{base_name}_R1_trimmed.fastq.gz")
        output_r2 = os.path.join(args.outputDir, f"{base_name}_R2_trimmed.fastq.gz")
        html_report = os.path.join(args.outputDir, f"{base_name}_fastp.html")
        json_report = os.path.join(args.outputDir, f"{base_name}_fastp.json")

        # Construct and run fastp command
        cmd = [
            "fastp",
            "-i", r1_path,
            "-I", r2_path,
            "-o", output_r1,
            "-O", output_r2,
            "-h", html_report,
            "-j", json_report,
            "-w", "16"  # Adjust the number of threads as needed
        ]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
    
    # Generate MultiQC report
    multiqc_cmd = ["multiqc", args.outputDir, "-o", args.outputDir]
    subprocess.run(multiqc_cmd, check=True)

if __name__ == "__main__":
    main()
