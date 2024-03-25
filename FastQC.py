#FastQC Script for NGS Pipeline (Michelle)
# Weilan testing with RNA-seq data @ BCM Feb 02, 2024 

#sys.argv[1] = input directory where fastq.gz files are in the main directory (fastq_data/)
inputDir="/home/weilan/MECP2/0_Raw/"
#sys.argv[2] = output directory where FastQC reports for each sample and general MultiQC will be stored (FastQC/)
outputDir="/home/weilan/MECP2/0_Raw/FastQC/"

import os,sys, glob
  
print("Started FastQC analysis")

#grabbing paths of fastq.gz files
files=glob.glob(inputDir+"*.fastq.gz")
print(files)

#empty list for processed samples
processed=[]

#creating directory for FastQC and MultiQC reports
os.system("mkdir "+inputDir+"FastQC")
os.system("mkdir "+inputDir+"MultiQC")

#loop to merge lanes and generate FastQC reports
for file in files:
    #sample name
    #sample=file.split("/")[-1].split("_")[0]
    sample="_".join(file.split("/")[-1].split("_")[:2])
    #skip sample if already processed
    if sample in processed:
        pass
    #processing sample
    else:
        processed.append(sample)
        print("\nProcessing:"+sample+"\n")
        #creating directory for FastQC reports for each sample
        os.system("mkdir "+inputDir+"FastQC/"+sample)
        #grabbing paths of R1/R2 files and sorting by ascending lanes
        r1=glob.glob(inputDir+sample+"_*R1_001.fastq.gz")
        r2=glob.glob(inputDir+sample+"_*R2_001.fastq.gz")
        ##r1=glob.glob(inputDir+"*/"+sample+"_*R1*.fastq.gz")
        ##r2=glob.glob(inputDir+"*/"+sample+"_*R2*.fastq.gz")
        r1.sort();r2.sort();
        #multiple lanes for R1/R2
        if len(r1) >1:
            #merging lanes for R1
            cmd = "cat "+" ".join(r1)+" > "+inputDir+"FastQC/"+sample+"/"+sample+"_R1_001.fastq.gz"
            print(cmd+"\n")
            os.system(cmd)
            #merging lanes for R2
            cmd = "cat "+" ".join(r2)+" > "+inputDir+"FastQC/"+sample+"/"+sample+"_R2_001.fastq.gz"
            print(cmd+"\n")
            os.system(cmd)
           #creating FastQC reports
            cmd="fastqc"+ " -t 20 " + inputDir+"FastQC/"+sample+"/"+sample+"_R1_001.fastq.gz "+inputDir+"FastQC/"+sample+"/"+sample+"_R2_001.fastq.gz -o "+inputDir+"FastQC/"+sample+"/"
            print(cmd+"\n")
            os.system(cmd)
        #single lane for R1/R2    
        else:
            #grabbing paths for R1/R2 files
            r1=r1[0]
            r2=r2[0]
            #creating FastQC reports
            cmd="fastqc"+ " -t 20 " +r1+" "+r2+" -o "+inputDir+"FastQC/"+sample
            print(cmd+"\n")
            os.system(cmd)

#generating MultiQC reports
cmd = "multiqc "+inputDir+" -o "+inputDir+"MultiQC/"
os.system(cmd)
