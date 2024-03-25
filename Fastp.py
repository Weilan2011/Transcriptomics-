#sys.argv[1] = input directory where fastq.gz files are in the main directory (fastqDir/sampleId_*.fastq.gz)
inputDir="/home/weilan/ENCODE/0_Raw/"
#sys.argv[2] = output directory where  trimmed fastq.gz files in individual sample folders and MultiQC report will be stored (Fastp)
outputDir="/home/weilan/ENCODE/1_Fastp/"
#sys.argv[3] = number indicating how what number trim this is (1, 2, 3, etc)
trim=1

import glob, os, sys

print("Started trimming with Fastp")

#grabbing paths of fastq.gz files
files=glob.glob(inputDir+"*.fastq.gz")
##files=glob.glob(inputDir+"*/*.fastq.gz")

#empty list for processed samples
processed=[]

#making a directory for Trim number (change for each trim)
os.system("mkdir "+outputDir+"Trim_"+str(trim))

#loop to trim fastq.gz files and generate Fastp reports
for file in files:
    #sample name
    #sample=file.split("/")[-1].split("_")[0]
    sample="_".join(file.split("/")[-1].split("_")[:3])
    #skip sample if already processed
    if sample in processed:
        pass
    #trimming the sample
    else:
        processed.append(sample)
        print("Trimming:" +sample+"\n")
        #path to sample directory
        samplePath=outputDir+"Trim_"+str(trim)+"/"+sample+"/"
        #creating directory for trimmed output
        os.system("mkdir "+samplePath)
        #grabbing paths of R1/R2 for all lanes
        ##r1=glob.glob(inputDir+sample+"_*R1*.fastq.gz")
        ##r2=glob.glob(inputDir+sample+"_*R2*.fastq.gz")
        r1=glob.glob(inputDir+sample+"_*R1*.fastq.gz")
        r2=glob.glob(inputDir+sample+"_*R2*.fastq.gz")
        r1.sort();r2.sort();

        #multiple lanes
        if len(r1) >1:
            #loop to process each lane for each sample
            for i in range(0, len(r1)):
                r1Lane=r1[i]
                r2Lane=r2[i]
                #finding the lane number 
                idx1=r1Lane.find("00")
                laneR1=r1Lane[idx1:idx1+3]
                idx2=r2Lane.find("00")
                laneR2=r2Lane[idx2:idx2+3]
                #general fastp command 
                    #change as necessary for target trimming
                cmd = "fastp -i "+r1Lane+" -I "+r2Lane+" --detect_adapter_for_pe -o "+samplePath+sample+"_R1_"+laneR1+".fastq.gz -O "+samplePath+sample+"_R2_"+laneR2+".fastq.gz -h "+samplePath+sample+"_"+laneR1+".fastp.html -j "+samplePath+sample+"_"+laneR1+".fastp.json"
                print("\nTrimming for multiple lanes: "+sample+" lane "+laneR1+"\n")
                os.system(cmd)
            
        #single lanes    
        else:
            r1=r1[0]
            r2=r2[0]
            #Trim front (-f # -F #)
            cmd = "fastp -i "+r1+" -I "+r2+" -o "+samplePath+"/"+sample+"_R1_001.fastq.gz -O "+samplePath+"/"+sample+"_R2_001.fastq.gz -h "+samplePath+"/"+sample+".fastp.html -j "+samplePath+"/"+sample+".fastp.json -w 16"
            print("Trimming for single lane:"+sample)
            os.system(cmd)

#generating MultiQC reports
cmd = "multiqc "+outputDir+" -o "+outputDir+"/MultiQC/"
os.system(cmd)
