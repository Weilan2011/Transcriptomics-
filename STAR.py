#sys.argv[1] = input directory where trimmed fastq.gz files are in the Fastp directory (Fastp/Trim_#)
inputDir="/home/weilan/RBM17/RBM_8162/2_Fastp/Trim_1"
#sys.argv[2] = output directory where STAR alignment output will be stored in individual sample folders (STAR)
outputDir="/home/weilan/RBM17/RBM_8162/3_STAR/"
#sys.argv[3] = STAR genome directory for read length 
STARdir="/home/weilan/STAR/GRCh38_STAR_indexes"
#sys.argv[4] = number of cores to use (usually 10)
cores=20
#sys.argv[5] = overhang = readlength-1
overhang=100
##REMEMBER to change ohang (readLength-1) and ank to appropriate values 

import os, sys, glob, shutil

files=glob.glob(inputDir+"*/*.fastq.gz")
print(files)

processed=[]

#Generating STAR alignments
for file in files:

    #sample=file.split("/")[-1].split("_")[0]
    sample="_".join(file.split("/")[-1].split("_")[:1])

    
    if sample in processed:
        pass

    else:
        processed.append(sample)
        print("\nProcessing:"+sample+"\n")

        os.system("mkdir "+outputDir+sample)
    
        r1=glob.glob(inputDir+"*/"+sample+"*R1_001*.fastq.gz")
        r2=glob.glob(inputDir+"*/"+sample+"*R2_001*.fastq.gz")
        r1.sort();r2.sort();
        
        #multiple lanes for R1/R2
        if len(r1) >1:
            r1 = ",".join(r1)
            r2 = ",".join(r2)

            
        #single lane for R1/R2    
        else:
            r1=r1[0]
            r2=r2[0]

        #Alignment using STAR 
        genomeDir = STARdir
        ohang = overhang;
        # --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2
        cmd = 'STAR --runThreadN '+str(cores)+' --genomeDir '+ genomeDir + ' --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --sjdbOverhang '+str(ohang)+' --outFileNamePrefix '+outputDir+sample+"/"+sample+'_ --readFilesCommand gunzip -c --readFilesIn '+r1+" "+r2+' --outSAMstrandField intronMotif'
        print(cmd+"\n")
        os.system(cmd)

