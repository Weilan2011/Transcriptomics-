
#STARdir=/mnt/local_storage/michelle/data/Experiments/HudaZoghbi/Antonia_CB/STAR/
#STARdir=/mnt/local_storage/michelle/data/Experiments/HudaZoghbi/Carolyn_BS/STAR/
#STARdir=/mnt/local_storage/michelle/data/Experiments/HudaZoghbi/BS-CB_SCA_2/CB/STAR/
#STARdir=/mnt/local_storage/michelle/data/Training/SH-SY5Y_Test/STAR_outSAMstrand_intronMotif/
#STARdir=/mnt/local_storage/michelle/data/Training/SH-SY5Y_Test/STAR_twoPassMode_intronMotif/
STARdir=/mnt/local_storage/michelle/data/Training/PreProcessing_Overview/STAR/


for bamfile in `ls $STARdir*/*Aligned.sortedByCoord.out.bam`; do
    echo "indexing $bamfile"
    samtools index $bamfile
    echo "renaming $bamfile"
    newName=$(echo $bamfile | sed "s#_Aligned.sortedByCoord.out.bam#".sorted.bam"#")
    mv $bamfile $newName
done


for bamIndex in `ls $STARdir*/*.bam.bai`; do
    echo "renaming $bamIndex"
    newName=$(echo $bamIndex | sed "s#_Aligned.sortedByCoord.out.bam.bai#".sorted.bam.bai"#")
    mv $bamIndex $newName
done


for counts in `ls $STARdir*/*ReadsPerGene.out.tab`; do
    echo "renaming $counts"
    newName=$(echo $counts | sed "s#_ReadsPerGene.out.tab#".counts.txt"#")
    mv $counts $newName
done

 
