for infile in *_control_bams.txt; do
    base=$(basename "${infile}" "_control_bams.txt")
    target_file="${base}_target_bams.txt"
    
    # Read the contents of the control and target BAM lists
    control_bams=$(cat "${infile}")
    target_bams=$(cat "${target_file}")
    
    # Assuming the Python script can accept BAM paths from stdin or as a direct argument
    # This example shows passing them as arguments
    # Update the CrypSplice2.0_modular.py command as necessary based on its input requirements
    python3 /home/weilan/CrypSplice2.0/CrypSplice2.0/CrypSplice2.0_modular.py CrypticJunctions \
        -c1 "${control_bams}" \
        -c2 "${target_bams}" \
        -gtf /home/weilan/HISAT2/HISAT2_Indexes/Homo_sapiens.GRCh38.111.gtf \
        -fasta /home/weilan/HISAT2/HISAT2_Indexes/Homo_sapiens.GRCh38.dna.primary_assembly.fasta \
        -s 2 \
        -o "${base}_CrySplice_out"
done
