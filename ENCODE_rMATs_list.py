def create_id_to_bam_mapping(file_path):
    mapping = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            control_id, bam_file = parts[0], parts[1]
            mapping[control_id] = bam_file
    return mapping

def update_control_ids_with_bam_names(mapping_file, target_file, updated_file):
    mapping = create_id_to_bam_mapping(mapping_file)
    
    with open(target_file, 'r') as infile, open(updated_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')
            # Assuming the control IDs you need to replace are in the second column
            before_ids, control_ids = parts[0], parts[1]
            updated_ids = [mapping.get(control_id.strip(), 'ID_NOT_FOUND') for control_id in control_ids.split(',')]
            outfile.write(f"{before_ids}\t{','.join(updated_ids)}\n")

# Usage
mapping_file = '/home/weilan/ENCODE/2_HISAT2/ENCODE_rMATS_HepG2_Control_pairs.txt'  # This should be replaced with the path to your first file
target_file = '/home/weilan/ENCODE/2_HISAT2/ENCODE_rMATS_HepG2_pairs.txt'  # This should be replaced with the path to your second file
updated_file = '/home/weilan/ENCODE/2_HISAT2/ENCODE_rMATS_HepG2_pairs_eidted.txt'  # Output file

update_control_ids_with_bam_names(mapping_file, target_file, updated_file)





