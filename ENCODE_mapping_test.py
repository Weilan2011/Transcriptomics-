def create_mapping_from_file(mapping_file_path):
    """Create a mapping from control IDs to BAM file names from a given file."""
    mapping = {}
    with open(mapping_file_path, 'r') as file:
        next(file)  # Skip the header
        for line in file:
            parts = line.strip().split('\t')
            control_id, bam_file = parts[0].strip(), parts[1].strip()
            mapping[control_id] = bam_file
    return mapping

def replace_control_ids_in_file(mapping, file_path, output_path):
    """Replace control IDs with BAM file names in a target file and save to a new file."""
    with open(file_path, 'r') as file:
        content = file.read()

    # Remove quotes and replace control IDs with BAM file names
    content = content.replace('"', '')
    for control_id, bam_file in mapping.items():
        content = content.replace(control_id, bam_file)

    # Save the updated content to a new file
    with open(output_path, 'w') as new_file:
        new_file.write(content)

if __name__ == "__main__":
    # Paths to your files
    mapping_file_path = '/home/weilan/ENCODE/2_HISAT2/ENCODE_rMATS_HepG2_Control_pairs.txt'
    file_path_needs_replacement = '/home/weilan/ENCODE/2_HISAT2/ENCODE_rMATS_HepG2_pairs.txt'
    output_path = '/home/weilan/ENCODE/2_HISAT2/ENCODE_rMATS_HepG2_pairs_updated.txt'
    
    # Create the mapping and replace control IDs in the target file
    mapping = create_mapping_from_file(mapping_file_path)
    replace_control_ids_in_file(mapping, file_path_needs_replacement, output_path)

    print(f"Updated file has been saved to: {output_path}")
