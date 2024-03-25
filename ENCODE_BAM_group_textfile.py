import sys
import os

# Check for correct usage
if len(sys.argv) != 2:
    print("Usage: python script.py <input_file_path>")
    sys.exit(1)

input_file_path = sys.argv[1]
output_folder = "/home/weilan/ENCODE/3_rMATS"  # Set the output folder as specified

def process_file(input_file_path, output_folder):
    with open(input_file_path, 'r') as file:
        # Skip the header line
        next(file)
        for line in file:
            group, target_files, control_files = line.strip().split('\t')
            
            # Remove spaces after commas
            target_files_cleaned = ','.join([filename.strip() for filename in target_files.split(',')])
            control_files_cleaned = ','.join([filename.strip() for filename in control_files.split(',')])
            
            # Prepare filenames with the output folder
            target_output_path = os.path.join(output_folder, f'{group}_target_bams.txt')
            control_output_path = os.path.join(output_folder, f'{group}_control_bams.txt')
            
            # Ensure the output folder exists
            os.makedirs(output_folder, exist_ok=True)
            
            # Write target files
            with open(target_output_path, 'w') as target_file:
                target_file.write(target_files_cleaned)
            
            # Write control files
            with open(control_output_path, 'w') as control_file:
                control_file.write(control_files_cleaned)
            
            print(f"Generated files: {target_output_path}, {control_output_path}")

if __name__ == "__main__":
    process_file(input_file_path, output_folder)




