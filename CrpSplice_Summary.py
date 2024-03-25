import os
import shutil
import glob
import pandas as pd

# Set the base directory where your folders containing rMATS outputs are located
base_dir = "/home/weilan/ENCODE/5_CrypSplice/ALS_output"

# Define the path for the common folder to store all renamed and consolidated files
common_folder = "/home/weilan/ENCODE/5_CrypSplice/"

# Check if the common folder exists, if not, create it
if not os.path.exists(common_folder):
    os.makedirs(common_folder)
    print(f"Created the directory: {common_folder}")

# Exclude 'consolidated_outputs' from the folders to process
folders = [os.path.join(base_dir, f) for f in os.listdir(base_dir)
           if os.path.isdir(os.path.join(base_dir, f)) and f != "consolidated_outputs"]

print(f"Found folders: {folders}")

# Loop through each subdirectory to process the files within
for folder in folders:
    folder_name = os.path.basename(folder)
    # Ensure we only list files, avoiding potential directory matches
    files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and f.endswith("Novel_Junctions.txt")]

    print(f"Processing folder: {folder_name}, found files: {files}")

    for file in files:
        new_file_name = f"{folder_name}_{file}"
        new_file_path = os.path.join(common_folder, new_file_name)
        original_file_path = os.path.join(folder, file)
        
        shutil.copy(original_file_path, new_file_path)  # Changed from shutil.move to shutil.copy
        print(f"Copied: {original_file_path} to {new_file_path}")

print("Files have been renamed and copied to the common folder successfully.")

### update the Gene Smybols ####
# Load the simplified GTF reference into a dictionary
gtf_df = pd.read_csv('/home/weilan/simplified_GRCh38.111.gtf.tsv', header=None, sep=' ', names=['gene_id', 'gene_name'])
gtf_dict = dict(zip(gtf_df['gene_id'].str.replace('"','').str.replace(';',''), gtf_df['gene_name'].str.replace('"','').str.replace(';','')))

# Function to add gene names to a file
def add_gene_names(file_path, gtf_dict):
    # Read the file
    df = pd.read_csv(file_path, sep='\t')
    
    # Remove potential quotes and semicolons from gene IDs in your file if present
    df['gene'] = df['gene'].str.replace('"','').str.replace(';','')
    
    # Map gene_id to gene_name using the gtf_dict
    df['gene_name'] = df['gene'].map(gtf_dict)
    
    # Save the updated DataFrame back to a file
    df.to_csv(file_path.replace('.txt', '_with_gene_names.txt'), index=False, sep='\t')

# Process each text file in the folder
for file_path in glob.glob('/home/weilan/ENCODE/5_CrypSplice/ALS_consolidated/*.txt'):
    add_gene_names(file_path, gtf_dict)

# Summarize the updated results
# Define the directory where your files are located
directory_path = "/home/weilan/ENCODE/5_CrypSplice/ALS_consolidated/Gene_symbols_udated"

# Initialize an empty dataframe for aggregating significant junction IDs
significant_juncIDs_df = pd.DataFrame()

# List all files in the directory
file_paths = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith('.txt')]

# Process each file
for file_path in file_paths:
    df = pd.read_csv(file_path, sep='\t')
    # Filter for significant junction IDs (adj.pVal < 0.05)
    significant_df = df[(df['adj.pVal'] < 0.05)]
  
    # Extract cell line and gene target from the file name
    parts = file_path.split('/')[-1].split('_')
    cell_line = parts[0]
    gene_target = parts[1]
    
    # Rename JS_diff column to include cell line and gene target
    significant_df = significant_df.rename(columns={'JS_diff': f'JS_diff_{cell_line}_{gene_target}',
                                                    'gene_name': 'gene_name',
                                                    'annotation': 'annotation'})
    
    # Select relevant columns, now including "annotation"
    significant_df = significant_df[['juncID', f'JS_diff_{cell_line}_{gene_target}', 'gene_name', 'annotation']]
    
    # Merge with the aggregated dataframe
    if significant_juncIDs_df.empty:
        significant_juncIDs_df = significant_df
    else:
        # Ensuring "annotation" is carried over and merged properly
        significant_juncIDs_df = pd.merge(significant_juncIDs_df, significant_df, on=['juncID', 'gene_name', 'annotation'], how='outer')

# After processing all files, save the aggregated dataframe to a CSV
output_path = os.path.join(directory_path, 'aggregated_significant_juncIDs_with_annotation.csv')
significant_juncIDs_df.to_csv(output_path, index=False)

print(f"Aggregated data saved to '{output_path}'")



    




