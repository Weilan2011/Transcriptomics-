import pandas as pd

# Define the path to your input file
file_path = '/home/weilan/ENCODE/2_HISAT2/ENCODE_rMATS_HepG2_pairs_rMATS.txt'  # Update this to the actual path of your file

# Read the data, assuming tab-separated values; adjust if necessary
data = pd.read_csv(file_path, delimiter='\t')

# Initialize a dictionary to hold the regrouped data
regrouped_data = {}

# Loop over each row in the DataFrame
for index, row in data.iterrows():
    # Split target and control files into lists
    target_files = row['Target_file'].split(', ')
    controls = row['Controlled_by'].split(', ')

    # Process each target file
    for target_file in target_files:
        # Extract the group name (e.g., HepG2_AATF from HepG2_AATF_rep1.bam)
        group_name = "_".join(target_file.split('_')[:2])
        
        # Add target file to the group in the dictionary
        if group_name not in regrouped_data:
            regrouped_data[group_name] = {'target_files': [], 'control_files': controls}
        regrouped_data[group_name]['target_files'].append(target_file)
    
# Convert the regrouped data into a format suitable for DataFrame construction
processed_data = []
for group, files in regrouped_data.items():
    processed_data.append({
        'Group': group,
        'Target_Files': ", ".join(files['target_files']),
        'Control_Files': ", ".join(files['control_files'])
    })

# Convert the list of dictionaries into a DataFrame
processed_df = pd.DataFrame(processed_data)

# Save the processed DataFrame to a new file
processed_df.to_csv('/home/weilan/ENCODE/2_HISAT2/Processed_ENCODE_rMATS_HepG2_regrouped.txt', index=False, sep='\t')

print("Data has been processed and saved.")
