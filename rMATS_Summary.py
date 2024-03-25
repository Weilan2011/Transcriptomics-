import os
import shutil
import pandas as pd

# Set the base directory where your folders containing rMATS outputs are located
base_dir = "/home/weilan/ENCODE/4_ALS/rMATS_output/"

# Define the path for the common folder to store all renamed and consolidated files
common_folder = "/home/weilan/ENCODE/4_ALS/ALS_consolidated_outputs"

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
    files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and f.endswith("JCEC.txt")]

    print(f"Processing folder: {folder_name}, found files: {files}")

    for file in files:
        new_file_name = f"{folder_name}_{file}"
        new_file_path = os.path.join(common_folder, new_file_name)
        original_file_path = os.path.join(folder, file)
        
        shutil.copy(original_file_path, new_file_path)  # Changed from shutil.move to shutil.copy
        print(f"Copied: {original_file_path} to {new_file_path}")

print("Files have been renamed and copied to the common folder successfully.")

# Initialize a list to store summary data
summary_data = []

# Loop through each file in the common folder
for file in os.listdir(common_folder):
    if file.endswith("JCEC.txt"):
        # Construct the full file path
        file_path = os.path.join(common_folder, file)
        
        # Read the JCEC file into a DataFrame
        try:
            df = pd.read_csv(file_path, sep='\t')
        except Exception as e:
            print(f"Error reading {file}: {e}")
            continue
        
        # Calculate the total number of events
        total_events = len(df)
        
        # Calculate the number of significant events based on criteria
        significant_events = df[(df['FDR'] < 0.05) & (df['IncLevelDifference'].abs() > 0.2)].shape[0]
        
        # Calculate the percentage of significant events
        percent_significant = (significant_events / total_events) * 100 if total_events > 0 else 0

        # Calculate sum for IJC and SJC for both groups across all samples and events    
        df['IJC_SAMPLE_1_sum'] = df['IJC_SAMPLE_1'].apply(lambda x: sum(map(int, x.split(','))) if isinstance(x, str) else x)
        df['SJC_SAMPLE_1_sum'] = df['SJC_SAMPLE_1'].apply(lambda x: sum(map(int, x.split(','))) if isinstance(x, str) else x)
        df['IJC_SAMPLE_2_sum'] = df['IJC_SAMPLE_2'].apply(lambda x: sum(map(int, x.split(','))) if isinstance(x, str) else x)
        df['SJC_SAMPLE_2_sum'] = df['SJC_SAMPLE_2'].apply(lambda x: sum(map(int, x.split(','))) if isinstance(x, str) else x)


        # Calculate the total sum across all events for each group
        group_1_i_sum = df['IJC_SAMPLE_1_sum'].sum()
        group_1_s_sum = df['SJC_SAMPLE_1_sum'].sum()
        group_2_i_sum = df['IJC_SAMPLE_2_sum'].sum()
        group_2_s_sum = df['SJC_SAMPLE_2_sum'].sum()
        
        # Calculate total IJC and SJC for both groups
        total_i_sum = group_1_i_sum + group_2_i_sum
        total_s_sum = group_1_s_sum + group_2_s_sum

        # Calculate the percentages for each group
        group_1_i_percent = round((group_1_i_sum / (group_1_i_sum + group_1_s_sum)) * 100, 2) if (group_1_i_sum + group_1_s_sum) > 0 else 0
        group_1_s_percent = round((group_1_s_sum / (group_1_i_sum + group_1_s_sum)) * 100, 2) if (group_1_i_sum + group_1_s_sum) > 0 else 0
        group_2_i_percent = round((group_2_i_sum / (group_2_i_sum + group_2_s_sum)) * 100, 2) if (group_2_i_sum + group_2_s_sum) > 0 else 0
        group_2_s_percent = round((group_2_s_sum / (group_2_i_sum + group_2_s_sum)) * 100, 2) if (group_2_i_sum + group_2_s_sum) > 0 else 0
    
        # Add the summary information to the list
        summary_data.append({
            'File': file, 
            'Total Events': total_events, 
            'Significant Events': significant_events,
            'Percent Significant': round(percent_significant, 2),
            'Group 1 IJC Total': group_1_i_sum,
            'Group 1 SJC Total': group_1_s_sum,
            'Group 2 IJC Total': group_2_i_sum,
            'Group 2 SJC Total': group_2_s_sum,
            'Group 1 IJC Percent': group_1_i_percent,
            'Group 1 SJC Percent': group_1_s_percent,
            'Group 2 IJC Percent': group_2_i_percent,
            'Group 2 SJC Percent': group_2_s_percent
        })
        
# Convert the summary data to a DataFrame for easier viewing and manipulation
summary_df = pd.DataFrame(summary_data)

# Optional: Save the summary to a CSV file for further analysis or reporting
summary_csv_path = os.path.join(common_folder, "summary_of_JCEC_files.csv")
summary_df.to_csv(summary_csv_path, index=False)
print(f"Summary saved to {summary_csv_path}")
