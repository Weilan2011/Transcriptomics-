# ENCODE protein binding protein raw fastq filename re-formating by Weilan  
# Prepare your list: Create a text file (e.g., rename_list.txt) that contains two columns: the current file names and the new file names, separated by a delimiter (e.g., a space or a tab). 

import os

# Read the rename list from a file
rename_list_file = '/home/weilan/ENCODE/1_Fastp/K562_Trimed/Set7/fastp_file_rename2.txt'

# Specify the folder containing the files you want to rename
folder_path = '/home/weilan/ENCODE/1_Fastp/K562_Trimed/Set7'

# Read the rename list and store it in a dictionary
rename_dict = {}
with open(rename_list_file, 'r') as f:
    for line in f:
        current_name, new_name = line.strip().split()
        rename_dict[current_name] = new_name

# Iterate through the files in the folder
for filename in os.listdir(folder_path):
    if filename in rename_dict:
        current_path = os.path.join(folder_path, filename)
        new_filename = rename_dict[filename]
        new_path = os.path.join(folder_path, new_filename)

        # Rename the file
        os.rename(current_path, new_path)
        print(f'Renamed: {current_path} -> {new_path}')
        
# Make sure to replace 'rename_list.txt' with the actual path to your rename list file and 'path/to/your/folder' with the folder path containing the files you want to rename.

# This Python script reads the rename list, creates a dictionary mapping old names to new names, and then iterates through the files in the folder. If a file matches an entry in the rename list, it renames the file accordingly.


