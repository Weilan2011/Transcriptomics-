import pandas as pd
import os
import glob

# Step 1: Discover text files in a directory
dir_path = "/home/weilan/STAR/subset/Regtools/HISAT2_regtools/BAM/"  # Replace with the path to your directory containing text files
file_paths = glob.glob(os.path.join(dir_path, "*.txt"))

# Step 2: Read and summarize counts from text files
def summarize_counts(file_path):
    # Read counts data from text file
    counts_data = pd.read_table(file_path)  # Adjust parameters based on file format
    
    # Calculate total counts for column named "Column4"
    total_counts = counts_data["name"].count()
    
    # Calculate counts for categories in two columns
    category_counts_col1 = counts_data["anchor"].value_counts()
    category_counts_col2 = counts_data["known_junction"].value_counts()
    
    # Return summarized counts data
    return pd.DataFrame({
        "File": [os.path.basename(file_path)],
        "Total_Counts": [total_counts],
        "Category_Counts_Column1": [category_counts_col1.to_dict()],
        "Category_Counts_Column2": [category_counts_col2.to_dict()]
    })

# Apply the summarize_counts function to each file path
summarized_counts = [summarize_counts(file_path) for file_path in file_paths]

# Step 3: Combine summarized counts data into one data frame
combined_counts = pd.concat(summarized_counts)

# Step 4: Write combined counts data to a new text file
output_file_path = os.path.join(dir_path, "STAR_combined_counts_summary.txt")
combined_counts.to_csv(output_file_path, sep="\t", index=False)

