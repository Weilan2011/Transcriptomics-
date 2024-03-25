import os
import pandas as pd

# Define the directory containing the text files
directory = "/home/weilan/RBM17/RBM_8162/4_regtools/"

# Initialize an empty DataFrame to store the combined data
combined_data = pd.DataFrame()

# Iterate over each text file in the directory
for filename in os.listdir(directory):
    if filename.endswith(".txt"):
        filepath = os.path.join(directory, filename)
        
        # Read the current text file into a DataFrame
        df = pd.read_csv(filepath, sep="\t", header=None, usecols=[1, 2, 4, 10, 13], names=["Start", "End", "Counts", "Junction_Type", "Known_junction"])
        
        # Add a column indicating the method or file name
        method_name = os.path.splitext(filename)[0]  # Extract method from file name
        df[method_name] = df["Counts"]  # Use the counts as values for the method column
        
        # Merge the current DataFrame with the combined data
        if combined_data.empty:
            combined_data = df.drop(columns=["Counts"])
        else:
            combined_data = pd.merge(combined_data, df.drop(columns=["Counts"]), how="outer", on=["Start", "End", "Junction_Type", "Known_junction"])

# Fill NaN values with 0
combined_data.fillna(0, inplace=True)

# Write the combined data to a new file
combined_file_path = "/home/weilan/RBM17/RBM_8162/4_regtools/TEST_combined_data.csv"
combined_data.to_csv(combined_file_path, index=False)
