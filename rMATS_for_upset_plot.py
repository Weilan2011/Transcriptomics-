import os
import pandas as pd

def get_columns_and_names_for_as_type(filename):
    """
    Determines the columns and names to use based on the alternative splicing type indicated in the filename,
    focusing on geneSymbol, FDR, and IncLevelDifference for significance evaluation.
    Correctly formats the filename to include the AS type without duplication.
    """
    as_types_info = {
        'SE': ([2, 19, 22],
               ["geneSymbol", "FDR", "IncLevelDifference"]),
        'MXE': ([2, 21, 24],
                ["geneSymbol", "FDR", "IncLevelDifference"]),
        'A3SS': ([2, 19, 22],
                 ["geneSymbol", "FDR", "IncLevelDifference"]),
        'A5SS': ([2, 19, 22],
                 ["geneSymbol", "FDR", "IncLevelDifference"]),
        'RI': ([2, 19, 22],
               ["geneSymbol", "FDR", "IncLevelDifference"]),
    }
    for as_type, (columns, names) in as_types_info.items():
        if f"_{as_type}." in filename:
            base_filename = filename.replace(f"_output_{as_type}_{as_type}", f"_{as_type}").replace(".MATS.JCEC.txt", "")
            return columns, names, as_type, base_filename
    return None, None, None, None

def aggregate_values(df):
    # Use the largest absolute 'IncLevelDifference' for each 'geneSymbol' within each 'condition'
    df['abs_IncLevelDifference'] = df['IncLevelDifference'].abs()
    agg_df = df.sort_values('abs_IncLevelDifference', ascending=False).drop_duplicates(['geneSymbol', 'condition'])
    agg_df.drop('abs_IncLevelDifference', axis=1, inplace=True)
    return agg_df

directory = "/home/weilan/ENCODE/4_ALS/ALS_consolidated_outputs"
all_data_original = []
all_data_binary = []

for filename in os.listdir(directory):
    if filename.endswith(".MATS.JCEC.txt"):
        filepath = os.path.join(directory, filename)
        columns, names, as_type, base_filename = get_columns_and_names_for_as_type(filename)
        if columns is None or names is None or as_type is None or base_filename is None:
            continue

        df = pd.read_csv(filepath, sep="\t", usecols=columns, skiprows=1)
        df.columns = names
        
        # Drop rows where 'geneSymbol' is empty and select only significant events
        df = df.dropna(subset=['geneSymbol'])
        df = df[(df['FDR'] < 0.05) & (df['IncLevelDifference'].abs() > 0.2)]

        # Set condition to modified filename
        condition = f"{base_filename}_{as_type}"
        df['condition'] = condition

        # Separate handling for original and binary versions
        df_original = df.copy()
        df_binary = df.copy()
        df_binary['IncLevelDifference'] = 1  # All significant entries are marked as 1
        
        all_data_original.append(df_original)
        all_data_binary.append(df_binary)

# Aggregate values for the original DataFrame to ensure unique geneSymbol-condition pairs with the largest abs IncLevelDifference
combined_original = pd.concat(all_data_original)
combined_original_aggregated = aggregate_values(combined_original)

combined_binary = pd.concat(all_data_binary)
combined_binary_aggregated = aggregate_values(combined_binary)

# Since combined_binary_aggregated now only marks significant changes, IncLevelDifference is already set to 1, we pivot directly
pivot_original = combined_original_aggregated.pivot(index='geneSymbol', columns='condition', values='IncLevelDifference').fillna(0).reset_index()
pivot_binary = combined_binary_aggregated.pivot(index='geneSymbol', columns='condition', values='IncLevelDifference').fillna(0).reset_index()

# Optionally, save the pivoted tables to files
output_filepath_original = os.path.join(directory, "significant_events_pivoted_original.csv")
output_filepath_binary = os.path.join(directory, "significant_events_pivoted_binary.csv")

pivot_original.to_csv(output_filepath_original, index=False)
pivot_binary.to_csv(output_filepath_binary, index=False)

print(f"Pivoted table with largest absolute IncLevelDifference for significant events saved to {output_filepath_original}")
print(f"Pivoted table with binary IncLevelDifference (1 for significant) saved to {output_filepath_binary}")
