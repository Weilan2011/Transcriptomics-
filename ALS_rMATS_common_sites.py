import os
import pandas as pd

def get_columns_and_names_for_as_type(filename):
    """
    Determines the columns and names to use based on the alternative splicing type indicated in the filename.
    """
    as_types_info = {
        'SE': ([2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 22],
               ["geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd",
                "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "FDR", "IncLevelDifference"]),
        'MXE': ([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 21, 24],
                ["geneSymbol", "chr", "strand", "1stExonStart_0base", "1stExonEnd",
                 "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE",
                 "downstreamES", "downstreamEE", "FDR", "IncLevelDifference"]),
        'A3SS': ([2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 22],
                 ["geneSymbol", "chr", "strand", "longExonStart_0base", "longExonEnd",
                  "shortES", "shortEE", "flankingES", "flankingEE", "FDR", "IncLevelDifference"]),
        'A5SS': ([2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 22],
                 ["geneSymbol", "chr", "strand", "longExonStart_0base", "longExonEnd",
                  "shortES", "shortEE", "flankingES", "flankingEE", "FDR", "IncLevelDifference"]),
        'RI': ([2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 22],
               ["geneSymbol", "chr", "strand", "riExonStart_0base", "riExonEnd", "upstreamES",
                "upstreamEE", "downstreamES", "downstreamEE", "FDR", "IncLevelDifference"]),
    }
    for as_type, (columns, names) in as_types_info.items():
        if f"_{as_type}." in filename:
            return columns, names, as_type
    return None, None, None

directory = "/home/weilan/ENCODE/4_ALS/ALS_consolidated_outputs"
combined_data = {}

for filename in os.listdir(directory):
    if filename.endswith(".MATS.JCEC.txt"):
        filepath = os.path.join(directory, filename)
        columns, names, as_type = get_columns_and_names_for_as_type(filename)
        if columns is None or names is None or as_type is None:
            continue

        df = pd.read_csv(filepath, sep="\t", usecols=columns, skiprows=1)
        df.columns = names

        # Filter based on FDR and IncLevelDifference
        df_filtered = df[(df['FDR'] < 0.05) & (df['IncLevelDifference'].abs() > 0.2)].copy()

        # Rename IncLevelDifference columns to include the filename or a unique identifier
        base_filename = filename.replace(".MATS.JCEC.txt", "")
        df_filtered.rename(columns={
            "IncLevelDifference": f"{base_filename}_IncLevelDifference"
        }, inplace=True)

        # Merge with combined data
        if as_type not in combined_data:
            combined_data[as_type] = df_filtered
        else:
            combined_data[as_type] = pd.merge(
                combined_data[as_type], df_filtered, 
                on=names[:-1],  # Assuming the last name is IncLevelDifference
                how='outer'
            )

# After merging all files, fill missing values with 0, then save the consolidated data for each AS type
for as_type, df in combined_data.items():
    # Fill missing values with 0
    df.fillna(0, inplace=True)
    
    output_filename = f"Common_sites_{as_type}.csv"
    output_filepath = os.path.join(directory, output_filename)
    df.to_csv(output_filepath, index=False)
    print(f"Filled and combined {as_type} data saved to {output_filepath}")


