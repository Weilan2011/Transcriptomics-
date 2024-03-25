# Adjust the script to use the correct column name 'symbol' for gene identifiers

significant_genes_df = pd.DataFrame()

for condition, file_path in file_paths.items():
    # Load dataset
    df = pd.read_csv(file_path, sep="\t")
    
    # Filter for significant genes based on criteria (padj < 0.05, |log2FoldChange| > 0.6)
    significant_df = df[(df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0.6)]
    
    # If significant genes exist, proceed
    if not significant_df.empty:
        # Extract gene symbols and log2FoldChange for significant genes
        temp_df = significant_df[['symbol', 'log2FoldChange']].copy()
        temp_df.rename(columns={'log2FoldChange': condition}, inplace=True)
        temp_df.set_index('symbol', inplace=True)
        
        # Merge with the accumulating DataFrame
        if significant_genes_df.empty:
            significant
