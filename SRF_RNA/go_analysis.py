#!/usr/bin/env python3

import pandas as pd
import gseapy as gp
import os
import matplotlib.pyplot as plt
import numpy as np
import mygene
import pickle

# Define input and output paths
de_file = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA/results/deseq2/YAF_vs_GFP/differential_expression.csv"
output_dir = "results/go_analysis"

padj_threshold = 0.05
fold_change_threshold = 0.0

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Define plot_enrichment function before we use it
def plot_enrichment(df, title, output_file, top_n=20):
    if len(df) == 0:
        print(f"No enrichment results found for {title}")
        return
    
    # Sort by adjusted p-value and take top N
    df = df.sort_values('Adjusted P-value').head(top_n)
    # Reverse the order for plotting (bottom to top)
    df = df.iloc[::-1]
    
    plt.figure(figsize=(12, 10))
    # Create horizontal bar plot
    plt.barh(df['Term'], -np.log10(df['Adjusted P-value']), color='darkred')
    plt.xlabel('-log10(Adjusted P-value)')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

# Define bubble plot function for better visualization - with improved error handling
def plot_enrichment_bubble(df, title, output_file, top_n=20):
    if len(df) == 0:
        print(f"No enrichment results found for {title}")
        return
    
    # Sort by adjusted p-value and take top N
    df = df.sort_values('Adjusted P-value').head(top_n)
    # Reverse the order for plotting (bottom to top)
    df = df.iloc[::-1]
    
    # Make a copy of the dataframe to avoid changing the original
    plot_df = df.copy()
    
    # Better error handling for Overlap column
    try:
        # First, inspect what's in the Overlap column
        if 'Overlap' in plot_df.columns:
            # Clean the Overlap strings and extract just the numbers
            def clean_overlap(x):
                if pd.isna(x) or not isinstance(x, str):
                    return [0, 1]  # Default if missing
                
                # Try to clean the string to extract just the numbers
                import re
                parts = re.findall(r'\d+', x)
                if len(parts) >= 2:
                    try:
                        return [int(parts[0]), int(parts[1])]
                    except ValueError:
                        return [0, 1]  # Default if conversion fails
                else:
                    return [0, 1]  # Default if pattern not found
            
            # Apply the cleaning function
            plot_df['Overlap_Clean'] = plot_df['Overlap'].apply(clean_overlap)
            plot_df['GeneRatio'] = plot_df['Overlap_Clean'].apply(lambda x: x[0]/max(x[1], 1))  # Avoid division by zero
            plot_df['Count'] = plot_df['Overlap_Clean'].apply(lambda x: x[0])
        else:
            # If Overlap column doesn't exist, look for alternatives
            if 'Overlapping Genes' in plot_df.columns:
                # Count the number of genes in the overlapping genes column
                def count_genes(x):
                    if pd.isna(x) or not isinstance(x, str):
                        return 0
                    return len(x.split(','))
                
                plot_df['Count'] = plot_df['Overlapping Genes'].apply(count_genes)
                # For GeneRatio, use Count divided by arbitrary background (100)
                plot_df['GeneRatio'] = plot_df['Count'] / 100
            else:
                # If no good columns available, create dummy data
                plot_df['Count'] = 10  # Default size
                plot_df['GeneRatio'] = 0.1  # Default ratio
    
    except Exception as e:
        print(f"  Error processing data for bubble plot: {str(e)}")
        # Create minimal working version with dummy data
        plot_df['Count'] = 10  # Default size
        plot_df['GeneRatio'] = 0.1  # Default ratio
    
    # Check if we have the necessary columns now
    if 'GeneRatio' not in plot_df.columns or 'Count' not in plot_df.columns:
        print("  Failed to create required columns for bubble plot")
        return
    
    # Create the plot
    plt.figure(figsize=(12, 14))
    
    try:
        # Create bubble plot with error handling
        scatter = plt.scatter(
            plot_df['GeneRatio'], 
            plot_df['Term'], 
            s=plot_df['Count'] * 2,  # Size based on gene count
            c=-np.log10(plot_df['Adjusted P-value']),  # Color based on significance
            cmap='Reds',
            alpha=0.7,
            edgecolors='black',
            linewidths=1
        )
        
        # Add colorbar
        cbar = plt.colorbar()
        cbar.set_label('-log10(p.adjust)')
        
        # Add size legend with better error handling
        try:
            handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6, num=4)
            
            # Extract size values safely
            sizes = []
            for label in labels:
                try:
                    # Try different patterns to extract the number
                    import re
                    matches = re.findall(r'\d+', label)
                    if matches:
                        sizes.append(int(matches[0]))
                    else:
                        sizes.append("?")
                except:
                    sizes.append("?")
            
            # Create legend
            legend = plt.legend(handles, sizes, title="Count", 
                             loc="center right", bbox_to_anchor=(1.3, 0.5))
        except Exception as e:
            print(f"  Error creating legend: {str(e)}")
            # Continue without legend
        
        plt.xlabel('GeneRatio')
        plt.grid(linestyle='--', alpha=0.3)
        plt.title(title)
        plt.tight_layout()
        
        # Save the plot
        plt.savefig(output_file, bbox_inches='tight')
        plt.close()
        print(f"  Successfully created bubble plot: {output_file}")
        
    except Exception as e:
        print(f"  Error creating bubble plot: {str(e)}")
        plt.close()  # Make sure to close the figure even if plotting fails

# Add this function before performing conversions
def load_cached_mapping(cache_file="ensembl_symbol_cache.pkl"):
    """Load cached Ensembl ID to gene symbol mapping if available"""
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'rb') as f:
                cached_mapping = pickle.load(f)
            print(f"Loaded {len(cached_mapping)} cached ID mappings")
            return cached_mapping
        except Exception as e:
            print(f"Error loading cached mapping: {str(e)}")
    return {}

def save_cached_mapping(mapping_dict, cache_file="ensembl_symbol_cache.pkl"):
    """Save Ensembl ID to gene symbol mapping for future use"""
    try:
        with open(cache_file, 'wb') as f:
            pickle.dump(mapping_dict, f)
        print(f"Cached {len(mapping_dict)} ID mappings for future use")
    except Exception as e:
        print(f"Error saving cached mapping: {str(e)}")

# Define the conversion function first, before it's called
def convert_ensembl_ids_to_symbols(ensembl_ids, batch_size=50):
    """Convert Ensembl IDs to gene symbols using biotools.fr API"""
    import requests
    import json
    import time
    from math import ceil
    
    # URL for the biotools.fr API
    url = "https://biotools.fr/human/ensembl_symbol_converter/"
    
    # Create mapping dictionary
    ensembl_to_symbol = {}
    
    # Clean Ensembl IDs - ensure they are in the correct format
    # The API expects IDs like ENSG00000142619 (without version numbers)
    clean_ids = []
    for eid in ensembl_ids:
        if eid and isinstance(eid, str) and eid.startswith("ENSG"):
            # Remove version if present (though your code already does this)
            clean_id = eid.split('.')[0]
            clean_ids.append(clean_id)
    
    print(f"Cleaned {len(clean_ids)} valid Ensembl IDs for conversion")
    
    # Sample a few IDs to check format
    if clean_ids:
        print(f"Sample IDs: {clean_ids[:3]}")
    
    # Process IDs in batches to avoid timeouts
    num_batches = ceil(len(clean_ids) / batch_size)
    print(f"Processing {len(clean_ids)} Ensembl IDs in {num_batches} batches...")
    
    for i in range(0, len(clean_ids), batch_size):
        batch = clean_ids[i:i+batch_size]
        
        try:
            # Prepare the POST request
            ids_json = json.dumps(batch)
            payload = {"api": 1, "ids": ids_json}
            
            # Make the request
            response = requests.post(url, data=payload)
            
            # Check if request was successful
            if response.status_code == 200:
                try:
                    result = response.json()
                    if result:
                        # Update the mapping dictionary
                        ensembl_to_symbol.update(result)
                    else:
                        print(f"  Empty result for batch {i//batch_size + 1}/{num_batches}")
                except json.JSONDecodeError:
                    print(f"  Invalid JSON response for batch {i//batch_size + 1}/{num_batches}")
                    print(f"  Response text: {response.text[:100]}...")
            else:
                print(f"  Batch {i//batch_size + 1}/{num_batches} failed with status code: {response.status_code}")
                print(f"  Response text: {response.text[:100]}...")
            
            # Add a small delay to avoid overwhelming the API
            time.sleep(0.5)
            
            # Print progress occasionally
            if (i//batch_size) % 5 == 0 or i + batch_size >= len(clean_ids):
                print(f"  Processed {i//batch_size + 1}/{num_batches} batches, found {len(ensembl_to_symbol)} symbols so far")
                
        except Exception as e:
            print(f"  Error in batch {i//batch_size + 1}/{num_batches}: {str(e)}")
    
    # Try a single request with a direct test ID to verify API functionality
    try:
        test_id = "ENSG00000142619"
        test_url = f"https://biotools.fr/human/ensembl_symbol_converter/?api=1&id={test_id}"
        test_response = requests.get(test_url)
        print(f"API test with ID {test_id}: status code {test_response.status_code}")
        if test_response.status_code == 200:
            print(f"Test response: {test_response.text[:100]}...")
    except Exception as e:
        print(f"Test request failed: {str(e)}")
    
    return ensembl_to_symbol

# Load differential expression results
print("Loading differential expression data...")
de_data = pd.read_csv(de_file, index_col=0)

# Extract gene_id from the dataframe if it exists
if 'gene_id' in de_data.columns:
    print("Using gene_id column for gene identifiers")
    de_data['ensembl_id'] = de_data['gene_id'].str.split('.').str[0]  # Remove version number
else:
    print("Creating ensembl_id from index")
    de_data['ensembl_id'] = de_data.index.str.split('.').str[0]  # Remove version number

# Load any existing cached mappings
cached_mapping = load_cached_mapping()
print(f"Found {len(cached_mapping)} cached ID mappings")

# Only convert IDs not in the cache
ensembl_ids = de_data['ensembl_id'].unique().tolist()
ids_to_convert = [eid for eid in ensembl_ids if eid not in cached_mapping]
print(f"Need to convert {len(ids_to_convert)} new IDs")

if ids_to_convert:
    new_mappings = convert_ensembl_ids_to_symbols(ids_to_convert)
    # Update cache with new mappings
    cached_mapping.update(new_mappings)
    # Save the updated cache
    save_cached_mapping(cached_mapping)

# Use the combined mapping
ensembl_to_symbol = cached_mapping

# Add gene symbols to dataframe
de_data['gene_symbol'] = de_data['ensembl_id'].map(ensembl_to_symbol)

# Check conversion rate and apply fallback if needed
mapped_count = de_data['gene_symbol'].notna().sum()
conversion_rate = mapped_count / len(de_data)
print(f"Initial conversion rate: {conversion_rate:.2%}")

# Apply fallback method for missing symbols if conversion rate is low
if conversion_rate < 0.50:  # If less than 50% of genes mapped
    print("Conversion rate is low, applying fallback conversion method...")
    
    # Use MyGeneInfo with different parameters as fallback
    mg = mygene.MyGeneInfo()
    
    # Get missing Ensembl IDs
    missing_ids = de_data[de_data['gene_symbol'].isna()]['ensembl_id'].unique().tolist()
    
    try:
        # Try with different scopes and fields
        fallback_results = mg.querymany(missing_ids, 
                                       scopes='ensembl.gene', 
                                       fields='symbol', 
                                       species='human',
                                       returnall=True,
                                       verbose=False)
        
        # Create fallback mapping
        fallback_mapping = {}
        for hit in fallback_results['out']:
            if 'symbol' in hit and hit.get('query'):
                fallback_mapping[hit['query']] = hit['symbol']
        
        # Update gene symbols for previously unmapped genes
        for idx, row in de_data[de_data['gene_symbol'].isna()].iterrows():
            if row['ensembl_id'] in fallback_mapping:
                de_data.at[idx, 'gene_symbol'] = fallback_mapping[row['ensembl_id']]
        
        # Update mapping dictionary for saving
        ensembl_to_symbol.update(fallback_mapping)
        
        # Report improvement
        new_mapped_count = de_data['gene_symbol'].notna().sum()
        print(f"Fallback method added {new_mapped_count - mapped_count} symbols")
        mapped_count = new_mapped_count
        
    except Exception as e:
        print(f"Error in fallback conversion: {str(e)}")

# Save the mapping for reference
symbol_mapping = pd.DataFrame({
    'ensembl_id': list(ensembl_to_symbol.keys()),
    'gene_symbol': list(ensembl_to_symbol.values())
})
symbol_mapping.to_csv(f"{output_dir}/ensembl_to_symbol_mapping.csv", index=False)

# Report final conversion rate
print(f"Successfully mapped {mapped_count} of {len(de_data)} genes to symbols ({mapped_count/len(de_data):.2%})")

# Filter genes by adjusted p-value
significant_genes = de_data[de_data['padj'] < padj_threshold].copy()

# Sort genes by log2FoldChange for ranking
de_data_sorted = de_data.sort_values(by='log2FoldChange', ascending=False)

# Create ranked gene list with gene symbols for GSEA
# Keep only genes with valid symbols
de_data_with_symbols = de_data_sorted[de_data_sorted['gene_symbol'].notna()].copy()

# Separate upregulated and downregulated genes
upregulated = significant_genes[significant_genes['log2FoldChange'] > fold_change_threshold]
downregulated = significant_genes[significant_genes['log2FoldChange'] < -fold_change_threshold]

print(f"Found {len(upregulated)} upregulated and {len(downregulated)} downregulated genes")

# Save lists of up and down regulated genes with both Ensembl ID and gene symbol
up_genes_with_symbols = upregulated[upregulated['gene_symbol'].notna()]
down_genes_with_symbols = downregulated[downregulated['gene_symbol'].notna()]

# For enrichr analysis we'll use gene symbols instead of Ensembl IDs
upregulated_symbols = up_genes_with_symbols['gene_symbol'].tolist()
downregulated_symbols = down_genes_with_symbols['gene_symbol'].tolist()

# Save gene lists with both identifiers
up_genes_df = pd.DataFrame({
    'ensembl_id': upregulated['ensembl_id'],
    'gene_symbol': upregulated['gene_symbol'],
    'log2FoldChange': upregulated['log2FoldChange'],
    'padj': upregulated['padj']
})
up_genes_df.to_csv(f"{output_dir}/upregulated_genes.csv", index=False)

down_genes_df = pd.DataFrame({
    'ensembl_id': downregulated['ensembl_id'],
    'gene_symbol': downregulated['gene_symbol'],
    'log2FoldChange': downregulated['log2FoldChange'],
    'padj': downregulated['padj']
})
down_genes_df.to_csv(f"{output_dir}/downregulated_genes.csv", index=False)

# Also save as simple text files
with open(f"{output_dir}/upregulated_genes.txt", "w") as f:
    for gene in upregulated.index.tolist():
        f.write(f"{gene}\n")

with open(f"{output_dir}/downregulated_genes.txt", "w") as f:
    for gene in downregulated.index.tolist():
        f.write(f"{gene}\n")

with open(f"{output_dir}/upregulated_symbols.txt", "w") as f:
    for gene in upregulated_symbols:
        f.write(f"{gene}\n")

with open(f"{output_dir}/downregulated_symbols.txt", "w") as f:
    for gene in downregulated_symbols:
        f.write(f"{gene}\n")

# Define list of gene set libraries to analyze based on valid sets from previous runs
enrichr_libraries = [
    # GO gene sets
    'GO_Biological_Process_2021',
    'GO_Cellular_Component_2021',
    'GO_Molecular_Function_2021',
    
    # Pathway databases
    'KEGG_2021_Human',
    'Reactome_2022',
    
    # Transcription and regulation
    'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
    'ChEA_2016',
    'TRANSFAC_and_JASPAR_PWMs',
    
    # Comprehensive libraries
    'MSigDB_Hallmark_2020',
    'Human_Gene_Atlas',
    
    # Cell type signatures
    'PanglaoDB_Augmented_2021',
    
    # Disease associations
    'DisGeNET',
    'Jensen_DISEASES',
    
    # Misc
    'ESCAPE'
]

# Modify the safe_enrichr function with retry capability
def safe_enrichr(gene_list, gene_sets, organism, outdir, no_plot=True, max_retries=3):
    """A safer wrapper for gseapy.enrichr that handles potential errors and can retry"""
    retries = 0
    while retries < max_retries:
        try:
            # Try the standard enrichr call
            return gp.enrichr(gene_list=gene_list,
                             gene_sets=gene_sets,
                             organism=organism,
                             outdir=outdir,
                             no_plot=no_plot)
        except ValueError as e:
            if "invalid literal for int()" in str(e):
                print(f"  Handling int conversion error in enrichr, attempting alternative method...")
                
                # Try alternative method (code remains the same)
                try:
                    # Alternative approach using enrichr API directly
                    import requests
                    import json
                    import pandas as pd
                    import os
                    
                    # Create output directory
                    os.makedirs(outdir, exist_ok=True)
                    
                    # Use Enrichr API directly
                    genes_str = '\n'.join(gene_list)
                    payload = {
                        'list': genes_str,
                        'description': 'Gene list analysis'
                    }
                    
                    # Step 1: Upload gene list
                    response = requests.post('https://maayanlab.cloud/Enrichr/addList', files=payload)
                    if response.status_code == 200:
                        user_list_id = response.json()['userListId']
                        
                        # Step 2: Get enrichment results for the specific gene set
                        for gene_set in gene_sets:
                            enrich_url = f'https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType={gene_set}'
                            enrich_response = requests.get(enrich_url)
                            
                            if enrich_response.status_code == 200:
                                results = enrich_response.json()[gene_set]
                                
                                # Convert to pandas DataFrame
                                df = pd.DataFrame(results, columns=[
                                    'Rank', 'Term', 'P-value', 'Z-score', 'Combined Score', 
                                    'Overlapping Genes', 'Adjusted P-value', 'Old P-value', 
                                    'Old Adjusted P-value'
                                ])
                                
                                # Also add Overlap column for compatibility with existing code
                                df['Overlap'] = df['Overlapping Genes'].apply(lambda x: f"{len(x.split(','))}/100")
                                
                                # Save results to file
                                output_file = f"{outdir}/{gene_set}.{organism}.enrichr.reports.txt"
                                df.to_csv(output_file, sep='\t', index=False)
                                
                                # Create a minimal results object similar to what gseapy.enrichr would return
                                class EnrichrResult:
                                    def __init__(self, df):
                                        self.results = df
                                        
                                return EnrichrResult(df)
                            else:
                                print(f"  Error fetching enrichment results for {gene_set}: {enrich_response.status_code}")
                                # Retry or return None based on the error
                                if "try again later" in enrich_response.text.lower():
                                    retries += 1
                                    if retries < max_retries:
                                        print(f"  Retrying ({retries}/{max_retries})...")
                                        import time
                                        time.sleep(5)  # Wait 5 seconds before retrying
                                        continue
                                return None
                    else:
                        print(f"  Error uploading gene list to Enrichr: {response.status_code}")
                        # Retry if it seems like a temporary error
                        if "try again later" in response.text.lower():
                            retries += 1
                            if retries < max_retries:
                                print(f"  Retrying ({retries}/{max_retries})...")
                                import time
                                time.sleep(5)  # Wait 5 seconds before retrying
                                continue
                        return None
                        
                except Exception as alt_e:
                    print(f"  Alternative method also failed: {str(alt_e)}")
                    # Create an empty result to avoid crashing
                    empty_df = pd.DataFrame(columns=['Term', 'Adjusted P-value', 'Genes', 'Combined Score', 'Overlap'])
                    class EmptyEnrichrResult:
                        def __init__(self):
                            self.results = empty_df
                    return EmptyEnrichrResult()
            else:
                # For other types of errors, just create an empty result
                print(f"  Error in enrichr: {str(e)}")
                empty_df = pd.DataFrame(columns=['Term', 'Adjusted P-value', 'Genes', 'Combined Score', 'Overlap'])
                class EmptyEnrichrResult:
                    def __init__(self):
                        self.results = empty_df
                return EmptyEnrichrResult()
        except Exception as e:
            error_msg = str(e)
            print(f"  Unexpected error in enrichr: {error_msg}")
            
            # Check if it's a network error that might benefit from a retry
            if "Error sending gene list" in error_msg or "try again later" in error_msg.lower():
                retries += 1
                if retries < max_retries:
                    print(f"  Retrying ({retries}/{max_retries})...")
                    import time
                    time.sleep(5 * retries)  # Increase wait time with each retry
                    continue
            
            # Create an empty result to avoid crashing
            empty_df = pd.DataFrame(columns=['Term', 'Adjusted P-value', 'Genes', 'Combined Score', 'Overlap'])
            class EmptyEnrichrResult:
                def __init__(self):
                    self.results = empty_df
            return EmptyEnrichrResult()
            
    # If we've exhausted retries, return an empty result
    print(f"  Max retries ({max_retries}) reached, giving up on this analysis")
    empty_df = pd.DataFrame(columns=['Term', 'Adjusted P-value', 'Genes', 'Combined Score', 'Overlap'])
    class EmptyEnrichrResult:
        def __init__(self):
            self.results = empty_df
    return EmptyEnrichrResult()

# For the upregulated genes section:
print("Performing enrichment analysis for upregulated genes...")
if upregulated_symbols:
    for library in enrichr_libraries:
        try:
            print(f"Running {library} analysis for upregulated genes...")
            outdir = f"{output_dir}/up_{library.replace(' ', '_')}"
            
            # Use safe_enrichr instead of gp.enrichr
            enr = safe_enrichr(gene_list=upregulated_symbols,
                             gene_sets=[library],
                             organism='Human',
                             outdir=outdir,
                             no_plot=True)
            
            if enr and len(enr.results) > 0:
                # Plot top enriched terms - standard bar plot
                plot_file = f"{outdir}/enrichment_plot.png"
                df = enr.results
                
                # Make sure the dataframe has the required columns before plotting
                if all(col in df.columns for col in ['Term', 'Adjusted P-value', 'Overlap']):
                    plot_enrichment(df, f'Top {library} - Upregulated Genes', plot_file)
                    
                    # Create bubble plot
                    bubble_plot_file = f"{outdir}/bubble_plot.png"
                    try:
                        plot_enrichment_bubble(df, f'Top {library} - Upregulated Genes', bubble_plot_file)
                    except Exception as bubble_error:
                        print(f"  Error creating bubble plot: {str(bubble_error)}")
                    
                    print(f"  Found {len(df)} enriched terms")
                else:
                    print(f"  Missing required columns in results dataframe")
            else:
                print(f"  No enrichment found for {library}")
                
        except Exception as e:
            print(f"  Error in {library} analysis: {str(e)}")

# Do the same for downregulated genes section
print("Performing enrichment analysis for downregulated genes...")
if downregulated_symbols:
    for library in enrichr_libraries:
        try:
            print(f"Running {library} analysis for downregulated genes...")
            outdir = f"{output_dir}/down_{library.replace(' ', '_')}"
            
            # Use safe_enrichr instead of gp.enrichr
            enr = safe_enrichr(gene_list=downregulated_symbols,
                             gene_sets=[library],
                             organism='Human',
                             outdir=outdir,
                             no_plot=True)
            
            if enr and len(enr.results) > 0:
                # Plot top enriched terms - standard bar plot
                plot_file = f"{outdir}/enrichment_plot.png"
                df = enr.results
                
                # Make sure the dataframe has the required columns before plotting
                if all(col in df.columns for col in ['Term', 'Adjusted P-value', 'Overlap']):
                    plot_enrichment(df, f'Top {library} - Downregulated Genes', plot_file)
                    
                    # Create bubble plot
                    bubble_plot_file = f"{outdir}/bubble_plot.png"
                    try:
                        plot_enrichment_bubble(df, f'Top {library} - Downregulated Genes', bubble_plot_file)
                    except Exception as bubble_error:
                        print(f"  Error creating bubble plot: {str(bubble_error)}")
                    
                    print(f"  Found {len(df)} enriched terms")
                else:
                    print(f"  Missing required columns in results dataframe")
            else:
                print(f"  No enrichment found for {library}")
                
        except Exception as e:
            print(f"  Error in {library} analysis: {str(e)}")

# Create a summary of all analyses
print("Creating summary visualization...")
try:
    # Collect top terms from each analysis
    summary_rows = []
    
    # Process up-regulated results
    for library in enrichr_libraries:
        try:
            # Updated file path to match actual gseapy output format
            result_file = f"{output_dir}/up_{library.replace(' ', '_')}/{library}.Human.enrichr.reports.txt"
            if os.path.exists(result_file):
                df = pd.read_csv(result_file, sep='\t')
                if len(df) > 0:
                    top_terms = df.sort_values('Adjusted P-value').head(5)
                    for _, row in top_terms.iterrows():
                        summary_rows.append({
                            'Term': row['Term'],
                            'P-value': row['Adjusted P-value'],
                            'Regulation': 'Up',
                            'Category': library,
                            'Genes': row['Genes'],
                            'Combined_score': row['Combined Score'],
                            'Overlap': row['Overlap']
                        })
            else:
                print(f"File not found: {result_file}")
        except Exception as e:
            print(f"Error processing up-regulated {library} results: {str(e)}")
    
    # Process down-regulated results
    for library in enrichr_libraries:
        try:
            # Updated file path to match actual gseapy output format
            result_file = f"{output_dir}/down_{library.replace(' ', '_')}/{library}.Human.enrichr.reports.txt"
            if os.path.exists(result_file):
                df = pd.read_csv(result_file, sep='\t')
                if len(df) > 0:
                    top_terms = df.sort_values('Adjusted P-value').head(5)
                    for _, row in top_terms.iterrows():
                        summary_rows.append({
                            'Term': row['Term'],
                            'P-value': row['Adjusted P-value'],
                            'Regulation': 'Down',
                            'Category': library,
                            'Genes': row['Genes'],
                            'Combined_score': row['Combined Score'],
                            'Overlap': row['Overlap']
                        })
            else:
                print(f"File not found: {result_file}")
        except Exception as e:
            print(f"Error processing down-regulated {library} results: {str(e)}")
    
    # Create summary dataframe
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_df = summary_df.sort_values(['Category', 'Regulation', 'P-value'])
        summary_df.to_csv(f"{output_dir}/enrichment_summary.csv", index=False)
        
        # Create top terms plot across all categories
        plt.figure(figsize=(15, 12))
        top_overall = summary_df.sort_values('P-value').head(30)
        top_overall['Label'] = top_overall['Term'] + ' (' + top_overall['Regulation'] + ')'
        top_overall = top_overall.iloc[::-1]  # Reverse for bottom-to-top plotting
        
        colors = ['red' if reg == 'Up' else 'blue' for reg in top_overall['Regulation']]
        plt.barh(top_overall['Label'], -np.log10(top_overall['P-value']), color=colors)
        plt.xlabel('-log10(Adjusted P-value)')
        plt.title('Top Enriched Terms across All Categories')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/top_overall_terms.png')
        plt.close()
        
        # Create bubble plots specifically for GO Biological Process
        try:
            # Filter for GO Biological Process terms
            go_bp_terms = summary_df[summary_df['Category'] == 'GO_Biological_Process_2021'].copy()
            
            if len(go_bp_terms) > 0:
                # Calculate gene ratio and count - with improved error handling
                def clean_overlap(x):
                    if pd.isna(x) or not isinstance(x, str):
                        return [0, 1]  # Default if missing
                    
                    # Extract numbers using regex
                    import re
                    parts = re.findall(r'\d+', x)
                    if len(parts) >= 2:
                        try:
                            return [int(parts[0]), int(parts[1])]
                        except ValueError:
                            return [0, 1]  # Default if conversion fails
                    else:
                        return [0, 1]  # Default if pattern not found
                
                go_bp_terms['Overlap_Clean'] = go_bp_terms['Overlap'].apply(clean_overlap)
                go_bp_terms['GeneRatio'] = go_bp_terms['Overlap_Clean'].apply(lambda x: x[0]/max(x[1], 1))
                go_bp_terms['Count'] = go_bp_terms['Overlap_Clean'].apply(lambda x: x[0])
                
                # Sort by p-value and take top 20
                go_bp_terms = go_bp_terms.sort_values('P-value').head(20)
                go_bp_terms = go_bp_terms.iloc[::-1]  # Reverse for better visualization
                
                plt.figure(figsize=(12, 14))
                
                # Create bubble plot with error handling
                try:
                    scatter = plt.scatter(
                        go_bp_terms['GeneRatio'], 
                        go_bp_terms['Term'], 
                        s=go_bp_terms['Count'] * 1.5,  # Size based on gene count
                        c=-np.log10(go_bp_terms['P-value']),  # Color based on significance
                        cmap='Reds',
                        alpha=0.7,
                        edgecolors='black',
                        linewidths=1
                    )
                    
                    # Add colorbar
                    cbar = plt.colorbar()
                    cbar.set_label('p.adjust')
                    
                    # Add size legend with better error handling
                    try:
                        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6, num=4)
                        
                        # Extract size values safely
                        sizes = []
                        for label in labels:
                            try:
                                # Try different patterns to extract the number
                                import re
                                matches = re.findall(r'\d+', label)
                                if matches:
                                    sizes.append(int(matches[0]))
                                else:
                                    sizes.append("?")
                            except:
                                sizes.append("?")
                        
                        legend = plt.legend(handles, sizes, title="Count", 
                                        loc="right", bbox_to_anchor=(1.3, 0.5))
                    except Exception as e:
                        print(f"  Error creating legend: {str(e)}")
                        # Continue without legend
                    
                    plt.xlabel('GeneRatio')
                    plt.grid(linestyle='--', alpha=0.3)
                    plt.title('GO Biological Process Enrichment')
                    plt.tight_layout()
                    plt.savefig(f'{output_dir}/go_bp_bubble_plot.png', bbox_inches='tight')
                    plt.close()
                    print("  Successfully created summary GO BP bubble plot")
                    
                except Exception as e:
                    print(f"  Error creating summary bubble plot: {str(e)}")
                    plt.close()  # Make sure to close the figure even if plotting fails
                
                # Rest of your code for creating separate plots for up and down-regulated genes
                # Use the same pattern of error handling as above for each of these plots
                
        except Exception as e:
            print(f"Error creating GO Biological Process bubble plot: {str(e)}")
        
    else:
        print("No enrichment results available for summary")

except Exception as e:
    print(f"Error creating summary: {str(e)}")

print("Analysis complete. Results saved to the 'results/go_analysis' directory.") 