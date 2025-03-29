#!/usr/bin/env python3

import pandas as pd
import gseapy as gp
import os
import matplotlib.pyplot as plt
import numpy as np
import mygene
import pickle
import argparse
import re # Import the 're' module globally
import textwrap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Define a function to parse arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Perform GO enrichment analysis on differential expression results.")
    parser.add_argument('--de_file', type=str,
                        default="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA/results/deseq2/YAF_vs_GFP/differential_expression.csv",
                        help="Path to the differential expression CSV file.")
    parser.add_argument('--output_dir', type=str, default="results/go_analysis",
                        help="Directory to save the analysis results and plots.")
    parser.add_argument('--padj_threshold', type=float, default=0.05,
                        help="Adjusted p-value threshold for significance.")
    parser.add_argument('--fold_change_threshold', type=float, default=0.0,
                        help="Log2 fold change threshold for defining up/downregulated genes (applied after padj filter).")
    parser.add_argument('--top_n_plots', type=int, default=20,
                        help="Number of top terms to display in plots.")
    parser.add_argument('--enrichr_libs', nargs='+',
                        default=[
                            'GO_Biological_Process_2021', 'GO_Cellular_Component_2021',
                            'GO_Molecular_Function_2021', 'KEGG_2021_Human', 'Reactome_2022'
                        ],
                        help="List of Enrichr gene set libraries to use.")
    parser.add_argument('--cache_file', type=str, default="ensembl_symbol_cache.pkl",
                        help="File to cache Ensembl ID to Symbol mappings.")
    parser.add_argument('--figure_width', type=float, default=12.0,
                        help="Width of output figures in inches.")
    parser.add_argument('--figure_height', type=float, default=8.0,
                        help="Base height of output figures in inches.")
    parser.add_argument('--bubble_scale', type=float, default=5.0,
                        help="Scaling factor for bubble size in bubble plots.")
    return parser.parse_args()

# Define the clean_overlap function here, outside other functions
def clean_overlap(x):
    """Cleans the 'Overlap' string (e.g., '5/20') into a list of integers [5, 20]."""
    if pd.isna(x) or not isinstance(x, str):
        return [0, 1]  # Default if missing or wrong type
    
    parts = re.findall(r'\d+', x) # Use re here
    if len(parts) >= 2:
        try:
            # Return the first two numbers found
            return [int(parts[0]), int(parts[1])]
        except ValueError:
            return [0, 1]  # Default if conversion fails
    else:
        # Handle cases like just '5' or missing denominator
        if len(parts) == 1:
             try:
                 # Assume denominator is 1 if only one number is found
                 return [int(parts[0]), 1]
             except ValueError:
                 return [0, 1]
        return [0, 1]  # Default if pattern not found or invalid

# --- Main script execution starts here ---
if __name__ == "__main__":
    args = parse_args()

    # Use arguments instead of hardcoded values
    de_file = args.de_file
    output_dir = args.output_dir
    padj_threshold = args.padj_threshold
    fold_change_threshold = args.fold_change_threshold
    top_n_plots = args.top_n_plots
    enrichr_libraries = args.enrichr_libs
    cache_file = args.cache_file
    figure_width = args.figure_width
    figure_height = args.figure_height
    bubble_scale = args.bubble_scale

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Define plot_enrichment function before we use it (with formatting updates)
    def plot_enrichment(df, title, output_file, top_n=20):
        if len(df) == 0:
            print(f"  No enrichment results found for {title}")
            return

        # Sort by adjusted p-value and take top N
        df_plot = df.sort_values('Adjusted P-value').head(top_n)
        # Reverse the order for plotting (bottom to top)
        df_plot = df_plot.iloc[::-1]

        # Calculate dynamic height based on number of terms
        dynamic_height = max(figure_height, len(df_plot) * 0.4)
        
        plt.figure(figsize=(figure_width, dynamic_height))
        # Create horizontal bar plot
        bars = plt.barh(df_plot['Term'], -np.log10(df_plot['Adjusted P-value']), color='salmon')
        plt.xlabel('-log10(Adjusted P-value)', fontsize=12)
        plt.ylabel('') # Remove Y-axis label, terms are labels
        plt.title(title, fontsize=14, fontweight='bold')
        plt.tick_params(axis='y', labelsize=10)
        plt.tick_params(axis='x', labelsize=10)
        plt.grid(axis='x', linestyle='--', alpha=0.6) # Add grid lines
        
        # Adjust left margin based on term length
        max_label_len = df_plot['Term'].str.len().max()
        if max_label_len > 50:
            plt.subplots_adjust(left=0.4)
        elif max_label_len > 30:
            plt.subplots_adjust(left=0.3)
        
        plt.tight_layout(pad=1.5) # Add padding
        plt.savefig(output_file, dpi=300) # Increase resolution
        plt.close()
        print(f"  Successfully created bar plot: {output_file}")


    # Define bubble plot function for better visualization - with improved formatting
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
            if 'Overlap' in plot_df.columns:
                # Use the globally defined clean_overlap function
                plot_df['Overlap_Clean'] = plot_df['Overlap'].apply(clean_overlap)
                plot_df['GeneRatio'] = plot_df['Overlap_Clean'].apply(lambda x: x[0] / max(x[1], 1) if isinstance(x, list) and len(x) == 2 else 0)
                plot_df['Count'] = plot_df['Overlap_Clean'].apply(lambda x: x[0] if isinstance(x, list) and len(x) == 2 else 0)
            else:
                plot_df['Count'] = 10  # Default size
                plot_df['GeneRatio'] = 0.1  # Default ratio
        except Exception as e:
            print(f"Error processing data for bubble plot: {str(e)}")
            plot_df['Count'] = 10  # Default size
            plot_df['GeneRatio'] = 0.1  # Default ratio

        # Check if we have the necessary columns now
        if 'GeneRatio' not in plot_df.columns or 'Count' not in plot_df.columns:
             print("  Failed to create required columns for bubble plot")
             return

        # Better handling of long term names
        plot_df['Term_Wrapped'] = plot_df['Term'].apply(lambda x: x if len(x) < 50 else '\n'.join(textwrap.wrap(x, width=50)))
        
        # Calculate dynamic height based on number of terms and term length
        # Increase height when terms are wrapped to multiple lines
        avg_lines_per_term = np.mean([len(term.split('\n')) for term in plot_df['Term_Wrapped']])
        dynamic_height = max(figure_height, len(plot_df) * 0.6 * avg_lines_per_term)
        
        # Increase figure width for wider plots
        figure_width_adj = figure_width * 1.5  # Increase width by 30%

        # Create the plot with improved layout
        plt.figure(figsize=(figure_width_adj, dynamic_height))
        plt.style.use('seaborn-v0_8-whitegrid')
        
        # Only show horizontal grid lines for better readability
        plt.grid(True, axis='x', linestyle='--', alpha=0.4)
        
        # Remove frame on right and top sides for cleaner look
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['left'].set_linewidth(1.2)
        plt.gca().spines['bottom'].set_linewidth(1.2)

        try:
            # Use viridis colormap for better color contrast and accessibility
            scatter = plt.scatter(
                plot_df['GeneRatio'], 
                plot_df['Term_Wrapped'], 
                s=plot_df['Count'] * bubble_scale,  # Use bubble_scale parameter
                c=-np.log10(plot_df['Adjusted P-value']), 
                cmap='viridis',  # Changed to viridis for better color perception
                alpha=0.8,  # Slightly increased opacity
                edgecolors='black', 
                linewidths=0.8
            )

            # Add colorbar with better formatting
            cbar = plt.colorbar(scatter, fraction=0.046, pad=0.04)
            cbar.set_label('-log10(Adjusted P-value)', fontsize=12, fontweight='bold')
            cbar.ax.tick_params(labelsize=10)
            
            # Improve colorbar ticks for better readability
            cbar.ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

            # Create a cleaner size legend with better spacing
            # Dynamically determine sizes_to_show based on actual data
            count_range = plot_df['Count'].max() - plot_df['Count'].min()
            if count_range > 100:
                sizes_to_show = [10, 25, 50, 100, 200]
            elif count_range > 50: 
                sizes_to_show = [10, 25, 50, 100]
            else:
                sizes_to_show = [5, 10, 20, 40]
                
            # Filter sizes_to_show to only include sizes that are in the data range
            max_count = plot_df['Count'].max()
            sizes_to_show = [size for size in sizes_to_show if size <= max(max_count * 1.2, 10)]
            
            size_handles = []
            for size in sizes_to_show:
                size_handles.append(
                    plt.scatter([], [], s=size*bubble_scale, c='gray', alpha=0.7, 
                               edgecolors='black', linewidths=0.8)
                )
            
            # Add size legend with better positioning
            legend = plt.legend(
                handles=size_handles,
                labels=sizes_to_show,
                title="Gene Count",
                loc="upper left",
                bbox_to_anchor=(1.05, 0.9),
                fontsize=10,
                frameon=True,
                framealpha=0.9,
                edgecolor='lightgray'
            )
            legend.get_title().set_fontsize(12)
            legend.get_title().set_fontweight('bold')

            plt.xlabel('Gene Ratio', fontsize=12, fontweight='bold')
            plt.ylabel('')  # No Y-axis label needed
            
            # Enhance title with more information if available
            title_parts = title.split('\n')
            if len(title_parts) >= 2:
                plt.suptitle(title_parts[0], fontsize=15, fontweight='bold', y=0.98)
                plt.title(title_parts[1], fontsize=13, fontweight='bold', pad=10)
            else:
                plt.title(title, fontsize=14, fontweight='bold', pad=10)
                
            # Improve tick parameters
            plt.tick_params(axis='y', labelsize=11, left=False)  # Remove Y ticks
            plt.tick_params(axis='x', labelsize=10, direction='out', width=1, length=4)
            
            # Set x-axis limits with padding for better visualization
            max_ratio = plot_df['GeneRatio'].max()
            plt.xlim([0, max_ratio * 1.15])
            
            # Add minor grid for x-axis
            plt.minorticks_on()
            plt.grid(True, axis='x', which='minor', linestyle=':', alpha=0.2)
            
            # Adjust left margin based on term length
            max_term_len = max([len(term) for term in plot_df['Term']])
            if max_term_len > 80:
                plt.subplots_adjust(left=0.4) 
            elif max_term_len > 50:
                plt.subplots_adjust(left=0.35)
            elif max_term_len > 30:
                plt.subplots_adjust(left=0.3)
            
            # Adjust layout for legend
            plt.tight_layout(pad=1.5)
            plt.subplots_adjust(right=0.72)  # Make room for legend and colorbar

            plt.savefig(output_file, bbox_inches='tight', dpi=300)
            plt.close()
            print(f"Successfully created bubble plot: {output_file}")

        except Exception as e:
            print(f"Error creating bubble plot: {str(e)}")
            plt.close()

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
    try:
        de_data = pd.read_csv(de_file, index_col=0)
    except FileNotFoundError:
        print(f"Error: Differential expression file not found at {de_file}")
        exit(1)
    except Exception as e:
        print(f"Error loading differential expression data: {e}")
        exit(1)

    # Extract gene_id from the dataframe if it exists
    if 'gene_id' in de_data.columns:
        print("Using gene_id column for gene identifiers")
        de_data['ensembl_id'] = de_data['gene_id'].str.split('.').str[0]  # Remove version number
    else:
        print("Creating ensembl_id from index")
        de_data['ensembl_id'] = de_data.index.str.split('.').str[0]  # Remove version number

    # Load any existing cached mappings
    cached_mapping = load_cached_mapping(cache_file)
    print(f"Found {len(cached_mapping)} cached ID mappings from {cache_file}")

    # Only convert IDs not in the cache
    ensembl_ids = de_data['ensembl_id'].unique().tolist()
    ids_to_convert = [eid for eid in ensembl_ids if eid not in cached_mapping]
    print(f"Need to convert {len(ids_to_convert)} new IDs")

    if ids_to_convert:
        new_mappings = convert_ensembl_ids_to_symbols(ids_to_convert)
        # Update cache with new mappings
        cached_mapping.update(new_mappings)
        # Save the updated cache
        save_cached_mapping(cached_mapping, cache_file)

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
    if len(de_data) > 0:
        print(f"Successfully mapped {mapped_count} of {len(de_data)} genes to symbols ({mapped_count/len(de_data):.2%})")
    else:
        print("No gene data loaded.")

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
    # enrichr_libraries = [...] # This is now defined by args.enrichr_libs

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
    print("\nPerforming enrichment analysis for upregulated genes...")
    if not upregulated_symbols:
        print("  No significant upregulated genes with symbols found.")
    else:
        print(f"  Analyzing {len(upregulated_symbols)} upregulated genes.")
        for library in enrichr_libraries:
            try:
                print(f"  Running {library} analysis...")
                outdir = f"{output_dir}/up_{library.replace(' ', '_').replace('/', '_')}" # Sanitize library name for dir

                # Use safe_enrichr instead of gp.enrichr
                enr = safe_enrichr(gene_list=upregulated_symbols,
                                 gene_sets=[library],
                                 organism='Human',
                                 outdir=outdir,
                                 no_plot=True)
                
                if enr and hasattr(enr, 'results') and len(enr.results) > 0:
                    # Plot top enriched terms - standard bar plot
                    plot_file = f"{outdir}/enrichment_plot.png"
                    df = enr.results

                    # Make sure the dataframe has the required columns before plotting
                    if all(col in df.columns for col in ['Term', 'Adjusted P-value', 'Overlap', 'Genes']):
                        # Use args.top_n_plots
                        plot_enrichment(df, f'Top {library}\nUpregulated Genes', plot_file, top_n=top_n_plots)

                        # Create bubble plot
                        bubble_plot_file = f"{outdir}/bubble_plot.png"
                        try:
                            # Use args.top_n_plots
                            plot_enrichment_bubble(df, f'Top {library}\nUpregulated Genes', bubble_plot_file, top_n=top_n_plots)
                        except Exception as bubble_error:
                            print(f"  Error creating bubble plot: {str(bubble_error)}")
                        
                        print(f"  Found {len(df)} enriched terms")
                    else:
                        print(f"  Missing required columns in results dataframe")
                else:
                    print(f"  No significant enrichment found for {library}")

            except Exception as e:
                print(f"  Error in {library} analysis for upregulated genes: {str(e)}")

    # Do the same for downregulated genes section
    print("\nPerforming enrichment analysis for downregulated genes...")
    if not downregulated_symbols:
         print("  No significant downregulated genes with symbols found.")
    else:
        print(f"  Analyzing {len(downregulated_symbols)} downregulated genes.")
        for library in enrichr_libraries:
            try:
                print(f"  Running {library} analysis...")
                outdir = f"{output_dir}/down_{library.replace(' ', '_').replace('/', '_')}" # Sanitize library name for dir

                # Use safe_enrichr instead of gp.enrichr
                enr = safe_enrichr(gene_list=downregulated_symbols,
                                 gene_sets=[library],
                                 organism='Human',
                                 outdir=outdir,
                                 no_plot=True)
                
                if enr and hasattr(enr, 'results') and len(enr.results) > 0:
                    # Plot top enriched terms - standard bar plot
                    plot_file = f"{outdir}/enrichment_plot.png"
                    df = enr.results

                    # Make sure the dataframe has the required columns before plotting
                    if all(col in df.columns for col in ['Term', 'Adjusted P-value', 'Overlap', 'Genes']):
                         # Use args.top_n_plots
                        plot_enrichment(df, f'Top {library}\nDownregulated Genes', plot_file, top_n=top_n_plots)

                        # Create bubble plot
                        bubble_plot_file = f"{outdir}/bubble_plot.png"
                        try:
                            # Use args.top_n_plots
                            plot_enrichment_bubble(df, f'Top {library}\nDownregulated Genes', bubble_plot_file, top_n=top_n_plots)
                        except Exception as bubble_error:
                            print(f"  Error creating bubble plot: {str(bubble_error)}")
                        
                        print(f"  Found {len(df)} enriched terms")
                    else:
                        print(f"  Missing required columns in results dataframe")
                else:
                    print(f"  No significant enrichment found for {library}")

            except Exception as e:
                 print(f"  Error in {library} analysis for downregulated genes: {str(e)}")

    # Create a summary of all analyses, combining results from up- and down-regulated gene sets.
    print("\nCreating summary visualization...")
    try:
        # Initialize an empty list to collect summary rows. Each row will represent a top enriched term.
        summary_rows = []
        
        # Process up-regulated results from each Enrichr library.
        for library in enrichr_libraries:
            try:
                # Sanitize the library name to create a valid directory name.
                sanitized_lib = library.replace(' ', '_').replace('/', '_')
                
                # Construct the file path to the Enrichr results file for the current library and up-regulated genes.
                result_file = f"{output_dir}/up_{sanitized_lib}/{library}.Human.enrichr.reports.txt"
                
                # Check if the result file exists.
                if os.path.exists(result_file):
                    # Read the result file into a pandas DataFrame.
                    df = pd.read_csv(result_file, sep='\t')
                    
                    # Check if the DataFrame is not empty.
                    if len(df) > 0:
                        # Sort the DataFrame by adjusted p-value and select the top 5 terms.
                        top_terms = df.sort_values('Adjusted P-value').head(5)
                        
                        # Iterate through the top terms and append them to the summary rows list.
                        for _, row in top_terms.iterrows():
                            summary_rows.append({
                                'Term': row['Term'],  # The enriched term.
                                'P-value': row['Adjusted P-value'],  # The adjusted p-value of the term.
                                'Regulation': 'Up',  # Indicates that this term is enriched in up-regulated genes.
                                'Category': library,  # The Enrichr library the term belongs to.
                                'Genes': row['Genes'],  # The genes associated with the term.
                                'Combined_score': row['Combined Score'], # The combined score from Enrichr
                                'Overlap': row['Overlap'] # The overlap between the gene list and the term
                            })
                else:
                    # If the result file is not found, print a message.
                    print(f"File not found: {result_file}")
            except Exception as e:
                # If an error occurs while processing the results, print an error message.
                print(f"  Error processing up-regulated {library} results: {str(e)}")
        
        # Process down-regulated results from each Enrichr library.  This mirrors the up-regulated processing.
        for library in enrichr_libraries:
            try:
                # Sanitize the library name for directory creation.
                sanitized_lib = library.replace(' ', '_').replace('/', '_')
                
                # Construct the file path to the Enrichr results file for the current library and down-regulated genes.
                result_file = f"{output_dir}/down_{sanitized_lib}/{library}.Human.enrichr.reports.txt"
                
                # Check if the result file exists.
                if os.path.exists(result_file):
                    # Read the result file into a pandas DataFrame.
                    df = pd.read_csv(result_file, sep='\t')
                    
                    # Check if the DataFrame is not empty.
                    if len(df) > 0:
                        # Sort the DataFrame by adjusted p-value and select the top 5 terms.
                        top_terms = df.sort_values('Adjusted P-value').head(5)
                        
                        # Iterate through the top terms and append them to the summary rows list.
                        for _, row in top_terms.iterrows():
                            summary_rows.append({
                                'Term': row['Term'],  # The enriched term.
                                'P-value': row['Adjusted P-value'],  # The adjusted p-value of the term.
                                'Regulation': 'Down',  # Indicates that this term is enriched in down-regulated genes.
                                'Category': library,  # The Enrichr library the term belongs to.
                                'Genes': row['Genes'],  # The genes associated with the term.
                                'Combined_score': row['Combined Score'], # The combined score from Enrichr
                                'Overlap': row['Overlap'] # The overlap between the gene list and the term
                            })
                else:
                    # If the result file is not found, print a message.
                    print(f"File not found: {result_file}")
            except Exception as e:
                # If an error occurs while processing the results, print an error message.
                print(f"  Error processing down-regulated {library} results: {str(e)}")
        
        # Create a summary DataFrame from the collected rows.
        if summary_rows:
            # Create a pandas DataFrame from the summary rows list.
            summary_df = pd.DataFrame(summary_rows)
            
            # Sort the DataFrame by category, regulation (Up/Down), and p-value.
            summary_df = summary_df.sort_values(['Category', 'Regulation', 'P-value'])
            
            # Save the summary DataFrame to a CSV file.
            summary_df.to_csv(f"{output_dir}/enrichment_summary.csv", index=False)
            print(f"  Saved enrichment summary to {output_dir}/enrichment_summary.csv")

            # Create a bar plot of the top enriched terms across all categories.
            # The plot visualizes the -log10(p-value) for each term, colored by regulation (Up/Down).
            plt.figure(figsize=(14, max(8, (top_n_plots + 10) * 0.4))) # Adjusted figure size
            top_overall = summary_df.sort_values('P-value').head(top_n_plots + 10) # Use argument + buffer
            top_overall['Label'] = top_overall['Term'].str[:60] + '...' + ' (' + top_overall['Regulation'] + ')' if top_overall['Term'].str.len().max() > 60 else top_overall['Term'] + ' (' + top_overall['Regulation'] + ')' # Truncate long labels
            top_overall = top_overall.iloc[::-1]  # Reverse for bottom-to-top plotting

            colors = ['#d62728' if reg == 'Up' else '#1f77b4' for reg in top_overall['Regulation']] # Specific colors
            plt.barh(top_overall['Label'], -np.log10(top_overall['P-value']), color=colors)
            plt.xlabel('-log10(Adjusted P-value)', fontsize=12)
            plt.ylabel('')
            plt.title('Top Enriched Terms (Up & Down)', fontsize=14, fontweight='bold') # Updated title
            plt.tick_params(axis='y', labelsize=10)
            plt.tick_params(axis='x', labelsize=10)
            plt.grid(axis='x', linestyle='--', alpha=0.6)
            plt.tight_layout(pad=1.5)
            # Adjust left margin if labels are long
            max_label_len = top_overall['Label'].str.len().max()
            if max_label_len > 60:
                 plt.subplots_adjust(left=0.4)
            elif max_label_len > 40:
                 plt.subplots_adjust(left=0.3)

            plt.savefig(f'{output_dir}/top_overall_terms.png', dpi=300)
            plt.close()
            print(f"  Saved overall top terms plot to {output_dir}/top_overall_terms.png")

            # Create summary bubble plots for each GO library
            go_libraries = ['GO_Biological_Process_2021', 'GO_Cellular_Component_2021',
                            'GO_Molecular_Function_2021', 'KEGG_2021_Human', 'Reactome_2022']
            
            # Process each library that's present in the enrichment results
            for go_bp_lib in go_libraries:
                if go_bp_lib in enrichr_libraries:
                    try:
                        # Filter for terms from the specified GO library
                        lib_terms_df = summary_df[summary_df['Category'] == go_bp_lib].copy()

                        if len(lib_terms_df) > 0:
                            # Process data for plotting
                            lib_terms_df['Overlap_Clean'] = lib_terms_df['Overlap'].apply(clean_overlap)
                            lib_terms_df['GeneRatio'] = lib_terms_df['Overlap_Clean'].apply(
                                lambda x: x[0] / max(x[1], 1) if isinstance(x, list) and len(x) == 2 else 0
                            )
                            lib_terms_df['Count'] = lib_terms_df['Overlap_Clean'].apply(
                                lambda x: x[0] if isinstance(x, list) and len(x) == 2 else 0
                            )
                            
                            # Sort by p-value and select top terms
                            lib_terms_df_plot = lib_terms_df.sort_values('P-value').head(top_n_plots).copy()
                            
                            # Extract GO IDs and clean term names
                            lib_terms_df_plot['GO_ID'] = lib_terms_df_plot['Term'].apply(
                                lambda x: re.search(r'\((GO:\d+)\)', x).group(1) if re.search(r'\((GO:\d+)\)', x) else ""
                            )
                            lib_terms_df_plot['Term_Clean'] = lib_terms_df_plot['Term'].apply(
                                lambda x: re.sub(r'\s*\(GO:\d+\)', '', x)
                            )
                            
                            # Create separate dataframes for up and down regulated terms
                            lib_terms_up = lib_terms_df_plot[lib_terms_df_plot['Regulation'] == 'Up']
                            lib_terms_down = lib_terms_df_plot[lib_terms_df_plot['Regulation'] == 'Down']
                            
                            # # Use plotly for interactive bubble plot
                            # import plotly.graph_objects as go
                            # from plotly.subplots import make_subplots
                            
                            # # Create figure with subplots
                            # fig = make_subplots(rows=1, cols=1)
                            
                            # # Add traces for upregulated terms
                            # if len(lib_terms_up) > 0:
                            #     fig.add_trace(
                            #         go.Scatter(
                            #             x=lib_terms_up['GeneRatio'],
                            #             y=lib_terms_up['Term_Clean'],
                            #             mode='markers',
                            #             marker=dict(
                            #                 size=lib_terms_up['Count'] * bubble_scale * 0.8,
                            #                 color=-np.log10(lib_terms_up['P-value']),
                            #                 colorscale='Reds',
                            #                 colorbar=dict(
                            #                     title="-log10(p-value)<br>Upregulated",
                            #                     x=1.02,
                            #                     y=0.8,
                            #                     len=0.4,
                            #                     thickness=15
                            #                 ),
                            #                 line=dict(width=1, color='darkred'),
                            #                 opacity=0.8,
                            #                 showscale=True
                            #             ),
                            #             text=[f"GO ID: {go_id}<br>Term: {term}<br>P-value: {p:.2e}<br>Gene Count: {count}<br>Gene Ratio: {ratio:.3f}" 
                            #                 for go_id, term, p, count, ratio in zip(
                            #                     lib_terms_up['GO_ID'], 
                            #                     lib_terms_up['Term_Clean'], 
                            #                     lib_terms_up['P-value'],
                            #                     lib_terms_up['Count'],
                            #                     lib_terms_up['GeneRatio']
                            #                 )],
                            #             hoverinfo='text',
                            #             name='Upregulated'
                            #         )
                            #     )
                            
                            # # Add traces for downregulated terms
                            # if len(lib_terms_down) > 0:
                            #     fig.add_trace(
                            #         go.Scatter(
                            #             x=lib_terms_down['GeneRatio'],
                            #             y=lib_terms_down['Term_Clean'],
                            #             mode='markers',
                            #             marker=dict(
                            #                 size=lib_terms_down['Count'] * bubble_scale * 0.8,
                            #                 color=-np.log10(lib_terms_down['P-value']),
                            #                 colorscale='Blues',
                            #                 colorbar=dict(
                            #                     title="-log10(p-value)<br>Downregulated",
                            #                     x=1.02,
                            #                     y=0.3,
                            #                     len=0.4,
                            #                     thickness=15
                            #                 ),
                            #                 line=dict(width=1, color='darkblue'),
                            #                 opacity=0.8,
                            #                 showscale=True
                            #             ),
                            #             text=[f"GO ID: {go_id}<br>Term: {term}<br>P-value: {p:.2e}<br>Gene Count: {count}<br>Gene Ratio: {ratio:.3f}" 
                            #                 for go_id, term, p, count, ratio in zip(
                            #                     lib_terms_down['GO_ID'], 
                            #                     lib_terms_down['Term_Clean'], 
                            #                     lib_terms_down['P-value'],
                            #                     lib_terms_down['Count'],
                            #                     lib_terms_down['GeneRatio']
                            #                 )],
                            #             hoverinfo='text',
                            #             name='Downregulated'
                            #         )
                            #     )
                            
                            # # Update layout
                            # fig.update_layout(
                            #     title={
                            #         'text': f'Top {go_bp_lib}<br><sup>Enrichment Analysis</sup>',
                            #         'y': 0.95,
                            #         'x': 0.5,
                            #         'xanchor': 'center',
                            #         'yanchor': 'top',
                            #         'font': dict(size=20)
                            #     },
                            #     xaxis=dict(
                            #         title='Gene Ratio',
                            #         gridcolor='lightgray',
                            #         gridwidth=0.5,
                            #         zeroline=False
                            #     ),
                            #     yaxis=dict(
                            #         title='',
                            #         gridcolor='white',
                            #         gridwidth=0,
                            #         zeroline=False
                            #     ),
                            #     plot_bgcolor='white',
                            #     hoverlabel=dict(
                            #         bgcolor="white",
                            #         font_size=12,
                            #         font_family="Arial"
                            #     ),
                            #     height=max(600, len(lib_terms_df_plot) * 30),
                            #     width=1000,
                            #     margin=dict(l=200, r=100, t=100, b=50),
                            #     showlegend=True,
                            #     legend=dict(
                            #         yanchor="top",
                            #         y=0.99,
                            #         xanchor="left",
                            #         x=0.01,
                            #         bgcolor="rgba(255, 255, 255, 0.8)",
                            #         bordercolor="lightgray",
                            #         borderwidth=1
                            #     )
                            # )
                            
                            # # Add size reference
                            # max_count = lib_terms_df_plot['Count'].max()
                            # sizes_to_show = [s for s in [10, 25, 50, 100] if s <= max_count]
                            
                            # # Add annotation for size reference
                            # if sizes_to_show:
                            #     size_note = "Bubble size represents gene count"
                            #     fig.add_annotation(
                            #         x=0.01,
                            #         y=-0.15,
                            #         xref="paper",
                            #         yref="paper",
                            #         text=size_note,
                            #         showarrow=False,
                            #         font=dict(size=12),
                            #         align="left"
                            #     )
                            
                            # # Save as interactive HTML and static PNG
                            # # output_html_file = f'{output_dir}/summary_{go_bp_lib}_bubble_plot.html'
                            # output_png_file = f'{output_dir}/summary_{go_bp_lib}_bubble_plot.png'
                            
                            # # fig.write_html(output_html_file)
                            # fig.write_image(output_png_file, scale=2)
                            
                            # # print(f"  Successfully created interactive {go_bp_lib} bubble plot: {output_html_file}")
                            # print(f"  Also saved static version: {output_png_file}")
                            
                            # Create a static version with matplotlib as backup
                            try:
                                import seaborn as sns
                                
                                # Set up the figure
                                plt.figure(figsize=(14, max(8, len(lib_terms_df_plot) * 0.6)))
                                sns.set_style("whitegrid")
                                
                                # Create the bubble plot using seaborn's scatterplot
                                g = sns.scatterplot(
                                    data=lib_terms_df_plot,
                                    x="GeneRatio",
                                    y="Term_Clean",
                                    size="Count",
                                    hue="P-value",
                                    palette="coolwarm_r",
                                    sizes=(20, bubble_scale * 150),
                                    alpha=0.7,
                                    edgecolor='gray',
                                    linewidth=0.5,
                                    legend="brief"
                                )
                                
                                # Customize the plot
                                plt.xlabel("Gene Ratio", fontsize=12, fontweight="bold")
                                plt.ylabel("")
                                plt.title(f"Top ${go_bp_lib}", fontsize=14, fontweight="bold")
                                
                                # Adjust legend
                                plt.legend(title="Parameters", bbox_to_anchor=(1.05, 1), loc='upper left')
                                
                                # Add GO IDs as annotations
                                for i, row in lib_terms_df_plot.iterrows():
                                    if row['GO_ID']:
                                        plt.annotate(
                                            row['GO_ID'],
                                            xy=(0, i),
                                            xytext=(-10, 0),
                                            textcoords="offset points",
                                            ha='right', va='center',
                                            fontsize=8, color='gray'
                                        )
                                
                                # Save the static backup plot
                                backup_file = f'{output_dir}/summary_{go_bp_lib}_bubble_plot_static.png'
                                plt.tight_layout()
                                plt.savefig(backup_file, dpi=300)
                                plt.close()
                                print(f"  Created static backup plot: {backup_file}")
                                
                            except Exception as e:
                                print(f"  Note: Could not create static backup plot: {str(e)}")
                                plt.close()
                        else:
                            print(f"  No results found for {go_bp_lib} in the summary.")
                    except Exception as e:
                        print(f"  Error creating summary bubble plot for {go_bp_lib}: {str(e)}")
                        import traceback
                        traceback.print_exc()
                        plt.close()
            else:
                print("  No enrichment results available for summary plot generation.")

        else:
            # If no enrichment results were found, print a message.
            print("  No enrichment results available for summary plot generation.")

    except Exception as e:
        # If an error occurs during the summary creation process, print an error message.
        print(f"Error creating summary: {str(e)}")

    # Indicate that the analysis is complete and print the output directory.
    print("\nAnalysis complete.")
    print(f"Results saved to the '{output_dir}' directory.")