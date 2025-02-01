import pandas as pd
import gseapy as gp
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse

def perform_enrichr(gene_list, description, output_dir):
    """Perform enrichment analysis for a given gene list"""
    print(f"\nPerforming enrichment analysis for {description}...")
    
    genes = list(gene_list)
    
    if len(genes) < 10:
        print(f"Warning: Too few genes ({len(genes)}) for reliable enrichment analysis. Need at least 10 genes.")
        return None
    
    # Define the gene set libraries to use
    databases = [
        'GO_Biological_Process_2021',
        'GO_Molecular_Function_2021',
        'GO_Cellular_Component_2021',
        'KEGG_2021_Human',
        'WikiPathways_2021_Human',
        'Reactome_2022'
    ]
    
    # If none of the above work, fall back to these simpler versions
    fallback_databases = [
        'GO_Biological_Process_2018',
        'GO_Molecular_Function_2018',
        'GO_Cellular_Component_2018',
        'KEGG_2019_Human',
        'WikiPathways_2019_Human',
        'Reactome_2016'
    ]
    
    results_dict = {}
    # Try each database
    for db in databases + fallback_databases:
        if db in results_dict:  # Skip if we already have results for this database type
            continue
            
        try:
            print(f"Trying enrichment with {db}...")
            enr = gp.enrichr(gene_list=genes,
                            gene_sets=db,
                            organism='Human',
                            outdir=None,
                            cutoff=0.05,
                            no_plot=True)
        
            # Filter results to only include terms with at least 5 genes
            if not enr.results.empty:
                enr.results['Gene_Count'] = enr.results['Genes'].str.count(',') + 1
                enr.results = enr.results[enr.results['Gene_Count'] >= 5]
            
                if not enr.results.empty:  # Only save if we have results after filtering
                    results_dict[db] = enr.results
                    
                    # Save filtered results
                    safe_db_name = db.replace('/', '_').replace(' ', '_')
                    enr.results.to_csv(os.path.join(output_dir, f'enrichment_{description.lower()}_{safe_db_name.lower()}.csv'))
                    print(f"Successfully processed {db} with {len(enr.results)} enriched terms")
        except Exception as e:
            print(f"Error processing {db}: {str(e)}")
            continue

    def plot_top_terms(results, title, output_file):
        if len(results) == 0:
            print(f"No significant terms found for {title}")
            return
            
        # Sort by adjusted p-value and take top 15
        results = results.sort_values('Adjusted P-value').head(15)
        
        if len(results) == 0:
            print(f"No significant terms found for {title}")
            return
        
        # Create figure
        plt.figure(figsize=(12, 8))
        
        # Create color palette based on gene counts
        gene_counts = results['Genes'].str.count(',') + 1
        colors = plt.cm.viridis(np.linspace(0, 1, len(gene_counts)))
        
        # Create bar plot
        bars = plt.barh(range(len(results)), 
                       -np.log10(results['Adjusted P-value']),
                       color=colors)
        
        # Customize plot
        plt.yticks(range(len(results)), 
                  [term[:50] + '...' if len(term) > 50 else term for term in results['Term']])
        plt.xlabel('-log10(Adjusted P-value)')
        plt.title(f'{title}\nColor intensity indicates number of genes')
        
        # Add value labels for gene counts
        for i, bar in enumerate(bars):
            width = bar.get_width()
            gene_count = results.iloc[i]['Genes'].count(',') + 1
            plt.text(width, bar.get_y() + bar.get_height()/2,
                    f'{gene_count}',
                    ha='left', va='center', fontweight='bold')
        
        # Add colorbar to show gene count scale
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, 
                                  norm=plt.Normalize(vmin=min(gene_counts), 
                                                   vmax=max(gene_counts)))
        cbar = plt.colorbar(sm)
        cbar.set_label('Number of genes', rotation=270, labelpad=15)
        
        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()

    # Plot results for each database
    for db_name, results in results_dict.items():
        if not results.empty:
            plot_top_terms(
                results,
                f'Top GO Terms - {db_name.replace("_", " ")} - {description}',
                os.path.join(output_dir, f'go_terms_{description.lower()}_{db_name.lower()}.png')
            )
    
    # Create combined heatmap for biological process and molecular function
    def create_combined_heatmap(bp_results, mf_results, description):
        bp_top = bp_results.head(10) if not bp_results.empty else pd.DataFrame()
        mf_top = mf_results.head(10) if not mf_results.empty else pd.DataFrame()
        
        if bp_top.empty and mf_top.empty:
            print(f"No significant terms found for combined heatmap - {description}")
            return
            
        bp_top['Category'] = 'Biological Process'
        mf_top['Category'] = 'Molecular Function'
        combined = pd.concat([bp_top, mf_top])
        
        if not combined.empty:
            plt.figure(figsize=(12, len(combined) * 0.4))
            
            sns.scatterplot(data=combined,
                          y='Term',
                          x='Category',
                          size=-np.log10(combined['Adjusted P-value']),
                          hue=-np.log10(combined['Adjusted P-value']),
                          palette='viridis')
            
            plt.title(f'GO Term Enrichment - {description}')
            plt.xticks(rotation=45)
            plt.legend(title='-log10(Adj P-value)', bbox_to_anchor=(1.05, 1), loc='upper left')
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'go_terms_heatmap_{description.lower()}.png'),
                       bbox_inches='tight', dpi=300)
            plt.close()
    
    if 'GO_Biological_Process_2021' in results_dict and 'GO_Molecular_Function_2021' in results_dict:
        create_combined_heatmap(
            results_dict['GO_Biological_Process_2021'],
            results_dict['GO_Molecular_Function_2021'],
            description
        )
    
    return results_dict

def main():
    parser = argparse.ArgumentParser(description='Perform GO enrichment analysis on gene lists')
    parser.add_argument('--input', required=True, help='Input file containing gene list')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--description', required=True, help='Description for the analysis')
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read gene list
    with open(args.input, 'r') as f:
        genes = set(line.strip() for line in f if line.strip())
    
    print(f"Read {len(genes)} genes from {args.input}")
    print("Sample of genes:", list(genes)[:5])
    
    # Perform enrichment analysis
    perform_enrichr(genes, args.description, args.output_dir)

if __name__ == "__main__":
    main()
