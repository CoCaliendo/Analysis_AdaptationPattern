#!/usr/bin/env python
# coding: utf-8

"""
Processing InterProScan Annotations for Chironomus riparius

This script processes InterProScan output files for different treatment groups
(Blau, Gold, Rot, Grün) to extract and analyze protein domain annotations.
It identifies unique domains, calculates their frequencies, and compares
annotations between treatment groups.

Author: Cosima Caliendo
"""

import pandas as pd
import matplotlib.pyplot as plt

# Define file paths and colors for visualization
BASE_PATH = "/path/to/annotation/files"
OUTPUT_PATH = "/path/to/output"
COLORS = {
    'rot': '#e76f51',   # Red/orange
    'blau': '#7e7f9a',  # Blue
    'gold': '#e9c46a',  # Yellow/gold
    'grun': '#81b29a'   # Green
}

# Function to load and process domain data
def process_domain_data(file_path, treatment, column_indices=(5, 11)):
    """
    Loads InterProScan output and extracts domain information.
    
    Args:
        file_path (str): Path to the InterProScan output file
        treatment (str): Treatment name (blau, gold, rot, grun)
        column_indices (tuple): Columns to extract (default: PFAM=5, InterPro=11)
    
    Returns:
        tuple: DataFrames with processed domain data
    """
    # Load data
    domains = pd.read_csv(file_path, sep="\t", header=None)
    
    # Process PFAM domains (column 5)
    pfam_counts = domains[column_indices[0]].value_counts().reset_index()
    pfam_counts.columns = ['Entry', 'Count']
    
    # Process InterPro domains (column 11)
    ip_domains = pd.DataFrame(domains[column_indices[1]])
    
    # Save outputs
    pfam_output = f"{OUTPUT_PATH}/domains_pfam_uniq_{treatment}_all_0.csv"
    ip_output = f"{OUTPUT_PATH}/domains_ip_uniq_{treatment}_all.csv"
    
    pfam_counts.to_csv(pfam_output, sep='\t', index=False)
    ip_domains.to_csv(ip_output, sep='\t', index=False)
    
    return pfam_counts, ip_domains

# Process domains for each treatment group
treatments = ['blau', 'gold', 'rot', 'grun']
domain_data = {}

for treatment in treatments:
    # Adjust file path as needed
    file_path = f"{BASE_PATH}/interproscan_output_{treatment}_all_0.tsv"
    
    # Handle the special character in 'grün'
    safe_treatment = treatment.replace('ü', 'u')
    
    # Process data
    pfam_counts, ip_domains = process_domain_data(file_path, safe_treatment)
    
    # Store results
    domain_data[safe_treatment] = {
        'pfam': pfam_counts,
        'ip': ip_domains
    }

# Create visualizations for the top 15 domains in each treatment
def create_domain_comparison_plot(data_dict, domain_type='pfam', n_top=15):
    """
    Creates a 2x2 subplot comparing top domains across treatments.
    
    Args:
        data_dict (dict): Dictionary containing domain data for each treatment
        domain_type (str): Type of domain to compare (pfam or ip)
        n_top (int): Number of top domains to include
    """
    # Set up the plot with 2x2 subplots
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    flattened_axes = fig.axes
    
    # Prepare data for plotting
    plot_data = []
    for treatment in ['rot', 'blau', 'gold', 'grun']:
        df = data_dict[treatment][domain_type].head(n_top)
        plot_data.append(df)
    
    # Plot histograms
    for i, ax in enumerate(axs.flatten()):
        treatment = treatments[i]
        ax.bar(plot_data[i].iloc[:, 0], plot_data[i].iloc[:, 1], color=COLORS[treatment])
        ax.set_title(treatment.capitalize(), fontsize=16)
        ax.set_xlabel('Protein Families', fontsize=16)
        ax.set_ylabel('Counts', fontsize=16)
        ax.set_xticklabels(plot_data[i].iloc[:, 0], rotation=45, ha='right', fontsize=11)
    
    # Adjust layout
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_PATH}/histograms_{domain_type}_all.pdf")
    plt.close()

# Create plots for PFAM domains
create_domain_comparison_plot(domain_data, domain_type='pfam')

# Process tail data (0.0001% threshold)
def process_tail_data():
    """
    Processes data from the 0.0001% tail of the distribution.
    This represents the most significant domains based on statistical thresholds.
    """
    tail_data = {}
    
    for treatment in treatments:
        # Adjust path as needed for tail data
        file_path = f"{OUTPUT_PATH}/domains_pfam_{treatment}_0.tsv"
        
        # Load tail data
        domains = pd.read_csv(file_path, sep="\t", header=None)
        
        # Extract InterPro domains
        ip_domains = pd.DataFrame(domains[11])
        ip_domains.to_csv(f"{OUTPUT_PATH}/domains_ip_uniq_{treatment}.csv", sep='\t', index=False)
        
        # Store for analysis
        tail_data[treatment] = ip_domains
    
    return tail_data

# Process tail data and create visualizations
tail_data = process_tail_data()

# Find overlapping domains between treatment groups
def find_domain_overlaps(domain_data):
    """
    Identifies domains that overlap between different treatment groups.
    
    Args:
        domain_data (dict): Dictionary containing domain data for each treatment
        
    Returns:
        dict: Dictionary of DataFrames showing overlapping domains
    """
    overlaps = {}
    treatments = list(domain_data.keys())
    
    # Pairwise overlaps
    for i in range(len(treatments)):
        for j in range(i+1, len(treatments)):
            t1, t2 = treatments[i], treatments[j]
            key = f"{t1}_{t2}"
            overlaps[key] = pd.merge(domain_data[t1]['pfam'], domain_data[t2]['pfam'], on='Entry')
    
    # Three-way overlaps
    for i in range(len(treatments)):
        for j in range(i+1, len(treatments)):
            for k in range(j+1, len(treatments)):
                t1, t2, t3 = treatments[i], treatments[j], treatments[k]
                key = f"{t1}_{t2}_{t3}"
                pair_key = f"{t1}_{t2}"
                overlaps[key] = pd.merge(overlaps[pair_key], domain_data[t3]['pfam'], on='Entry')
    
    # Four-way overlap
    if len(treatments) >= 4:
        three_way_key = f"{treatments[0]}_{treatments[1]}_{treatments[2]}"
        all_key = "_".join(treatments)
        overlaps[all_key] = pd.merge(overlaps[three_way_key], domain_data[treatments[3]]['pfam'], on='Entry')
    
    return overlaps

# Find domain overlaps
overlaps = find_domain_overlaps(domain_data)

# Save all overlapping domains
overlaps["all"].to_csv(f"{OUTPUT_PATH}/overlapping_genes_allFour.csv", sep='\t')

# Extract cytochrome P450 genes for each treatment
def extract_cytochrome_genes():
    """
    Extracts Cytochrome P450 genes from each treatment group for further analysis.
    These genes are often involved in detoxification and adaptation.
    """
    cytochrome_genes = {}
    
    for treatment in treatments:
        safe_treatment = treatment.replace('ü', 'u')
        
        # Get raw domains data with all columns
        file_path = f"{BASE_PATH}/interproscan_output_{safe_treatment}_all_0.tsv"
        domains = pd.read_csv(file_path, sep="\t", header=None)
        
        # Add column names for clarity
        domains.columns=['marker', 'smth1','count','db', 'db_nr','Entry',
                         '6','7','8','9','10','11','12','13']
        
        # Filter for cytochrome P450 genes
        cytochrome = domains[domains['Entry'].str.contains("Cytochrome P450")]
        cytochrome_genes[safe_treatment] = cytochrome
    
    return cytochrome_genes

# Extract cytochrome genes
cytochrome_genes = extract_cytochrome_genes()

# Save cytochrome genes data
for treatment, data in cytochrome_genes.items():
    data.to_csv(f"{OUTPUT_PATH}/cytochrome_p450_{treatment}.csv", sep='\t', index=False)

print("Analysis complete. Output files saved to", OUTPUT_PATH)
