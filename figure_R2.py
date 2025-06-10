import pandas as pd
import plotnine as p9
from plotnine import element_text
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from itertools import product
import requests
from plotnine import ggplot, aes, theme_minimal, geom_bar, geom_boxplot, element_text, annotate, geom_jitter, theme_bw, scale_color_manual, labs, ggtitle, theme, scale_fill_manual, geom_text, geom_segment
from collections import Counter
import os
import pandas as pd
import itertools
import warnings
warnings.filterwarnings("ignore")

def rotate_string(s):
    if not s:
        return []
    rotations = [s[i:] + s[:i] for i in range(len(s))]
    return rotations
	
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])
    
def unit(seq):
    all_strings = rotate_string(seq)
    all_strings.extend([x[::-1] for x in all_strings])
    all_strings.sort()
    return all_strings[0]

def is_self_reverse_complement(motif):
    """Check if an motif is its own reverse complement."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_comp = "".join(complement[base] for base in reversed(motif))
    return motif == rev_comp

def generate_motifs_by_complementarity(max_length):
    """Generate all RNA motifs and separate them into self-reverse-complementary and non-complementary."""
    bases = ['A', 'T', 'G', 'C']
    self_complementary = set()
    non_complementary = set()
    for length in range(1, max_length + 1):
        for motif in product(bases, repeat=length):
            motif_str = "".join(motif)
            if is_self_reverse_complement(motif_str):
                self_complementary.add(motif_str)
            else:
                non_complementary.add(motif_str)
    return sorted(self_complementary), sorted(non_complementary)

#############################################################################
## read in the datasets and generate rev and non-rev complementary motifs
#############################################################################

## generate reverse complementary and other sequence motifs
reverse_complementary, non_reverse_complementary = generate_motifs_by_complementarity(4)

## read in the datasets
allele_info = pd.read_csv('data/allele_data.tsv', delimiter='\t')
result_info = pd.read_csv('data/result_data.tsv', delimiter='\t')
merged_data = pd.merge(allele_info, result_info, on='location', how='left')
merged_data = merged_data[merged_data['fSTR']=='YES']

# Initialize list to collect results
results = []

# Loop through each location
for location, group in merged_data.groupby('location'):
    within_diffs = []
    between_diffs = []
    
    # Group by cluster within this location
    clusters = group.groupby('cluster')
    
    # 1. Calculate within-cluster differences
    for _, cluster_data in clusters:
        alleles = cluster_data['allele'].tolist()
        # Get all pairwise differences within the cluster
        for a1, a2 in itertools.combinations(alleles, 2):
            diff = abs(a1 - a2)
            within_diffs.append(diff)
    
    # 2. Calculate between-cluster differences
    # Get all pairs between clusters
    cluster_items = {name: data['allele'].tolist() for name, data in clusters}
    cluster_names = list(cluster_items.keys())
    
    for i in range(len(cluster_names)):
        for j in range(i+1, len(cluster_names)):
            alleles1 = cluster_items[cluster_names[i]]
            alleles2 = cluster_items[cluster_names[j]]
            for a1 in alleles1:
                for a2 in alleles2:
                    diff = abs(a1 - a2)
                    between_diffs.append(diff)
    
    # Compute averages (handle empty cases)
    within_avg = np.mean(within_diffs) if within_diffs else np.nan
    between_avg = np.mean(between_diffs) if between_diffs else np.nan
    
    results.append({
        'location': location,
        'within_cluster_avg_diff': within_avg,
        'between_cluster_avg_diff': between_avg
    })

# Make a new DataFrame
result_df = pd.DataFrame(results)
result_df = result_df.fillna(0)
print(result_df)
count = (result_df['within_cluster_avg_diff'] < result_df['between_cluster_avg_diff']).sum()
print(f"Number of times within-cluster avg diff < between-cluster avg diff: {count}")

data = {
    'average_allele_difference': [
        result_df['within_cluster_avg_diff'].mean(),
        result_df['between_cluster_avg_diff'].mean()
    ],
    'type': ['within_clusters', 'between_clusters']
}

summary_df = pd.DataFrame(data)
print(summary_df)

bar_plot = (
    ggplot(summary_df, aes(x='type', y='average_allele_difference', fill='type')) +
    geom_bar(stat='identity') +
    theme_bw() +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=22),legend_title=element_text(size=24),legend_text=element_text(size=22),plot_title=element_text(size=28)) +
    labs(y='Mean Allele Difference', x='Comparison Type')
)

bar_plot.save('figure_R2.png', dpi=200, width=20, height=5, limitsize=False)