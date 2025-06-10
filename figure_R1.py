import pandas as pd
import plotnine as p9
from plotnine import element_text
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from itertools import product
import requests
from plotnine import ggplot, aes, geom_boxplot, element_text, annotate, geom_jitter, theme_bw, scale_color_manual, labs, ggtitle, theme, scale_fill_manual, geom_text, geom_segment
from collections import Counter
from scipy.stats import chisquare
import os
import pandas as pd

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

def calc_corr(group):
    return group['allele'].corr(group['percent'])

#############################################################################
## read in the datasets and generate rev and non-rev complementary motifs
#############################################################################

## generate reverse complementary and other sequence motifs
reverse_complementary, non_reverse_complementary = generate_motifs_by_complementarity(3)
reverse_complementary, _ = generate_motifs_by_complementarity(4)

## read in the datasets
allele_info = pd.read_csv('data/allele_data.tsv', delimiter='\t')
result_info = pd.read_csv('data/result_data.tsv', delimiter='\t')
region_info = pd.read_csv('data/region_data.tsv', delimiter='\t')
merged_data = pd.merge(allele_info, result_info, on='location', how='left')
merged_data = merged_data[merged_data['fSTR']=='YES']
merged_data = pd.merge(allele_info, region_info, on='location', how='left')
merged_data['unit'] = merged_data['unit'].apply(unit)

#############################################################################
## figure S4a
#############################################################################
merged_data_filtered = merged_data[merged_data['unit'].isin(non_reverse_complementary + reverse_complementary)]

accessibility = { 'location':[], 'allele':[], 'unit':[], 'accessibility':[], 'length':[] }
for i,row in merged_data_filtered.iterrows():
    accessibility['location'].append(row['location'])
    accessibility['allele'].append(row['allele'])
    accessibility['accessibility'].append(row['structure'].count('.'))
    accessibility['unit'].append(row['unit'])
    accessibility['length'].append(len(row['structure']))
accessibility = pd.DataFrame(accessibility)
accessibility['percent'] = accessibility['accessibility'] / accessibility['length']
accessibility = accessibility[accessibility['length'] - 100 < 30]
accessibility = accessibility[accessibility['accessibility'] < 100]

# remove units with a low number of occurances in the experimental data
unit_counts = accessibility['unit'].value_counts()
valid_units = unit_counts[unit_counts >= 10].index
accessibility = accessibility[accessibility['unit'].isin(valid_units)]

# Group by 'unit' and apply the correlation function
result = accessibility.groupby('unit').apply(calc_corr).reset_index()
result.columns = ['unit', 'correlation']
result['type'] = result['unit'].apply(lambda x: 'reverse complementary' if x in reverse_complementary else 'non-reverse complementary')
result = result.sort_values('correlation', ascending=True)
result['unit'] = pd.Categorical(result['unit'], categories=result['unit'], ordered=True)

# Create the bar chart
plot = (
    p9.ggplot() +
    p9.geom_bar(data=result, mapping=p9.aes(x='unit', y='correlation', fill='type'), stat='identity') +
    p9.labs(title='array length vs accessibility correlation coefficient',x='motif',y='correlation',fill='type') +
    p9.theme_bw() +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=22),legend_title=element_text(size=24),legend_text=element_text(size=22),plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=90, hjust=.5))
)

# Save the plot as a 300 dpi PNG file
plot.save('figure_R1.png', dpi=200, width=20, height=5, limitsize=False)
