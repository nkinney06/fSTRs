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
import os
import pandas as pd
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

def make_plot(accessibility,unit,panel):
    
    # Downsample: point density in figures 4c and 4d were downsampled for clarity
    # downsampled_data = accessibility.sample(frac=0.15, random_state=42)
    
    plot = ( # creating boxplot of rpkm values vs secondary structure assignments
        ggplot(accessibility, aes(x='allele', y='percent', group='allele', color='allele'))
        + geom_jitter(aes(fill='allele'), alpha=0.5, width=0.35, size=5, show_legend=True, color="black")  # Fixed fill color for jitter points
        # + geom_jitter(downsampled_data, aes(fill='allele'), alpha=0.5, width=0.35, size=5, show_legend=True, color="black")  
        + geom_boxplot(show_legend=False, color="black", alpha=0.7, outlier_shape='', fill="#ffffff")  # Fill color mapped to combinedAlleleCluster for boxplot
        + theme_bw()  # White background theme
        + theme(
            text=element_text(size=24),         # Increase all text
            axis_title=element_text(size=26),   # Adjust axis titles
            axis_text=element_text(size=22),    # Adjust axis tick labels
            legend_title=element_text(size=24), # Adjust legend title
            legend_text=element_text(size=22),  # Adjust legend items
            plot_title=element_text(size=28)    # Adjust plot title
        )
        + ggtitle(unit)  # Title with beta symbol , β={???}
        + labs(y="accessibility", x="allele length")  # Axis labels
    )
    plot.save(f'./figure_3{panel}.png', width=12, height=6, dpi=300)

#############################################################################
## read in the datasets and generate rev and non-rev complementary motifs
#############################################################################

## generate reverse complementary and other sequence motifs
reverse_complementary, non_reverse_complementary = generate_motifs_by_complementarity(4)

## read in the datasets
allele_info = pd.read_csv('data/allele_data.tsv', delimiter='\t')
result_info = pd.read_csv('data/result_data.tsv', delimiter='\t')
region_info = pd.read_csv('data/region_data.tsv', delimiter='\t')
merged_data = pd.merge(allele_info, result_info, on='location', how='left')
merged_data = merged_data[merged_data['fSTR']=='YES']
merged_data = pd.merge(allele_info, region_info, on='location', how='left')
merged_data['unit'] = merged_data['unit'].apply(unit)

#############################################################################
## figure 4a
#############################################################################
merged_data_filtered = merged_data[merged_data['unit'] == 'A']
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
corr_coeff, p_value = stats.pearsonr(accessibility['allele'], accessibility['percent'])
print(f'poly A: correlation coefficient = {corr_coeff}, p value = {p_value}')
make_plot(accessibility,'poly-A accessibility','a')

#############################################################################
## figure 4b
#############################################################################
merged_data_filtered = merged_data[merged_data['unit'] == 'AT']
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
corr_coeff, p_value = stats.pearsonr(accessibility['allele'], accessibility['percent'])
print(f'poly AT: correlation coefficient = {corr_coeff}, p value = {p_value}')
make_plot(accessibility,'poly-AT accessibility','b')

#############################################################################
## figure 4c
#############################################################################
merged_data_filtered = merged_data[merged_data['unit'].isin(non_reverse_complementary)]
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
corr_coeff, p_value = stats.pearsonr(accessibility['allele'], accessibility['percent'])
print(f'non-reverse complementary: correlation coefficient = {corr_coeff}, p value = {p_value}')
make_plot(accessibility,'non reverse complementary fSTRs','c')

#############################################################################
## figure 4d
#############################################################################
merged_data_filtered = merged_data[merged_data['unit'].isin(reverse_complementary)]
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
corr_coeff, p_value = stats.pearsonr(accessibility['allele'], accessibility['percent'])
print(f'reverse complementary: correlation coefficient = {corr_coeff}, p value = {p_value}')
make_plot(accessibility,'reverse complementary fSTRs','d')