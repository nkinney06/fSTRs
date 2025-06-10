import pandas as pd
import plotnine as p9
from plotnine import element_text
import scipy.stats as stats
import seaborn as sns
from itertools import product
import matplotlib.pyplot as plt
import numpy as np

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

def tally_differences(string1, string2): # tally secondary structure changes due to insertion.
    if (string1 == string2):
        return
    length_string_1 = sum(1 for char in string1 if char.isalpha())
    length_string_2 = sum(1 for char in string2 if char.isalpha())
    for char1, char2 in zip(string1, string2):
        if char1 == '-' or char2 == '-':
            continue
        if (length_string_1 < length_string_2):
            pair = (char1, char2)
            difference_tally[pair] = difference_tally.get(pair, 0) + 1
    return difference_tally
    
def tallies_to_matrix(tallies):
    letters = set(char for pair in tallies.keys() for char in pair)
    letters = sorted(letters)  # Sort letters to create a consistent matrix
    matrix = pd.DataFrame(0, index=letters, columns=letters)
    for (char1, char2), count in tallies.items():
        matrix.loc[char1, char2] = count
    return matrix

def save_heatmap(matrix, title, filename="heatmap.png"):
    # left-handed stem (L),
	# right-handed stem (R),
    # internal loops (I), 
    # bulges (B), 
    # hairpins (H), 
    # multiloops (M), 
    # external loops (X), 
    # ends (E)
    plt.figure(figsize=(7, 5.25))  # Adjust the size as needed
    ax = sns.heatmap(matrix, annot=True, fmt=".2f", cmap="coolwarm", cbar=True, linewidths=0.5, vmin=0, vmax=.05)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=14)
    plt.title(title, fontsize=20)
    plt.xlabel("longer variant", fontsize=20)
    plt.ylabel("shorter variant", fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

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

## generate reverse complementary and other sequence motifs
reverse_complementary, non_reverse_complementary = generate_motifs_by_complementarity(4)

# Read the data
matrix_info = pd.read_csv('data/matrix_data.tsv', delimiter='\t')
result_info = pd.read_csv('data/result_data.tsv', delimiter='\t')
merged_data = pd.merge(matrix_info, result_info, on='location', how='left')

#############################################################################
## figure 3a
#############################################################################
df = merged_data
difference_tally = {}
for i,row in df.iterrows():
    tally_differences(row['first_alignment'],row['second_alignment'])
matrix = tallies_to_matrix(difference_tally)
row_sums = matrix.sum(axis=1)
normalized_matrix = matrix.div(row_sums, axis=0)
if not matrix.empty:
    save_heatmap(normalized_matrix, "all fSTRs", f"figure_4a.png")

#############################################################################
## figure 3b
#############################################################################
merged_data['seq'] = merged_data['unit'].apply(unit)
df = merged_data[merged_data['seq'] == 'A']
difference_tally = {}
for i,row in df.iterrows():
    tally_differences(row['first_alignment'],row['second_alignment'])
matrix = tallies_to_matrix(difference_tally)
row_sums = matrix.sum(axis=1)
normalized_matrix = matrix.div(row_sums, axis=0)
if not matrix.empty:
    save_heatmap(normalized_matrix, "poly-A fSTRs", f"figure_4b.png")

#############################################################################
## figure 3c
#############################################################################
merged_data['seq'] = merged_data['unit'].apply(unit)
df = merged_data[merged_data['seq'] == 'AT']
difference_tally = {}
for i,row in df.iterrows():
    tally_differences(row['first_alignment'],row['second_alignment'])
matrix = tallies_to_matrix(difference_tally)
row_sums = matrix.sum(axis=1)
normalized_matrix = matrix.div(row_sums, axis=0)
if not matrix.empty:
    save_heatmap(normalized_matrix, "poly-AT fSTRs", f"figure_4c.png")
    
#############################################################################
## figure 3e (I changed the order during manuscript prep)
#############################################################################
merged_data['seq'] = merged_data['unit'].apply(unit)
df = merged_data[merged_data['seq'].isin(non_reverse_complementary)]
difference_tally = {}
for i,row in df.iterrows():
    tally_differences(row['first_alignment'],row['second_alignment'])
matrix = tallies_to_matrix(difference_tally)
row_sums = matrix.sum(axis=1)
normalized_matrix = matrix.div(row_sums, axis=0)
if not matrix.empty:
    save_heatmap(normalized_matrix, "non reverse complementary fSTRs", f"figure_4e.png")
    
#############################################################################
## figure 3f
#############################################################################
merged_data['seq'] = merged_data['unit'].apply(unit)
df = merged_data[merged_data['seq'].isin(reverse_complementary)]
difference_tally = {}
for i,row in df.iterrows():
    tally_differences(row['first_alignment'],row['second_alignment'])
matrix = tallies_to_matrix(difference_tally)
row_sums = matrix.sum(axis=1)
normalized_matrix = matrix.div(row_sums, axis=0)
if not matrix.empty:
    save_heatmap(normalized_matrix, "reverse complementary fSTRs", f"figure_4f.png")
    
#############################################################################
## figure 3d
#############################################################################
merged_data['seq'] = merged_data['unit'].apply(unit)
df = merged_data[merged_data['seq'] == 'ACG']
difference_tally = {}
for i,row in df.iterrows():
    tally_differences(row['first_alignment'],row['second_alignment'])
matrix = tallies_to_matrix(difference_tally)
row_sums = matrix.sum(axis=1)
normalized_matrix = matrix.div(row_sums, axis=0)
if not matrix.empty:
    save_heatmap(normalized_matrix, "CAG fSTRs", f"figure_4d.png")
    
