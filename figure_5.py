import random
import matplotlib.pyplot as plt
from itertools import product
from plotnine import ggplot, aes, geom_boxplot, element_text, annotate, geom_jitter, theme_bw, scale_color_manual, labs, ggtitle, theme, scale_fill_manual, geom_text, geom_segment
import RNA
import seaborn as sns
import os
import pandas as pd

def sequence_context(N):
    """Generate a random RNA sequence of length N."""
    nucleotides = ['A', 'U', 'C', 'G']
    return ''.join(random.choices(nucleotides, k=N))

def generate_rna(left_context,motif,variant,right_context):
    sequence = motif * 100
    my_sequnce = left_context + sequence[0:variant] + right_context
    return my_sequnce

def fold_RNA(rna_sequence):
    fc = RNA.fold_compound(rna_sequence)
    mfe_structure, mfe = fc.mfe()
    return mfe_structure

def is_self_reverse_complement(motif):
    """Check if an RNA motif is its own reverse complement."""
    complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    rev_comp = "".join(complement[base] for base in reversed(motif))
    return motif == rev_comp

def generate_motifs_by_complementarity(max_length):
    """Generate all RNA motifs and separate them into self-reverse-complementary and non-complementary."""
    bases = ['A', 'U', 'G', 'C']
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
    plot = ( # creating boxplot of rpkm values vs secondary structure assignments
        ggplot(accessibility, aes(x='allele', y='percent', group='allele', color='allele'))
        + geom_jitter(aes(fill='allele'), alpha=0.5, width=0.35, size=5, show_legend=True, color="black")  # Fixed fill color for jitter points
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
    plot.save(f'./figure_5{panel}.png', width=12, height=6, dpi=300)

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
    plt.figure(figsize=(8, 6))  # Adjust the size as needed
    sns.heatmap(matrix, annot=True, fmt=".3f", cmap="coolwarm", cbar=True, linewidths=0.5, vmin=0, vmax=.05)
    plt.title(title, fontsize=20)
    plt.xlabel("longer variant", fontsize=20)
    plt.ylabel("shorter variant", fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

## need access to bpRNA-align.py
if os.path.exists("bpRNA_align.py"):
    print("found bpRNA_align.py, continuing")
else:
    print("need access to the following scripts available at https://github.com/BLasher113/bpRNA_align/tree/main")
    print("bpRNA.pl")
    print("bpRNA_align.py")
    print("bpRNA_align_module.py")
    print("bpRNA_align_module.pyc")
    exit()

## generate reverse complementary and other sequence motifs
reverse_complementary, non_reverse_complementary = generate_motifs_by_complementarity(4)

############################################################
## figure 5a
############################################################
accessibility = { 'accessibility':[], 'allele':[] }
for sequence in range(1000):
    allele = random.randint(0,30)
    unit = random.choice(reverse_complementary)
    myRNA = generate_rna(sequence_context(50),unit,allele,sequence_context(50))
    structure = fold_RNA(myRNA)
    accessibility['accessibility'].append(structure.count("."))
    accessibility['allele'].append(allele)
accessibility = pd.DataFrame(accessibility)
accessibility['percent'] = accessibility['accessibility'] / ( 100 + accessibility['allele'] )
make_plot(accessibility,'reverse complementary simulation','a')

############################################################
## figure 5b
############################################################
matrix_info = { 'first_allele':[], 'second_allele':[], 'bpRNA_align_score':[], 'first_alignment':[], 'second_alignment':[] }

for i in range(100):
    reference_allele = random.randint(0,30)
    left_context = sequence_context(50)
    right_context = sequence_context(50)
    unit = random.choice(reverse_complementary)

    # write .dbn files for bpRNA
    for variant in range(5):    
        RNA_sequence = generate_rna(left_context,unit,reference_allele + variant,right_context)
        fc = RNA.fold_compound(RNA_sequence)
        mfe_structure, mfe = fc.mfe()
        with open(f"structure_{reference_allele + variant}.dbn", "w") as f:
            f.write(RNA_sequence + "\n")
            f.write(mfe_structure + "\n")
            
    # Run bpRNA_align on the .dbn files
    with open("input_files.txt", "w") as f:
        for variant in range(5):
            f.write(f"structure_{reference_allele + variant}.dbn\n")
    os.system("python3 bpRNA_align.py -a True -f input_files.txt -w 5 -o similarity_matrix.txt")
    bpRNA_result = pd.read_csv("similarity_matrix.txt", delim_whitespace=True)

    # save bpRNA-align results
    for i,comparison in bpRNA_result.iterrows():
        matrix_info['first_allele'].append(int(comparison['name_1'].replace('structure_','')))
        matrix_info['second_allele'].append(int(comparison['name_2'].replace('structure_','')))
        matrix_info['bpRNA_align_score'].append(comparison['score'])
        matrix_info['first_alignment'].append(comparison['alignment_1'])
        matrix_info['second_alignment'].append(comparison['alignment_2'])    

df = pd.DataFrame(matrix_info)
difference_tally = {}
for i,row in df.iterrows():
    tally_differences(row['first_alignment'],row['second_alignment'])
matrix = tallies_to_matrix(difference_tally)
row_sums = matrix.sum(axis=1)
normalized_matrix = matrix.div(row_sums, axis=0)
if not matrix.empty:
    save_heatmap(normalized_matrix, "reverse complementary simulation", f"figure_5b.png")

############################################################
## figure 5c
############################################################
accessibility = { 'accessibility':[], 'allele':[] }
for sequence in range(1000):
    allele = random.randint(0,30)
    unit = random.choice(non_reverse_complementary)
    myRNA = generate_rna(sequence_context(50),unit,allele,sequence_context(50))
    structure = fold_RNA(myRNA)
    accessibility['accessibility'].append(structure.count("."))
    accessibility['allele'].append(allele)
accessibility = pd.DataFrame(accessibility)
accessibility['percent'] = accessibility['accessibility'] / ( 100 + accessibility['allele'] )
make_plot(accessibility,'non reverse complementary simulation','c')

############################################################
## figure 5d
############################################################
matrix_info = { 'first_allele':[], 'second_allele':[], 'bpRNA_align_score':[], 'first_alignment':[], 'second_alignment':[] }

for i in range(100):
    reference_allele = random.randint(0,30)
    left_context = sequence_context(50)
    right_context = sequence_context(50)
    unit = random.choice(non_reverse_complementary)

    # write .dbn files for bpRNA
    for variant in range(5):    
        RNA_sequence = generate_rna(left_context,unit,reference_allele + variant,right_context)
        fc = RNA.fold_compound(RNA_sequence)
        mfe_structure, mfe = fc.mfe()
        with open(f"structure_{reference_allele + variant}.dbn", "w") as f:
            f.write(RNA_sequence + "\n")
            f.write(mfe_structure + "\n")
            
    # Run bpRNA_align on the .dbn files
    with open("input_files.txt", "w") as f:
        for variant in range(5):
            f.write(f"structure_{reference_allele + variant}.dbn\n")
    os.system("python3 bpRNA_align.py -a True -f input_files.txt -w 5 -o similarity_matrix.txt")    
    bpRNA_result = pd.read_csv("similarity_matrix.txt", delim_whitespace=True)

    # save bpRNA-align results
    for i,comparison in bpRNA_result.iterrows():
        matrix_info['first_allele'].append(int(comparison['name_1'].replace('structure_','')))
        matrix_info['second_allele'].append(int(comparison['name_2'].replace('structure_','')))
        matrix_info['bpRNA_align_score'].append(comparison['score'])
        matrix_info['first_alignment'].append(comparison['alignment_1'])
        matrix_info['second_alignment'].append(comparison['alignment_2'])    

df = pd.DataFrame(matrix_info)
difference_tally = {}
for i,row in df.iterrows():
    tally_differences(row['first_alignment'],row['second_alignment'])
matrix = tallies_to_matrix(difference_tally)
row_sums = matrix.sum(axis=1)
normalized_matrix = matrix.div(row_sums, axis=0)
if not matrix.empty:
    save_heatmap(normalized_matrix, "non reverse complementary simulation", f"figure_5d.png")