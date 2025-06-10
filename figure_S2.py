import pandas as pd
import plotnine as p9
from plotnine import element_text
import scipy.stats as stats
import numpy as np

def rotate_string(s):
    if not s:
        return []
    rotations = [s[i:] + s[:i] for i in range(len(s))]
    return rotations
	
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])
    
def unit(unit):
    all_strings = rotate_string(unit)
    all_strings.extend([x[::-1] for x in all_strings])
    all_strings.sort()
    return all_strings[0]

figure_2_units = ['T','A','GT','AC','CT','AG','AT','GTT','AAC','CCT','AGG','C','G','ACG','CGT','ACT']

#######################################################################
## characterize sequence motifs for the entire set of STRs
#######################################################################

# Read the data
df = pd.read_csv('data/region_data.tsv', delimiter='\t')

# generalize the sequence motifs
df['seq'] = df['unit'].apply(unit)

# count number of fSTRs for each unit
seq_fstr_counts = df.groupby(['seq']).size().reset_index(name='count')

# use same motifs as figure 2
seq_fstr_counts = seq_fstr_counts[seq_fstr_counts['seq'].isin(figure_2_units)]

# Create the stacked bar chart
plot = (
    p9.ggplot() +
    p9.geom_bar(data=seq_fstr_counts, mapping=p9.aes(x='seq', y='count', fill='seq'), stat='identity') +
    p9.geom_text(data=seq_fstr_counts, mapping=p9.aes(x='seq', y='count', label='count'), color='black', size=16, nudge_y=6000) +
    p9.labs(title='STRs',x='sequence motif',y='Count',fill='seq') +
    p9.theme_bw() +
    p9.theme(legend_position='none') +
    p9.coord_flip() +
    p9.scale_y_continuous(limits=(0, 54000),breaks=range(0, 50000, 20000)) +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=18),plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=0, hjust=.5))
)

# Save the plot as a 300 dpi PNG file
plot.save('figure_S2a.png', dpi=200, width=5, height=10, limitsize=False)

#######################################################################
## characterize sequence motifs for fSTRs
#######################################################################

# Read the data
rd = pd.read_csv('data/result_data.tsv', delimiter='\t')
df = pd.read_csv('data/region_data.tsv', delimiter='\t')

# generalize the sequence motifs
df['seq'] = df['unit'].apply(unit)

# filter the fSTRs
df = pd.merge(df, rd[['location', 'fSTR']], on='location', how='left')
df = df[df['fSTR'] == 'YES']

# count number of fSTRs for each unit
seq_fstr_counts = df.groupby(['seq']).size().reset_index(name='count')

# use same motifs as figure 2
seq_fstr_counts = seq_fstr_counts[seq_fstr_counts['seq'].isin(figure_2_units)]

# Create the stacked bar chart
plot = (
    p9.ggplot() +
    p9.geom_bar(data=seq_fstr_counts, mapping=p9.aes(x='seq', y='count', fill='seq'), stat='identity') +
    p9.geom_text(data=seq_fstr_counts, mapping=p9.aes(x='seq', y='count', label='count'), color='black', size=16, nudge_y=300) +
    p9.labs(title='fSTRs',x='sequence motif',y='Count',fill='seq') +
    p9.theme_bw() +
    p9.theme(legend_position='none') +
    p9.coord_flip() +
    p9.scale_y_continuous(limits=(0, 3300)) +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=18),plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=0, hjust=.5))
)

# Save the plot as a 300 dpi PNG file
plot.save('figure_S2b.png', dpi=200, width=5, height=10, limitsize=False)

#######################################################################
## characterize sequence motifs for efSTRs
#######################################################################

# Read the data
rd = pd.read_csv('data/result_data.tsv', delimiter='\t')
df = pd.read_csv('data/region_data.tsv', delimiter='\t')

# generalize the sequence motifs
df['seq'] = df['unit'].apply(unit)

# filter the fSTRs
df = pd.merge(df, rd[['location', 'efSTR']], on='location', how='left')
df = df[df['efSTR'] == 'YES']

# count number of fSTRs for each unit
seq_fstr_counts = df.groupby(['seq']).size().reset_index(name='count')

# use same motifs as figure 2
seq_fstr_counts = seq_fstr_counts[seq_fstr_counts['seq'].isin(figure_2_units)]

# Create the stacked bar chart
plot = (
    p9.ggplot() +
    p9.geom_bar(data=seq_fstr_counts, mapping=p9.aes(x='seq', y='count', fill='seq'), stat='identity') +
    p9.geom_text(data=seq_fstr_counts, mapping=p9.aes(x='seq', y='count', label='count'), color='black', size=16, nudge_y=4) +
    p9.labs(title='efSTRs',x='sequence motif',y='Count',fill='seq') +
    p9.theme_bw() +
    p9.theme(legend_position='none') +
    p9.coord_flip() +
    p9.scale_y_continuous(limits=(0, 72)) +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=18),plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=0, hjust=.5))
)

# Save the plot as a 300 dpi PNG file
plot.save('figure_S2c.png', dpi=200, width=5, height=10, limitsize=False)