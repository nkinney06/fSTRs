import pandas as pd
import plotnine as p9
from plotnine import element_text
import scipy.stats as stats
import numpy as np

# Read the data
df = pd.read_csv('data/result_data.tsv', delimiter='\t')
rg = pd.read_csv('data/region_data.tsv', delimiter='\t')

# need to reduce the number of annotation categories so place in priority order.
priority_order = {'cds':'cds', 'TSS':'cds', 'utr3':'utr3', 'utr5':'utr5', 'exon':'exon', 'nc_exon':'exon', 'intron':'intron', 'nc_intron':'intron'}
def get_highest_priority(fstr):
    fstr_components = fstr.split(",")
    for k,v in priority_order.items():
        if k in fstr_components:
            return v
    return fstr  # In case no priority match is found, return as is
df['feature'] = df['feature'].apply(get_highest_priority)

# make a unit length column
df['unit_length'] = df['unit'].str.len()

# add the amino acid column for coding STRs
df = pd.merge(df, rg[['location', 'amino']], on='location', how='left')

#######################################################################
## figure S3a characterize all coding STRs
#######################################################################

region_fstr_counts = df[(df['amino'] != 'NaN') & (df['fSTR'] == 'YES')].groupby(['amino']).size().reset_index(name='count')

# Define a function to extract the motif
def get_motif(seq):
    return ''.join(sorted(set(seq))) if len(set(seq)) > 1 else seq[0]  # Single-letter motif if all same

# Apply motif extraction
region_fstr_counts['motif'] = region_fstr_counts['amino'].apply(get_motif)

# Group by motif and sum counts
motif_counts = region_fstr_counts.groupby('motif', as_index=False)['count'].sum()

# save the motifs
valid_motifs = motif_counts['motif'].to_list()

plot = (
    p9.ggplot() +
    p9.geom_bar(data=motif_counts, mapping=p9.aes(x='motif', y='count', fill='motif'), stat='identity') +
    p9.geom_text(data=motif_counts, mapping=p9.aes(x='motif', y='count', label='count'), color='black', size=16, nudge_y=10) +
    p9.labs(title=str(sum(motif_counts['count'])) + ' coding fSTRs') +
    p9.theme_bw() +
    p9.theme(legend_position='none') +
    p9.coord_flip() +
    p9.labs(x='amino acid motif') + 
    p9.scale_y_continuous(limits=(0, 125)) +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=16),plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=90, hjust=.5))
)

plot.save('figure_S3a.png', dpi=300, width=5, height=20)

#######################################################################
## figure S3b characterize all coding STRs
#######################################################################

region_fstr_counts = df[df['amino'] != 'NaN'].groupby(['amino']).size().reset_index(name='count')

# Define a function to extract the motif
def get_motif(seq):
    return ''.join(sorted(set(seq))) if len(set(seq)) > 1 else seq[0]  # Single-letter motif if all same

# Apply motif extraction
region_fstr_counts['motif'] = region_fstr_counts['amino'].apply(get_motif)

# Group by motif and sum counts
motif_counts = region_fstr_counts.groupby('motif', as_index=False)['count'].sum()

# save the motifs
motif_counts = motif_counts[motif_counts['motif'].isin(valid_motifs)]

plot = (
    p9.ggplot() +
    p9.geom_bar(data=motif_counts, mapping=p9.aes(x='motif', y='count', fill='motif'), stat='identity') +
    p9.geom_text(data=motif_counts, mapping=p9.aes(x='motif', y='count', label='count'), color='black', size=16, nudge_y=20) +
    p9.labs(title=str(sum(motif_counts['count'])) + ' coding STRs') +
    p9.theme_bw() +
    p9.theme(legend_position='none') +
    p9.coord_flip() +
    p9.labs(x='amino acid motif') + 
    p9.scale_y_continuous(limits=(0, 255)) +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=16),plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=90, hjust=.5))
)

plot.save('figure_S3b.png', dpi=300, width=5, height=20)

#######################################################################
## figure S1f characterize all coding efSTRs
#######################################################################

region_fstr_counts = df[(df['amino'] != 'NaN') & (df['efSTR'] == 'YES')].groupby(['amino']).size().reset_index(name='count')

# Define a function to extract the motif
def get_motif(seq):
    return ''.join(sorted(set(seq))) if len(set(seq)) > 1 else seq[0]  # Single-letter motif if all same

# Apply motif extraction
region_fstr_counts['motif'] = region_fstr_counts['amino'].apply(get_motif)

# Group by motif and sum counts
motif_counts = region_fstr_counts.groupby('motif', as_index=False)['count'].sum()

# save the motifs
for motif in valid_motifs:
    if ( motif not in motif_counts['motif'].to_list() ):
        motif_counts = pd.concat([motif_counts, pd.DataFrame([{'motif': motif, 'count': 0}])], ignore_index=True)

plot = (
    p9.ggplot() +
    p9.geom_bar(data=motif_counts, mapping=p9.aes(x='motif', y='count', fill='motif'), stat='identity') +
    p9.geom_text(data=motif_counts, mapping=p9.aes(x='motif', y='count', label='count'), color='black', size=16, nudge_y=.5) +
    p9.labs(title=str(sum(motif_counts['count'])) + ' coding efSTRs') +
    p9.theme_bw() +
    p9.theme(legend_position='none') +
    p9.coord_flip() +
    p9.labs(x='amino acid motif') + 
    p9.scale_y_continuous(limits=(0, 10)) +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=16),plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=90, hjust=.5))
)

plot.save('figure_S3c.png', dpi=300, width=5, height=20)
exit()