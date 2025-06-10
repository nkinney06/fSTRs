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
## figure S1a characterize all transcribed STRs
#######################################################################

region_fstr_counts = df.groupby(['feature']).size().reset_index(name='count')

plot = (
    p9.ggplot() +
    p9.geom_bar(data=region_fstr_counts, mapping=p9.aes(x='feature', y='count', fill='feature'), stat='identity') +
    p9.geom_text(data=region_fstr_counts, mapping=p9.aes(x='feature', y='count', label='count'), color='black', size=16, nudge_y=2000) +
    p9.labs(title=str(sum(region_fstr_counts['count'])) + ' STRs') +
    p9.theme_bw() +
    p9.theme(legend_position='none') +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=14),plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=0, hjust=.5))
)

plot.save('figure_S1a.png', dpi=300, width=4, height=3.5)

#######################################################################
## figure S1b characterize fSTRs
#######################################################################

region_fstr_counts = df[df['fSTR'] == 'YES'].groupby(['feature']).size().reset_index(name='count')

plot = (
    p9.ggplot() +
    p9.geom_bar(data=region_fstr_counts, mapping=p9.aes(x='feature', y='count', fill='feature'), stat='identity') +
    p9.geom_text(data=region_fstr_counts, mapping=p9.aes(x='feature', y='count', label='count'), color='black', size=16, nudge_y=550) +
    p9.labs(title=str(sum(region_fstr_counts['count'])) + ' fSTRs') +
    p9.theme_bw() +
    p9.theme(legend_position='none') +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=14),plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=0, hjust=.5))
)

plot.save('figure_S1b.png', dpi=300, width=4, height=3.5)

#######################################################################
## figure S1c characterize fSTRs
#######################################################################

region_fstr_counts = df[df['efSTR'] == 'YES'].groupby(['feature']).size().reset_index(name='count')
exons = pd.DataFrame([{'feature': 'exon', 'count': 0}])
region_fstr_counts = pd.concat([region_fstr_counts, exons], ignore_index=True)

plot = (
    p9.ggplot() +
    p9.geom_bar(data=region_fstr_counts, mapping=p9.aes(x='feature', y='count', fill='feature'), stat='identity') +
    p9.geom_text(data=region_fstr_counts, mapping=p9.aes(x='feature', y='count', label='count'), color='black', size=16, nudge_y=12) +
    p9.labs(title=str(sum(region_fstr_counts['count'])) + ' efSTRs') +
    p9.theme_bw() +
    p9.theme(legend_position='none') +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=14),plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=0, hjust=.5))
)

plot.save('figure_S1c.png', dpi=300, width=4, height=3.5)

