import pandas as pd
import plotnine as p9
from plotnine import element_text
import scipy.stats as stats
import numpy as np

# Read the data
df = pd.read_csv('data/result_data.tsv', delimiter='\t')

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

# Calculate overall fSTR proportion
print("\noverall fSTR counts")
summary = df['fSTR'].value_counts()
print(summary)

print("\noverall efSTR counts")
summary = df['efSTR'].value_counts()
print(summary)

# proportions used for statistical testing
overall_fstr_counts = df['fSTR'].value_counts(normalize=True)
print("\noverall fSTR Proportion:")
print(overall_fstr_counts)

##############################################################################################################
## check significance of STRs by feature, unit, and length againt overal proportion of YES and NO fSTR entries
##############################################################################################################

def check_significance(row):
    observed = row.values
    total = sum(observed)
    expected_yes = total * overall_fstr_counts['YES']
    expected_no = total * overall_fstr_counts['NO']
    expected = [expected_yes, expected_no]
    chi2, p_value = stats.chisquare(observed, f_exp=expected)
    if ( observed[0] / sum(observed) > overall_fstr_counts['YES'] ):
        return ("+",p_value)
    else:
        return ("-",p_value)

################################################################################################
## figure 2a
################################################################################################

# Count YES and NO entries for each feature, save p-values
region_fstr_counts = df.groupby(['feature', 'fSTR']).size().unstack(fill_value=0)
region_p_values = {}
for feature in region_fstr_counts.index:
    row = region_fstr_counts.loc[feature, ['YES', 'NO']]
    region_p_values[feature] = check_significance(row)

# these lines sort the bars from largest to smallest, It was difficult to figure out and I dont even understand.
region_fstr_counts['total'] = region_fstr_counts.sum(axis=1)
sort_order = region_fstr_counts['total'].sort_values(ascending=False).index
region_fstr_counts_reset = region_fstr_counts.reset_index()
region_fstr_counts_melted = region_fstr_counts_reset.melt( id_vars='feature', var_name='fSTR_status', value_name='count' )
region_fstr_counts_melted = region_fstr_counts_melted[region_fstr_counts_melted['fSTR_status'] != 'total'] # Remove the 'total' column from plotting
region_fstr_counts_melted['count'] = pd.to_numeric(region_fstr_counts_melted['count'])
region_fstr_counts_melted = region_fstr_counts_melted.sort_values(['feature', 'count'])
region_fstr_counts_melted['feature'] = pd.Categorical( region_fstr_counts_melted['feature'], categories=sort_order, ordered=True )

# Prepare significance data for annotations
significance_data = []
for feature in sort_order:
    max_height = region_fstr_counts.loc[feature, ['YES', 'NO']].sum()
    p_value = region_p_values[feature][1]
    significance_data.append({'feature': feature,'y_pos': max_height + 2000,'significance': region_p_values[feature][0] if p_value < 0.05 else ''})
significance_df = pd.DataFrame(significance_data)

# Create the stacked bar chart
plot = (
    p9.ggplot() +
    p9.geom_bar(data=region_fstr_counts_melted, mapping=p9.aes(x='feature', y='count', fill='fSTR_status'), stat='identity', position='stack') +
    p9.geom_text(data=significance_df, mapping=p9.aes(x='feature', y='y_pos', label='significance'), color='black', size=20, nudge_y=50) +
    p9.labs(title='fSTRs by feature',x='Feature',y='Count',fill='fSTR') +
    p9.theme_bw() +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=22),legend_title=element_text(size=24),legend_text=element_text(size=22),plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=45, hjust=1))
)

# Save the plot as a 300 dpi PNG file
plot.save('figure_2a.png', dpi=300, width=7, height=3.5)

# Print p-values
print("\Feature P-values:")
for feature, p_value in region_p_values.items():
    print(f"{feature}: {p_value}")

#######################################################################
## figure 2b
#######################################################################

# Count YES and NO entries for each unit_length
unit_length_fstr_counts = df.groupby(['unit_length', 'fSTR']).size().unstack(fill_value=0)

# Calculate p-values for each unit_length
unit_length_p_values = {}
for unit_length in unit_length_fstr_counts.index:
    row = unit_length_fstr_counts.loc[unit_length, ['YES', 'NO']]
    unit_length_p_values[unit_length] = check_significance(row)

# these lines sort the bars from largest to smallest
unit_length_fstr_counts['total'] = unit_length_fstr_counts.sum(axis=1)
sort_order = [1,2,3,4,5,6,7,8] # sort_order = unit_length_fstr_counts['total'].sort_values(ascending=False).index
unit_length_fstr_counts_reset = unit_length_fstr_counts.reset_index()
unit_length_fstr_counts_melted = unit_length_fstr_counts_reset.melt(id_vars='unit_length', var_name='fSTR_status', value_name='count')
unit_length_fstr_counts_melted = unit_length_fstr_counts_melted[unit_length_fstr_counts_melted['fSTR_status'] != 'total']
unit_length_fstr_counts_melted['count'] = pd.to_numeric(unit_length_fstr_counts_melted['count'])
unit_length_fstr_counts_melted = unit_length_fstr_counts_melted.sort_values(['unit_length', 'count'])
unit_length_fstr_counts_melted['unit_length'] = pd.Categorical(unit_length_fstr_counts_melted['unit_length'],categories=sort_order,ordered=True)

# Prepare significance data for annotations
significance_data = []
for unit_length in sort_order:
    max_height = unit_length_fstr_counts.loc[unit_length, ['YES', 'NO']].sum()
    p_value = unit_length_p_values[unit_length][1]
    significance_data.append({'unit_length': unit_length,'y_pos': max_height + 1200,'significance': unit_length_p_values[unit_length][0] if p_value < 0.05 else ''})
significance_df = pd.DataFrame(significance_data)

# Create the stacked bar chart
plot = (
    p9.ggplot() +
    p9.geom_bar(data=unit_length_fstr_counts_melted, mapping=p9.aes(x='unit_length', y='count', fill='fSTR_status'), stat='identity', position='stack') +
    p9.geom_text(data=significance_df,mapping=p9.aes(x='unit_length', y='y_pos', label='significance'),color='black',size=20,nudge_y=42) +
    p9.labs(title='fSTRs by motif length', x='motif length', y='Count', fill='fSTR') +
    p9.theme_bw() +
    p9.theme(text=element_text(size=24), axis_title=element_text(size=26), axis_text=element_text(size=22), legend_title=element_text(size=24), legend_text=element_text(size=22), plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=0, hjust=.5))
)

# Save the plot as a 300 dpi PNG file
plot.save('figure_2b.png', dpi=300, width=7, height=3.5)

# Print p-values
print("\nunit_length P-values:")
for unit_length, p_value in unit_length_p_values.items():
    print(f"{unit_length}: {p_value}")

#######################################################################
## figure 2c
#######################################################################  

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
    
# generalize the sequence motifs
df['seq'] = df['unit'].apply(unit)

# count number of fSTRs for each unit
seq_fstr_counts = df.groupby(['seq', 'fSTR']).size().unstack(fill_value=0)

# Calculate p-values for each seq
seq_p_values = {}
for seq in seq_fstr_counts.index:
    row = seq_fstr_counts.loc[seq, ['YES', 'NO']]
    seq_p_values[seq] = check_significance(row)

# Calculate total count for sorting
seq_fstr_counts['total'] = seq_fstr_counts.sum(axis=1)
sort_order = seq_fstr_counts['total'].sort_values(ascending=False).index
seq_fstr_counts_reset = seq_fstr_counts.reset_index()
seq_fstr_counts_melted = seq_fstr_counts_reset.melt(id_vars='seq',var_name='fSTR_status',value_name='count')
seq_fstr_counts_melted = seq_fstr_counts_melted[seq_fstr_counts_melted['fSTR_status'] != 'total']
seq_fstr_counts_melted['count'] = pd.to_numeric(seq_fstr_counts_melted['count'])
seq_fstr_counts_melted = seq_fstr_counts_melted.sort_values(['seq', 'count'])
seq_fstr_counts_melted = seq_fstr_counts_melted[seq_fstr_counts_melted['seq'].str.len() < 4]

# remove low frequency units from the plot
threshold = 200
over_threshold_categories = seq_fstr_counts_melted.loc[seq_fstr_counts_melted['count'] > threshold, 'seq'].unique()
seq_fstr_counts_melted = seq_fstr_counts_melted[seq_fstr_counts_melted['count'] > threshold]
sort_order = [ x for x in sort_order if x in over_threshold_categories]
seq_fstr_counts_melted['seq'] = pd.Categorical(seq_fstr_counts_melted['seq'], categories=sort_order, ordered=True)

# Prepare significance data for annotations
significance_data = []
for seq in sort_order:
    max_height = seq_fstr_counts.loc[seq, ['YES', 'NO']].sum()
    p_value = seq_p_values[seq][1]
    significance_data.append({'seq': seq, 'y_pos': max_height + 500, 'significance': seq_p_values[seq][0] if p_value < 0.05 else ''})
significance_df = pd.DataFrame(significance_data)

# Create the stacked bar chart
plot = (
    p9.ggplot() +
    p9.geom_bar(data=seq_fstr_counts_melted, mapping=p9.aes(x='seq', y='count', fill='fSTR_status'), stat='identity', position='stack') +
    p9.geom_text(data=significance_df, mapping=p9.aes(x='seq', y='y_pos', label='significance'), color='black', size=20, nudge_y=4) +
    p9.labs(title='fSTRs by motif',x='motif',y='Count',fill='fSTR') +
    p9.theme_bw() +
    p9.theme(text=element_text(size=24),axis_title=element_text(size=26),axis_text=element_text(size=22),legend_title=element_text(size=24),legend_text=element_text(size=22),plot_title=element_text(size=28)) +
    p9.theme(axis_text_x=p9.element_text(rotation=90, hjust=.5))
)

# Save the plot as a 300 dpi PNG file
plot.save('figure_2c.png', dpi=200, width=28, height=5, limitsize=False)

# Print p-values
print("\nseq P-values:")
for seq, p_value in seq_p_values.items():
    if ( len(seq) < 4 ):
        print(f"{seq}: {p_value}")
