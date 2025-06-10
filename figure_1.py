import pandas as pd
from PIL import Image, ImageDraw
import os
import RNA
import time
from plotnine import ggplot, aes, geom_boxplot, element_text, annotate, geom_jitter, theme_bw, scale_color_manual, labs, ggtitle, theme, scale_fill_manual, geom_text, geom_segment
import glob
import argparse
from scipy.stats import f_oneway, kruskal
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import warnings
warnings.filterwarnings("ignore")

def ss_image(sequence, mfe_structure, filename, exemplar):
    RNA.PS_rna_plot(sequence, mfe_structure, filename + ".ps")
    try:
        image = Image.open(filename + ".ps")
        image = image.convert("RGB")
        width, height = image.size
        if ( exemplar > 0 ): # highlight the exemplar
            mask = Image.new("L", (width, height), 0)
            draw = ImageDraw.Draw(mask)
            outer_rectangle = (0, 0, width, height)  # Full image dimensions
            inner_rectangle = (10, 10, width - 10, height - 10)  # Inner border dimensions
            corner_radius = 30
            draw.rounded_rectangle(outer_rectangle, radius=corner_radius, fill=255)
            draw.rounded_rectangle(inner_rectangle, radius=corner_radius - 10, fill=0)
            black_image = Image.new("RGB", (width, height), "red")  # Solid black
            bordered_image = Image.composite(black_image, image, mask)
            bordered_image.save(filename + ".png")
        else:
            image.save(filename + ".png")
    except Exception as e:
        pass

def concatenate_images_vertically(image_paths, output_path="concatenated_image.png"):
    images = [Image.open(img) for img in image_paths]
    total_width = max(img.width for img in images)
    total_height = sum(img.height for img in images)
    new_image = Image.new("RGB", (total_width, total_height), "white")
    y_offset = 0
    for img in images:
        new_image.paste(img, (0, y_offset))
        y_offset += img.height
    new_image.save(output_path)

def concatenate_images_horizontally(image_paths, output_path):
    images = [Image.open(image_path) for image_path in image_paths]
    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    max_height = max(heights)
    new_image = Image.new('RGB', (total_width, max_height), "white")
    x_offset = 0
    for img in images:
        new_image.paste(img, (x_offset, 0))
        x_offset += img.width
    new_image.save(output_path)

####################################################################################################################
## read in the command line arguments (genomic location)
####################################################################################################################
parser = argparse.ArgumentParser(description='Process a genomic location.')
parser.add_argument(
    'location',
    nargs='?',  # makes the argument optional
    default='chr12:111447548-111447572',  # default value if not provided
    help='Genomic location (e.g., chr11:7672623-7672642)'
)

# Parse the arguments
args = parser.parse_args()

# Use the location
print(f"Using genomic location: {args.location}")

####################################################################################################################
## read in the datasets
####################################################################################################################
allele_info = pd.read_csv('data/allele_data.tsv', sep='\t')
RNAseq_info = pd.read_csv('data/RNAseq_data.tsv', sep='\t')
sample_info = pd.read_csv('data/sample_data.tsv', sep='\t')
result_info = pd.read_csv('data/result_data.tsv', sep='\t')
result_info = result_info[result_info['efSTR'] == "YES"]
result_info = result_info[result_info['location'] == args.location]
if result_info.empty:
    print(f"No data found for location: {args.location}")

####################################################################################################################
## figure 1a
####################################################################################################################
for i,STR in result_info.iterrows():
    location = STR['location'].replace(":","_").replace("-","_") #
    print(f"making cluster plot for STR at {location}")
    
    alleles = allele_info[allele_info['location'] == STR['location']]
    cluster_list = list(set(alleles['cluster'].to_list()))
    for cluster in cluster_list:
        clusterDataFrame = alleles[(alleles['cluster'] == cluster)]
        for i, structure in clusterDataFrame.iterrows():
            ss_image(structure['sequence'], structure['structure'], f"cluster_{cluster}_structure_{i}", structure['exemplar'])
        concatenate_images_vertically([f"cluster_{cluster}_structure_{i}.png" for i,_ in clusterDataFrame.iterrows()],f"cluster_{cluster}.png")
    concatenate_images_horizontally([f"cluster_{cluster}.png" for cluster in cluster_list], f"./figure_1a.png")
    for file in glob.glob("cluster*.p*"):
        os.remove(file)
        
####################################################################################################################
## figure 1b
####################################################################################################################
for i,STR in result_info.iterrows():
    location = STR['location'].replace(":","_").replace("-","_")
    print(f"making expression plot for STR at {location}")
    
    alleles = allele_info[allele_info['location'] == STR['location']]
    cluster_dict = dict(zip(alleles['allele'], alleles['cluster'])) 
    
    # map cluster labels to sample genotypes
    sampleAlleles = sample_info[sample_info['location'] == STR['location']]
    sampleAlleles['firstAlleleCluster'] = sampleAlleles['first_allele'].map(cluster_dict)
    sampleAlleles['secondAlleleCluster'] = sampleAlleles['second_allele'].map(cluster_dict)
    
    # attach rpkm data using gene name
    expression = RNAseq_info[RNAseq_info['gene'] == STR['gene']]   
    sampleExpression = pd.merge(sampleAlleles, expression, on='Sample', how='inner')
    sampleExpression['combinedAlleleCluster'] = sampleExpression['firstAlleleCluster'].astype(str) + "|" + sampleExpression['secondAlleleCluster'].astype(str)
    sampleExpression['combinedAlleleCluster'] = pd.Categorical(sampleExpression['combinedAlleleCluster'])
    
    # remove outliers
    Q1 = sampleExpression['rpkm'].quantile(0.05)
    Q3 = sampleExpression['rpkm'].quantile(0.95)
    IQR = Q3 - Q1
    filteredSampleExpression = sampleExpression[(sampleExpression['rpkm'] >= (Q1 - 1.5 * IQR)) & (sampleExpression['rpkm'] <= (Q3 + 1.5 * IQR))]
    filteredSampleExpression['combinedAlleleCluster'] = filteredSampleExpression['combinedAlleleCluster'].cat.remove_unused_categories()
    value_counts = filteredSampleExpression['combinedAlleleCluster'].value_counts()
    filteredSampleExpression = filteredSampleExpression.groupby('combinedAlleleCluster').filter(lambda x: len(x) > 1)
    filteredSampleExpression['combinedAlleleCluster'] = filteredSampleExpression['combinedAlleleCluster'].cat.remove_unused_categories()
    
    # statistical association with rpkm?
    anova_result = f_oneway(*[filteredSampleExpression[filteredSampleExpression['combinedAlleleCluster'] == group]['rpkm'] for group in filteredSampleExpression['combinedAlleleCluster'].unique()])
    group_mapping = {group: idx for idx, group in enumerate(filteredSampleExpression['combinedAlleleCluster'].cat.categories)}
    tukey = pairwise_tukeyhsd(endog=filteredSampleExpression['rpkm'], groups=filteredSampleExpression['combinedAlleleCluster'], alpha=0.05)
    significant_pairs = tukey.summary().data[1:]  # Tukey results excluding header
    significant_pairs = [row for row in significant_pairs if row[3] < 0.05]  # p-value < 0.05
    
    plot = ( # creating boxplot of rpkm values vs secondary structure assignments
        ggplot(filteredSampleExpression, aes(x='combinedAlleleCluster', y='rpkm', group='combinedAlleleCluster', color='Superpopulation'))
        + geom_jitter(aes(fill='Superpopulation'), alpha=0.5, width=0.35, size=5, show_legend=True, color="black")  # Fixed fill color for jitter points
        + geom_boxplot(show_legend=False, color="black", alpha=0.7, outlier_shape='', fill="#ffffff")  # Fill color mapped to combinedAlleleCluster for boxplot
        + scale_color_manual(values=['#ee6677', '#ccbb44', '#aa3377', '#4477aa', '#66ccee'])  # Custom color scale for 'Superpopulation'
        + scale_fill_manual(values=['#EE6677','#4477AA'])  # Custom fill colors for jitterplot
        + theme_bw()  # White background theme
        + theme(
            text=element_text(size=32),         # Increase all text
            axis_title=element_text(size=34),   # Adjust axis titles
            axis_text=element_text(size=30),    # Adjust axis tick labels
            legend_title=element_text(size=32), # Adjust legend title
            legend_text=element_text(size=30),  # Adjust legend items
            plot_title=element_text(size=36)    # Adjust plot title
        )
        + ggtitle(f"{STR['gene']}, {STR['location']}")
        + labs(y="rpkm", x="structure assignments", color="Ethnicity", fill="Ethnicity")  # Axis labels
    )

    # Annotate significant pairs
    y_max = filteredSampleExpression['rpkm'].max()
    y_step = (y_max * 0.1)  # Space between annotation lines
    for i, (group1, group2, meandiff, p_adj, lower, upper, reject) in enumerate(significant_pairs):
        if group1 in group_mapping and group2 in group_mapping:
            x1 = group_mapping[group1] + 1  # +1 because ggplot index starts at 1
            x2 = group_mapping[group2] + 1  # +1 because ggplot index starts at 1
            y = y_max + (i + 1) * y_step  # Incremental y-coordinate for each pair
            if ( p_adj < .01 ):
                plot += annotate('segment', x=x1, xend=x2, y=y, yend=y, color='black')  # Horizontal line
                plot += annotate('text', x=(x1 + x2) / 2, y=y + y_step * 0.1, label='***', size=22, ha='center')  # Significance mark
            elif ( p_adj < .05 ):
                plot += annotate('segment', x=x1, xend=x2, y=y, yend=y, color='black')  # Horizontal line
                plot += annotate('text', x=(x1 + x2) / 2, y=y + y_step * 0.1, label='*', size=22, ha='center')  # Significance mark
            else:
                plot += annotate('segment', x=x1, xend=x2, y=y, yend=y, color='white')  # Horizontal line
            
    plot.save(f'./figure_1b.png', width=15, height=10, dpi=300)
