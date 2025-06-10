import pandas as pd
import RNA
from PIL import Image, ImageDraw
import numpy as np
from sklearn.cluster import AffinityPropagation
import os
from plotnine import ggplot, aes, geom_boxplot, element_text, annotate, geom_jitter, theme_bw, scale_color_manual, labs, ggtitle, theme, scale_fill_manual, geom_text, geom_segment
from scipy.stats import f_oneway, kruskal
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import time
import glob
import warnings
warnings.filterwarnings("ignore")

def transcribe(dna):
    bases = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'U'}
    return ''.join([bases[base] for base in dna])

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def dataframe_to_symmetric_matrix(data):
    structures = sorted(set(data['name_1']).union(set(data['name_2'])))
    n = len(structures)
    structure_to_index = {name: i for i, name in enumerate(structures)}
    matrix = np.zeros((n, n))
    for _, row in data.iterrows():
        i = structure_to_index[row['name_1']]
        j = structure_to_index[row['name_2']]
        score = row['score']
        matrix[i, j] = score
        matrix[j, i] = score  # Fill symmetric entry
    return matrix

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

#########################################################################################################
## main programming section
######################################################################################################### 

# read the RNAseq data
RNAseq_info = pd.read_csv('RNAseq_data.tsv', sep='\t')

# read the STRs
region_info = pd.read_csv('region_data.tsv', sep='\t')

# uncomment the following three lines if you just want efSTRs
# fSTR_info = pd.read_csv('result_data.tsv', sep='\t')
# fSTR_loci = fSTR_info[fSTR_info['efSTR']=='YES']['location'].to_list()
# region_info = region_info[region_info['location'].isin(fSTR_loci)]

# read the samples
sample_info = pd.read_csv('sample_data.tsv', sep='\t')

# this will save the ViennaRNA results
allele_info = { 'location':[], 'allele':[], 'sequence':[], 'structure':[], 'cluster':[], 'exemplar':[] }

# this will save the bpRNA-align results
matrix_info = { 'location':[], 'first_allele':[], 'second_allele':[], 'bpRNA_align_score':[], 'first_alignment':[], 'second_alignment':[] }

# this will save the results of this paper
result_info = { 'location':[], 'unit':[], 'gene':[], 'feature':[], 'fSTR':[], 'eSTR':[] }

#########################################################################################################
## main analysis loop for each STR
######################################################################################################### 
    
for i,STR in region_info.iterrows(): # loop through all the region_info
    print(f"\nworking on STR {i+1} out of {len(region_info)} ({STR['location']})")

    # get the samples that have genotype for this STR
    samples = sample_info[sample_info['location'] == STR['location']]
    if ( len(samples) == 0 ):
        print("\tno samples with STR genotypes available, skipping")
        continue
        
    # get the variants in the collection of samples
    allele_lengths = sorted([ int(x) for x in set(samples['first_allele'].to_list() + samples['second_allele'].to_list()) ])
    if ( len(allele_lengths) == 1 ):
        print("\tno variation in STR, skipping")
        continue
    
    #########################################################################################################
    ## fold variants with ViennaRNA, then compare secondary structure with bpRNA-align
    ######################################################################################################### 
    
    for variant in allele_lengths:
        print(f"\tfolding {variant}bp variant")
        DNA_sequence = ( STR['leftFlanking'] + (STR['unit'] * 100)[:variant] + STR['rightFlanking'] )
        RNA_sequence = transcribe(DNA_sequence) if STR['strand'] == "pos" else transcribe(reverse_complement(DNA_sequence))
        
        # predict the variant structure with ViennaRNA
        fc = RNA.fold_compound(RNA_sequence)
        mfe_structure, mfe = fc.mfe()

        # append entry to the allele_info
        allele_info['location'].append(STR['location'])
        allele_info['allele'].append(variant)
        allele_info['sequence'].append(RNA_sequence)
        allele_info['structure'].append(mfe_structure)
        
        # Writing structures to dbn files
        with open(f"structure_{variant}.dbn", "w") as f:
            f.write(RNA_sequence + "\n")
            f.write(mfe_structure + "\n")
            
    # Run bpRNA_align on the .dbn files
    with open("input_files.txt", "w") as f:
        for variant in allele_lengths:
            f.write(f"structure_{variant}.dbn\n")
    os.system("python3 bpRNA_align.py -a True -f input_files.txt -w 5 -o similarity_matrix.txt")    
    bpRNA_result = pd.read_csv("similarity_matrix.txt", delim_whitespace=True)
    similarity_matrix = dataframe_to_symmetric_matrix(bpRNA_result)
    
    # save bpRNA-align results
    for i,comparison in bpRNA_result.iterrows():
        matrix_info['location'].append(STR['location'])
        matrix_info['first_allele'].append(int(comparison['name_1'].replace('structure_','')))
        matrix_info['second_allele'].append(int(comparison['name_2'].replace('structure_','')))
        matrix_info['bpRNA_align_score'].append(comparison['score'])
        matrix_info['first_alignment'].append(comparison['alignment_1'])
        matrix_info['second_alignment'].append(comparison['alignment_2'])    
    
    # check if the STR alters secondary structure (fSTR)
    if ( np.max(similarity_matrix) - np.min(similarity_matrix) < 100 ):
        result_info['location'].append(STR['location'])
        result_info['unit'].append(STR['unit'])
        result_info['gene'].append(STR['gene'])
        result_info['feature'].append(STR['feature'])
        result_info['fSTR'].append('NO')
        result_info['eSTR'].append('NA')
        for variant in allele_lengths: # we are not going to try to cluster the RNA structures
            allele_info['cluster'].append(0)
            allele_info['exemplar'].append(0)
        print("\tSTR does not affect RNA secondary structure, moving on")
        continue
        
    # if there is a difference run affinity propagation clustering
    aff_prop = AffinityPropagation(affinity='precomputed', damping=.5)
    aff_prop.fit(similarity_matrix)
    labels = [ x+1 for x in aff_prop.labels_ ]
    exemplars = aff_prop.cluster_centers_indices_
    allele_info['cluster'].extend(labels)
    allele_info['exemplar'].extend([x if i in exemplars else 0 for i,x in enumerate(labels)])
    cluster_dict = { variant:labels[i] for i,variant in enumerate(allele_lengths) }
    cluster_list = sorted(labels)
    if ( max(cluster_list) == 0 ): # only one cluster   ---> important update on 4/27/25 <--- 
        result_info['location'].append(STR['location'])
        result_info['unit'].append(STR['unit'])
        result_info['gene'].append(STR['gene'])
        result_info['feature'].append(STR['feature'])
        result_info['fSTR'].append('NO')
        result_info['eSTR'].append('NA')
        print("\tSTR does not affect RNA secondary structure, moving on")
        continue
    else:
        print(f"\tfound {len(cluster_list)} RNA clusters (fSTR found)")
    
    # Retrieve sample genotypes for the STR
    expression = RNAseq_info[RNAseq_info['gene'] == STR['gene']]
    if ( expression.empty ): # no samples with rpkm data available
        result_info['location'].append(STR['location'])
        result_info['unit'].append(STR['unit'])
        result_info['gene'].append(STR['gene'])
        result_info['feature'].append(STR['feature'])
        result_info['fSTR'].append('YES')
        result_info['eSTR'].append('NA')
        print("\tNo samples with RPKM values found, moving on")
        continue
    else:
        print(f"\tfound {expression.shape[0]} samples rpkm value for {STR['gene']}")

    # map cluster labels to sample genotypes and attach rpkm data using gene name
    sampleAlleles = sample_info[sample_info['location'] == STR['location']]
    sampleAlleles['firstAlleleCluster'] = sampleAlleles['first_allele'].map(cluster_dict)
    sampleAlleles['secondAlleleCluster'] = sampleAlleles['second_allele'].map(cluster_dict)
    sampleExpression = pd.merge(sampleAlleles, expression, on='Sample', how='inner')
    sampleExpression['combinedAlleleCluster'] = sampleExpression['firstAlleleCluster'].astype(str) + "|" + sampleExpression['secondAlleleCluster'].astype(str)
    sampleExpression['combinedAlleleCluster'] = pd.Categorical(sampleExpression['combinedAlleleCluster'])
    
    # remove outliers
    Q1 = sampleExpression['rpkm'].quantile(0.05)
    Q3 = sampleExpression['rpkm'].quantile(0.95)
    IQR = Q3 - Q1
    filteredSampleExpression = sampleExpression[(sampleExpression['rpkm'] >= (Q1 - 1.5 * IQR)) & (sampleExpression['rpkm'] <= (Q3 + 1.5 * IQR))]

    # remove singleton categories
    filteredSampleExpression = filteredSampleExpression.groupby('combinedAlleleCluster').filter(lambda x: len(x) > 1)
    if ( filteredSampleExpression.empty ):                                   # no samples, so we cant
        result_info['location'].append(STR['location'])                      # check for association with
        result_info['unit'].append(STR['unit'])                              # gene expression
        result_info['gene'].append(STR['gene'])                     
        result_info['feature'].append(STR['feature'])
        result_info['fSTR'].append('YES')
        result_info['eSTR'].append('NA')
        print("\tNot enough samples in any clustered genotype category, moving on")
        continue

    # do we have multiple categories for association testing?
    if ( filteredSampleExpression['combinedAlleleCluster'].nunique() == 1 ): # only one combinedAlleleCluster
        result_info['location'].append(STR['location'])                      # so we cant check for association with
        result_info['unit'].append(STR['unit'])                              # gene expression
        result_info['gene'].append(STR['gene'])
        result_info['feature'].append(STR['feature'])
        result_info['fSTR'].append('YES')
        result_info['eSTR'].append('NA') 
        print("\tOnly one clustered genotype category, moving on")
        continue
    else:
        print(f"\tadequate samples available for association testing")
        
    filteredSampleExpression['combinedAlleleCluster'] = filteredSampleExpression['combinedAlleleCluster'].cat.remove_unused_categories()
    value_counts = filteredSampleExpression['combinedAlleleCluster'].value_counts()
    
    #########################################################################################################
    ## test for association with gene expression
    #########################################################################################################  
        
    # statistical association with rpkm?
    anova_result = f_oneway(*[filteredSampleExpression[filteredSampleExpression['combinedAlleleCluster'] == group]['rpkm'] for group in filteredSampleExpression['combinedAlleleCluster'].unique()])
    group_mapping = {group: idx for idx, group in enumerate(filteredSampleExpression['combinedAlleleCluster'].cat.categories)}
    tukey = pairwise_tukeyhsd(endog=filteredSampleExpression['rpkm'], groups=filteredSampleExpression['combinedAlleleCluster'], alpha=0.05)
    significant_pairs = tukey.summary().data[1:]  # Tukey results excluding header
    significant_pairs = [row for row in significant_pairs if row[3] < 0.05]  # p-value < 0.05

    # if we found an association make a figure and continue to the next STR
    if ( len(significant_pairs) > 0 ):
        print(f"\tmaking plot of rpkm association for {STR['location']} variants with {STR['gene']}")
        
        plot = ( # creating boxplot of rpkm values vs secondary structure assignments
            ggplot(filteredSampleExpression, aes(x='combinedAlleleCluster', y='rpkm', group='combinedAlleleCluster', color='Superpopulation'))
            + geom_jitter(aes(fill='Superpopulation'), alpha=0.5, width=0.35, size=5, show_legend=True, color="black")  # Fixed fill color for jitter points
            + geom_boxplot(show_legend=False, color="black", alpha=0.7, outlier_shape='', fill="#ffffff")  # Fill color mapped to combinedAlleleCluster for boxplot
            + scale_color_manual(values=['#ee6677', '#ccbb44', '#aa3377', '#4477aa', '#66ccee'])  # Custom color scale for 'Superpopulation'
            + scale_fill_manual(values=['#EE6677','#4477AA'])  # Custom fill colors for jitterplot
            + theme_bw()  # White background theme
            + theme(
                text=element_text(size=24),         # Increase all text
                axis_title=element_text(size=26),   # Adjust axis titles
                axis_text=element_text(size=22),    # Adjust axis tick labels
                legend_title=element_text(size=24), # Adjust legend title
                legend_text=element_text(size=22),  # Adjust legend items
                plot_title=element_text(size=28)    # Adjust plot title
            )
            + ggtitle(f"{STR['gene']}, {STR['location']}")  # Title with beta symbol , β={???}
            + labs(y="rpkm", x="structure assignments")  # Axis labels
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
                    plot += annotate('text', x=(x1 + x2) / 2, y=y + y_step * 0.1, label='***', size=14, ha='center')  # Significance mark
                elif ( p_adj < .05 ):
                    plot += annotate('segment', x=x1, xend=x2, y=y, yend=y, color='black')  # Horizontal line
                    plot += annotate('text', x=(x1 + x2) / 2, y=y + y_step * 0.1, label='*', size=14, ha='center')  # Significance mark
                else:
                    plot += annotate('segment', x=x1, xend=x2, y=y, yend=y, color='white')  # Horizontal line
                
        plot.save(f"./plots/boxplot_{STR['id']}.png", width=10, height=6, dpi=300)
                
        # make 2ndary structure images
        alleles = { key:value[-1 * len(allele_lengths):] for key,value in allele_info.items() }
        alleles = pd.DataFrame(alleles)
        for cluster in list(set(cluster_list)):
            clusterDataFrame = alleles[(alleles['cluster'] == cluster)]
            for i, structure in clusterDataFrame.iterrows():
                ss_image(structure['sequence'], structure['structure'], f"cluster_{cluster}_structure_{i}", structure['exemplar'])
            concatenate_images_vertically([f"cluster_{cluster}_structure_{i}.png" for i,_ in clusterDataFrame.iterrows()],f"cluster_{cluster}.png")
        concatenate_images_horizontally([f"cluster_{cluster}.png" for cluster in list(set(cluster_list))], f"./plots/cluster_{STR['id']}.png")
        
        # save the results
        result_info['location'].append(STR['location'])                 
        result_info['unit'].append(STR['unit'])                         
        result_info['gene'].append(STR['gene'])
        result_info['feature'].append(STR['feature'])
        result_info['fSTR'].append('YES')
        result_info['eSTR'].append('YES') 
        print(f"\tStatistical different between clustered genotypes found! finished with {STR['location']}")
        
    else: # no statistical differnece between genotype categories
        result_info['location'].append(STR['location'])                 
        result_info['unit'].append(STR['unit'])                         
        result_info['gene'].append(STR['gene'])
        result_info['feature'].append(STR['feature'])
        result_info['fSTR'].append('YES')
        result_info['eSTR'].append('NO') 
        print("\tNo differnce found in gene expression between clustered genotypes, moving on")

    # clean up intermediate files
    os.remove("similarity_matrix.txt")
    os.remove("input_files.txt")
    for variant in allele_lengths:
        os.remove(f"structure_{variant}.dbn")
        if os.path.exists(f"structure_{variant}.st"):
            os.remove(f"structure_{variant}.st")
    for file in glob.glob("cluster*.p*"):
        os.remove(file)

#########################################################################################################
## save the results
#########################################################################################################  
# uncomment the following lines if you want to overwrite
# allele_info = pd.DataFrame(allele_info)
# allele_info.to_csv('allele_data.tsv', sep='\t', index=False)
# matrix_info = pd.DataFrame(matrix_info)
# matrix_info.to_csv('matrix_data.tsv', sep='\t', index=False)
# result_info = pd.DataFrame(result_info)
# result_info.to_csv('result_data.tsv', sep='\t', index=False)

