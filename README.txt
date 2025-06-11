This is supplementary material for the paper:

SHORT TANDEM REPEAT VARIANTS ARE POSSIBLY ASSOCIATED 
WITH RNA SECONDARY STRUCTURE AND GENE EXPRESSION

Publication is open access at Plos One
contact the corresponding author for any questions

WHAT'S IN THIS DIRECTORY?
.
├── README.txt				# limited instructions for running our pipeline
├── envs
│   └── environment_linux.yaml		# anaconda environment
├── data
│   ├── make_datasets.py		# generates allele_data.tsv, matrix_data.tsv, and result_data.tsv
│   ├── make_datasets_stdout.txt	# audit trail for the make_datasets.py script
│   ├── region_data.tsv			# list of transcribed STRs investigated in this study
│   ├── sample_data.tsv			# sample genotypes for transcribed STRs
│   ├── RNAseq_data.tsv			# RPKM values each sample with gene names
│   ├── allele_data.tsv			# allele structures predicted by ViennaRNA
│   ├── matrix_data.tsv			# similarity matrix used for clustering and 2° structure assignments
│   ├── result_data.tsv			# summary of final results
│   └── plots			# folder for plots produced by make_datasets.py
├── figure_1.py				# makes a cluster plot and barplot for any region in the dataset
├── figure_2.py				# characterization of fSTRs
├── figure_3.py				# accessibility plots
├── figure_4.py				# secondary structure plots
├── figure_5.py				# accessibility based on simulated sequences
├── figure_S1.py			# extended characterization of fSTRs
├── figure_S2.py			# extended characterization of fSTRs
├── figure_S3.py			# extended characterization of fSTRs
├── figure_R1.py			# plot of accessibility correlation coefficients for STR motifs
├── figure_R2.py			# characterization of within and between cluster allele lengths
└── figures
    ├── figure_1.png
    ├── figure_2.png
    ├── figure_3.png
    ├── figure_4.png
    ├── figure_5.png
    ├── figure_6.png
    ├── figure_S1.png
    ├── figure_S2.png
    └── figure_S3.png

WHAT'S NOT IN THIS DIRECTORY OR THE DATA DIRECTORY:
1. bpRNA.pl
2. bpRNA_align.py
3. bpRNA_align_module.py
4. bpRNA_align_module.pyc

THE bpRNA_align SCRIPTS ARE AVAILABLE AT: https://github.com/BLasher113/bpRNA_align/tree/main
   Two scripts need access to these scripts 
   (make_datasets.py and figure_1.py). I did 
   not want to  redistribute these codes, users 
   interested reproducing our results need to
   copy and paste these scripts to this direc-
   tory and the data directory. If you have any
   problems please contact the corresponding
   author.

INSTALLATION (unfortunately not one size fits all, I'm using windows subsystem for linux
              with conda 24.9.2. Contact the corresponding author with issues and I will
			  try to help)

	CREATE THE ANACONDA ENVIRONMENT
		conda env create -f ./envs/environment_linux.yaml

	ACTIVATE THE ENVIRONMENT
		conda activate fSTRs
	
	ADDITIONAL THINGS YOU MIGHT NEED TO TRY TO GET THE PIPELINE WORKING:
		conda install -c bioconda viennarna
		conda install -c conda-forge gcc
		conda install -c conda-forge libxcrypt

	PERL MODULES THAT ARE USED
		cpan install Text::Levenshtein
		cpan install Set::Object
		cpan install Graph
		
	YOU NEED TO DECOMPRESS SIX TSV FILES IN THE DATA DIRECTORY
		7z e RNAseq_data.7z
		7z e allele_data.7z
		7z e matrix_data.7z
		7z e region_data.7z
		7z e result_data.7z
		7z e sample_data.7z
		
NOTES
1. This environment may have more packages than are actually needed. 
2. In particular, pytorch is not needed but could help extend the work.
3. All the scripts run on the command line. For example: python figure_1.py
4. figure_1.py takes an optional argument: python figure_1.py chr11:7672623-7672642
5. No need to run make_datasets.py. Only if you want to remake the datasets
6. figure_R1.py outputs figure_R1.png which only appears in the response letter
7. figure_R2.png does not appear anywhere but the data is cited in the response letter
8. There is no script for figure_6.png. I made that at https://www.lucidchart.com/
9. If anything seems broken please contact the corresponding author.
10.These scripts are somewhat commented but the use of pandas can make things hard to read.
11.Contact the authors with any questions.
