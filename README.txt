README
------
Unsupervised discovery of microbial population discovery within metagenomes using nucleotide base composition

Isaam Saeed
05 Aug 2011
Revised: 14 Mar 2012

To run the two-tiered binning procedure (perl scripts):

Please ensure the following perl packages are installed on your system (in addition to BioPerl):

Term::ProgressBar;
Data::Dumper;
Bio::Tools::SeqWords;
Statistics::Descriptive;

1. Ensure that you change the directory and filenames in "generate_Features.sh" and "Two-Tiered_Binning.R" in accordance with your data set.

	- The test data set "simBG" is included in the "sample_data/" directory and the current parameter settings have been set in accordance with the data set

2. Run the script "generate_Features.sh"
	- This will generate the OFDEG/GC and tetranucleotide frequency vectors for clustering
	
3. Run the R script "Two-Tiered_Binning.R" (either interactively or directly from the terminal) 
	- Note that the following R-packages are required: Matrix, mclust, covRobust
	- Please ensure that the parameters are set accordingly
