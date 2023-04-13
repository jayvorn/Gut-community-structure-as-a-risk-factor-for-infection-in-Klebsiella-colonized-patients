# Gut-community-structure-as-a-risk-factor-for-infection-in-Klebsiella-colonized-patients
Repository for data an code associated with the project Gut community structure as a risk factor for infection in Klebsiella-colonized patients

Hi there! Thanks for taking a look at this. The intention of this repository is to store R files and other files associated with the above project.
The goal of this repository is to provide open-access to the analysis performed for the above project, such that it can be replicated by any end-user.
If you decide to run this analysis yourself, make sure to change all of the directories! 

The data processing files are:
mothur_processing_case_control_16S.batch (data processing with mothur)
mothur_processing_case_control_16S_ASV.batch (data processing with mothur)
asv_en.R (elastic net processing with mikropml using ASVs as input data)
combine_asv_en.R (assembly of elastic net data using ASVs as input data)
asv_en.batch (batch file for running asv_en.R)
combine_asv_en.batch (batch file for running combine_asv_en.R)
Makefile (makefile rule structure for making .Rds files from mikropml)

*Note* For machine learning models in mikropml, I am providing a single data input example; however, I am happy to provide any additional set of processing files.

The analysis files are:
case_control_ASV_analysis.R (analysis and visualization of mothur output)
case_control_otu_analysis.R (analysis and visualization of mothur output)
gut_community_analysis.R (analysis and visualization of mothur output)
ml_analysis.R (analysis and visualization mikropml output)

I have also included a fasta files with the sequences for ASV1 and ASV19
