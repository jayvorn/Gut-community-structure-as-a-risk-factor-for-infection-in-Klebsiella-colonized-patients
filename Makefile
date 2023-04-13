#RULE STRUCTURE
#Rule
#target:prerequisites
#(tab)recipe

processed_data/asv_en_results_%.Rds	:	raw_data/final_case_control.asv.ASV.cons.taxonomy\
												raw_data/final_case_control.asv.ASV.subsample.shared\
												raw_data/clinical_metadata.xlsx
	./code/asv_en.R $*
	
SEEDS = $(shell seq 1 1 100)
EN_RDS = $(patsubst %,processed_data/asv_en_results_%.Rds,$(SEEDS))

processed_data/asv_en_%.tsv : code/combine_asv_en.R $(EN_RDS)
	$^
