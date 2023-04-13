#!/usr/bin/env Rscript

library(readxl)
library(tidyverse)
library(mikropml)
library(glue)
library(dplyr)

seed<-as.numeric(commandArgs(trailingOnly=TRUE))

#READIN IN CLINICAL DATA##############
clinical_metadata<-read_excel(path="raw_data/clinical_metadata.xlsx") 

age_quantile<-quantile(clinical_metadata$age, probs = seq(0, 1, 1/5), na.rm = TRUE)
elixhauser_quantile<-quantile(clinical_metadata$elixhauser_weighted, probs = seq(0, 1, 1/5), na.rm = TRUE)
hgb_quantile<-quantile(clinical_metadata$baseline_hgb_g_per_dl, probs = seq(0, 1, 1/5), na.rm = TRUE)
creat_quantile<-quantile(clinical_metadata$baseline_creat_mg_per_dl, probs = seq(0, 1, 1/5), na.rm = TRUE)
alb_quantile<-quantile(clinical_metadata$baseline_alb_g_per_dl, probs = seq(0, 1, 1/5), na.rm = TRUE)
prot_quantile<-quantile(clinical_metadata$baseline_prot_g_per_dl, probs = seq(0, 1, 1/5), na.rm = TRUE)

final_clinical_metadata<-clinical_metadata %>%
  mutate(age = case_when(clinical_metadata$age > age_quantile[5] ~ "q5",
                         clinical_metadata$age >= age_quantile[4] ~ "q4",
                         clinical_metadata$age >= age_quantile[3] ~ "q3",
                         clinical_metadata$age >= age_quantile[2] ~ "q2",
                         clinical_metadata$age >= age_quantile[1] ~ "q1",
                         is.na(clinical_metadata$age) ~ "missing",
                         TRUE ~ as.character(clinical_metadata$age))) %>%
  mutate(elixhauser_weighted = case_when(clinical_metadata$elixhauser_weighted > elixhauser_quantile[5] ~ "q5",
                                         clinical_metadata$elixhauser_weighted >= elixhauser_quantile[4] ~ "q4",
                                         clinical_metadata$elixhauser_weighted >= elixhauser_quantile[3] ~ "q3",
                                         clinical_metadata$elixhauser_weighted >= elixhauser_quantile[2] ~ "q2",
                                         clinical_metadata$elixhauser_weighted >= elixhauser_quantile[1] ~ "q1",
                                         is.na(clinical_metadata$elixhauser_weighted) ~ "missing",
                                         TRUE ~ as.character(clinical_metadata$elixhauser_weighted))) %>%
  mutate(baseline_hgb_g_per_dl = case_when(clinical_metadata$baseline_hgb_g_per_dl > hgb_quantile[5] ~ "q5",
                                           clinical_metadata$baseline_hgb_g_per_dl >= hgb_quantile[4] ~ "q4",
                                           clinical_metadata$baseline_hgb_g_per_dl >= hgb_quantile[3] ~ "q3",
                                           clinical_metadata$baseline_hgb_g_per_dl >= hgb_quantile[2] ~ "q2",
                                           clinical_metadata$baseline_hgb_g_per_dl >= hgb_quantile[1] ~ "q1",
                                           is.na(clinical_metadata$baseline_hgb_g_per_dl) ~ "missing",
                                           TRUE ~ as.character(clinical_metadata$baseline_hgb_g_per_dl))) %>%
  mutate(baseline_creat_mg_per_dl = case_when(clinical_metadata$baseline_creat_mg_per_dl > creat_quantile[5] ~ "q5",
                                              clinical_metadata$baseline_creat_mg_per_dl >= creat_quantile[4] ~ "q4",
                                              clinical_metadata$baseline_creat_mg_per_dl >= creat_quantile[3] ~ "q3",
                                              clinical_metadata$baseline_creat_mg_per_dl >= creat_quantile[2] ~ "q2",
                                              clinical_metadata$baseline_creat_mg_per_dl >= creat_quantile[1] ~ "q1",
                                              is.na(clinical_metadata$baseline_creat_mg_per_dl) ~ "missing",
                                              TRUE ~ as.character(clinical_metadata$baseline_creat_mg_per_dl))) %>%
  mutate(baseline_alb_g_per_dl = case_when(clinical_metadata$baseline_alb_g_per_dl > alb_quantile[5] ~ "q5",
                                           clinical_metadata$baseline_alb_g_per_dl >= alb_quantile[4] ~ "q4",
                                           clinical_metadata$baseline_alb_g_per_dl >= alb_quantile[3] ~ "q3",
                                           clinical_metadata$baseline_alb_g_per_dl >= alb_quantile[2] ~ "q2",
                                           clinical_metadata$baseline_alb_g_per_dl >= alb_quantile[1] ~ "q1",
                                           is.na(clinical_metadata$baseline_alb_g_per_dl) ~ "missing",
                                           TRUE ~ as.character(clinical_metadata$baseline_alb_g_per_dl))) %>%
  mutate(baseline_prot_g_per_dl = case_when(clinical_metadata$baseline_prot_g_per_dl > prot_quantile[5] ~ "q5",
                                            clinical_metadata$baseline_prot_g_per_dl >= prot_quantile[4] ~ "q4",
                                            clinical_metadata$baseline_prot_g_per_dl >= prot_quantile[3] ~ "q3",
                                            clinical_metadata$baseline_prot_g_per_dl >= prot_quantile[2] ~ "q2",
                                            clinical_metadata$baseline_prot_g_per_dl >= prot_quantile[1] ~ "q1",
                                            is.na(clinical_metadata$baseline_prot_g_per_dl) ~ "missing",
                                            TRUE ~ as.character(clinical_metadata$baseline_prot_g_per_dl))) %>%
  select(group, case_control) %>%  
  select(case_control, everything()) %>%
  arrange(group)

asv_counts<-read_tsv(file="raw_data/final_case_control.asv.ASV.subsample.shared") %>%
  select(-label, -numOtus) %>%
  rename(group = Group) %>%
  pivot_longer(-group, names_to="otu", values_to = "count")
complete_taxonomy<-read_tsv(file=glue("raw_data/final_case_control.asv.ASV.cons.taxonomy")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";") %>%
  mutate(genus = str_replace(genus,
                             "(.*)_unclassified", "Unclassified *\\1*"),
         genus = str_replace(genus,
                             "^(\\S*)$", "*\\1*")) %>%
  mutate(pretty_name = paste(.$otu,.$genus))
asv_rel_abun_taxon_level<-inner_join(asv_counts, complete_taxonomy, by = "otu") %>%
  group_by(group) %>%
  mutate(rel_abund = 100*(count / sum(count))) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu", "pretty_name"),
               names_to="level",
               values_to="taxon") 
asv_rel_abund<-asv_rel_abun_taxon_level %>%
  filter(level == "pretty_name") %>%
  group_by(group, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop")


asv_data<-inner_join(final_clinical_metadata, asv_rel_abund, by="group") %>% 
  pivot_wider(names_from = taxon, values_from = rel_abund) %>%
  select(-group) %>%
  select(case_control, everything())

#PROCESSING DATA####################
processed_data<-preprocess_data(asv_data, outcome_colname = "case_control")$dat_transformed

#WRITING MODEL FORMULA###############
hp<-list(alpha = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), 
               lambda = c(1e-4, 2.5e-4, 5e-4, 7.5e-4,
                          1e-3, 2.5e-3, 5e-3, 7.5e-3,
                          1e-2, 2.5e-2, 5e-2, 7.5e-2,
                          1e-1, 2.5e-1, 5e-1, 7.5e-1,
                          1, 2.5, 5, 7.5, 10))

results<-run_ml(processed_data, 
                method = "glmnet",
                outcome_colname = "case_control",
                hyperparameters = hp,
                find_feature_importance = TRUE,
                seed = seed)

#WRITING RESULTS##########
saveRDS(results, file=glue("processed_data/asv_en_results_{seed}.Rds"))

