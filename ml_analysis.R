#Load required packages
library(readxl)
library(ggplot2)
library(tidyverse)
library(glue)
library(ggtext)
library(RColorBrewer)
library(purrr)

#SET ENVIRONMENT########################
#Set working directory
setwd("/Volumes/Seagate Backup Plus Drive/mothur/Case_control/analysis/ml_models/final_analysis")
#Set seed for all graphing
set.seed(20220212)
#Set path for all files
path<-"/Volumes/Seagate Backup Plus Drive/mothur/Case_control/analysis/ml_models/final_analysis"
#Create graph output directory
dir.create(final_ml_graphs)

#TAXON PERFORMANCE##############
#read in elastic net performance 
ASV_en_performance<-read_tsv(glue("{path}/tax_levels_elastic_net/ASV/processed_data/taxon_data_ASV_performance.tsv")) %>%
  mutate(taxon = "ASV")
otu_en_performance<-read_tsv(glue("{path}/tax_levels_elastic_net/otu/processed_data/taxon_data_otu_performance.tsv")) %>%
  mutate(taxon = "OTU")
genus_en_performance<-read_tsv(glue("{path}/tax_levels_elastic_net/genus/processed_data/taxon_data_genus_performance.tsv")) %>%
  mutate(taxon = "Genus")
family_en_performance<-read_tsv(glue("{path}/tax_levels_elastic_net/family/processed_data/taxon_data_family_performance.tsv")) %>%
  mutate(taxon = "Family")
order_en_performance<-read_tsv(glue("{path}/tax_levels_elastic_net/order/processed_data/taxon_data_order_performance.tsv")) %>%
  mutate(taxon = "Order")
class_en_performance<-read_tsv(glue("{path}/tax_levels_elastic_net/class/processed_data/taxon_data_class_performance.tsv")) %>%
  mutate(taxon = "Class")
phylum_en_performance<-read_tsv(glue("{path}/tax_levels_elastic_net/phylum/processed_data/taxon_data_phylum_performance.tsv")) %>%
  mutate(taxon = "Phylum")

#read in random forest performance 
ASV_rf_performance<-read_tsv(glue("{path}/tax_levels_random_forest/ASV/processed_data/taxon_data_ASV_performance.tsv")) %>%
  mutate(taxon = "ASV")
otu_rf_performance<-read_tsv(glue("{path}/tax_levels_random_forest/otu/processed_data/taxon_data_otu_performance.tsv")) %>%
  mutate(taxon = "OTU")
genus_rf_performance<-read_tsv(glue("{path}/tax_levels_random_forest/genus/processed_data/taxon_data_genus_performance.tsv")) %>%
  mutate(taxon = "Genus")
family_rf_performance<-read_tsv(glue("{path}/tax_levels_random_forest/family/processed_data/taxon_data_family_performance.tsv")) %>%
  mutate(taxon = "Family")
order_rf_performance<-read_tsv(glue("{path}/tax_levels_random_forest/order/processed_data/taxon_data_order_performance.tsv")) %>%
  mutate(taxon = "Order")
class_rf_performance<-read_tsv(glue("{path}/tax_levels_random_forest/class/processed_data/taxon_data_class_performance.tsv")) %>%
  mutate(taxon = "Class")
phylum_rf_performance<-read_tsv(glue("{path}/tax_levels_random_forest/phylum/processed_data/taxon_data_phylum_performance.tsv")) %>%
  mutate(taxon = "Phylum")

#Combine en performance data
taxon_en_performance_data<-bind_rows(ASV_en_performance, otu_en_performance, genus_en_performance, family_en_performance, order_en_performance, class_en_performance, phylum_en_performance) %>%
  mutate(log_AUC = log10(AUC))

#Taxon en performance ANOVA
taxon_en_performance_aov<-TukeyHSD(aov(log_AUC ~ taxon, data = taxon_en_performance_data))

#Graph en performance data
taxon_en_performance_data %>%
  select(seed, taxon, AUC) %>%
  mutate(taxon = factor(taxon, levels=c("Phylum","Class","Order","Family","Genus","OTU","ASV"))) %>%
  ggplot(aes(y=AUC, x=taxon, color = taxon)) +
  geom_point(show.legend = FALSE, position = position_jitter(), alpha = 0.5, size = 0.75) +
  stat_summary(fun = median,
               geom = "pointrange",
               fun.max = function(y) median(y) + sd(y),
               fun.min = function(y) median(y) - sd(y),
               orientation = "x",
               color = "black", 
               size = 0.5, 
               inherit.aes = TRUE) +
  scale_y_continuous(limits = c(0.15, 1.05), breaks = c(0.3, 0.5, 0.7, 0.9)) +
  labs(x = NULL) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank(), 
        panel.grid.major =  element_line(color = "grey", size = 0.5, linetype = 2)) +
  annotate("segment", x = c(1,2,3,4,5,6), xend = c(7,7,7,7,7,7), y = c(0.2, 0.25, 0.3, 1, 0.95 ,0.9), yend = c(0.2, 0.25, 0.3, 1, 0.95 ,0.9), color = "black") +
  geom_richtext(data=tibble(x=4, y=0.17), fill = NA, label.color = NA, label="*p* = 5.1e-4" ,aes(x=x, y=y), inherit.aes=FALSE, size=2) +
  geom_richtext(data=tibble(x=4.5, y=0.22), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=2) +
  geom_richtext(data=tibble(x=5, y=0.27), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=2) +
  geom_richtext(data=tibble(x=5.5, y=1.02), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=2) +
  geom_richtext(data=tibble(x=6, y=0.97), fill = NA, label.color = NA, label="*p* = 3.5e-4" ,aes(x=x, y=y), inherit.aes=FALSE, size=2) +
  geom_richtext(data=tibble(x=6.5, y=0.92), fill = NA, label.color = NA, label="*p* = 1.3e-3" ,aes(x=x, y=y), inherit.aes=FALSE, size=2)
ggsave(glue("final_graphs/elastic_net_taxon_performance_data.pdf"), height = 8, width = 10, unit = "cm")

#Combine rf performance data
taxon_rf_performance_data<-bind_rows(ASV_rf_performance, otu_rf_performance, genus_rf_performance, family_rf_performance, order_rf_performance, class_rf_performance, phylum_rf_performance) %>%
  mutate(log_AUC = log10(AUC))

#Taxon rf performance ANOVA
taxon_rf_performance_aov<-TukeyHSD(aov(log_AUC ~ taxon, data = taxon_rf_performance_data))

#Graph rf performance data
taxon_rf_performance_data %>%
  select(seed, taxon, AUC) %>%
  mutate(taxon = factor(taxon, levels=c("Phylum","Class","Order","Family","Genus","OTU","ASV"))) %>%
  ggplot(aes(y=AUC, x=taxon, color = taxon)) +
  geom_point(show.legend = FALSE, position = position_jitter(), alpha = 0.5, size = 0.75) +
  stat_summary(fun = median,
               geom = "pointrange",
               fun.max = function(y) median(y) + sd(y),
               fun.min = function(y) median(y) - sd(y),
               orientation = "x",
               color = "black", 
               size = 0.5, 
               inherit.aes = TRUE) +
  scale_y_continuous(limits = c(0.15, 1.05), breaks = c(0.3, 0.5, 0.7, 0.9)) +
  labs(x = NULL) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank(), 
        panel.grid.major =  element_line(color = "grey", size = 0.5, linetype = 2)) +
  annotate("segment", x = c(1,2,3,4,5,6), xend = c(7,7,7,7,7,7), y = c(0.2, 0.25, 0.3, 1, 0.95 ,0.9), yend = c(0.2, 0.25, 0.3, 1, 0.95 ,0.9), color = "black") +
  geom_richtext(data=tibble(x=4, y=0.17), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=2) +
  geom_richtext(data=tibble(x=4.5, y=0.22), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=2) +
  geom_richtext(data=tibble(x=5, y=0.27), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=2) +
  geom_richtext(data=tibble(x=5.5, y=1.02), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=2) +
  geom_richtext(data=tibble(x=6, y=0.97), fill = NA, label.color = NA, label="*p* = 6.7e-3" ,aes(x=x, y=y), inherit.aes=FALSE, size=2) +
  geom_richtext(data=tibble(x=6.5, y=0.92), fill = NA, label.color = NA, label="*p* = 3.2e-3" ,aes(x=x, y=y), inherit.aes=FALSE, size=2)
ggsave(glue("final_graphs/random_forest_taxon_performance_data.pdf"), height = 8, width = 10, unit = "cm")

#ITERATED DATA PERFORMANCE#########
#Read in en performance data
clin_en_performance<-read_tsv(glue("{path}/elastic_net_performance/clinical_elastic_net/processed_data/clin_en_performance.tsv")) %>%
  mutate(data = "clinical")
geno_en_performance<-read_tsv(glue("{path}/elastic_net_performance/kp_genotype_elastic_net/processed_data/kp_genotype_en_performance.tsv")) %>%
  mutate(data = "genotype")
asv_en_performance<-read_tsv(glue("{path}/elastic_net_performance/asv_elastic_net/processed_data/asv_en_performance.tsv")) %>%
  mutate(data = "ASV")
clin_geno_en_performance<-read_tsv(glue("{path}/elastic_net_performance/clinical_genotype_elastic_net/processed_data/clin_geno_en_performance.tsv")) %>%
  mutate(data = "clinical_genotype")
asv_geno_en_performance<-read_tsv(glue("{path}/elastic_net_performance/asv_genotype_elastic_net/processed_data/asv_geno_en_performance.tsv")) %>%
  mutate(data = "ASV_genotype")
asv_clin_en_performance<-read_tsv(glue("{path}/elastic_net_performance/asv_clinical_elastic_net/processed_data/asv_clin_en_performance.tsv")) %>%
  mutate(data = "ASV_clinical")
all_en_performance<-read_tsv(glue("{path}/elastic_net_performance/all_data_elastic_net/processed_data/all_en_performance.tsv")) %>%
  mutate(data = "ASV_clinical_genotype")

#Read in rf performance data
clin_rf_performance<-read_tsv(glue("{path}/random_forest_performance/clinical_random_forest/processed_data/clin_rf_performance.tsv")) %>%
  mutate(data = "clinical")
geno_rf_performance<-read_tsv(glue("{path}/random_forest_performance/kp_genotype_random_forest/processed_data/kp_genotype_rf_performance.tsv")) %>%
  mutate(data = "genotype")
asv_rf_performance<-read_tsv(glue("{path}/random_forest_performance/asv_random_forest/processed_data/asv_rf_performance.tsv")) %>%
  mutate(data = "ASV")
clin_geno_rf_performance<-read_tsv(glue("{path}/random_forest_performance/clinical_genotype_random_forest/processed_data/clin_geno_rf_performance.tsv")) %>%
  mutate(data = "clinical_genotype")
asv_geno_rf_performance<-read_tsv(glue("{path}/random_forest_performance/asv_genotype_random_forest/processed_data/asv_geno_rf_performance.tsv")) %>%
  mutate(data = "ASV_genotype")
asv_clin_rf_performance<-read_tsv(glue("{path}/random_forest_performance/asv_clinical_random_forest/processed_data/asv_clin_rf_performance.tsv")) %>%
  mutate(data = "ASV_clinical")
all_rf_performance<-read_tsv(glue("{path}/random_forest_performance/all_data_random_forest/processed_data/all_rf_performance.tsv")) %>%
  mutate(data = "ASV_clinical_genotype")

#Combine en performance data
iter_en_performance_data<-bind_rows(clin_en_performance, geno_en_performance, asv_en_performance, clin_geno_en_performance, asv_geno_en_performance, asv_clin_en_performance, all_en_performance) %>%
  mutate(log_AUC = log10(AUC))

#Taxon en performance ANOVA
iter_en_performance_aov<-TukeyHSD(aov(log_AUC ~ data, data = iter_en_performance_data))

#Graph en performance data
iter_en_performance_data %>%
  select(seed, data, AUC) %>%
  filter(data != "clinical_genotype") %>%
  mutate(data = factor(data, levels=c("ASV","clinical","ASV_clinical","genotype","ASV_genotype","ASV_clinical_genotype"))) %>%
  ggplot(aes(y=AUC, x=data, color = data)) +
  geom_point(show.legend = FALSE, position = position_jitter(), alpha = 0.5, size = 0.75) +
  scale_x_discrete(breaks = c("ASV","clinical","ASV_clinical","genotype","ASV_genotype","ASV_clinical_genotype"),
                   labels = c("ASV","clinical","ASV+<br>clinical","genotype","ASV+<br>genotype","ASV+<br>clinical+<br>genotype")) +
  stat_summary(fun = median,
               geom = "pointrange",
               fun.max = function(y) median(y) + sd(y),
               fun.min = function(y) median(y) - sd(y),
               orientation = "x",
               color = "black", 
               size = 0.5, 
               inherit.aes = TRUE) +
  scale_y_continuous(limits = c(0.2, 1.0), breaks = c(0.3, 0.5, 0.7, 0.9)) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_blank(), 
        panel.grid.major =  element_line(color = "grey", size = 0.5, linetype = 2)) + 
  annotate("segment", x = c(2,4,4,5), xend = c(3,5,6,6), y = c(0.9,0.4, 0.3 ,0.9), yend = c(0.9,0.4, 0.3 ,0.9), color = "black") +
  geom_richtext(data=tibble(x=2.5, y=0.94), fill = NA, label.color = NA, label="*p* = 7e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=4.5, y=0.36), fill = NA, label.color = NA, label="*p* = 9.2e-4" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=5, y=0.26), fill = NA, label.color = NA, label="*p* = 1.4e-4" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=5.5, y=0.94), fill = NA, label.color = NA, label="*p* = 0.99" ,aes(x=x, y=y), inherit.aes=FALSE, size=3)
ggsave(glue("final_graphs/elastic_net_iterative_performance_data.pdf"), height = 8, width = 9, unit = "cm")

#Combine rf performance data
iter_rf_performance_data<-bind_rows(clin_rf_performance, geno_rf_performance, asv_rf_performance, clin_geno_rf_performance, asv_geno_rf_performance, asv_clin_rf_performance, all_rf_performance) %>%
  mutate(log_AUC = log10(AUC))

#Taxon en performance ANOVA
iter_rf_performance_aov<-TukeyHSD(aov(log_AUC ~ data, data = iter_rf_performance_data))

#Graph rf performance data
iter_rf_performance_data %>%
  select(seed, data, AUC) %>%
  filter(data != "clinical_genotype") %>%
  mutate(data = factor(data, levels=c("ASV","clinical","ASV_clinical","genotype","ASV_genotype","ASV_clinical_genotype"))) %>%
  ggplot(aes(y=AUC, x=data, color = data)) +
  geom_point(show.legend = FALSE, position = position_jitter(), alpha = 0.5, size = 0.75) +
  scale_x_discrete(breaks = c("ASV","clinical","ASV_clinical","genotype","ASV_genotype","ASV_clinical_genotype"),
                   labels = c("ASV","clinical","ASV+<br>clinical","genotype","ASV+<br>genotype","ASV+<br>clinical+<br>genotype")) +
  stat_summary(fun = median,
               geom = "pointrange",
               fun.max = function(y) median(y) + sd(y),
               fun.min = function(y) median(y) - sd(y),
               orientation = "x",
               color = "black", 
               size = 0.5, 
               inherit.aes = TRUE) +
  scale_y_continuous(limits = c(0.2, 1), breaks = c(0.3, 0.5, 0.7, 0.9)) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_blank(), 
        panel.grid.major =  element_line(color = "grey", size = 0.5, linetype = 2)) + 
  annotate("segment", x = c(2,4,4,5), xend = c(3,5,6,6), y = c(0.9,0.4, 0.3 ,0.45), yend = c(0.9,0.4, 0.3 ,0.45), color = "black") +
  geom_richtext(data=tibble(x=2.5, y=0.94), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=4.5, y=0.36), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=5, y=0.26), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=5.5, y=0.41), fill = NA, label.color = NA, label="*p* = 0.7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3)
ggsave(glue("final_graphs/random_forest_iterative_performance_data.pdf"), height = 8, width = 9, unit = "cm")

#FEATURE WEIGHTS########
#Write function to extract en asv feature weights
get_asv_weights<-function(file_name){
  model<-readRDS(file_name) %>%
    pluck("trained_model")
  coef(model$finalModel, model$bestTune$lambda) %>%
    as.matrix() %>%
    as.tibble(rownames = "feature") %>%
    rename(weight = s1) %>%
    mutate(seed = str_replace(file_name, "elastic_net_performance/asv_elastic_net/processed_data/asv_en_results_(\\d*).Rds", "\\1"))
}

#Create en file list
asv_files<-list.files(path = "elastic_net_performance/asv_elastic_net/processed_data",
                                pattern = "asv_en_results_(\\d*).Rds",
                                full.names = TRUE)

#Run en feature weight extraction function
asv_weights<-map_dfr(asv_files, get_asv_weights)

#Remove intercept variable and categorize case control status
asv_weight_summary<-asv_weights %>%
  group_by(feature) %>%
  summarize(mean_weight = mean(weight)) %>%
  filter(feature != "(Intercept)") %>%
  mutate(direction = case_when(mean_weight > 0 ~ "control",
                               mean_weight < 0 ~ "case")) %>%
  mutate(feature = str_remove(feature, "^`")) %>%
  mutate(feature = str_remove(feature, "`$"))

#Write function to extract en clin feature weights
get_clin_weights<-function(file_name){
  model<-readRDS(file_name) %>%
    pluck("trained_model")
  coef(model$finalModel, model$bestTune$lambda) %>%
    as.matrix() %>%
    as.tibble(rownames = "feature") %>%
    rename(weight = s1) %>%
    mutate(seed = str_replace(file_name, "elastic_net_performance/clinical_elastic_net/processed_data/clin_en_results_(\\d*).Rds", "\\1"))
}

#Create en file list
clin_files<-list.files(path = "elastic_net_performance/clinical_elastic_net/processed_data",
                      pattern = "clin_en_results_(\\d*).Rds",
                      full.names = TRUE)

#Run en feature weight extraction function
clin_weights<-map_dfr(clin_files, get_clin_weights)

#Remove intercept variable and categorize case control status
clin_weight_summary<-clin_weights %>%
  group_by(feature) %>%
  summarize(mean_weight = mean(weight)) %>%
  filter(feature != "(Intercept)") %>%
  mutate(direction = case_when(mean_weight > 0 ~ "control",
                               mean_weight < 0 ~ "case")) %>%
  mutate(feature = str_remove(feature, "^`")) %>%
  mutate(feature = str_remove(feature, "`$"))


#Write function to extract en asv_clin feature weights
get_asv_clin_weights<-function(file_name){
  model<-readRDS(file_name) %>%
    pluck("trained_model")
  coef(model$finalModel, model$bestTune$lambda) %>%
    as.matrix() %>%
    as.tibble(rownames = "feature") %>%
    rename(weight = s1) %>%
    mutate(seed = str_replace(file_name, "elastic_net_performance/asv_clinical_elastic_net/processed_data/asv_clin_en_results_(\\d*).Rds", "\\1"))
}

#Create en file list
asv_clin_files<-list.files(path = "elastic_net_performance/asv_clinical_elastic_net/processed_data",
                       pattern = "asv_clin_en_results_(\\d*).Rds",
                       full.names = TRUE)

#Run en feature weight extraction function
asv_clin_weights<-map_dfr(asv_clin_files, get_asv_clin_weights)

#Remove intercept variable and categorize case control status
asv_clin_weight_summary<-asv_clin_weights %>%
  group_by(feature) %>%
  summarize(mean_weight = mean(weight)) %>%
  filter(feature != "(Intercept)") %>%
  mutate(direction = case_when(mean_weight > 0 ~ "control",
                               mean_weight < 0 ~ "case")) %>%
  mutate(feature = str_remove(feature, "^`")) %>%
  mutate(feature = str_remove(feature, "`$"))

#ITERATED IN FEATURE IMPORTANCE DATA#########
#Read in en feature importance
clin_en_feature_importance<-read_tsv(glue("{path}/elastic_net_performance/clinical_elastic_net/processed_data/clin_en_feature_importance.tsv")) %>%
  mutate(data = "clinical")
geno_en_feature_importance<-read_tsv(glue("{path}/elastic_net_performance/kp_genotype_elastic_net/processed_data/kp_genotype_en_feature_importance.tsv")) %>%
  mutate(data = "genotype")
asv_en_feature_importance<-read_tsv(glue("{path}/elastic_net_performance/asv_elastic_net/processed_data/asv_en_feature_importance.tsv")) %>%
  mutate(data = "ASV")
clin_geno_en_feature_importance<-read_tsv(glue("{path}/elastic_net_performance/clinical_genotype_elastic_net/processed_data/clin_geno_en_feature_importance.tsv")) %>%
  mutate(data = "clinical_genotype")
asv_geno_en_feature_importance<-read_tsv(glue("{path}/elastic_net_performance/asv_genotype_elastic_net/processed_data/asv_geno_en_feature_importance.tsv")) %>%
  mutate(data = "ASV_genotype")
asv_clin_en_feature_importance<-read_tsv(glue("{path}/elastic_net_performance/asv_clinical_elastic_net/processed_data/asv_clin_en_feature_importance.tsv")) %>%
  mutate(data = "ASV_clinical")
all_en_feature_importance<-read_tsv(glue("{path}/elastic_net_performance/all_data_elastic_net/processed_data/all_en_feature_importance.tsv")) %>%
  mutate(data = "ASV_clinical_genotype")

#Read in rf feature importance
clin_rf_feature_importance<-read_tsv(glue("{path}/random_forest_performance/clinical_random_forest/processed_data/clin_rf_feature_importance.tsv")) %>%
  mutate(data = "clinical")
geno_rf_feature_importance<-read_tsv(glue("{path}/random_forest_performance/kp_genotype_random_forest/processed_data/kp_genotype_rf_feature_importance.tsv")) %>%
  mutate(data = "genotype")
asv_rf_feature_importance<-read_tsv(glue("{path}/random_forest_performance/asv_random_forest/processed_data/asv_rf_feature_importance.tsv")) %>%
  mutate(data = "ASV")
clin_geno_rf_feature_importance<-read_tsv(glue("{path}/random_forest_performance/clinical_genotype_random_forest/processed_data/clin_geno_rf_feature_importance.tsv")) %>%
  mutate(data = "clinical_genotype")
asv_geno_rf_feature_importance<-read_tsv(glue("{path}/random_forest_performance/asv_genotype_random_forest/processed_data/asv_geno_rf_feature_importance.tsv")) %>%
  mutate(data = "ASV_genotype")
asv_clin_rf_feature_importance<-read_tsv(glue("{path}/random_forest_performance/asv_clinical_random_forest/processed_data/asv_clin_rf_feature_importance.tsv")) %>%
  mutate(data = "ASV_clinical")
all_rf_feature_importance<-read_tsv(glue("{path}/random_forest_performance/all_data_random_forest/processed_data/all_rf_feature_importance.tsv")) %>%
  mutate(data = "ASV_clinical_genotype")

#Graph asv en feature importance 
asv_en_feature_importance %>%
  rename(feature = names) %>%
  group_by(feature) %>%
  summarize(mean = mean(perf_metric_diff * 1),
            l_quartile = quantile(perf_metric_diff * 1, prob = 0.25),
            u_quartile = quantile(perf_metric_diff * 1, prob = 0.75)) %>%
  filter(mean > 0.001) %>%
  inner_join(., asv_weight_summary, by = "feature") %>%
  mutate(feature = fct_reorder(feature, mean)) %>%
  ggplot(aes(x=mean, y=feature, xmin=l_quartile, xmax = u_quartile, color = direction)) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey") +
  geom_point(show.legend = FALSE, size = 3) +
  geom_linerange(show.legend = FALSE) +
  scale_color_manual(breaks = c("case", "control"),
                     values = c("red", "black")) +
  xlim(-0.002, 0.2) +
  labs(x = "Reduction in AUC when removed",
       subtitle = "Cutoff: AUC reduction > 0.001") +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
ggsave(glue("final_graphs/elastic_net_asv_feature_importance.pdf"), height = 8, width = 12, unit = "cm")

#Graph asv rf feature importance 
asv_rf_feature_importance %>%
  rename(feature = names) %>%
  group_by(feature) %>%
  summarize(mean = mean(perf_metric_diff * 1),
            l_quartile = quantile(perf_metric_diff * 1, prob = 0.25),
            u_quartile = quantile(perf_metric_diff * 1, prob = 0.75)) %>%
  filter(mean > 0.005) %>%
  mutate(feature = fct_reorder(feature, mean)) %>%
  ggplot(aes(x=mean, y=feature, xmin=l_quartile, xmax = u_quartile)) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey") +
  geom_point(show.legend = FALSE, size = 3) +
  geom_linerange(show.legend = FALSE) +
  xlim(-0.002, 0.2) +
  labs(x = "Reduction in AUC when removed",
       subtitle = "Cutoff: AUC reduction > 0.005") +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
ggsave(glue("final_graphs/random_forest_asv_feature_importance.pdf"), height = 8, width = 12, unit = "cm")

#Convert clinical variable names into user-friendly names
clin_pretty_names<-read_excel("{path}/clin_pretty_names.xlsx")
clin_weight_summary<-inner_join(clin_weight_summary, clin_pretty_names, by = c("feature" = "names")) %>%
                                  rename(names = feature) %>%
                                  rename(feature = feature.y)

#Graph top clin en feature importance 
clin_en_feature_importance %>%
  inner_join(., clin_pretty_names, by = "names") %>%
  group_by(feature) %>%
  summarize(mean = mean(perf_metric_diff * 1),
            l_quartile = quantile(perf_metric_diff * 1, prob = 0.25),
            u_quartile = quantile(perf_metric_diff * 1, prob = 0.75)) %>%
  filter(mean > 0.005) %>%
  inner_join(., clin_weight_summary, by = "feature") %>%
  mutate(feature = fct_reorder(feature, mean)) %>%
  ggplot(aes(x=mean, y=feature, xmin=l_quartile, xmax = u_quartile, color = direction)) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey") +
  geom_point(show.legend = FALSE, size = 3) +
  geom_linerange(show.legend = FALSE) +
  scale_color_manual(breaks = c("case", "control"),
                     values = c("red", "black")) +
  xlim(-0.002, 0.2) +
  labs(x = "Reduction in AUC when removed",
       subtitle = "Cutoff: AUC reduction > 0.005") +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
ggsave(glue("final_graphs/elastic_net_clin_feature_importance.pdf"), height = 8, width = 12, unit = "cm")

#Graph abx only clin en feature importance 
clin_en_feature_importance %>%
  inner_join(., clin_pretty_names, by = "names") %>%
  group_by(feature) %>%
  summarize(mean = mean(perf_metric_diff * 1),
            l_quartile = quantile(perf_metric_diff * 1, prob = 0.25),
            u_quartile = quantile(perf_metric_diff * 1, prob = 0.75)) %>%
  inner_join(., clin_weight_summary, by = "feature") %>%
  filter(feature == "prior aminoglycoside" | 
           feature == "prior combo betalactam" |
           feature == "prior carbapenem" |
           feature == "prior 3rd or 4th<br>generation cephalosporin" |
           feature == "prior cephalosporin" |
           feature == "high risk antibiotics" |
           feature == "prior penicillin" |
           feature == "prior penicillin") %>%
  mutate(feature = fct_reorder(feature, mean)) %>%
  ggplot(aes(x=mean, y=feature, xmin=l_quartile, xmax = u_quartile, color = direction)) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey") +
  geom_point(show.legend = FALSE, size = 3) +
  geom_linerange(show.legend = FALSE) +
  scale_color_manual(breaks = c("case", "control"),
                     values = c("red", "black")) +
  xlim(-0.01, 0.2) +
  labs(x = "Reduction in AUC when removed") +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
ggsave(glue("final_graphs/elastic_net_clin_abx_feature_importance.pdf"), height = 8, width = 12, unit = "cm")

#Convert clinical variable names into user-friendly names
clin_asv_pretty_names<-read_excel("{path}/clin_asv_pretty_names.xlsx")
asv_clin_weight_summary<-inner_join(asv_clin_weight_summary, clin_asv_pretty_names, by = c("feature" = "names")) %>%
  rename(names = feature) %>%
  rename(feature = feature.y)

#Graph top asv_clin en feature importance 
asv_clin_en_feature_importance %>%
  inner_join(., clin_asv_pretty_names, by = "names") %>%
  group_by(feature) %>%
  summarize(mean = mean(perf_metric_diff * 1),
            l_quartile = quantile(perf_metric_diff * 1, prob = 0.25),
            u_quartile = quantile(perf_metric_diff * 1, prob = 0.75)) %>%
  filter(mean > 0.001) %>%
  inner_join(., asv_clin_weight_summary, by = "feature") %>%
  mutate(feature = fct_reorder(feature, mean)) %>%
  ggplot(aes(x=mean, y=feature, xmin=l_quartile, xmax = u_quartile, color = direction)) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey") +
  geom_point(show.legend = FALSE, size = 3) +
  geom_linerange(show.legend = FALSE) +
  scale_color_manual(breaks = c("case", "control"),
                     values = c("red", "black")) +
  xlim(-0.002, 0.2) +
  labs(x = "Reduction in AUC when removed",
       subtitle = "Cutoff: AUC reduction > 0.001") +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
ggsave(glue("final_graphs/elastic_net_clin_asv_feature_importance.pdf"), height = 10, width = 12, unit = "cm")

#Graph abx only asv_clin en feature importance 
asv_clin_en_feature_importance %>%
  inner_join(., clin_asv_pretty_names, by = "names") %>%
  group_by(feature) %>%
  summarize(mean = mean(perf_metric_diff * 1),
            l_quartile = quantile(perf_metric_diff * 1, prob = 0.25),
            u_quartile = quantile(perf_metric_diff * 1, prob = 0.75)) %>%
  inner_join(., asv_clin_weight_summary, by = "feature") %>%
  filter(feature == "prior aminoglycoside" | 
           feature == "prior combo betalactam" |
           feature == "prior carbapenem" |
           feature == "prior 3rd or 4th<br>generation cephalosporin" |
           feature == "prior cephalosporin" |
           feature == "high risk antibiotics" |
           feature == "prior penicillin" |
           feature == "prior penicillin") %>%
  mutate(feature = fct_reorder(feature, mean)) %>%
  ggplot(aes(x=mean, y=feature, xmin=l_quartile, xmax = u_quartile, color = direction)) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey") +
  geom_point(show.legend = FALSE, size = 3) +
  geom_linerange(show.legend = FALSE) +
  scale_color_manual(breaks = c("case", "control"),
                     values = c("red", "black")) +
  xlim(-0.002, 0.2) +
  labs(x = "Reduction in AUC when removed") +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
ggsave(glue("final_graphs/elastic_net_clin_asv_abx_feature_importance.pdf"), height = 8, width = 12, unit = "cm")

#REMOVAL OF TOP ASV FEATURES PERFORMANCE#########
#Read in en performance data
all_ASV_en_performance<-read_tsv(glue("{path}/elastic_net_performance/asv_elastic_net/processed_data/asv_en_performance.tsv")) %>%
  mutate(data = "All ASVs")
cutoff_en_performance<-read_tsv(glue("{path}/elastic_net_performance/asv_cutoff_important_features_elastic_net/processed_data/asv_cutoff_en_performance.tsv")) %>%
  mutate(data = "-Top ASVs")

#Combine en performance
cutoff_en_performance_data<-bind_rows(all_ASV_en_performance, cutoff_en_feature_importance)

#Test difference between model performance
AUC_pairwise_performance <- pairwise.t.test(cutoff_en_performance_data$AUC, g=cutoff_en_performance_data$data, p.adjust.method = "BH", paired=FALSE)

#Graph en performance data
cutoff_en_performance_data %>%
  select(seed, data, AUC) %>%
  filter(data != "clinical_genotype") %>%
  mutate(data = factor(data, levels=c("All ASVs","-Top ASVs"))) %>%
  ggplot(aes(y=AUC, x=data, color = data)) +
  geom_point(show.legend = FALSE, position = position_jitter(), alpha = 0.5, size = 0.75) +
  geom_line(data=tibble(x=c(1,2), y=c(0.91, 0.91)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=0.95), fill = NA, label.color = NA, label="*p* = 4.2e-9" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  stat_summary(fun = median,
               geom = "pointrange",
               fun.max = function(y) median(y) + sd(y),
               fun.min = function(y) median(y) - sd(y),
               orientation = "x",
               color = "black", 
               size = 0.5, 
               inherit.aes = TRUE) +
  scale_y_continuous(limits = c(0.2, 1), breaks = c(0.3, 0.5, 0.7, 0.9)) +
  scale_color_brewer(palette = "Paired") +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_blank())
ggsave(glue("final_graphs/elastic_net_cutoff_performance_data.pdf"), height = 8, width = 4, unit = "cm")

#Read in en feature importance data
cutoff_en_feature_importance<-read_tsv(glue("{path}/elastic_net_performance/asv_cutoff_important_features_elastic_net/processed_data/asv_cutoff_en_feature_importance.tsv"))

#Graph en feature importance data
cutoff_en_feature_importance %>%
  rename(feature = names) %>%
  group_by(feature) %>%
  summarize(mean = mean(perf_metric_diff * 1),
            l_quartile = quantile(perf_metric_diff * 1, prob = 0.25),
            u_quartile = quantile(perf_metric_diff * 1, prob = 0.75)) %>%
  mutate(feature = fct_reorder(feature, mean)) %>%
  filter(mean > 0.005) %>%
  ggplot(aes(x=mean, y=feature, xmin=l_quartile, xmax = u_quartile)) +
  geom_point(show.legend = FALSE, size = 3, color = "darkgrey") +
  geom_linerange(show.legend = FALSE, color = "darkgrey") +
  xlim(-0.001, 0.05) +
  labs(x = "Reduction in AUC when removed",
       subtitle = "Cutoff: AUC reduction > 0.005") +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank())

#TRAINED MODEL HP EXAMPLE#########
#Read in trained model data
ASV_trained_model<-read_tsv(glue("{path}/elastic_net_performance/asv_elastic_net/processed_data/asv_en_trained_model.tsv")) %>%
  mutate(log_lambda = log(lambda))
#Convert lambda values to characters to permit graphing
ASV_trained_model$lambda<-as.character(ASV_trained_model$lambda) 
ASV_trained_model$lambda<-fct_reorder(ASV_trained_model$lambda, ASV_trained_model$log_lambda)

#Graph trained model performance data
ASV_trained_model %>%
  ggplot(aes(x=alpha, y=lambda, fill=AUC)) +
  geom_tile(color = "black", lwd = 0.25, linetype = 1) + 
  scale_fill_gradientn(colors=topo.colors(7), limits = c(0.4,0.8)) +
  scale_x_continuous(limits = c(-0.05,1.05), expand = c(0,0), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  theme_classic() +
  theme(title = element_text()) +
  theme(axis.line = element_blank())
ggsave(glue("final_graphs/ASV_trained_elastic_net.pdf"), height = 8, width = 12, unit = "cm")