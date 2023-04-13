#LOAD REQUIRED PACKAGES##############
library(readxl)
library(tidyverse)
library(glue)
library(dplyr)
library(ggtext)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(broom)
library(purrr)

#SET ENVIRONMENT#####################################################
#Set working directory
setwd("/Volumes/Seagate Backup Plus Drive/mothur/Case_control/ASV/analysis")
#Set seed for all graphing
set.seed(20220212)
#Create graph output directory
dir.create(final_graphs)

#READ IN METADATA#####################################################
metadata<-read_excel(path="raw_data/case_control_metadata.xlsx")

#ASSIGNING TAG PREFIXES#####################################################
tag_prefix<-"final_case_control.asv"
graph_prefix<-"case_control_ASV_status"

#DEFINING AESTETICS#####################################################
group_name<-"Group"
group_breaks<-c("case", "control")
group_colors<-c("red", "black")
group_shapes<-c(16,15)
group_labels<-c("Case", "Control")

#DEFINING SAMPLE SIZE#####################################################
#Write sample size function
run_samplesize<-function(tag){
  samplesize<-glue('./mothur "#count.groups(shared={tag}.shared, inputdir=raw_data)"')
  system(samplesize)
}

#Run samples size function
run_samplesize(tag_prefix)

#CREATING SUBSAMPLE SHARE FILE#####################################################
#Assign sample size REMEMBER TO CUSTOMIZE FOR EACH PROJECT
subsample_size<-4438

#Create subsample function
run_subsample<-function(tag, size){
  subsample<-glue('./mothur "#sub.sample(shared={tag}.shared, size={size}, inputdir=raw_data, outputdir=raw_data)"')
  system(subsample)
}

#Running subsample
run_subsample(tag_prefix, subsample_size)

#ALPHA DIVERSITY#####################################################
#Write alpha diversity function
run_alpha<-function(tag){
  alpha<-glue('./mothur "#summary.single(shared={tag}.ASV.subsample.shared, subsample=F, inputdir=raw_data, outputdir=raw_data)"')
  system(alpha)
}

#Running alpha diversity function
run_alpha(tag_prefix)

#Create alpha diversity metadata file
alpha_diversity <- read_tsv(glue("raw_data/{tag_prefix}.ASV.subsample.groups.summary")) %>%
  mutate(invsimpson = 1/simpson)
metadata_alpha <- inner_join(metadata, alpha_diversity, by="group")

#Test for differences in alpha diversity
alpha_kw <- kruskal.test(invsimpson~case_control, data=metadata_alpha )
alpha_pairwise_invsimpson <- pairwise.t.test(metadata_alpha$invsimpson, g=metadata_alpha$case_control, p.adjust.method = "BH")
alpha_pairwise_shannon <- pairwise.t.test(metadata_alpha$shannon, g=metadata_alpha$case_control, p.adjust.method = "BH")
alpha_pairwise_chao <- pairwise.t.test(metadata_alpha$chao, g=metadata_alpha$case_control, p.adjust.method = "BH")

#Graph inverse simpson
metadata_alpha %>%
  ggplot(aes(x=case_control, y=invsimpson, shape=case_control, color=case_control, fill=case_control)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0) + 
  geom_jitter(show.legend = FALSE, width=0.25, size = 0.75, alpha = 0.5) +
  geom_line(data=tibble(x=c(1,2), y=c(47, 47)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=48.5), fill = NA, label.color = NA, label="*p* = 0.24" ,aes(x=x, y=y), inherit.aes=FALSE, size=4) +
  labs(x=NULL, 
       y="Inverse Simpson") +
  scale_x_discrete(breaks=group_breaks,
                   labels=group_labels) +
  scale_color_manual(breaks=group_breaks,
                     values=group_colors,
                     labels=group_labels) +
  scale_shape_manual(breaks=group_breaks,
                     values=group_shapes,
                     labels=group_labels) +
  scale_fill_manual(breaks=group_breaks,
                    values=group_colors,
                    labels=group_labels) + 
  theme_classic() +
  theme(axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank())
ggsave(glue("graphs/{graph_prefix}_inverse_simpson.pdf"), height = 10, width = 6, unit = "cm")

#Graph shannon
metadata_alpha %>%
  ggplot(aes(x=case_control, y=shannon, shape=case_control, color=case_control, fill=case_control)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0) + 
  geom_jitter(show.legend = FALSE, width=0.25, size = 0.75, alpha = 0.5) +
  geom_line(data=tibble(x=c(1,2), y=c(5.5, 5.5)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=5.75), fill = NA, label.color = NA, label="*p* = 0.039" ,aes(x=x, y=y), inherit.aes=FALSE, size=4) +
  labs(x=NULL, 
       y="Shannon") +
  scale_x_discrete(breaks=group_breaks,
                   labels=group_labels) +
  scale_color_manual(breaks=group_breaks,
                     values=group_colors,
                     labels=group_labels) +
  scale_shape_manual(breaks=group_breaks,
                     values=group_shapes,
                     labels=group_labels) +
  scale_fill_manual(breaks=group_breaks,
                    values=group_colors,
                    labels=group_labels) + 
  theme_classic() +
  theme(axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank())
ggsave(glue("graphs/{graph_prefix}_shannon.pdf"), height = 10, width = 6, unit = "cm")

#Graph chao
metadata_alpha %>%
  ggplot(aes(x=case_control, y=chao, shape=case_control, color=case_control, fill=case_control)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0) + 
  geom_jitter(show.legend = FALSE, width=0.25, size = 0.75, alpha = 0.5) +
  geom_line(data=tibble(x=c(1,2), y=c(1700, 1700)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=1800), fill = NA, label.color = NA, label="*p* = 0.053" ,aes(x=x, y=y), inherit.aes=FALSE, size=4) +
  labs(x=NULL, 
       y="Chao") +
  scale_x_discrete(breaks=group_breaks,
                   labels=group_labels) +
  scale_color_manual(breaks=group_breaks,
                     values=group_colors,
                     labels=group_labels) +
  scale_shape_manual(breaks=group_breaks,
                     values=group_shapes,
                     labels=group_labels) +
  scale_fill_manual(breaks=group_breaks,
                    values=group_colors,
                    labels=group_labels) + 
  theme_classic() +
  theme(axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank())
ggsave(glue("graphs/{graph_prefix}_chao.pdf"), height = 10, width = 6, unit = "cm")

#BETA DIVERSITY#####################################################
#Define variables
permutations<-10000
pairwise_p<-as.numeric()

#Create distance matrix function
run_beta<-function(tag){
  beta<-glue('./mothur "#dist.shared(shared={tag}.shared, output=square, calc=thetayc, inputdir=raw_data)"')
  system(beta)
}

#Run distance matrix function
run_beta(tag_prefix)

#Import distance matrix
distance <- read_tsv(file=glue("raw_data/{tag_prefix}.thetayc.ASV.square.dist"), skip=1, col_names=FALSE)
colnames(distance)<-c("group", distance$X1)
metadata_distance<-inner_join(metadata, distance, by="group")

#Write AMOVA function
run_amova<-function(x, tag){
  distance<-metadata_distance %>%
    select(all_of(.[["group"]])) %>%
    as.dist()
  test<-adonis(x, 
               data=metadata_distance, 
               permutations = permutations)
  capture.output(test, file=glue("processed_data/case_control_{tag}.amova_summary"))
  return(test$aov.tab$`Pr(>F)`[1])
}

#Run AMOVA function
case_control_var<-distance~case_control
pairwise_p['case_control']<-run_amova(case_control_var, 'case_control')

#p-value adjustment
adjust_pairwise_p<-p.adjust(pairwise_p, method="BH") %>% signif(digits=3)

#PCoA#####################################################
#Read in shared file
shared_file<-read_tsv(file=glue("raw_data/{tag_prefix}.ASV.subsample.shared"))

#Join share and metadata
shared_design<-inner_join(shared_file, metadata, by=c("Group" = "group"))

#PCoA function
run_pcoa<-function(x, tag){
  x <-shared_design %>%
    select(Group, x) %>%
    na.omit()
  groups<-capture.output(cat(x$Group, sep="-"))
  dist<-glue('./mothur "#dist.shared(shared=raw_data/{tag_prefix}.ASV.subsample.shared, output=square, calc=thetayc, inputdir=raw_data, groups={groups})"')
  system(dist)
  rename<-glue('./mothur "#rename.file(shared=raw_data/{tag_prefix}.ASV.subsample.thetayc.ASV.square.dist, prefix={tag})"')
  system(rename)
  pcoa<-glue('./mothur "#pcoa(phylip=raw_data/{tag}.dist, outputdir=processed_data)"')
  system(pcoa)
}

#Run PCoA function
run_pcoa("case_control", "case_control")

#Adjust shared for missing samples
adjusted_shared<-read_tsv(glue("raw_data/{tag_prefix}.ASV.subsample.shared")) %>%
  filter(Group != "PR14363") %>%
  filter(Group != "PR19028")
write_tsv(adjusted_shared, glue("processed_data/{tag_prefix}.ASV.subsample.adjust.shared"))

#Run function to create spearman corr axes file
run_axes<-glue('./mothur "#corr.axes(axes=processed_data/case_control.pcoa.axes, shared=processed_data/{tag_prefix}.ASV.subsample.adjust.shared, outputdir=processed_data, method=spearman, numaxes=2)"')
system(run_axes)

#Graph case control PCoA
tag <- "case_control"
axes <- read_tsv(file=glue("processed_data/{tag}.pcoa.axes"))
arrows<-read_tsv(glue("processed_data/{tag_prefix}.ASV.subsample.adjust.spearman.corr.axes")) %>%
  filter(OTU == "ASV000001" | OTU == "ASV000002"| OTU == "ASV000012"| OTU == "ASV000021"| OTU == "ASV000195")
loadings <- read_tsv(file=glue("processed_data/{tag}.pcoa.loadings"))
pcoa <- inner_join(metadata, axes, by="group")
loading1 <- loadings %>% filter(axis == 1) %>% pull(loading) %>% round(1)
loading2 <- loadings %>% filter(axis == 2) %>% pull(loading) %>% round(1)
adjusted_p <- adjust_pairwise_p['case_control']
ggplot(pcoa, aes(x=axis1, y=axis2, color=case_control, fill=case_control, shape=case_control)) +
  labs(x=glue("PCoA Axis 1 ({loading1}%)"), 
       y=glue("PCoA Axis 2 ({loading2}%)"),
       subtitle = glue("Adjusted *p*-value = {adjusted_p}"))+
  #ASV000001 Klebsiella arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = 0.847, yend = 0.182), size = 0.3,color = "darkgrey",
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_richtext(data=tibble(x = 1.1,y = 0.2), fill = NA, label.color = NA, label="ASV000001<br>*Klebsiella*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  #ASV000002 Enterococcus arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = -0.0422 , yend =-0.644), size = 0.3,color = "darkgrey", 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_richtext(data=tibble(x = -0.06,y = -0.8), fill = NA, label.color = NA, label="ASV000002<br>*Enterococcus*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  #ASV000012 Streptococcus arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = -0.0365 , yend =0.169),size = 0.3,color = "darkgrey", 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_segment(inherit.aes = FALSE, aes(x = -0.015, y = 0.08, xend = -0.75, yend = 0.5), size = 0.3,color = "darkgrey",linetype = "dashed") +
  geom_richtext(data=tibble(x = -0.75,y = 0.6), fill = NA, label.color = NA, label="ASV000012<br>*Streptococcus*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  #ASV000021 Akkermansia arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = -0.0864, yend = 0.203),size = 0.3, color = "darkgrey", 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_segment(inherit.aes = FALSE,aes(x = -0.04, y = 0.1, xend = 0.75, yend = 0.5),size = 0.3,color = "darkgrey",linetype = "dashed") +
  geom_richtext(data=tibble(x = 0.75,y = 0.6), fill = NA, label.color = NA, label="ASV000021<br>*Akkermansia*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  #ASV000195 Anaerostipes arrow
  geom_segment(inherit.aes = FALSE,aes(x = 0, y = 0, xend = 0.089, yend = -0.0401),size = 0.3,color = "darkgrey", 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_segment(inherit.aes = FALSE,aes(x = 0.04, y = -0.02, xend = 0.3, yend = -0.6),size = 0.3,color = "darkgrey",linetype = "dashed") +
  geom_richtext(data=tibble(x = 0.5,y = -0.725), fill = NA, label.color = NA, label="ASV000195<br>*Anaerostipes*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  stat_ellipse(level=0.95,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  geom_point(size = 0.75, alpha = 0.5) +
  coord_fixed(xlim=c(-1, 1.4), ylim = c(-1, 1)) +
  scale_color_manual(name=group_name,
                     breaks=group_breaks,
                     values=group_colors,
                     labels=group_labels) +
  scale_shape_manual(name=group_name,
                     breaks=group_breaks,
                     values=group_shapes,
                     labels=group_labels) +
  scale_fill_manual(name=group_name,
                    breaks=group_breaks,
                    values=group_colors,
                    labels=group_labels) +
  theme_classic() +
  theme(legend.key.height = unit(0.2, "cm"),
        legend.position = "bottom",
        legend.background = element_rect(color="black", size=0.4),
        legend.margin = margin(t=3, b=2, r=2, l=2),
        legend.title.align = 0.5,
        legend.title = element_text(face="bold", size=12),
        legend.text = element_markdown(size=10),
        axis.line = element_line(size = 0.5),
        axis.title = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 10),
        plot.subtitle = element_markdown(size = 8)
        )
ggsave(glue("graphs/{graph_prefix}_pcoa.pdf"), width = 12, height = 12, unit = "cm")

#READ AND JOIN ASV AND TAXONOMY FILES#####################################################
#Read in asv-specific metadata
asv_metadata<-read_excel("raw_data/case_control_metadata.xlsx") %>%
  select(group, case_control)

#Read in asv counts
asv_counts<-read_tsv(file=glue("raw_data/{tag_prefix}.ASV.subsample.shared")) %>%
  select(-label, -numOtus) %>%
  rename(group = Group) %>%
  pivot_longer(-group, names_to="otu", values_to = "count")

#Determine LOD
nseqs <- asv_counts %>%
  group_by(group) %>%
  summarize(N = sum(count), .group = "drop") %>%
  count(N) %>%
  pull(N)
LOD<-100* 1/nseqs

#Read in taxonomy and re-format taxon names
taxonomy<-read_tsv(file=glue("raw_data/{tag_prefix}.ASV.cons.taxonomy")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";") %>%
  mutate(pretty_otu = str_replace(string=otu,
                                  pattern="ASV",
                                  replacement = "ASV"),
         genus = str_replace(string=genus,
                             pattern="(.*)",
                             replacement="*\\1*"),
         genus = str_replace(string=genus,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified<br>*\\1*"),
         taxon = glue("{pretty_otu} {genus}")) %>%
  select(otu, taxon)

#Read in complete taxonomy
complete_taxonomy<-read_tsv(file=glue("raw_data/{tag_prefix}.ASV.cons.taxonomy")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";")

#Create relative abundance df
asv_rel_abun<-inner_join(asv_metadata, asv_counts, by="group") %>%
  inner_join(., taxonomy, by = "otu") %>%
  group_by(group) %>%
  mutate(rel_abund = 100*(count / sum(count))) %>%
  ungroup() %>%
  select(-count) %>%
  mutate(rel_abund = if_else(rel_abund ==0, 2/3 * LOD, rel_abund))

#Create relative abundance df sorted by taxonomic levels
asv_rel_abun_taxon_level<-inner_join(asv_metadata, asv_counts, by="group") %>%
  inner_join(., complete_taxonomy, by = "otu") %>%
  group_by(group) %>%
  mutate(rel_abund = 100*(count / sum(count))) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
               names_to="level",
               values_to="taxon") 

#ASV#####################################################
#Graph filtered asvs
asv_rel_abun %>%
  filter(taxon == "ASV000001 *Klebsiella*" | taxon == "ASV000002 *Enterococcus*" | taxon == "ASV000012 *Streptococcus*" | taxon == "ASV000021 *Akkermansia*"| taxon == "ASV000195 *Anaerostipes*") %>%
  mutate(taxon = factor(taxon, levels=c("ASV000195 *Anaerostipes*","ASV000021 *Akkermansia*","ASV000012 *Streptococcus*", "ASV000002 *Enterococcus*", "ASV000001 *Klebsiella*"))) %>%
  ggplot(aes(y=taxon, x=rel_abund, shape=case_control, fill=case_control, color = case_control)) +
  geom_boxplot(width = 0.8, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0, fatten = 0.75) + 
  geom_vline(xintercept=LOD, linetype = "dashed", color = "grey") +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 0.75, alpha=0.5) + 
  labs(y=NULL, 
       x="Relative Abundance (%)") +
  coord_trans(x="log10") +
  scale_y_discrete(breaks=c("ASV000195 *Anaerostipes*","ASV000021 *Akkermansia*","ASV000012 *Streptococcus*", "ASV000002 *Enterococcus*", "ASV000001 *Klebsiella*"),
                   labels =c("ASV000195<br>*Anaerostipes*","ASV000021<br>*Akkermansia*","ASV000012<br>*Streptococcus*", "ASV000002<br>*Enterococcus*", "ASV000001<br>*Klebsiella*")) +
  scale_x_continuous(limits=c(NA, 100),
                     breaks=c(0.01, 0.1, 1, 10, 100),
                     labels=c(0.01, 0.1, 1, 10, 100)) +
  scale_shape_manual(name=group_name,
                     breaks=group_breaks,
                     values=group_shapes,
                     labels=group_labels) +
  scale_color_manual(name=group_name,
                      breaks=group_breaks,
                      values=group_colors,
                      labels=group_labels) +
  scale_fill_manual(name=group_name,
                    breaks=group_breaks,
                    values=group_colors,
                    labels=group_labels) + 
  theme_classic() +
  theme(legend.key.height = unit(0.75, "cm"),
        legend.position = "bottom",
        legend.background = element_rect(color="black", size=0.6),
        legend.margin = margin(t=1, b=2, r=2, l=2),
        legend.title.align = 0.5,
        legend.title = element_markdown(face="bold", size=12),
        legend.text = element_markdown(size=10),
        axis.title = element_markdown(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text.x = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(color = "black",size = 10))
ggsave(glue("graphs/{graph_prefix}_significant_ASVs.pdf"), width = 12, height = 10, unit = "cm")

#Filter Klebsiella asvs
kleb_asv<-complete_taxonomy %>%
  filter(genus == "Klebsiella") %>%
  select(otu) %>%
  inner_join(., asv_rel_abun, by = "otu")

#Graph Klebsiella asvs
kleb_asv %>%
  filter(taxon == "ASV000001 *Klebsiella*" | taxon == "ASV000019 *Klebsiella*" | taxon == "ASV005157 *Klebsiella*" | taxon == "ASV016367 *Klebsiella*" ) %>%
  mutate(taxon = factor(taxon, levels=c("ASV016367 *Klebsiella*", "ASV005157 *Klebsiella*", "ASV000019 *Klebsiella*", "ASV000001 *Klebsiella*"))) %>%
  ggplot(aes(y=taxon, x=rel_abund, shape=case_control, fill=case_control, color=case_control)) +
  geom_vline(xintercept=LOD, linetype = "dashed", color = "grey") +
  geom_boxplot(width = 0.8, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0, fatten = 0.75) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 0.75, alpha=0.5) + 
  labs(y=NULL, 
       x="Relative Abundance (%)") +
  coord_trans(x="log10") +
  scale_x_continuous(limits=c(NA, 100),
                     breaks=c(0.01, 0.1, 1, 10, 100),
                     labels=c(0.01, 0.1, 1, 10, 100)) +
  scale_shape_manual(name=group_name,
                     breaks=group_breaks,
                     values=group_shapes,
                     labels=group_labels) +
  scale_color_manual(name=group_name,
                     breaks=group_breaks,
                     values=group_colors,
                     labels=group_labels) +
  scale_fill_manual(name=group_name,
                    breaks=group_breaks,
                    values=group_colors,
                    labels=group_labels) + 
  theme_classic() +
  theme(legend.key.height = unit(0.75, "cm"),
        legend.position = "bottom",
        legend.background = element_rect(color="black", size=0.6),
        legend.margin = margin(t=1, b=2, r=2, l=2),
        legend.title.align = 0.5,
        legend.title = element_markdown(face="bold", size=12),
        legend.text = element_markdown(size=10),
        axis.title = element_markdown(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text.x = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(color = "black", size = 10))
ggsave(glue("graphs/{graph_prefix}_kleb_ASVs.pdf"), width = 12, height = 8, unit = "cm")

#Read in unformatted complete taxonomy
unf_complete_taxonomy<-read_tsv(file=glue("raw_data/{tag_prefix}.ASV.cons.taxonomy")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";")

#Creating list of Klebsiella asvs
kleb_list<-as.vector(unf_complete_taxonomy %>%
                       filter(genus == "Klebsiella") %>%
                       select(otu))

#Create file of shared Klebsiella counts - ASVs are from kleb_list
kleb_total_counts<-read_tsv(file=glue("raw_data/{tag_prefix}.shared")) %>%
  select(Group, 
         ASV000001,
         ASV000019,
         ASV005157,
         ASV016367,
         ASV017041,
         ASV028116,
         ASV028673,
         ASV037014,
         ASV038693,
         ASV044066,
         ASV051092,
         ASV056351,
         ASV064026,
         ASV068895,
         ASV069190,
         ASV074238,
         ASV077560,
         ASV082227,
         ASV084575,
         ASV088272,
         ASV088355,
         ASV088424,
         ASV089180,
         ASV090069,
         ASV092533,
         ASV096162,
         ASV100038,
         ASV104480,
         ASV110379,
         ASV113813) %>%
  rename(group = Group)
write_csv(kleb_total_counts, file="processed_data/kleb_ASV_counts.csv")

#NOTE: at this point, the kleb_ASV_counts.csv needs to manually edited to add the "Count" column used below by counting the total unique Klebsiella ASVs per sample

#Graph Klebsiella counts
read_csv(file="processed_data/kleb_ASV_counts.csv") %>%
  ggplot(aes(x=ASV_count, fill = "red")) +
  geom_bar(show.legend = FALSE) +
  labs(y="Count", x="*Klebsiella* ASVs") +
  scale_fill_manual(values=c("darkorchid")) + 
  theme_classic() + 
  theme(axis.title = element_markdown(size = 12, face = "bold"),
        axis.title.x = element_markdown(),
        axis.line = element_line(size = 0.5),
        axis.text = element_markdown(color = "black", size = 10),
        plot.subtitle = element_markdown(size = 8)
  )
ggsave(glue("graphs/Kleb_ASV_count_barchart.pdf"), width = 6, height = 8, unit = "cm")

#Create ASV000019 +/- list
ASV19<-kleb_asv %>%
  filter(taxon == "ASV000019 *Klebsiella*") %>%
  mutate(ASV19 = case_when(rel_abund > 0.0151 ~ "yes",
                           rel_abund < 0.0151 ~ "no")) %>%
  select(group, ASV19)

#Graph ASV000001 stratified by the presence of ASV000019
kleb_asv %>%
  filter(taxon == "ASV000001 *Klebsiella*") %>%
  inner_join(., ASV19, by = "group") %>%
  select(group, rel_abund, ASV19) %>%
  ggplot(aes(x=ASV19, y=rel_abund, shape=ASV19, fill=ASV19, color=ASV19)) +
  geom_hline(yintercept=LOD, linetype = "dashed", color = "grey") +
  geom_boxplot(width = 0.8, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0, fatten = 0.75, show.legend = FALSE) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 0.75, alpha=0.5, show.legend = FALSE) + 
  labs(x=NULL, 
       y="Relative Abundance (%)") +
  coord_trans(y="log10") +
  scale_y_continuous(limits=c(NA, 100),
                     breaks=c(0.01, 0.1, 1, 10, 100),
                     labels=c(0.01, 0.1, 1, 10, 100)) +
  scale_x_discrete(breaks = c("yes", "no"),
                   labels = c("+ASV000019", "-ASV000019")) +
  scale_shape_manual(breaks = c("yes", "no"),
                     values=c(21, 22)) +
  scale_color_manual(breaks = c("yes", "no"),
                     values=c("cornflowerblue", "darkorange")) +
  scale_fill_manual(breaks = c("yes", "no"),
                    values=c("cornflowerblue", "darkorange")) + 
  theme_classic() +
  theme(legend.key.height = unit(0.75, "cm"),
        legend.position = "bottom",
        legend.background = element_rect(color="black", size=0.6),
        legend.margin = margin(t=1, b=2, r=2, l=2),
        legend.title.align = 0.5,
        legend.title = element_markdown(face="bold", size=12),
        legend.text = element_markdown(size=10),
        axis.title = element_markdown(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text.x = element_markdown(color = "black", size = 10, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_markdown(color = "black", size = 10))
ggsave(glue("graphs/ASV000019_ASV000001_rel_abund.pdf"), width = 4, height = 8, unit = "cm")

#LEFSE/LDA#####################################################
#Read in metadata, limited to design data
metadata<-read_excel("raw_data/case_control_metadata.xlsx") %>%
  select(group, case_control)

#Read in shared file
shared_file<-read_tsv(file=glue("raw_data/{tag_prefix}.ASV.subsample.shared"))

#Join share and metadata
shared_design<-inner_join(shared_file, metadata, by=c("Group" = "group"))

#Write LefSe function
run_lefse<-function(tag){
  #Select specific comparison groups
  x_y <-shared_design
  #Write shared file for pairwise comparison
  x_y %>%
    select(-case_control) %>%
    write_tsv(glue("processed_data/{tag}.shared"))
  #Write design file for pairwise comparison
  x_y %>%
    select(Group, case_control) %>%
    write_tsv(glue("processed_data/{tag}.design"))
  #Assign lefse function
  lefse<-glue('./mothur "#lefse(shared={tag}.shared, design={tag}.design, inputdir=processed_data)"')
  #Run lefse function
  system(lefse)
  return(glue("processed_data/{tag}.ASV.lefse_summary"))
}

#Run LefSe function
case_control_ASV_status<-run_lefse("case_control_ASV_status")

#Graph case control LDA
read_tsv(case_control_ASV_status) %>%
  drop_na(LDA) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(LDA >= 2.5) %>%
  mutate(LDA = if_else(Class == "control", -1*LDA, LDA),
         taxon = fct_reorder(taxon, LDA), 
         label_x = if_else(Class == "control", 0.1, -0.1),
         label_hjust = if_else(Class == "control", 0, 1)) %>%
  ggplot(aes(x=LDA, y=taxon, label=taxon, fill=Class)) +
  geom_col() +
  geom_richtext(aes(x=label_x, hjust=label_hjust), fill = NA, label.color = NA, size = 3) + 
  labs(y=NULL, x="LDA Score (log 10)") +
  scale_fill_manual(name=group_name,
                    breaks=group_breaks,
                    values=group_colors,
                    labels=group_labels) + 
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_markdown(size = 12, face = "bold"),
        axis.text.x = element_markdown(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10))
ggsave(glue("graphs/{graph_prefix}_ASV_LefSe.pdf"), width = 12, height = 6, unit = "cm")

#ter ASV000001 ANALYSIS#######################
#Create list of samples with terF
terF_list<-read.csv("raw_data/gene_presence_absence_pivot.csv", skip = 1) %>%
  select(Gene, group_4492) %>%
  filter(group_4492 == 1 | group_4492 == 0) %>%
  mutate(group_4492 = case_when(group_4492 == 1 ~ "yes",
                                group_4492 == 0 ~ "no"))
key<-read_excel("raw_data/WGS_PR_key.xlsx")

#Graph ASV000001 abundance stratified by terF
inner_join(key, terF_list, by = c("prokka_ID" = "Gene")) %>%
  select(group, group_4492) %>%
  inner_join(.,asv_rel_abun, by = "group") %>%
  filter(otu == "ASV000001") %>%
  ggplot(aes(y=group_4492, x=rel_abund, shape=group_4492, fill=group_4492, color=group_4492)) +
  geom_vline(xintercept=LOD, linetype = "dashed", color = "grey") +
  geom_boxplot(width = 0.8, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0, fatten = 0.75, show.legend = FALSE) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 0.75, alpha=0.5, show.legend = FALSE) + 
  labs(y = NULL, 
       x="Relative Abundance (%)") +
  coord_trans(x="log10") +
  scale_y_discrete(breaks=c("no", "yes"),
                   labels =c("No *terF*", "*terF*")) +
  scale_x_continuous(limits=c(NA, 100),
                     breaks=c(0.01, 0.1, 1, 10, 100),
                     labels=c(0.01, 0.1, 1, 10, 100)) +
  scale_shape_manual(values = c(15,16)) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") + 
  theme_classic() +
  theme(legend.key.height = unit(0.75, "cm"),
        legend.position = "bottom",
        legend.background = element_rect(color="black", size=0.6),
        legend.margin = margin(t=1, b=2, r=2, l=2),
        legend.title.align = 0.5,
        legend.title = element_markdown(face="bold", size=12),
        legend.text = element_markdown(size=10),
        axis.title = element_markdown(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text.x = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(color = "black",size = 10))
ggsave(glue("graphs/terF_ASV000001.pdf"), width = 12, height = 5, unit = "cm")

