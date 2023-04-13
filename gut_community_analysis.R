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
setwd("/Volumes/Seagate Backup Plus Drive/mothur/Case_control/analysis")
#Set seed for all graphing
set.seed(20220212)
#Create graph output directory
dir.create("final_graphs")

#READ IN METADATA#####################################################
metadata<-read_excel(path="raw_data/case_control_metadata.xlsx")

#ASSIGNING TAG PREFIXES#####################################################
tag_prefix<-"Final_case_control.opti_mcc"
graph_prefix <- "gut_community"

#DEFINING AESTETICS#####################################################
group_name<-"Partition"
group_breaks<-c("Partition_1", "Partition_2")
group_colors<-c("#00274C", "#FFCB05")
group_shapes<-c(16,15)
group_labels<-c("Partition 1", "Partition 2")

#DEFINING SAMPLE SIZE#####################################################
#Write sample size function
run_samplesize<-function(tag){
  samplesize<-glue('./mothur "#count.groups(shared={tag}.shared, inputdir=raw_data)"')
  system(samplesize)
}

#Run samples size function
run_samplesize(case_control_prefix)

#CREATING SUBSAMPLE SHARE FILE#####################################################
#Assign sample size REMEMBER TO CUSTOMIZE FOR EACH PROJECT
subsample_size<-4438

#Create subsample function
run_subsample<-function(tag, size){
  subsample<-glue('./mothur "#sub.sample(shared={tag}.shared, size={size}, inputdir=raw_data, outputdir=raw_data)"')
  system(subsample)
}

#Run subsample function
run_subsample(case_control_prefix, subsample_size)


#RUN COMMUNITY TYPE#####################################################
#Create community type function
run_ct<-function(tag){
  ct<-glue('./mothur "#get.communitytype(shared={tag}.shared, inputdir=raw_data, outputdir=processed_data)"')
  system(ct)
}

#Run community type function
run_ct(tag_prefix)

#READ AND JOIN OTU AND TAXONOMY FILES#####################################################
#Read in otu-specific metadata
otu_metadata<-read_excel("raw_data/case_control_metadata.xlsx") %>%
  select(group, case_control)

#Read in otu counts
otu_counts<-read_tsv(file=glue("raw_data/{tag_prefix}.0.03.subsample.shared")) %>%
  select(-label, -numOtus) %>%
  rename(group = Group) %>%
  pivot_longer(-group, names_to="otu", values_to = "count")

#Determine LOD
nseqs <- otu_counts %>%
  group_by(group) %>%
  summarize(N = sum(count), .group = "drop") %>%
  count(N) %>%
  pull(N)
LOD<-100* 1/nseqs

#Read in taxonomy and re-format taxon names
taxonomy<-read_tsv(file=glue("raw_data/{tag_prefix}.0.03.cons.taxonomy")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";") %>%
  mutate(pretty_otu = str_replace(string=otu,
                                  pattern="Otu",
                                  replacement = "OTU"),
         genus = str_replace(string=genus,
                             pattern="(.*)",
                             replacement="*\\1*"),
         genus = str_replace(string=genus,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified<br>*\\1*"),
         taxon = glue("{pretty_otu} {genus}")) %>%
  select(otu, taxon)

#Read in complete taxonomy
complete_taxonomy<-read_tsv(file=glue("raw_data/{tag_prefix}.0.03.cons.taxonomy")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";")

#Create relative abundance df
otu_rel_abun<-inner_join(otu_metadata, otu_counts, by="group") %>%
  inner_join(., taxonomy, by = "otu") %>%
  group_by(group) %>%
  mutate(rel_abund = 100*(count / sum(count))) %>%
  ungroup() %>%
  select(-count) %>%
  mutate(rel_abund = if_else(rel_abund ==0, 2/3 * LOD, rel_abund))

#Create relative abundance df sorted by taxonomic levels
otu_rel_abun_taxon_level<-inner_join(otu_metadata, otu_counts, by="group") %>%
  inner_join(., complete_taxonomy, by = "otu") %>%
  group_by(group) %>%
  mutate(rel_abund = 100*(count / sum(count))) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
               names_to="level",
               values_to="taxon") 

#OTU#####################################################
#Graph filtered otus
otu_rel_abun %>%
  filter(otu == "Otu00001" | otu == "Otu00002" | otu == "Otu00003" | otu == "Otu00004" | otu == "Otu00005" ) %>%
  mutate(taxon = factor(taxon, levels=c("OTU00005 *Peptoniphilus*","OTU00004 *Finegoldia*","OTU00003 *Escherichia/Shigella*","OTU00002 *Enterococcus*", "OTU00001 *Klebsiella*"))) %>%
  ggplot(aes(y=taxon, x=rel_abund, fill=taxon, color = taxon)) +
  geom_vline(xintercept=LOD, linetype = "dashed", color = "grey") +
  geom_boxplot(width = 0.8, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0, fatten = 0.75, show.legend = FALSE) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 2), size = 0.75, alpha=0.5, show.legend = FALSE) + 
  labs(y=NULL, 
       x="Relative Abundance (%)") +
  coord_trans(x="log10") +
  scale_y_discrete(breaks=c("OTU00005 *Peptoniphilus*","OTU00004 *Finegoldia*","OTU00003 *Escherichia/Shigella*","OTU00002 *Enterococcus*", "OTU00001 *Klebsiella*"),
                   labels =c("OTU00005<br>*Peptoniphilus*","OTU00004<br>*Finegoldia*","OTU00003<br>*Escherichia/Shigella*","OTU00002<br>*Enterococcus*", "OTU00001<br>*Klebsiella*")) +
  scale_x_continuous(limits=c(NA, 200),
                     breaks=c(0.01, 0.1, 1, 10, 100),
                     labels=c(0.01, 0.1, 1, 10, 100)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  theme_classic() +
  theme(axis.title = element_markdown(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text.x = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(color = "black",size = 10))
ggsave(glue("final_graphs/{graph_prefix}_top_OTUs.pdf"), width = 12, height =8, unit = "cm")

#2 COMMUNITY TYPE ALPHA###################
#Write alpha diversity function
run_alpha<-function(tag){
  alpha<-glue('./mothur "#summary.single(shared={tag}.0.03.subsample.shared, subsample=F, inputdir=raw_data, outputdir=raw_data)"')
  system(alpha)
}

#Run alpha diversity function
run_alpha(tag_prefix)

#Create alpha diversity metadata file
alpha_diversity <- read_tsv(glue("raw_data/{tag_prefix}.0.03.subsample.groups.summary")) %>%
  mutate(invsimpson = 1/simpson)
metadata_alpha <- inner_join(metadata_design, alpha_diversity, by="group")

#Test for differences in alpha diversity
alpha_pairwise_invsimpson <- pairwise.t.test(metadata_alpha$invsimpson, g=metadata_alpha$partition, p.adjust.method = "BH")
alpha_pairwise_shannon <- pairwise.t.test(metadata_alpha$shannon, g=metadata_alpha$partition, p.adjust.method = "BH")
alpha_pairwise_chao <- pairwise.t.test(metadata_alpha$chao, g=metadata_alpha$partition, p.adjust.method = "BH")

#Graph inverse simpson
metadata_alpha %>%
  ggplot(aes(x=partition, y=invsimpson, shape=partition, color=partition, fill=partition)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0) + 
  geom_jitter(show.legend = FALSE, width=0.25, alpha = 0.5) +
  geom_line(data=tibble(x=c(1,2), y=c(40, 40)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=41.5), fill = NA, label.color = NA, label="*p* < 2e-16" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
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
ggsave(glue("final_graphs/2_partition_inverse_simpson.pdf"), height = 10, width = 6, unit = "cm")

#Graph shannon
metadata_alpha %>%
  ggplot(aes(x=partition, y=shannon, shape=partition, color=partition, fill=partition)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0) + 
  geom_jitter(show.legend = FALSE, width=0.25, alpha = 0.5) +
  geom_line(data=tibble(x=c(1,2), y=c(5, 5)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=5.2), fill = NA, label.color = NA, label="*p* < 2e-16" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
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
ggsave(glue("final_graphs/2_partition_shannon.pdf"), height = 10, width = 6, unit = "cm")

#Graph chao
metadata_alpha %>%
  ggplot(aes(x=partition, y=chao, shape=partition, color=partition, fill=partition)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0) + 
  geom_jitter(show.legend = FALSE, width=0.25, alpha = 0.5) +
  geom_line(data=tibble(x=c(1,2), y=c(725, 725)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=750), fill = NA, label.color = NA, label="*p* < 2e-16" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
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
ggsave(glue("final_graphs/2_partition_chao.pdf"), height = 10, width = 6, unit = "cm")

#PCoA 2 COMMUNITY TYPE#####################################################
#Read in metadata, limited to design data
design<-read_tsv(file=glue("processed_data/{tag_prefix}.0.03.subsample.0.03.dmm.mix.design"), col_names = FALSE)
colnames(design) <- c("group", "partition")
shared_file<-read_tsv(file=glue("raw_data/{tag_prefix}.0.03.subsample.shared"))

#Join share and metadata
shared_design<-inner_join(shared_file, design, by=c("Group" = "group"))

#Write PCoA function
run_pcoa<-function(x, tag){
  x <-shared_design %>%
    select(Group, x) %>%
    na.omit()
  groups<-capture.output(cat(x$Group, sep="-"))
  dist<-glue('./mothur "#dist.shared(shared=raw_data/{tag_prefix}.0.03.subsample.shared, output=square, calc=thetayc, inputdir=raw_data, groups={groups})"')
  system(dist)
  rename<-glue('./mothur "#rename.file(shared=raw_data/{tag_prefix}.0.03.subsample.thetayc.0.03.square.dist, prefix={tag})"')
  system(rename)
  pcoa<-glue('./mothur "#pcoa(phylip=raw_data/{tag}.opti_mcc.dist, outputdir=processed_data)"')
  system(pcoa)
}

#Run PCoA function 
run_pcoa("partition", "ct_partition")

#Adjust shared for missing samples
adjusted_shared<-read_tsv(glue("raw_data/{tag_prefix}.0.03.subsample.shared")) %>%
  filter(Group != "PR14363") %>%
  filter(Group != "PR19028")
write_tsv(adjusted_shared, glue("processed_data/{tag_prefix}.0.03.subsample.adjust.shared"))

#Run function to create spearman corr axes file
run_axes<-glue('./mothur "#corr.axes(axes=processed_data/case_control.opti_mcc.pcoa.axes, shared=processed_data/{tag_prefix}.0.03.subsample.adjust.shared, outputdir=processed_data, method=spearman, numaxes=2)"')
system(run_axes)

#Graph 2 partition PCoA
tag <- "ct_partition"
axes <- read_tsv(file=glue("processed_data/{tag}.opti_mcc.pcoa.axes"))
arrows<-read_tsv(glue("processed_data/{tag_prefix}.0.03.subsample.adjust.spearman.corr.axes")) %>%
  filter(OTU == "Otu00001" | OTU == "Otu00002" | OTU == "Otu00003" | OTU == "Otu00004" | OTU == "Otu00005" 
         | OTU == "Otu00006"| OTU == "Otu00007" | OTU == "Otu00008" | OTU == "Otu00009" | OTU == "Otu00010")
loadings <- read_tsv(file=glue("processed_data/{tag}.opti_mcc.pcoa.loadings"))
pcoa <- inner_join(design, axes, by="group")
loading1 <- loadings %>% filter(axis == 1) %>% pull(loading) %>% round(1)
loading2 <- loadings %>% filter(axis == 2) %>% pull(loading) %>% round(1)
ggplot(pcoa, aes(x=axis1, y=axis2, color=partition, fill=partition, shape=partition)) +
  #Klebsiella arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = 0.886, yend = 0.3), size = 0.3,color = "darkgrey",
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_richtext(data=tibble(x = 1.1,y = 0.3), fill = NA, label.color = NA, label="OTU00001<br>*Klebsiella*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  #Enterococcus arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = 0.136, yend = -0.189), size = 0.3,color = "darkgrey", 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_segment(inherit.aes = FALSE, aes(x = 0.07, y = -0.11, xend = 0.07, yend = -0.6), size = 0.3,color = "darkgrey",linetype = "dashed") +
  geom_richtext(data=tibble(x = 0.07,y = -0.7), fill = NA, label.color = NA, label="OTU00002<br>*Enterococcus*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  #Escherichia/Shigella arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = -0.353, yend = 0.704),size = 0.3,color = "darkgrey", 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_richtext(data=tibble(x = -0.353,y = 0.8), fill = NA, label.color = NA, label="OTU00003<br>*Escherichia/Shigella*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  #Finegoldia arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = -0.528, yend = -0.334),size = 0.3,color = "darkgrey", 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_segment(inherit.aes = FALSE, aes(x = -0.25, y = -0.15, xend = -0.75, yend = -0.75), size = 0.3,color = "darkgrey",linetype = "dashed") +
  geom_richtext(data=tibble(x = -0.8,y = -.85), fill = NA, label.color = NA, label="OTU00004<br>*Fingoldia*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) + 
  #Peptoniphilus arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = -0.511, yend = -0.249),size = 0.3,color = "darkgrey", 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_segment(inherit.aes = FALSE, aes(x = -0.25, y = -0.1, xend = -0.75, yend = 0.25), size = 0.3,color = "darkgrey",linetype = "dashed") +
  geom_richtext(data=tibble(x = -0.8,y = 0.35), fill = NA, label.color = NA, label="OTU00005<br>*Peptoniphilus*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  stat_ellipse(level=0.95,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  geom_point(size = 0.75, alpha = 0.5) +
  coord_fixed(xlim=c(-1, 1.2), ylim = c(-1, 1)) +
  labs(x=glue("PCoA Axis 1 ({loading1}%)"), 
       y=glue("PCoA Axis 2 ({loading2}%)"))+
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
        plot.subtitle = element_markdown(size = 8))
ggsave(glue("final_graphs/2_partition_pcoa.pdf"), width = 12, height = 12, unit = "cm")

#3 COMMUNITY TYPE ALPHA###################
#Read in 3 partition design
partition_3_design<-read_tsv("raw_data/3_partition.design")

#Create alpha diversity metadata file
alpha_diversity <- read_tsv(glue("raw_data/{tag_prefix}.0.03.subsample.groups.summary")) %>%
  mutate(invsimpson = 1/simpson)
metadata_alpha <- inner_join(partition_3_design, alpha_diversity, by="group")

#Test for differences in alpha diversity
alpha_pairwise_invsimpson <- pairwise.t.test(metadata_alpha$invsimpson, g=metadata_alpha$partition, p.adjust.method = "BH")
alpha_pairwise_shannon <- pairwise.t.test(metadata_alpha$shannon, g=metadata_alpha$partition, p.adjust.method = "BH")
alpha_pairwise_chao <- pairwise.t.test(metadata_alpha$chao, g=metadata_alpha$partition, p.adjust.method = "BH")

#Graph inverse simpson
metadata_alpha %>%
  ggplot(aes(x=partition, y=invsimpson, shape=partition, color=partition, fill=partition)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0) + 
  geom_jitter(show.legend = FALSE, width=0.25, alpha = 0.5) +
  geom_line(data=tibble(x=c(1,2), y=c(40, 40)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=41.5), fill = NA, label.color = NA, label="*p* < 2e-16" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_line(data=tibble(x=c(2,3), y=c(41, 41)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=2.5, y=42.5), fill = NA, label.color = NA, label="*p* < 2e-16" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_line(data=tibble(x=c(1,3), y=c(45, 45)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=2, y=46.5), fill = NA, label.color = NA, label="*p* = 5.3e-5" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  labs(x=NULL, 
       y="Inverse Simpson") +
  scale_x_discrete(breaks=c("partition_1", "partition_2", "partition_3"),
                   labels=c("Partition 1", "Partition 2", "Partition 3")) +
  scale_color_manual(name="Partition",
                     breaks=c("partition_1", "partition_2", "partition_3"),
                     values=c("#e2e14c", "#8b324d", "#567243"),
                     labels=c("Partition 1", "Partition 2", "Partition 3")) +
  scale_shape_manual(name="Partition",
                     breaks=c("partition_1", "partition_2", "partition_3"),
                     values=c(16, 15, 17),
                     labels=c("Partition 1", "Partition 2", "Partition 3")) +
  scale_fill_manual(name="Partition",
                    breaks=c("partition_1", "partition_2", "partition_3"),
                    values=c("#e2e14c", "#8b324d", "#567243"),
                    labels=c("Partition 1", "Partition 2", "Partition 3")) +
  theme_classic() +
  theme(axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank())
ggsave(glue("final_graphs/3_partition_inverse_simpson.pdf"), height = 10, width = 8, unit = "cm")

#Graph shannon
metadata_alpha %>%
  ggplot(aes(x=partition, y=shannon, shape=partition, color=partition, fill=partition)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0) + 
  geom_jitter(show.legend = FALSE, width=0.25, alpha = 0.5) +
  geom_line(data=tibble(x=c(1,2), y=c(4.5, 4.5)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=4.65), fill = NA, label.color = NA, label="*p* = 3.1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_line(data=tibble(x=c(2,3), y=c(4.6, 4.6)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=2.5, y=4.75), fill = NA, label.color = NA, label="*p* < 2e-16" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_line(data=tibble(x=c(1,3), y=c(5, 5)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=2, y=5.15), fill = NA, label.color = NA, label="*p* < 2e-16" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  labs(x=NULL, 
       y="Shannon") +
  scale_x_discrete(breaks=c("partition_1", "partition_2", "partition_3"),
                   labels=c("Partition 1", "Partition 2", "Partition 3")) +
  scale_color_manual(name="Partition",
                     breaks=c("partition_1", "partition_2", "partition_3"),
                     values=c("#e2e14c", "#8b324d", "#567243"),
                     labels=c("Partition 1", "Partition 2", "Partition 3")) +
  scale_shape_manual(name="Partition",
                     breaks=c("partition_1", "partition_2", "partition_3"),
                     values=c(16, 15, 17),
                     labels=c("Partition 1", "Partition 2", "Partition 3")) +
  scale_fill_manual(name="Partition",
                    breaks=c("partition_1", "partition_2", "partition_3"),
                    values=c("#e2e14c", "#8b324d", "#567243"),
                    labels=c("Partition 1", "Partition 2", "Partition 3")) +
  theme_classic() +
  theme(axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank())
ggsave(glue("final_graphs/3_partition_shannon.pdf"), height = 10, width = 8, unit = "cm")

#Graph chao
metadata_alpha %>%
  ggplot(aes(x=partition, y=chao, shape=partition, color=partition, fill=partition)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0) + 
  geom_jitter(show.legend = FALSE, width=0.25, alpha = 0.5) +
  geom_line(data=tibble(x=c(1,2), y=c(750, 750)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=775), fill = NA, label.color = NA, label="*p* = 1.7e-6" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_line(data=tibble(x=c(2,3), y=c(775, 775)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=2.5, y=800), fill = NA, label.color = NA, label="*p* < 2e-16" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_line(data=tibble(x=c(1,3), y=c(850, 850)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=2, y=875), fill = NA, label.color = NA, label="*p* = 4.6e-11" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  labs(x=NULL, 
       y="Chao") +
  scale_x_discrete(breaks=c("partition_1", "partition_2", "partition_3"),
                   labels=c("Partition 1", "Partition 2", "Partition 3")) +
  scale_color_manual(name="Partition",
                     breaks=c("partition_1", "partition_2", "partition_3"),
                     values=c("#e2e14c", "#8b324d", "#567243"),
                     labels=c("Partition 1", "Partition 2", "Partition 3")) +
  scale_shape_manual(name="Partition",
                     breaks=c("partition_1", "partition_2", "partition_3"),
                     values=c(16, 15, 17),
                     labels=c("Partition 1", "Partition 2", "Partition 3")) +
  scale_fill_manual(name="Partition",
                    breaks=c("partition_1", "partition_2", "partition_3"),
                    values=c("#e2e14c", "#8b324d", "#567243"),
                    labels=c("Partition 1", "Partition 2", "Partition 3")) +
  theme_classic() +
  theme(axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank())
ggsave(glue("final_graphs/3_partition_chao.pdf"), height = 10, width = 8, unit = "cm")

#PCoA 3 COMMUNITY TYPE#####################################################
#Read in 3 partition design
partition_3_design<-read_tsv("raw_data/3_partition.design")
partition_3_shared_design<-inner_join(shared_file, partition_3_design, by=c("Group" = "group"))

#Write PCoA function
partition_3_run_pcoa<-function(x, tag){
  x <-partition_3_shared_design %>%
    select(Group, x) %>%
    na.omit()
  groups<-capture.output(cat(x$Group, sep="-"))
  dist<-glue('./mothur "#dist.shared(shared=raw_data/{tag_prefix}.0.03.subsample.shared, output=square, calc=thetayc, inputdir=raw_data, groups={groups})"')
  system(dist)
  rename<-glue('./mothur "#rename.file(shared=raw_data/{tag_prefix}.0.03.subsample.thetayc.0.03.square.dist, prefix={tag})"')
  system(rename)
  pcoa<-glue('./mothur "#pcoa(phylip=raw_data/{tag}.opti_mcc.dist, outputdir=processed_data)"')
  system(pcoa)
}

#Run PCoA
partition_3_run_pcoa("partition", "ct_partition_3")

#Graph 3 partition PCoA
tag <- "ct_partition_3"
axes <- read_tsv(file=glue("processed_data/{tag}.opti_mcc.pcoa.axes"))
arrows<-read_tsv(glue("processed_data/{tag_prefix}.0.03.subsample.adjust.spearman.corr.axes")) %>%
  filter(OTU == "Otu00001" | OTU == "Otu00002" | OTU == "Otu00003" | OTU == "Otu00004" | OTU == "Otu00005")
loadings <- read_tsv(file=glue("processed_data/{tag}.opti_mcc.pcoa.loadings"))
pcoa <- inner_join(partition_3_design, axes, by="group")
loading1 <- loadings %>% filter(axis == 1) %>% pull(loading) %>% round(1)
loading2 <- loadings %>% filter(axis == 2) %>% pull(loading) %>% round(1)
ggplot(pcoa, aes(x=axis1, y=axis2, color=partition, fill=partition, shape=partition)) +
  #Klebsiella arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = 0.886, yend = 0.3), size = 0.3,color = "darkgrey",
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_richtext(data=tibble(x = 1.1,y = 0.3), fill = NA, label.color = NA, label="OTU00001<br>*Klebsiella*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  #Enterococcus arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = 0.136, yend = -0.189), size = 0.3,color = "darkgrey", 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_segment(inherit.aes = FALSE, aes(x = 0.07, y = -0.11, xend = 0.07, yend = -0.6), size = 0.3,color = "darkgrey",linetype = "dashed") +
  geom_richtext(data=tibble(x = 0.07,y = -0.7), fill = NA, label.color = NA, label="OTU00002<br>*Enterococcus*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  #Escherichia/Shigella arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = -0.353, yend = 0.704),size = 0.3,color = "darkgrey", 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_richtext(data=tibble(x = -0.353,y = 0.8), fill = NA, label.color = NA, label="OTU00003<br>*Escherichia/Shigella*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  #Finegoldia arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = -0.528, yend = -0.334),size = 0.3,color = "darkgrey", 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_segment(inherit.aes = FALSE, aes(x = -0.25, y = -0.15, xend = -0.75, yend = -0.75), size = 0.3,color = "darkgrey",linetype = "dashed") +
  geom_richtext(data=tibble(x = -0.8,y = -.85), fill = NA, label.color = NA, label="OTU00004<br>*Fingoldia*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) + 
  #Peptoniphilus arrow
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0, xend = -0.511, yend = -0.249),size = 0.3,color = "darkgrey", 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_segment(inherit.aes = FALSE, aes(x = -0.25, y = -0.1, xend = -0.75, yend = 0.25), size = 0.3,color = "darkgrey",linetype = "dashed") +
  geom_richtext(data=tibble(x = -0.9,y = 0.35), fill = NA, label.color = NA, label="OTU00005<br>*Peptoniphilus*" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  stat_ellipse(level=0.95,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  geom_point(size = 0.75) +
  coord_fixed(xlim=c(-1, 1.2), ylim = c(-1, 1)) +
  labs(x=glue("PCoA Axis 1 ({loading1}%)"), 
       y=glue("PCoA Axis 2 ({loading2}%)"))+
  scale_color_manual(name="Partition",
                     breaks=c("partition_1", "partition_2", "partition_3"),
                     values=c("#e2e14c", "#8b324d", "#567243"),
                     labels=c("Partition 1", "Partition 2", "Partition 3")) +
  scale_shape_manual(name="Partition",
                     breaks=c("partition_1", "partition_2", "partition_3"),
                     values=c(16, 15, 17),
                     labels=c("Partition 1", "Partition 2", "Partition 3")) +
  scale_fill_manual(name="Partition",
                    breaks=c("partition_1", "partition_2", "partition_3"),
                    values=c("#e2e14c", "#8b324d", "#567243"),
                    labels=c("Partition 1", "Partition 2", "Partition 3")) +
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
        plot.subtitle = element_markdown(size = 8))
ggsave(glue("final_graphs/3_partition_pcoa.pdf"), width = 12, height = 12, unit = "cm")