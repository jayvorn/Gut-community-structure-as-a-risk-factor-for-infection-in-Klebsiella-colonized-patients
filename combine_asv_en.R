#!/usr/bin/env Rscript

library(tidyverse)
library(mikropml)
library(plyr)

rds_files <- commandArgs(trailingOnly = TRUE)

taxon_data_results <- map(rds_files, readRDS)

taxon_data_results %>%
  map_dfr(pluck, "feature_importance") %>%
  write_tsv("processed_data/asv_en_feature_importance.tsv")

taxon_data_results %>%
  map(pluck, "performance") %>%
  rbind.fill() %>%
  write_tsv("processed_data/asv_en_performance.tsv")

taxon_data_results %>%
  map(pluck, "trained_model") %>%
  map_dfr(pluck, "results") %>%
  write_tsv("processed_data/asv_en_trained_model.tsv")