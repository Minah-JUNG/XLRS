## UpSet Diagram for Eye-Related Genes 

rm(list = ls())

##################################################

library(dplyr)
library(tidyverse)
library(UpSetR)

##################################################
## Data
##################################################

dir <- "D:/Documents"
setwd(dir)

file_in <- "EyeRelatedGeneList_perGroups.txt"
data <- read.table(file_in, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

##################################################
## Preprocessing
##################################################

## character to vector
gene_list <- data %>%
  mutate(Genes = strsplit(Genes, ",\\s*")) %>%
  deframe()  # list 형태로 변환

## geneset
common_genes <- gene_list$Common
gene_list <- lapply(gene_list, function(x) c(x, common_genes))
gene_list <- gene_list[!names(gene_list) %in% "Common"]

## binary
gene_data <- fromList(gene_list) %>% as.data.frame()

##################################################
## UpSet Diagram
##################################################

## file
file_png <- "UpSet_Diagram_EyeRelatedGenes.png"

## PNG
png(file_png, width = 3200, height = 1600, res = 150)
UpSetR::upset(
  gene_data,
  nsets = length(names(gene_list)),
  order.by = "degree",
  line.size = 2,
  main.bar.color = "gray40",
  sets = rev(names(gene_list)),
  sets.bar.color = "steelblue",
  point.size = 5,
  text.scale = 3,
  keep.order = TRUE,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  mb.ratio = c(0.45, 0.55)
)
dev.off()

cat("UpSet diagram saved to:\n", file_png, "\n")

##################################################
## Done
##################################################
