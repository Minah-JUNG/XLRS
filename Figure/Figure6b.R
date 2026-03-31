## Heatmap of Eye-Related Gene Amino Acid Changes

rm(list = ls())

##################################################

library(pheatmap)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

##################################################
## data
##################################################

dir <- "D:/Documents"
setwd(dir)

file_path <- "Data_EyeRelated_AA.txt"
data <- read_tsv(file_path, show_col_types = FALSE)

## Row order: ROM1, PRPH2, RP1L1
data <- data[c(4, 1:3, 5:nrow(data)),]

## matrix
data_dummy <- data %>%
  select(-c(Gene, Count)) %>%
  column_to_rownames("AA Change") %>%
  as.matrix()

##################################################
## Heatmap
##################################################

##------------------------------------------------
## Annotation 
##------------------------------------------------

## gene information 
annotation_row_base <- data %>%
  select(`AA Change`, Gene) %>%
  column_to_rownames("AA Change")

## color for gene
gene_colors <- c("ROM1" = "#FFDF65", "PRPH2" = "#EA5148", "RP1L1" = "#00BCE1")

## gap for gene
gaps_row <- which(annotation_row_base$Gene[-1] != annotation_row_base$Gene[-nrow(annotation_row_base)])

##------------------------------------------------
## Plot 
##------------------------------------------------

# AA 상태 계산: 한 줄이라도 1 이상이면 “AA”, 아니면 “none”
annotation_row_extended <- annotation_row_base %>%
  mutate(`AA Status` = ifelse(rowSums(data_dummy) > 0, "AA", "none"))

## color 
aa_colors <- c("AA" = "navy", "none" = "#fdfdff")

## file
file_out <- paste0("Heatmap_EyeRelated_AA_base_legend_", Sys.Date(), ".png")

## plot
pheatmap(
  data_dummy,
  filename = file_out,
  cellwidth = 35,
  cellheight = 20,
  width = 8,
  height = 6,
  color = c("#fdfdff", "navy"),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 12,
  fontsize_col = 12,
  legend = FALSE,
  annotation_legend = TRUE,
  angle_col = 90,
  border_color = "grey40",
  main = "Amino Acid Changes in Genes from Eye-Related Clusters",
  display_numbers = FALSE,
  annotation_row = annotation_row_extended,
  annotation_colors = list(
    Gene = gene_colors,
    `AA Status` = aa_colors
  ),
  gaps_row = gaps_row
)

cat("Heatmap saved:", file_out, "\n")

##################################################
## Done
##################################################
