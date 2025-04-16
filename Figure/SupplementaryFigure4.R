## Network Data Plot (Eye-Related Highlight for subgroup)

rm(list = ls())

##################################################

library(readr)
library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(ggplot2)

##################################################
## Set
##################################################

dir <- "D:/Documents"
file_in <- "string_interactions_merged.tsv"
eye_gene_dir <- "D:/Documents/EyeGeneLists"
out_dir <- file.path(dir, "Figure_Network_EyeRelated_Highlight")
dir.create(out_dir, showWarnings = FALSE)

groups <- paste0("Group", LETTERS[1:9])  # GroupA ~ GroupI

##################################################
## Network
##################################################

for (group in groups) {
  
  cat("Processing:", group, "\n")
  
  ##-----------------------------
  ## Load Data
  ##-----------------------------
  
  file_path <- file.path(dir, group, "Total", file_in)
  if (!file.exists(file_path)) {
    warning("File not found for ", group, " → Skipping.")
    next
  }
  
  data <- read_tsv(file_path, show_col_types = FALSE) %>%
    rename(node1 = `#node1`) %>%
    select(Source = node1, Target = node2, CombinedScore = combined_score)
  
  nodes <- data %>%
    pivot_longer(cols = c(Source, Target), values_to = "Node") %>%
    distinct(Node)
  
  ##-----------------------------
  ## Eye-Related Genes
  ##-----------------------------
  
  eye_gene_file <- file.path(eye_gene_dir, paste0("ClusterEyeGenes_", group, ".txt"))
  if (!file.exists(eye_gene_file)) {
    warning("Eye gene list not found for ", group, " → Skipping.")
    next
  }
  
  genes_eye <- read.table(eye_gene_file, header = FALSE)$V1
  
  ##-----------------------------
  ## Build Graph
  ##-----------------------------
  
  graph <- graph_from_data_frame(d = data, vertices = nodes, directed = FALSE)
  
  V(graph)$degree_centrality <- degree(graph)
  V(graph)$is_eye_gene <- V(graph)$name %in% genes_eye
  
  ##-----------------------------
  ## Plot
  ##-----------------------------
  
  set.seed(2468)
  
  graph_plot <- ggraph(graph, layout = 'graphopt') +
    geom_edge_link(aes(alpha = CombinedScore, colour = CombinedScore, width = CombinedScore)) +
    scale_edge_colour_gradient(low = "grey50", high = "black") +
    scale_edge_width(range = c(0.2, 0.6)) +
    
    geom_node_point(aes(size = degree_centrality,
                        fill = ifelse(is_eye_gene, "red", "grey80")),
                    shape = 21, stroke = 0.2, colour = "grey60", alpha = 0.9) +
    scale_size_continuous(range = c(2, 6)) +
    scale_fill_identity() +
    
    geom_node_text(aes(label = ifelse(is_eye_gene, name, "")),
                   size = 3, fontface = "bold", repel = TRUE, max.overlaps = Inf) +
    
    theme_void() +
    theme(legend.position = "right") +
    guides(
      colour = guide_colourbar(title = "CombinedScore"),
      size = guide_legend(title = "Degree Centrality"),
      edge_width = "none",
      edge_color = "none",
      alpha = "none"
    ) +
    labs(title = paste("Network -", group))
  
  ##-----------------------------
  ## Save
  ##-----------------------------
  
  file_base <- paste0("Figure_Network_EyeRelated_Highlight_", group, "_", Sys.Date())
  ggsave(file.path(out_dir, paste0(file_base, ".svg")), graph_plot, width = 6, height = 5, unit = "in")
  ggsave(file.path(out_dir, paste0(file_base, ".png")), graph_plot, width = 6, height = 5, unit = "in", dpi = 300)
  ggsave(file.path(out_dir, paste0(file_base, ".jpg")), graph_plot, width = 6, height = 5, unit = "in", dpi = 300)
  
  cat("Saved:", file_base, "\n")
}

##################################################
## Done
##################################################