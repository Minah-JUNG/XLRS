## Network + Clustering(Louvain) + ORA per Group 

##################################################

rm(list = ls())

##################################################

library(readr)
library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)

##################################################

dir <- "D:/Documents"
file_in <- "string_interactions_short.tsv"
groups <- paste0("Group", LETTERS[1:9])

#--------------------------------------------------

for (grp in groups) {
  
  cat("== Processing", grp, "\n")
  
  #----- Load edges -----
  file_path <- file.path(dir, grp, file_in)
  if (!file.exists(file_path)) {
    cat("== File missing:", file_path, "\n")
    next
  }
  
  data <- read_tsv(file_path, show_col_types = FALSE) %>%
    rename(node1 = `#node1`) %>%
    select(Source = node1, Target = node2, Weight = combined_score)
  
  #----- Extract nodes -----
  nodes <- data %>%
    pivot_longer(cols = c(Source, Target), values_to = "Node") %>%
    distinct(Node)
  
  #----- Build graph -----
  graph <- graph_from_data_frame(d = data, vertices = nodes, directed = FALSE)
  V(graph)$degree_centrality <- degree(graph)
  
  #----- clustering -----
  set.seed(1234)  
  louvain <- cluster_louvain(graph)
  V(graph)$cluster <- louvain$membership
  cluster_info <- data.frame(gene = V(graph)$name, cluster = louvain$membership)
  
  #----- Save cluster info -----
  out_dir <- file.path(dir, grp, "ClusterGO")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  cluster_file <- paste0("network_graph_cluster_", grp, ".txt")
  write.table(cluster_info, file = file.path(out_dir, cluster_file), sep = "\t", row.names = FALSE)
  
  #----- Save cluster-colored network -----
  cluster_colors <- rainbow(length(unique(louvain$membership)))
  V(graph)$color_cluster <- cluster_colors[louvain$membership]
  
  graph_plot <- ggraph(graph, layout = 'fr') +
    geom_edge_link(aes(alpha = Weight, colour = Weight, width = Weight)) +
    scale_edge_colour_gradient(low = "grey50", high = "grey40") +
    scale_edge_width(range = c(0.3, 1)) +
    geom_node_point(aes(size = degree_centrality, fill = color_cluster),
                    shape = 21, stroke = 0.3, colour = "grey60") +
    scale_size_continuous(range = c(2, 8)) +
    scale_fill_identity() +
    theme_void() +
    theme(legend.position = "none") +
    labs(title = paste(grp, "clustering"))
  
  file_graph <- paste0("network_graph_", grp, "_cluster.svg")
  ggsave(filename = file.path(out_dir, file_graph), plot = graph_plot, width = 8, height = 6)
  
  cat("== Network graph saved:", file_graph, "\n")
  
  #----- ORA for each cluster (only if > 10 genes) -----
  wb <- createWorkbook()
  cluster_ids <- names(which(table(louvain$membership) > 10))
  
  for (cid in cluster_ids) {
    genes <- V(graph)$name[louvain$membership == as.integer(cid)]
    
    ora <- enrichGO(gene = genes,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db,
                    ont = "ALL",
                    pvalueCutoff = 0.5,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)
    
    result <- ora@result
    if (nrow(result) == 0) {
      cat("== No GO terms for cluster", cid, "\n")
      next
    }
    
    sheet_name <- paste0("Cluster", cid)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, result)
    freezePane(wb, sheet_name, firstRow = TRUE)
  }
  
  file_out_xlsx <- paste0("network_graph_multi_cluster_function_", grp, ".xlsx")
  saveWorkbook(wb, file = file.path(out_dir, file_out_xlsx), overwrite = TRUE)
  
  cat("== ORA results saved:", file_out_xlsx, "\n\n")
}

cat("=== Done with all groups ===\n")

##################################################
## Done
##################################################
