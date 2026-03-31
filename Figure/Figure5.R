## network data with clusters for Group B

##################################################

rm(list = ls())

##################################################

library(readr)
library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(stringr)
library(ggtext)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)

##################################################

dir <- "D:/Documents"
group <- "GroupC"
file_in <- "string_interactions_short.tsv"
cluster_file <- "network_graph_cluster_GroupC.txt"

##################################################

selected_clusters <- c(1, 6, 11, 10, 3, 13) # top 6 clusters
cluster_names <- paste("Cluster", letters[1:6])
cluster_colors <- c("#9e0142", "#f46d43", "#fee08b", "#66c2a5", "#3288bd", "#5e4fa2")
wrapped_cluster_names <- str_wrap(cluster_names, width = 35)
names(wrapped_cluster_names) <- selected_clusters
names(cluster_colors) <- selected_clusters

##------------------------------------------------

data_path <- file.path(dir, group, file_in)
cluster_info_path <- file.path(dir, group, cluster_file)
out_dir <- file.path(dir, group, "ClusterSelected")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

##------------------------------------------------

data <- read_tsv(data_path, show_col_types = FALSE) %>%
  rename(node1 = `#node1`) %>%
  select(Source = node1, Target = node2, CombinedScore = combined_score)

##------------------------------------------------

nodes_all <- data %>%
  pivot_longer(cols = c(Source, Target), values_to = "Node") %>%
  distinct(Node)

##------------------------------------------------

cluster_info <- read.table(cluster_info_path, header = TRUE)
filtered_nodes <- cluster_info %>% filter(cluster %in% selected_clusters)

##------------------------------------------------

nodes <- nodes_all %>%
  filter(Node %in% filtered_nodes$gene) %>%
  left_join(filtered_nodes, by = c("Node" = "gene"))

##------------------------------------------------

filtered_edges <- data %>% filter(Source %in% nodes$Node & Target %in% nodes$Node)

##################################################
## network graph
##################################################

graph <- graph_from_data_frame(d = filtered_edges, vertices = nodes, directed = FALSE)
V(graph)$degree_centrality <- degree(graph)

##------------------------------------------------

set.seed(1357)
graph_plot <- ggraph(graph, layout = 'fr') +
  geom_edge_link(aes(alpha = CombinedScore, colour = CombinedScore, width = CombinedScore)) +
  scale_edge_colour_gradient(low = "grey50", high = "grey40") +
  scale_edge_width(range = c(0.1, 0.5)) +
  geom_node_point(aes(size = degree_centrality, fill = factor(cluster, levels = selected_clusters)),
                  shape = 21, stroke = 0.2, colour = "grey40") +
  scale_size_continuous(range = c(1, 6)) +
  scale_fill_manual(values = cluster_colors, labels = wrapped_cluster_names, name = "Clusters") +
  theme_void() +
  theme(legend.position = 'right') +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 4)),
    size = guide_legend(order = 2, title = "DegreeCentrality"),
    edge_alpha = guide_legend(order = 3, title = "CombinedScore")
  ) +
  ggtitle(paste0("Network - ", group))

##------------------------------------------------

ggsave(file.path(out_dir, paste0("network_graph_", group, "_cluster_only_", Sys.Date(), ".png")),
       plot = graph_plot, width = 8, height = 5, dpi = 300)

cat("== Network image saved.\n\n")

##################################################
## ORA for each selected cluster
##################################################

dotplot_list <- list()

for (k in seq_along(selected_clusters)) {
  cluster_id <- selected_clusters[k]
  cluster_genes <- V(graph)$name[cluster_info$cluster == cluster_id]
  
  if (length(cluster_genes) == 0) {
    next
  }
  
  ora <- enrichGO(gene = cluster_genes,
                  keyType = "SYMBOL",
                  OrgDb = org.Hs.eg.db,
                  ont = "ALL",
                  pvalueCutoff = 0.5,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  ora_filtered <- ora@result %>%
    filter(ONTOLOGY %in% c("BP", "MF")) %>%
    arrange(p.adjust) %>%
    head(7)
  
  if (nrow(ora_filtered) == 0) next
  
  ora_filtered <- ora_filtered %>%
    mutate(GeneRatio = as.numeric(sapply(GeneRatio, function(x) eval(parse(text = x)))),
           GeneRatio_Group = case_when(
             GeneRatio < 0.1 ~ "< 0.1",
             GeneRatio <= 0.2 ~ "0.1 - 0.2",
             GeneRatio <= 0.3 ~ "0.2 - 0.3",
             TRUE ~ "> 0.3"
           ),
           GeneRatio_Group = factor(GeneRatio_Group, levels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.3", "> 0.3")),
           Description = str_to_sentence(Description),
           label = paste0("<b>", ID, "</b><br>", str_replace_all(str_wrap(Description, width = 40), "\n", "<br>")),
           logP = -log10(p.adjust))
  
  x_max <- max(ora_filtered$logP) * 1.1
  
  p <- ggplot(ora_filtered, aes(x = logP, y = reorder(label, logP), color = GeneRatio_Group)) +
    geom_segment(aes(x = 0, xend = logP, yend = label), size = 1.8, color = "gray50") +
    geom_point(size = 8) +
    scale_color_manual(values = c("< 0.1" = "#3288bd", 
                                  "0.1 - 0.2" = "#99d594", 
                                  "0.2 - 0.3" = "#fee08b", 
                                  "> 0.3" = "#d53e4f")) +
    labs(title = paste("Gene Ontology - Cluster", letters[k]), x = "-log10(p.adjust)", y = NULL) +
    xlim(0, x_max) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 20, face = "bold"),
      axis.text.y = ggtext::element_markdown(size = 16),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16)
    ) +
    geom_vline(xintercept = seq(1, x_max, 1), linetype = "dotted", color = "gray80")
  
  dotplot_list[[length(dotplot_list) + 1]] <- p
}

##------------------------------------------------

if (length(dotplot_list) > 0) {
  final_plot <- wrap_plots(dotplot_list, ncol = 3, nrow = 2) &
    theme(legend.position = "bottom")
  
  ggsave(file.path(out_dir, paste0("Dotplot_GO_2x3_", group, "_", Sys.Date(), ".png")),
         final_plot, width = 25, height = 18, units = "in")
  
  cat("== ORA plot saved.\n\n")
  
}

##################################################
## Done
##################################################
