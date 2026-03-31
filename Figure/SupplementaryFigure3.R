## Network Plot: Cluster Highlight (GroupC & GroupF)

rm(list = ls())

##################################################

library(readr)
library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(patchwork)

##################################################

dir <- "D:/Documents"

##################################################
## Setting
##################################################

group <- "GroupD"
#group <- "GroupG"
cluster_file <- paste0("network_graph_", group, "_cluster.txt")

## cluster name and color
if (group == "GroupD") {
  cluster_table <- tribble(
    ~cluster, ~Description, ~cluster_colors,
    6, "Extracellular structure and organization", "#9e0142",
    9, "Olfactory sensory processes and signal adaptation mechanisms", "#f46d43",
    1, "Sensory perception of light and photoreceptor cell regulation", "#5e4fa2",
    5, "Regulation of Nuclear transport and RNA localization", "#155a8a",
    2, "Regulation of Morphogenesis in Muscle and Neurons", "#33a02c",
    13, "RNA processing and cytoskeletal-based Motility", "#fb9a99"
  )
} else if (group == "GroupG") {
  cluster_table <- tribble(
    ~cluster, ~Description, ~cluster_colors,
    2, "Extracellular structure and organization", "#9e0142",
    6, "Mechanisms of olfactory sensory perception", "#f46d43",
    9, "Retina development, cilium assembly and light detection", "#5e4fa2",
    7, "Blood coagulation and protein activation cascades", "#ffb300",
    11, "Protein glycosylation and polyamine metabolism", "#a6cee3",
    16, "DNA repair mechanism and cell cycle regulation", "#b2df8a"
  )
}

cluster_table <- cluster_table %>%
  mutate(wrapped_name = paste0(str_wrap(Description, width = 35), "\n"))

##################################################
## Load Interaction Data
##################################################

file_in <- "string_interactions_short.tsv"

path <- file.path(dir, group, file_in)
data <- read_tsv(path, show_col_types = FALSE) %>%
  rename(node1 = `#node1`) %>%
  select(Source = node1, Target = node2, CombinedScore = combined_score)

edges <- data
nodes <- data %>%
  pivot_longer(cols = c(Source, Target), values_to = "Node") %>%
  distinct(Node)

##################################################
## Graph & Cluster Info
##################################################

graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
V(graph)$degree_centrality <- degree(graph)

## cluster information
cluster_info_path <- file.path(dir, group, cluster_file)
cluster_info <- read.table(cluster_info_path, header = FALSE, col.names = "cluster")

## node information
selected_clusters <- cluster_table$cluster
wrapped_cluster_names <- cluster_table$wrapped_name
names(wrapped_cluster_names) <- selected_clusters
cluster_colors <- cluster_table$cluster_colors
names(cluster_colors) <- selected_clusters

##################################################
## Full Network
##################################################

graph_plot_cluster <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(alpha = CombinedScore, colour = CombinedScore, width = CombinedScore)) +
  scale_edge_colour_gradient(low = "grey50", high = "grey40") +
  scale_edge_width(range = c(0.1, 0.5)) +
  geom_node_point(aes(
    size = degree_centrality,
    fill = factor(cluster_info$cluster, levels = selected_clusters)
  ), shape = 21, stroke = 0.2, colour = "grey40") +
  scale_size_continuous(range = c(1, 6)) +
  scale_fill_manual(
    values = cluster_colors,
    breaks = selected_clusters,
    labels = wrapped_cluster_names,
    name = "Clusters"
  ) +
  theme_void() +
  theme(legend.position = "right") +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 4)),
    size = guide_legend(order = 2, title = "DegreeCentrality"),
    edge_alpha = guide_legend(order = 3, title = "CombinedScore")
  ) +
  labs(size = "DegreeCentrality") +
  ggtitle(paste0("Network - ", group, " (Cluster)"))

## save
out_dir <- file.path(dir, group, "ClusterSelected")
dir.create(out_dir, showWarnings = FALSE)

ggsave(file.path(out_dir, paste0("network_graph_", group, "_cluster_", Sys.Date(), ".png")),
       graph_plot_cluster, width = 8, height = 5, dpi = 300)

##################################################
## Subgraph per Cluster
##################################################

plot_list <- list()

for (i in seq_len(nrow(cluster_table))) {
  cluster_id <- cluster_table$cluster[i]
  cluster_color <- cluster_table$cluster_colors[i]
  cluster_title <- cluster_table$wrapped_name[i]
  
  cluster_nodes <- which(cluster_info$cluster == cluster_id)
  subgraph <- induced_subgraph(graph, cluster_nodes)
  
  if (vcount(subgraph) == 0) {
    cat("Cluster", cluster_id, "is empty. Skipping...\n")
    next
  }
  
  V(subgraph)$degree_centrality <- degree(subgraph)
  
  p <- ggraph(subgraph, layout = 'fr') +
    geom_edge_link(aes(alpha = CombinedScore, width = CombinedScore), color = "grey20") +
    scale_edge_width(range = c(0.1, 0.5)) +
    geom_node_point(aes(size = degree_centrality), shape = 21, fill = cluster_color,
                    stroke = 0.2, colour = "grey40") +
    geom_node_text(aes(label = name), size = 2, repel = TRUE, max.overlaps = Inf) +
    scale_size_continuous(range = c(1, 6)) +
    theme_void() +
    ggtitle(cluster_title) +
    theme(plot.title = element_text(size = 12, face = "bold"),
          legend.position = "none")
  
  plot_list[[i]] <- p
  ggsave(file.path(out_dir, paste0("subgraph_cluster_", cluster_id, "_", Sys.Date(), ".png")),
         p, width = 4, height = 3, dpi = 300)
}

# 결합 저장
combined_plot_3col <- wrap_plots(plot_list, ncol = 3)

ggsave(file.path(out_dir, paste0("all_clusters_subgraphs_wide_", Sys.Date(), ".png")),
       combined_plot_3col, width = 10, height = 5, dpi = 300)

cat("Done for ", group, "\n")

##################################################
## Done
##################################################
