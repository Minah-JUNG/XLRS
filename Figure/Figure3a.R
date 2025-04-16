## network and GO enrichment for common and Group C

##################################################

rm(list = ls())

##################################################

library(readr)
library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggtext)
library(clusterProfiler)
library(org.Hs.eg.db)

##################################################

dir <- "D:/Documents"
file_in <- "string_interactions_short.tsv"
group_list <- c("Common", "GroupC")
date_tag <- Sys.Date()

##################################################
## Functions
##################################################

##------------------------------------------------
## Function: Load interaction data
##------------------------------------------------

load_network_data <- function(group) {
  path <- if (group == "Common") {
    file.path(dir, group, shell, file_in)
  } else {
    file.path(dir, group, "Total", shell, file_in)
  }
  
  data <- read_tsv(path, show_col_types = FALSE)
  colnames(data)[colnames(data) == "#node1"] <- "node1"
  
  edges <- data %>% select(Source = node1, Target = node2, CombinedScore = combined_score)
  nodes <- edges %>% pivot_longer(cols = c(Source, Target), values_to = "Node") %>% distinct(Node)
  
  list(edges = edges, nodes = nodes)
}

##------------------------------------------------
## Function: Create graph plot
##------------------------------------------------

plot_network <- function(edges, nodes, color_group = NULL, group_name = "Group") {
  graph <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  V(graph)$degree_centrality <- degree(graph)
  
  set.seed(2468)
  g <- ggraph(graph, layout = 'graphopt') +
    geom_edge_link(aes(alpha = CombinedScore, colour = CombinedScore, width = CombinedScore)) +
    scale_edge_colour_gradient(low = "grey50", high = "grey40") +
    scale_edge_width(range = c(0.1, 0.5)) +
    geom_node_point(
      aes(size = degree_centrality, fill = color_group, alpha = 0.95),
      shape = 21, stroke = 0.2, colour = "grey70"
    ) +
    scale_size_continuous(range = c(1, 6)) +
    scale_fill_manual(values = c("Common" = "grey30", "Group C" = "salmon")) +
    theme_void() +
    theme(legend.position = "right") +
    guides(
      fill = guide_legend(order = 1, title = "Category", override.aes = list(size = 3)),
      size = guide_legend(order = 2, title = "Degree Centrality"),
      colour = guide_colourbar(order = 3, title = "Combined Score"),
      edge_width = "none",
      edge_color = "none",
      alpha = "none"
    ) +
    ggtitle(paste0("Network - ", gsub("oup", "oup ", group_name)))
  
  list(graph = graph, plot = g)
}

##------------------------------------------------
## Function: Perform ORA
##------------------------------------------------

perform_ora <- function(gene_list, group_name) {
  ora <- enrichGO(
    gene = gene_list,
    keyType = "SYMBOL",
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pvalueCutoff = 0.5,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.5,
    readable = TRUE
  )
  
  res <- ora@result %>%
    filter(ONTOLOGY %in% c("BP", "MF")) %>%
    arrange(p.adjust) %>%
    head(15)
  
  if (nrow(res) == 0) return(NULL)
  
  res <- res %>%
    mutate(
      GeneRatio = as.numeric(sapply(GeneRatio, function(x) eval(parse(text = x)))),
      GeneRatio_Group = case_when(
        GeneRatio < 0.1 ~ "< 0.1",
        GeneRatio <= 0.2 ~ "0.1 - 0.2",
        GeneRatio <= 0.3 ~ "0.2 - 0.3",
        TRUE ~ "> 0.3"
      ),
      Description = stringr::str_to_sentence(Description),
      label = paste0("<b>", ID, "</b><br>", str_replace_all(str_wrap(Description, width = 45), "\n", "<br>")),
      logP = -log10(p.adjust)
    )
  
  return(res)
}

##------------------------------------------------
# Function: Draw lollipop dotplot
##------------------------------------------------

plot_dotplot <- function(ora_df, group_name) {
  x_max <- max(ora_df$logP, na.rm = TRUE) * 1.1
  base_plot <- ggplot(ora_df, aes(x = logP, y = reorder(label, logP), color = GeneRatio_Group)) +
    geom_segment(aes(x = 0, xend = logP, yend = label), size = 1.8, color = "gray50") +
    geom_point(size = 8) +
    scale_color_manual(values = c(
      "< 0.1" = "#3288bd", 
      "0.1 - 0.2" = "#99d594",
      "0.2 - 0.3" = "#fee08b", 
      "> 0.3" = "#d53e4f"
    )) +
    scale_x_continuous(limits = c(0, x_max), expand = c(0, 0)) +
    labs(
      title = paste("Gene Ontology -", gsub("oup", "oup ", group_name)),
      x = "-log10(p.adjust)", y = NULL, color = "GeneRatio"
    ) +
    theme_classic() +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      plot.title = element_text(size = 22, face = "bold"),
      axis.text.y = ggtext::element_markdown(size = 21),
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 20),
      panel.grid = element_blank(),
      axis.line.x = element_line(linewidth = 2.0, color = "grey80"),
      axis.line.y = element_line(linewidth = 2.0, color = "grey80")
    )
  
  base_plot + geom_vline(xintercept = ggplot_build(base_plot)$layout$panel_params[[1]]$x$breaks,
                         linetype = "dotted", color = "gray80")
}

##################################################
## Run
##################################################

for (group in group_list) {
  cat("Processing:", group, "\n")
  
  data <- load_network_data(group)
  if (group == "GroupC") {
    # compare to Common
    data_common <- load_network_data("Common")
    new_nodes <- setdiff(data$nodes$Node, data_common$nodes$Node)
    data$nodes <- data$nodes %>%
      mutate(color_group = factor(ifelse(Node %in% new_nodes, "Group C", "Common"),
                                  levels = c("Common", "Group C")))
  } else {
    data$nodes <- data$nodes %>% mutate(color_group = "Common")
  }
  
  # Plot network
  plot_data <- plot_network(data$edges, data$nodes, color_group = data$nodes$color_group, group_name = group)
  
  # Save network
  dir_plot <- file.path(dir, group, "Total", "ClusterGO")
  dir.create(dir_plot, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(dir_plot, paste0("network_graph_", group, "_", date_tag, ".png")),
         plot_data$plot, width = 6, height = 4)
  ggsave(file.path(dir_plot, paste0("network_graph_", group, "_", date_tag, ".svg")),
         plot_data$plot, width = 6, height = 4)
  
  # ORA
  geneset <- plot_data$graph %>% igraph::V() %>% names()
  ora_df <- perform_ora(geneset, group)
  if (!is.null(ora_df)) {
    plot_ora <- plot_dotplot(ora_df, group)
    
    ora_dir <- file.path(dir_plot, "ORA_Results")
    dir.create(ora_dir, recursive = TRUE, showWarnings = FALSE)
    
    write.table(ora_df, file.path(ora_dir, paste0("ORA_result_", group, "_", date_tag, ".txt")),
                sep = "\t", row.names = FALSE)
    
    ggsave(file.path(ora_dir, paste0("Lollipop_", group, "_", date_tag, ".svg")),
           plot_ora, width = 10.5, height = 13)
    ggsave(file.path(ora_dir, paste0("Lollipop_", group, "_", date_tag, ".png")),
           plot_ora, width = 10.5, height = 13, dpi = 300)
    ggsave(file.path(ora_dir, paste0("Lollipop_", group, "_", date_tag, ".jpg")),
           plot_ora, width = 10.5, height = 13, dpi = 300)
  } else {
    cat("No GO terms found for", group, "\n")
  }
}

##################################################
## Done
##################################################
