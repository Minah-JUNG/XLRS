## Network data plot & ORA (Eye-related Cluster, Group C & F)

rm(list = ls())

##################################################

library(dplyr)
library(tidyr)
library(tidyverse)
library(igraph)
library(ggraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)
library(colorspace)
library(officer)
library(rvg)
library(ggtext)

##################################################
## Data
##################################################

dir <- "D:/Documents/"
file_in <- "string_interactions_short.tsv"
groups <- c("GroupC", "GroupF")  

##################################################
## Loop per Group
##################################################

for (grp in groups) {
  cat("Processing:", grp, "\n")
  
  ##-----------------------------
  ## edge & node data
  ##-----------------------------
  
  edge_path <- file.path(dir, grp, file_in)
  data <- read_tsv(edge_path) %>% rename(node1 = `#node1`)
  edges <- data %>% select(Source = node1, Target = node2, CombinedScore = combined_score)
  
  # node
  nodes <- edges %>%
    pivot_longer(cols = c(Source, Target), values_to = "Node") %>%
    distinct(Node)
  
  # edge
  gene_file <- file.path(dir, grp, paste0("ClusterEyeGenes_", grp, ".txt"))
  genelist <- read.table(gene_file)$V1
  all_nodes <- bind_rows(nodes, tibble(Node = genelist)) %>% distinct(Node)
  
  ##-----------------------------
  ## network
  ##-----------------------------
  
  graph <- graph_from_data_frame(edges, vertices = all_nodes, directed = FALSE)
  V(graph)$degree_centrality <- degree(graph)
  V(graph)$fill_color <- rep(brewer.pal(12, "Set3"), length.out = gorder(graph))
  
  graph_plot <- ggraph(graph, layout = 'fr') +
    geom_edge_link(aes(alpha = CombinedScore, colour = CombinedScore, width = CombinedScore)) +
    scale_edge_colour_gradient(low = "grey20", high = "black") +
    scale_edge_width(range = c(0.9, 1.3)) +
    geom_node_point(aes(size = V(graph)$degree_centrality,
                        fill = V(graph)$fill_color),
                    shape = 21, stroke = 0.5, colour = "grey20") +
    geom_node_text(aes(label = names(V(graph))),
                   fontface = "bold", size = 5, repel = TRUE, max.overlaps = Inf) +
    scale_size_continuous(range = c(5, 10)) +
    scale_fill_identity() +
    theme_void() +
    theme(legend.position = 'right') +
    guides(
      colour = guide_colourbar(order = 1, title = "CombinedScore"),
      size = guide_legend(order = 2, title = "DegreeCentrality"),
      edge_width = "none", edge_color = "none", alpha = "none"
    ) +
    labs(size = "DegreeCentrality") +
    ggtitle(gsub("oup", "oup ", grp))
  
  ##-----------------------------
  ## Save network graph
  ##-----------------------------
  
  out_dir <- file.path(dir, grp, "NetworkImage")
  dir.create(out_dir, showWarnings = FALSE)
  
  ggsave(file.path(out_dir, paste0("network_graph_", grp, "_", Sys.Date(), ".svg")),
         graph_plot, width = 5, height = 4, unit = "in")
  ggsave(file.path(out_dir, paste0("network_graph_", grp, "_", Sys.Date(), ".png")),
         graph_plot, width = 5, height = 4, unit = "in", dpi = 300)
  ggsave(file.path(out_dir, paste0("network_graph_", grp, "_", Sys.Date(), ".jpg")),
         graph_plot, width = 5, height = 4, unit = "in", dpi = 300)
  
  # Save to PowerPoint
  editable_plot <- dml(ggobj = graph_plot)
  ppt_path <- file.path(out_dir, "network_plot.pptx")
  
  read_pptx() %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(editable_plot, location = ph_location(left = 1, top = 1, width = 12, height = 7)) %>%
    print(target = ppt_path)
  
  ##-----------------------------
  ## ORA
  ##-----------------------------
  
  ora <- enrichGO(gene = genelist,
                  keyType = "SYMBOL",
                  OrgDb = org.Hs.eg.db,
                  ont = "ALL",
                  pvalueCutoff = 0.5,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  ora_result <- ora@result %>%
    filter(ONTOLOGY %in% c("BP", "MF")) %>%
    arrange(p.adjust) %>%
    head(10)
  
  # Save ORA table
  ora_dir <- file.path(dir, grp, "ORA")
  dir.create(ora_dir, showWarnings = FALSE)
  write.xlsx(ora_result, file.path(ora_dir, paste0("ORAresult_", grp, ".xlsx")), rowNames = FALSE)
  
  if (nrow(ora_result) == 0) {
    message("No significant GO terms found for ", grp)
    next
  }
  
  # GeneRatio 및 logP 계산
  ora_result <- ora_result %>%
    mutate(GeneRatio = as.numeric(sapply(GeneRatio, function(x) eval(parse(text = x)))),
           logP = -log10(p.adjust),
           GeneRatio_Group = case_when(
             GeneRatio <= 0.1 ~ "≤ 0.1",
             GeneRatio <= 0.2 ~ "0.1 - 0.2",
             GeneRatio <= 0.3 ~ "0.2 - 0.3",
             TRUE ~ "> 0.3"
           ),
           GeneRatio_Group = factor(GeneRatio_Group, levels = c("≤ 0.1", "0.1 - 0.2", "0.2 - 0.3", "> 0.3")),
           Description = str_to_sentence(Description),
           label = paste0("<b>", ID, "</b><br>", str_replace_all(str_wrap(Description, 40), "\n", "<br>")),
           label = factor(label, levels = rev(unique(label)))
    )
  
  ##-----------------------------
  ## Barplot
  ##-----------------------------
  
  barplot_data <- ora_result %>%
    head(10) %>%
    mutate(label = paste0("<b style='font-size:14pt;'>", ID, "</b><br><span style='font-size:10pt;'>", Description, "</span>"),
           label = factor(label, levels = rev(unique(label))))
  
  barplot <- ggplot(barplot_data, aes(x = logP, y = label, fill = logP)) +
    geom_vline(xintercept = 1:5, size = 1, color = "grey", linetype = "dotted") +
    geom_col() +
    scale_fill_continuous_sequential(palette = "Heat 2", rev = TRUE) +
    labs(title = paste("Gene Ontology -", grp),
         x = "-log10(p.adjust)", y = NULL, fill = "-log10(p.adjust)") +
    theme_classic() +
    theme(
      axis.text.y = element_markdown(hjust = 0),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      plot.title = element_text(size = 20, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(file.path(out_dir, paste0("ORA_barplot_", grp, "_", Sys.Date(), ".png")),
         barplot, width = 15, height = 10, unit = "in", dpi = 300)
  
  pptx_barplot <- dml(ggobj = barplot)
  read_pptx() %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(pptx_barplot, location = ph_location(left = 1, top = 1, width = 12, height = 7)) %>%
    print(target = file.path(out_dir, paste0("ORA_barplot_", grp, "_", Sys.Date(), ".pptx")))
  
}
