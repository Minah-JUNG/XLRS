## Eye-related Cluster GO terms (top 5)

rm(list = ls())

##################################################

library(dplyr)
library(tidyr)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(stringr)

##################################################
## data
##################################################

dir <- "D:/Documents"
groups <- paste0("Group", LETTERS[1:9])
file_in <- "network_graph_multi_Cluster.xlsx"

##################################################
## Extract Eye-related GO Terms (Top 5 per cluster)
##################################################

combined_data <- data.frame()

for (grp in groups) {
  cluster_dir <- file.path(dir, "Cluster")
  setwd(cluster_dir)
  
  for (sheet in excel_sheets(file_in)) {
    data <- read_excel(file_in, sheet = sheet)
    
    if (any(grepl("visual", data, ignore.case = TRUE))) {
      eye_terms <- data %>%
        filter(ONTOLOGY != "CC") %>%
        mutate(Group = grp, Cluster = sheet) %>%
        arrange(p.adjust) %>%
        slice_head(n = 5)
      
      message("Detected eye-related cluster: ", grp, " / ", sheet)
      
      combined_data <- combined_data %>%
        bind_rows(eye_terms) %>%
        distinct(ONTOLOGY, ID, Description, Group, Cluster, .keep_all = TRUE)
    }
  }
}

##################################################
## Heatmap data
##################################################

heatmap_data <- combined_data %>%
  filter(ONTOLOGY != "CC") %>%
  mutate(log_p_adjust = -log10(p.adjust),
         Group = str_replace(Group, "oup", "oup ")) %>%  # GroupA → Group A
  select(Description, Group, log_p_adjust) %>%
  pivot_wider(names_from = Group, values_from = log_p_adjust, values_fill = NA)

heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) <- str_to_sentence(heatmap_data$Description)

##################################################
## Heatmap
##################################################

output_path <- file.path(dir, paste0("Heatmap_EyeRelatedCluster_padjust_", Sys.Date(), ".png"))

png(output_path, width = 8, height = 7, units = "in", res = 300)

pheatmap(heatmap_matrix,
         cellwidth = 25,
         cellheight = 20,
         color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 11,
         fontsize_col = 11,
         legend = TRUE,
         angle_col = 90,
         border_color = "grey",
         main = "GO Terms of Eye-related Gene Clusters",
         na_col = "white")

grid.text("-log10(p.adjust)", x = 0.97, y = 0.6, rot = 90,
          gp = gpar(fontsize = 10, fontface = "bold"))

dev.off()

message("Heatmap saved at: ", output_path)

##################################################
## Done
##################################################

