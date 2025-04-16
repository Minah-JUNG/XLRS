## RS1 mutation heatmap

rm(list=ls())

##################################################

library(pheatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(grid)

##################################################

dir <- "D:/Documents/"
setwd(dir)

##################################################

file.in <- "mutation.txt"
data <- read.table(file.in, check.names = FALSE, sep = "\t", header = T)
data$`Variant position` <- gsub("X_", "ChrX:", data$`Variant position`)

##------------------------------------------------

file.in <- "subgroup.txt"
subgroup <- read.table(file.in, check.names = FALSE, sep = "\t", header = T)

##################################################
## change gene names
##################################################

genename <- data$`Variant position`

##------------------------------------------------

data <- data %>% data.frame(row.names = genename)

##################################################
## dataset 
##################################################

###### heatmap ordering by group

data.plot <- data

sub.orders <- list()
unique.subgroups <- unique(subgroup$Subgroup)

for (sg in unique.subgroups) {
  samples <- subgroup %>%
    filter(Subgroup == sg) %>%
    rownames()
  if (length(samples) == 1) {
    sub.orders[[subgroup]] <- samples
  } else {
    data.sub <- data.plot %>% select(samples)
    hc.sub <- hclust(dist(data.sub %>% t))
    sub.orders[[subgroup]] <- samples[hc.sub$order]
  }
}

samples.order <- sub.orders %>% unlist %>% as.vector

##------------------------------------------------

data.plot <- data.plot[,samples.order]

##################################################

dev.off()

##------------------------------------------------

data.ann.col <- data %>% select(Type, Mutation) 

data.ann.col$'Mutation' <- as.factor(data.ann.col$'Mutation')
data.ann.col$'Type' <- as.factor(data.ann.col$'Type')

##
data.ann.col.color <- RColorBrewer::brewer.pal(n = length(unique(data.ann.col$Mutation)), name = "Accent") #Set1
names(data.ann.col.color) <- levels(as.factor(data.ann.col$Mutation))

##------------------------------------------------


file.out <- paste0("heatmap_mutation_", Sys.Date(), ".png")
png(file.out,
    width=1200,
    height=700)

##------------------------------------------------

color_palette <- colorRampPalette(c("#fdfdff", "#34495e"))(10)

##------------------------------------------------

pheatmap::pheatmap(data.plot,
                   annotation_col = data.ann.row,
                   annotation_colors = list(Subgroup = data.ann.row.color),
                   angle_col = 90,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   legend = FALSE,
                   fontsize = 16,        
                   fontsize_row = 16,    
                   fontsize_col = 18,
                   fontface_row = "bold",  # For bold row labels
                   fontface_col = "bold",  # For bold column labels
                   color = color_palette,
                   border_color = "grey70"
)

dev.off()

##################################################
## Done
##################################################



