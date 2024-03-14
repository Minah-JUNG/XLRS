## RS1 mutation heatmap

library(pheatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

##################################################

dir <- "D:/XLRS/Heatmap/"
setwd(dir)

##################################################

file.in <- "mutation.txt"
data <- read.table(file.in, check.names = FALSE, sep = "\t", header = T)

##------------------------------------------------

file.in <- "subgroup.txt"
subgroup <- read.table(file.in, check.names = FALSE, sep = "\t", header = T)

##################################################
## check gene name
##################################################

genename <- data$`Variant position` 

##------------------------------------------------

data <- data %>% data.frame(row.names = genename)

##################################################
## dataset 
##################################################

## Select column
data.ann.col <- data %>%
  select(Type, Mutation) 

data.ann.col$'Mutation' <- as.factor(data.ann.col$'Mutation')
data.ann.col$'Type' <- as.factor(data.ann.col$'Type')

##
data.ann.col.color <- RColorBrewer::brewer.pal(n = length(unique(data.ann.col$Mutation)), name = "Dark2")[c(4,2,3,5,6,1,7)] # color rearrangement
names(data.ann.col.color) <- levels(as.factor(data.ann.col$Mutation))

##
data.ann.row <- subgroup %>% 
  arrange()

## 
data.ann.row.color <- RColorBrewer::brewer.pal(n = length(subgroup[,1] %>% unique), name = "Set3")
names(data.ann.row.color) <- levels(as.factor(data.ann.row$Subgroup))

##
data.plot <- data %>%
  select(-matches(samples.removed)) %>%
  select(-c(Variant.position, HGVS, Type, Mutation))

##################################################
## heatmap ordering by group
##################################################

sub.orders <- list()
unique.subgroups <- unique(data.ann.row$Subgroup)

for (subgroup in unique.subgroups) {
  samples <- data.ann.row %>%
    filter(Subgroup == subgroup) %>%
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

## reordering for subgroup
samples.order.2 <- c("IRD0034",
                     "IRD0053",
                     "IRD0023",
                     "IRD0030",
                     "IRD0031",
                     "IRD0006",
                     "IRD0045",
                     "IRD0010",
                     "IRD0044",
                     "IRD0070",
                     "IRD0024",
                     "IRD0037",
                     "IRD0008",
                     "IRD0009",
                     "IRD0036",
                     "IRD0071",
                     "IRD0026",
                     "IRD0027",
                     "IRD0054",
                     "IRD0002",
                     "IRD0042",
                     "IRD0043",
                     "IRD0025",
                     "IRD0028",
                     "IRD0066",
                     "IRD0013",
                     "IRD0061",
                     "IRD0060",
                     "IRD0032",
                     "IRD0059",
                     "IRD0011",
                     "IRD0039",
                     "IRD0058",
                     "IRD0005",
                     "IRD0021",
                     "IRD0001",
                     "IRD0048",
                     "IRD0022",
                     "IRD0078",
                     "IRD0004")


##
data.plot <- data.plot[,samples.order.2]

##################################################
## plot 
##################################################

dev.off()

file.out <- "heatmap_mutation_indel_reorder_position.png"
png(file.out,
    width=800,
    height=700)

pheatmap::pheatmap(data.plot %>% t,
                   annotation_col = data.ann.col[,2, drop=FALSE],
                   annotation_colors = list(Mutation = data.ann.col.color, 
                                            Subgroup = data.ann.row.color),
                   annotation_row = data.ann.row,
                   angle_col = 90,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   legend = FALSE)

dev.off()

##################################################
##################################################
##################################################
##################################################
##################################################


