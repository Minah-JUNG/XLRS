## Subgroup Mutated Genes (Specific/Common barplot)

rm(list=ls())

##################################################

library(dplyr)
library(ggplot2)

##################################################

dir <- "D:/XLRS/SubgroupMutatedGenes"
setwd(dir)

##################################################

file.in <- "SubgroupMutation_intersection_groupJ.RData"
load(file.in) # result.intersection

##################################################

groups <- paste0("Group", LETTERS[1:9]) # GroupJ 제외

##################################################

group <- "GroupB"
names(result.intersection)

# count
DataCount <- data.frame(Group = character(), Specific = numeric(), Common = numeric(), stringsAsFactors = FALSE)
for (group in groups) {
  GroupSpecific <- result.intersection[[group]] %>% length
  # setdiff(names(result.intersection)[grepl(group, names(result.intersection))], group) GroupA를 포함한 이름 중  'GroupA' 제외
  GroupCommon <- result.intersection[setdiff(names(result.intersection)[grepl(group, names(result.intersection))], group)]
  cat("===", group, "DONE===", "\n", "\n", GroupCommon %>% names, "\n")
  GroupCommon <- result.intersection[setdiff(names(result.intersection)[grepl(group, names(result.intersection))], group)] %>% unlist %>% unique %>% length
  DataCount <- rbind(DataCount, data.frame(Group = group, Specific = GroupSpecific, Common = GroupCommon))
}

##################################################

DataCount$'All' <- result.intersection[paste(groups, collapse =":")] %>% unlist %>% unique %>% length
DataCount$'Multiple' <- DataCount$Common - DataCount$All

##################################################
##################################################

DataCount.melt <- reshape2::melt(DataCount %>% dplyr::select(-"Common"), # "Common" is excluded for higher category
                                 id.vars = "Group")

colnames(DataCount.melt)[2:3] <- c("Categories", "Number") 
DataCount.melt$Genes <- ifelse(DataCount.melt$Categories %in% c('All', 'Multiple'), 'Common', "Specific")

##################################################

color_palette <- rev(RColorBrewer::brewer.pal(3, "Pastel1"))

##################################################

DataCount.melt$Genes <- factor(DataCount.melt$Genes, levels=c("Specific", "Common"))
DataCount.melt$Categories <- factor(DataCount.melt$Categories, levels=c("Multiple", "All", "Specific"))

p <- ggplot(DataCount.melt, aes(x = Genes, y = Number, fill = Categories)) +
  geom_bar(stat = "identity", color = "gray") +
  geom_text(aes(label = Number), position = "stack", vjust = -0.5) + 
  facet_grid(. ~ Group) +
  theme_classic() +
  theme(strip.background = element_rect(colour = NA, fill = NA)) +
  scale_fill_manual(values = color_palette) 
p

##################################################

file.out <- paste0(dir, "SubgroupMutatedGenes_SpecificCommon_GroupJ_wide_stack.png")
ggsave(p, file = file.out, width = 14, height = 5)

##################################################
##################################################
