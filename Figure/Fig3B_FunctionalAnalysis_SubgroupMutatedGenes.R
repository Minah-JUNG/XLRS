## Functional Analysis
## Subgroup Mutated Genes

library(dplyr)
library(clusterProfiler)
library(GO.db)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr) ## ggarrange

##################################################

dir <- "D:/XLRS/FunctionalAnalysis/"
setwd(dir)

##################################################

list.files()

file.in <- "./TableSubgroupMutation/SubgroupMutation_intersection_groupJ.RData"
load(file.in) # result.intersection

##################################################

names(result.intersection)[[length(result.intersection)]]
geneset <- result.intersection[[length(result.intersection)]] # "GroupA:GroupB:GroupC:GroupD:GroupE:GroupF:GroupG:GroupH:GroupI"

##################################################

ego2 <- enrichGO(gene         = geneset,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = "SYMBOL",
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)

##################################################

result <- ego2@result

file.out <- "Result_ORA_except_GroupJ.txt"
write.table(result, file.out, sep="\t")

##################################################

a <- barplot(ego2) + theme_bw()
b <- dotplot(ego2) + theme_bw()

c <- ggpubr::ggarrange(a, b)

file.out <- "SubgroupMutatedGenes_function_plot_HighModerate_GroupJ.png"
ggsave(c, file = file.out, width = 10, height = 5, units = "in")

##################################################