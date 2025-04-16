## ORA 

rm(list = ls())

##################################################

library(clusterProfiler)
library(org.Hs.eg.db)

##################################################

ora <- enrichGO(gene = genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

##################################################

ora_result <- ora@result

##################################################
## Done
##################################################