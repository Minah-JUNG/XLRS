## Subgroup Mutated Genes

rm(list=ls())

##################################################

library(dplyr)
library(biomaRt)

##################################################

dir <- "D:/XLRS/Subgroup/"
setwd(dir)

##################################################

file.in <- list.files(pattern = "Final_")

##################################################
##  mutated genes according to subgroup (HIGH, MODERATED)
##################################################

## creating dataset

data <- list()
for (file in file.in) {
  temp <- read.table(file, check.names = FALSE) %>% dplyr::filter(V5 > 0 | V6 >0) # high, moderated filter 추가
  temp.genes <- unname(temp[, "V1"]) %>% unique # symbol
  data[[file]] <- temp.genes
}

names(data) <- paste0("Group", LETTERS[1:length(data)])
data <- data[!grepl("GroupJ", names(data))] # GroupJ 제외 

##################################################

## common (combination)
result.common <- list()
for (i in 2:length(data)) {
  for (comb in combn(names(data), i, simplify = FALSE)) {
    intersected.genes <- Reduce(intersect, data[comb])
    result.common[[paste(comb, collapse = ":")]] <- intersected.genes
  }
}

## specific
result.specific <- list()
for (group in names(data)) {
  unique.genes <- setdiff(data[[group]], data[!grepl(group, names(data))] %>% unlist %>% unique)
  result.specific[[group]] <- unique.genes
}

## combine
result.intersection <- c(result.specific, result.common)

## save
file.out <- "SubgroupMutation_intersection_groupJ.RData"
save(result.intersection, file = file.out)

##################################################

## statistics: count genes
result.intersection.number <- lapply(result.intersection, length) %>% as.matrix
rn <- rownames(result.intersection.number)
colnames(result.intersection.number) <- "n(genes)"

## save
file.out <- "SubgroupMutatedGenes_intersection_number.txt"
write.table(result.intersection.number, file.out, sep="\t")

##################################################

## statistics
## group-specific only

result.intersection.number <- as.data.frame(result.intersection.number)
result.intersection.number.specific <- result.intersection.number %>% 
  dplyr::filter(!grepl(":", rownames(.)))

##################################################

## genes per group number

subgroup.gene.number <- lapply(data, length) %>% as.matrix %>% data.frame
colnames(subgroup.gene.number) <- "n(Total genes)"

subgroup.gene.number.merge <- data.frame(
  TatalGenes = subgroup.gene.number$`n(Total genes)` %>% unlist %>% as.vector(),
  SpecificGenes = result.intersection.number.specific$`n(genes)` %>% unlist %>% as.vector(),
  row.names = rownames(subgroup.gene.number)) 

file.out <- "SubgroupMutatedGenes_specific_number.txt"
write.table(subgroup.gene.number.merge, file = file.out, sep = "\t")

##################################################
##################################################
##################################################
##################################################
