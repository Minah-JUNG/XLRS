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
## mutated genes according to subgroup (HIGH)
##################################################

data <- list()
for (file in file.in) {
  temp <- read.table(file, check.names = FALSE) %>% dplyr::filter(V5 > 0) # variant effect 'HIGH'
  temp.genes <- unname(temp[, "V1"]) %>% unique # symbol
  data[[file]] <- temp.genes
}

## Group naming
names(data) <- paste0("Group", LETTERS[1:length(data)])

## Except GroupJ
data <- data[!grepl("GroupJ", names(data))]

## one-hot encoding
result.dummy <- matrix(0, ncol = length(data), nrow = length(unique(unlist(data))))
rownames(result.dummy) <- unique(unlist(data))
colnames(result.dummy) <- names(data)

for (i in seq_along(data)) {
  result.dummy[match(data[[i]], rownames(result.dummy)), i] <- 1
}

result.dummy_df <- as.data.frame(result.dummy)

## result
gene.groups <- lapply(rownames(result.dummy_df), function(gene) {
  groups <- colnames(result.dummy_df)[result.dummy_df[gene, ] == 1]
  return(data.frame(Gene = gene, Count = length(groups), Groups = toString(groups)))
})

result <- do.call(rbind, gene.groups) %>% arrange(desc(Count))

file.out <- "SubgroupMutatedGenes_High_except_groupJ.txt" # Table 3
write.table(result, file.out, sep="\t")

##################################################
## Mutated genes with phenotype
##################################################

geneset <- result %>% dplyr::select(Gene) %>% unlist %>% as.vector # mutated gene list

##################################################

## DB
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

## Gene information
gene.info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description", "phenotype_description", "go_id"),
                   filters = "hgnc_symbol",
                   values = geneset,
                   mart = ensembl)

## multiple GO term
gene.info.go <- gene.info %>%
  group_by(ensembl_gene_id) %>%
  summarize(hgnc_symbol = unique(hgnc_symbol),
            description = unique(description), 
            phenotype_description = paste(unique(phenotype_description), collapse=","),
            GOterm = paste(go_id, collapse=","))

## select
gene.info.select <- gene.info.go[, colnames(gene.info.go) %in% c("hgnc_symbol", 
                                                                 "description", 
                                                                 "phenotype_description")]


## modification: removing after '[' 
gene.info.select$description <- gsub("\\[.*", "", gene.info.select$description)

## removing duplication (ensemblID)
gene.info.select <- gene.info.select[!duplicated(gene.info.select$hgnc_symbol),]
colnames(gene.info.select)[colnames(gene.info.select) == "hgnc_symbol"] <- "Gene"

## arrange
result.pn <- result
result.pn <- inner_join(result.pn, gene.info.select, by="Gene")

## save
file.sv <- "SubgroupMutatedGenes_description_phenotype_High_GroupJ.txt"
write.table(result.pn, file.sv, sep="\t")

##################################################
##################################################
##################################################
##################################################
