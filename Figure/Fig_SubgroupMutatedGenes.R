## Subgroup Mutated Genes
## 2023-12-26
## 2024-06-11: Heatmap - shared mutated genes
## 2024-06-12: Heatmap legend removed

rm(list=ls())

##################################################

library(dplyr)

##################################################

dir <- "D:/majung/작업/XLRS/03 Result/2023 12 26 SubgroupMore"
dir <- "D:/majung/작업/XLRS/03 Result/2023 12 26 SubgroupMore/TableSampleMutation"
setwd(dir)

##################################################

file.in <- list.files(pattern = "Final_")

##################################################
##################################################
##################################################
## mutated genes according to subgroup

data <- list()
for (file in file.in) {
  temp <- read.table(file, check.names = FALSE) %>% dplyr::filter(V5 > 0) # mutation list 정리용 high only
  #temp <- read.table(file, check.names = FALSE) %>% dplyr::filter(V5 > 0 | V6 >0) # high, moderated filter 추가
  temp.genes <- unname(temp[, "V1"]) %>% unique # symbol
  data[[file]] <- temp.genes
}

## Group naming
names(data) <- paste0("Group", LETTERS[1:length(data)])

## GroupJ(outgroup) 제외
data <- data[!grepl("GroupJ", names(data))]

# one-hot encoding
result.dummy <- matrix(0, ncol = length(data), nrow = length(unique(unlist(data))))
rownames(result.dummy) <- unique(unlist(data))
colnames(result.dummy) <- names(data)

for (i in seq_along(data)) {
  result.dummy[match(data[[i]], rownames(result.dummy)), i] <- 1
}

result.dummy_df <- as.data.frame(result.dummy)

##################################################

## Heatmap group-specific genes
## 2024-06-12: legend 삭제

##################################################
##################################################

result.dummy_df <- result.dummy_df %>%
  mutate(count = rowSums(result.dummy_df)) %>%
  arrange(-count)

# annotation (rowSum)
sum_vector <- result.dummy_df$count
annotation_col <- data.frame(
  count = sum_vector
)

#rownames(annotation_col) <- colnames(result.dummy_df)[1:9]

# for heatmap image
result.dummy_df <- result.dummy_df %>% dplyr::select(-count)
result.dummy_df <- as.matrix(result.dummy_df)

# image
dev.off()

png("../TableSubgroupMutation/Heatmap_SubgroupMutationHigh_no_legend.png", width=11, height=5, units = "in", res=1000)

pheatmap(
  t(result.dummy_df), 
  color = c("gray80", "navy"),    
  cluster_rows = FALSE,           
  cluster_cols = FALSE,           
  #legend_breaks = c(0, 1),        
  #legend_labels = c("0: none", "1: mutation"),  
  legend = FALSE,                  # legend 삭제
  angle_col = 90,                  # x축 label 회전
  border_color = "white",          
  main = "Mutated genes according to subgroup of XLRS patients",
  annotation_col = annotation_col  # row Sum 
)

dev.off()


##################################################
##################################################
##################################################
##################################################
