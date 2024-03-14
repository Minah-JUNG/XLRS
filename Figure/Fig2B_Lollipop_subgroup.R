## lollipop plot

##################################################

library(VariantAnnotation) ## load package for reading vcf file
library(TxDb.Hsapiens.UCSC.hg38.knownGene) ## load package for gene model
library(org.Hs.eg.db) ## load package for gene name
library(rtracklayer)
library(trackViewer)

##################################################

dir <- "D:/XLRS/lollipop/"
setwd(dir)

file.in <- "GenotypeGVCFs_XLRS_chrX.snp_filt_position_selected.vcf"

##################################################

## set the track range
gr <- GRanges("X", IRanges(18639688, 18672108, names="RS1"))

## read in vcf file
vcf <- readVcf(file.in, "hg38")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## get GRanges from VCF object 
mutation.frequency <- rowRanges(vcf)

## keep the metadata
mcols(mutation.frequency) <-
  cbind(mcols(mutation.frequency),
        VariantAnnotation::info(vcf))

## set colors
file.in <- "mutation_subgroup_color.txt"
mut.sub.col <- read.table(file.in, check.names = FALSE)

mutation.frequency$border <- "gray30"
mutation.frequency$subgroup <- mut.sub.col$V2
mutation.frequency$color <- mut.sub.col$V3

## plot Global Allele Frequency based on AC/AN
mutation.frequency$score <- round(mutation.frequency$AF*100)

## change the SNPs label rotation angle
mutation.frequency$label.parameter.rot <- 45

## keep sequence level style same
seqlevelsStyle(gr) <- seqlevelsStyle(mutation.frequency) <- "UCSC"
seqlevels(mutation.frequency) <- "chrX"

## extract transcripts in the range
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                         org.Hs.eg.db, gr=gr)

## check genes
lapply(trs, function(x) c(names(x), x$dat$symbol %>% unique)) %>% as.matrix # RS1: 1, 4

## subset the features to show the interested transcripts only
features <- range(trs[[1]]$dat) # RS1
names(features) <- trs[[1]]$name # define feature level
features$fill <- "lightsteelblue2" # feature color
features$height <- .04 # feature height

## set the legends
file.in <- "subgroup_color.txt"
sub.col <- read.table(file.in, check.names=FALSE)

legends <- list(labels = paste0("Group", LETTERS[1:10]),
                fill = sub.col$V2,
                color = "gray90")
              
##################################################           
## lollipop plot
##################################################

pdf(file = "lollipop_2.pdf", width = 13, height= 5)

lolliplot(mutation.frequency,
          features, 
          ranges=gr,
          #cex=0.7,
          label_on_feature = FALSE,
          legend=legends,
          type="circle")
dev.off()

##################################################
