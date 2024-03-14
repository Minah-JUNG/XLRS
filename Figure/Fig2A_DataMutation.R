# Data Description: Mutation

##################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)

##################################################

dir <- "D:/XLRS/BasicStatistics/"
setwd(dir)

##################################################

file.in <- "SampleMutation.txt"
data <- read.table(file.in, sep="\t")

##################################################

data <- data %>% 
  as.data.frame() %>%
  dplyr::select(Mutation) %>%
  dplyr::filter(!Mutation %in% "") %>% # NULL removed
  dplyr::filter(!(rownames(.) %in% c("IRD0029", "IRD0035", "IRD0057", "IRD0065", "IRD0077"))) # filter

data.melt <- reshape2::melt(data)

##------------------------------------------------

ggplot(data.melt, aes(x = Mutation)) +
  geom_bar(aes(fill = Mutation), color = "gray") +
  scale_fill_viridis_d() +  # viridis 
  labs(title = "Mutation", x = "Mutation", y = "Frequency") +
  theme_bw() +
  theme(legend.position = "none")

file.out <- "PlotMutation_filter_viridis.png"
ggsave(file.out, width=12, height=4, units = "in")


##################################################
##################################################
##################################################
##################################################
##################################################
