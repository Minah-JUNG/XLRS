# Data Description

##################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

##################################################

dir <- "D:/XLRS/BasicStatistics/"
setwd(dir)

##################################################

list.files()

file.in <- "DatasizeCoverage.txt"
data <- read.table(file.in, sep="\t")

##################################################

summary(data)

data <- data %>%
  as.data.frame() %>%
  dplyr::filter(!(rownames(.) %in% c("IRD0029", "IRD0035", "IRD0057", "IRD0065", "IRD0077")))


data.melt <- reshape2::melt(data) %>% as.data.frame

##------------------------------------------------

a <- ggplot(data, aes(x = DataSize)) +
  geom_density(color = "coral1", fill = "coral1", alpha = 0.5) +
  labs(title = "Density plot of data size (GB)",
       x = "Data size (GB)",
       y = "Density") +
  theme_bw()

##------------------------------------------------

b <- ggplot(data, aes(x = CoverageDepth)) +
  geom_histogram(binwidth = 5, color = "deepskyblue4", fill = "deepskyblue4", alpha = 0.5) +
  labs(title = "Frequency distribution of coverage depth",
       x = "Coverage depth (X)",
       y = "Frequency") +
  scale_y_continuous(breaks = seq(0, 13, 2)) +  # Set breaks to integer values
  theme_bw() 

##------------------------------------------------

c <- ggarrange(a, b, ncol = 1, nrow = 2)

##------------------------------------------------

file.out <- "PlotDensity_histogram.png"
ggsave(file.out, c, width = 7, height = 10, units = "in")


##------------------------------------------------
##################################################
