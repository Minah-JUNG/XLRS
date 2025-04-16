## Synonymous vs. Non-synonymous Mutations per Patient

rm(list = ls())

##################################################

library(dplyr)
library(ggplot2)

##################################################
## data
##################################################

result_df <- read.table("Nonsynonymous_Variant_Gene_Counts.txt", header = T)

##################################################
## plot
##################################################

##------------------------------------------------

synonymous_color <- "#E69A8D"
non_synonymous_color <- "#5F4B8B"

##------------------------------------------------

p <- ggplot(result_df, aes(x=Sample, y=Count, fill=MutationType)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("Synonymous"=synonymous_color, "Non-Synonymous"=non_synonymous_color)) +
  labs(title="Synonymous vs. Non-synonymous Mutations per Patient",
       y="n(Genes)", x="Sample", fill="Mutation Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=17),
        axis.text.y = element_text(size=17),
        axis.title.x = element_text(size=19, face="bold"),
        axis.title.y = element_text(size=19, face="bold"),
        plot.title = element_text(size=22, face="bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_rect(color="black", fill=NA, linewidth=1.5),
        legend.position = "bottom",  # 범례 위치 하단
        legend.title = element_text(size=17, face="bold"),
        legend.text = element_text(size=15))

##------------------------------------------------

ggsave("Synonymous_and_NonSynonymous.png", p, width=12, height=6, dpi=300)

##################################################
## Done
##################################################


