# -----------------------------------------------
# Scatter plot for gene mutation rate
# (by tumor size/stage)
# Last Modified: Apr.21, 2020
# -----------------------------------------------

library(ggplot2)
library(ggrepel)

load("cesanal_bysize.RData")  # size data
# load("cesanal_bystage.RData") # stage data

mut_rate <- data.frame(gene = names(analysis_size@mutrates_list[[1]]),
                       size_large = analysis_size@mutrates_list[[1]],
                       size_small = analysis_size@mutrates_list[[2]])

scatterplot_gene_mutrate <- ggplot(data = mut_rate, aes(x=size_large, y=size_small))+
  geom_point(size=1.5, color="steelblue", alpha=0.4) +
  ggtitle("KIRC - Gene-level mutation rates") +
  labs(subtitle = "covariate file: kirc_pca") +
  xlab("Mutation rate in tumors (size > 3)") +
  ylab("Mutation rate in tumors (size 0-3)") +
  geom_smooth(method="lm", color="navyblue") +
  geom_abline(slope=1, intercept=0, color = "darkred", linetype = "dashed") +
  theme_bw() + xlim(3e-7, 1.5e-6) + ylim(3e-7, 1.5e-6) +
  geom_text_repel(data=subset(mut_rate, size_large > 0.0000015 | size_small > 0.000005), aes(x=size_large, y=size_small, label=gene), size=2.5)

scatterplot_gene_mutrate
