# -----------------------------------------------------
# Linear-Scatter plot of mutation rate by stage/size
# (Highlighting genes of interest)
# Last Modified: May.15, 2020
# -----------------------------------------------------

library(cancereffectsizeR)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(scales)

# load data
load("cesanal_bystage.RData")  # stage data `kirc_02`
load("cesanal_bysize.RData")  # size data `analysis_size`
load("cesanal_size_med.RData")  # size data split by median value `analysis_size_med`

# Specify genes of interest
highlight_gene <- c("VHL", "PBRM1", "TTN", "SPEN")

# Change character string in scientific notation e.g. from e05 to 10^5
fancy_scientific <- function(l) {
  l <- format(l, scientific = TRUE)
  # turn the 'e+' into plotmath format
  l <- gsub("e","%*%10^", l)
  # return this as an expression
  parse(text=l)
}

# Generate linear scatter plot
jitter_scatter <- function(data, highlight_gene, palette, grp_name, ylim) {
  mutrates <- data
  highlight <- which(mutrates$gene %in% highlight_gene)
  mutrates_plot <- ggplot() +
    geom_jitter(data = mutrates[-highlight,], 
                aes(x=1, y= mutation_rate, color = mutation_rate), 
                shape=16, position=position_jitter(0.05), size = .75, alpha=0.7) +
    scale_colour_gradient(low=brewer.pal(9, palette)[2], 
                          high=brewer.pal(9, palette)[6]) +
    geom_jitter(data = mutrates[highlight,],
                aes(x=1, y= mutation_rate), shape=16,
                position=position_jitter(0.05, seed = 5), size = 3.5, alpha=1) +
    geom_text_repel(data = mutrates[highlight,], 
                    aes(x=1, y= mutation_rate, label = gene), 
                    position=position_jitter(0.05, seed = 5)) +
    ylab("Mutation Rate") + xlab("") + labs(title = grp_name) +
    theme(axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(color = "black"),
          legend.position="none") + 
    scale_y_continuous(labels = fancy_scientific, limits = ylim,
                       breaks = c(2.5e-7, 5e-7, 7.5e-7, 1e-6))+
    coord_flip()
  return(mutrates_plot)
}

# stage 1
mutrates_list <- kirc_02@mutrates_list$`Stage I`

# Build dataframe of genes and their mutation rates
mutrates <- data.frame(gene = names(mutrates_list), 
                       mutation_rate = mutrates_list)

jitter_scatter(mutrates, highlight_gene, "Greens", "Stage I", 
               ylim = c(1e-7, 1.2e-6))

ggsave("KIRC_mutrate_stage1_0429.png", device = "png", units = "in", height = 2, width = 8)

# stage 234
mutrates_list <- kirc_02@mutrates_list$Stages234

# Build dataframe of genes and their mutation rates
mutrates <- data.frame(gene = names(mutrates_list), 
                       mutation_rate = mutrates_list)

jitter_scatter(mutrates, highlight_gene, "Greens", "Stages II, III, IV", 
               ylim = c(1e-7, 1.2e-6))

# mutrates_plot
ggsave("KIRC_mutrate_stage234_0429.png", device = "png", units = "in", height = 2, width = 8)

# size large
mutrates_list <- analysis_size_med@mutrates_list$`Size >5`

# Build dataframe of genes and their mutation rates
mutrates <- data.frame(gene = names(mutrates_list), 
                       mutation_rate = mutrates_list)

jitter_scatter(mutrates, highlight_gene, "Blues", "Size > 5", 
               ylim = c(1e-7, 1.2e-6))

# mutrates_plot
ggsave("KIRC_mutrate_sizemed_large_0515.png", device = "png", units = "in", height = 2, width = 8)


# size small
mutrates_list <- analysis_size_med@mutrates_list$`Size 0-5`

# Build dataframe of genes and their mutation rates
mutrates <- data.frame(gene = names(mutrates_list),
                       mutation_rate = mutrates_list)

jitter_scatter(mutrates, highlight_gene, "Blues", "Size 0 ~ 5", 
               ylim = c(1e-7, 1.2e-6))

# mutrates_plot
ggsave("KIRC_mutrate_sizemed_small_0515.png", device = "png", units = "in", height = 2, width = 8)


