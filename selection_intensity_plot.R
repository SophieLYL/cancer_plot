# -----------------------------------------------
# Barplot for Selection Intensity (by stage/size) 
# Last Modified: Apr.21, 2020
# -----------------------------------------------

library(dplyr)
library(ggplot2)

# ---- Bar plot for selection intensity ----

# Data preparation
load("cesanal_bysize.RData")

gene_selection <- data.frame(varient = analysis_size@selection_results$variant,
                             gene = analysis_size@selection_results$gene,
                             selection_intensity = analysis_size@selection_results$selection_intensity,
                             size = analysis_size@selection_results$progression)
size_large <- gene_selection %>% filter(size == "Size >3")
size_small <- gene_selection %>% filter(size == "Size 0-3")

selection_filtered <- size_large %>% filter(selection_intensity > 1e4)

selection_filtered$varient <- 
  reorder(selection_filtered$varient, -selection_filtered$selection_intensity)
selection_filtered$gene <- 
  reorder(selection_filtered$gene, -selection_filtered$selection_intensity)

# target genes for highlight
target <- c(c("VHL", "HIF1A", "HIF3A"), 
            as.vector(head(arrange(selection_filtered, desc(selection_intensity))$gene, 10)))

selection_filtered <- selection_filtered %>%
  mutate(label = ifelse(gene %in% target, gene, NA)) %>%
  mutate(label = ifelse(is.na(label), label, levels(gene)[label])) %>%
  arrange(desc(selection_intensity))

# generate bar plot

gg_bar <- 
  ggplot() + 
  geom_bar(data = selection_filtered, 
           aes(x = varient, y = selection_intensity),
           stat = "identity", fill = "black") +
  geom_bar(data = filter(selection_filtered, !is.na(label)), 
           aes(x = varient, y = selection_intensity, fill = gene),
           stat = "identity") +
  geom_point(data = filter(selection_filtered, !is.na(label)), 
             aes(x = varient, y = selection_intensity), shape = 20) +
  theme_classic() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = unit(c(1, -1, 1, 3.7), "mm"),
    legend.position = c(0.6, 0.85),
    legend.title = element_blank()
  ) +
  geom_text(
    data = filter(selection_filtered, !is.na(label)), 
    aes(x = varient, y = selection_intensity, label = as.character(varient), 
        color = gene),
    nudge_y = 0.3, nudge_x = 1.5,
    angle = 60,
    size = 2.25
  ) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  labs(y = bquote("Selection intensity"))

gg_bar  # display


# ----------------------

selection_filtered <- size_small %>% filter(selection_intensity > 10^5)

selection_filtered$varient <- 
  reorder(selection_filtered$varient, -selection_filtered$selection_intensity)
selection_filtered$gene <- 
  reorder(selection_filtered$gene, -selection_filtered$selection_intensity)

# target genes for highlight
target <- c(c("VHL", "HIF1A", "HIF3A"), 
            as.vector(head(arrange(selection_filtered, desc(selection_intensity))$gene, 10)))

selection_filtered <- selection_filtered %>%
  mutate(label = ifelse(gene %in% target, gene, NA)) %>%
  mutate(label = ifelse(is.na(label), label, levels(gene)[label])) %>%
  arrange(desc(selection_intensity))

gg_bar <- 
  ggplot() + 
  geom_bar(data = selection_filtered, 
           aes(x = varient, y = selection_intensity/1e4),
           stat = "identity", fill = "black") +
  geom_bar(data = filter(selection_filtered, !is.na(label)), 
           aes(x = varient, y = selection_intensity/1e4, fill = gene),
           stat = "identity") +
  geom_point(data = filter(selection_filtered, !is.na(label)), 
             aes(x = varient, y = selection_intensity/1e4), shape = 20) +
  theme_classic() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = unit(c(1, -1, 1, 3.7), "mm"),
    legend.position = c(0.6, 0.85),
    legend.title = element_blank()
  ) +
  geom_text(
    data = filter(selection_filtered, !is.na(label)), 
    aes(x = varient, y = selection_intensity/1e4, label = as.character(varient), 
        color = gene),
    nudge_y = 0.18, nudge_x = 1.5,
    angle = 60,
    size = 2.25
  ) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  labs(y = bquote("Selection intensity ("~10^4~" )"))

gg_bar  # display



# plot with stage data

library(dplyr)
gene_selection_stage <- data.frame(varient = kirc_02@selection_results$variant,
                                   gene = kirc_02@selection_results$gene,
                                   selection_intensity = kirc_02@selection_results$selection_intensity,
                                   stage = kirc_02@selection_results$progression)
stage1 <- gene_selection_stage %>% filter(stage == "Stage I")
stage234 <- gene_selection_stage %>% filter(stage == "Stages234")

selection_filtered <- stage1 %>% filter(selection_intensity > 2.5e4)
# selection_filtered <- stage234 %>% filter(selection_intensity > 7e4)

selection_filtered$varient <- 
  reorder(selection_filtered$varient, -selection_filtered$selection_intensity)
selection_filtered$gene <- 
  reorder(selection_filtered$gene, -selection_filtered$selection_intensity)

# target genes for highlight
target <- c(c("VHL", "HIF1A", "HIF3A"), 
            as.vector(head(arrange(selection_filtered, desc(selection_intensity))$gene, 10)))

selection_filtered <- selection_filtered %>%
  mutate(label = ifelse(gene %in% target, gene, NA)) %>%
  mutate(label = ifelse(is.na(label), label, levels(gene)[label])) %>%
  arrange(desc(selection_intensity))

# generate bar plot

gg_bar <- 
  ggplot() + 
  geom_bar(data = selection_filtered, 
           aes(x = varient, y = selection_intensity/1e4),
           stat = "identity", fill = "black") +
  geom_bar(data = filter(selection_filtered, !is.na(label)), 
           aes(x = varient, y = selection_intensity/1e4, fill = gene),
           stat = "identity") +
  geom_point(data = filter(selection_filtered, !is.na(label)), 
             aes(x = varient, y = selection_intensity/1e4), shape = 20) +
  theme_classic() + 
  theme(
    axis.text.x = element_blank(),
    plot.margin = unit(c(1, -1, 1, 3.7), "mm"),
    legend.position = c(0.6, 0.85),
    legend.title = element_blank()
  ) +
  geom_text(
    data = filter(selection_filtered, !is.na(label)), 
    aes(x = varient, y = selection_intensity/1e4, label = as.character(varient), 
        color = gene),
    nudge_y = 0.5, nudge_x = 1.5,
    angle = 60,
    size = 2.25
  ) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(selection_filtered$selection_intensity)/1e4 + 1)) +
  labs(y = bquote("Selection intensity ("~10^4~" )"), x = "Varients",
       title = "Selection intensity of stage I tumor")

gg_bar  # display
