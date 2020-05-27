# -----------------------------------------------
# Box plot for gene selection intensity
# (by tumor size/stage)
# Last Modified: May.15, 2020
# -----------------------------------------------

library(dplyr)
library(ggplot2)
library(RColorBrewer)

load("cesanal_bysize.RData")  # size data (`analysis_size`)
load("cesanal_bystage.RData") # stage data (`kirc_02`)
load("cesanal_size_med.RData")  # size data break by median (`analysis_size_med`)

# size break by median (size >5 vs. size0-5)
size_data <- 
  data.frame(variant = analysis_size_med@selection_results$variant,
             gene = analysis_size_med@selection_results$gene,
             selection_intensity = analysis_size_med@selection_results$selection_intensity,
             group = analysis_size_med@selection_results$progression) %>%
  filter(selection_intensity > 1)

stage_data <- 
  data.frame(variant = kirc_02@selection_results$variant,
             gene = kirc_02@selection_results$gene,
             selection_intensity = kirc_02@selection_results$selection_intensity,
             group = kirc_02@selection_results$progression) %>%
  filter(selection_intensity > 1)

# top genes with highest total SI (TTN excluded)
top_gene_stage <- c("VHL", "PBRM1", "SETD2", "BAP1", 
                    "DNAH9", "DST", "HMCN1", "CSMD3", "MTOR", "MUC16")

# size >3 vs size 0-3
# top_gene_size <- c("VHL", "PBRM1", "SPEN", "MUC16", "KIAA0125",
#                    "SETD2", "DNMT3A", "CELSR1", "NBPF12", "ATM")

# size >5 vs size0-5
top_gene_size <- c("VHL", "PBRM1", "SETD2", "SPEN", "MUC16", 
                   "BAP1", "KIAA0125", "MTOR", "HMCN1", "DST")

stage_plot_data <- stage_data %>% 
  filter(gene %in% top_gene_stage) %>%
  mutate(gene = factor(gene, levels = top_gene_stage))

size_plot_data <- size_data %>% 
  filter(gene %in% top_gene_size) %>%
  mutate(gene = factor(gene, levels = top_gene_size),
         group = factor(group, levels = c("Size 0-5", "Size >5")))

# Generate box plot
si_boxplot <- function(data, group_names, genes, colormap, yticks) {
  palette <- brewer.pal(6, colormap)[c(2, 5)]
  myplt <-
    boxplot(selection_intensity ~ group*gene, data = data, boxwex=0.4,
            col = palette, xlab = "", ylab = "", 
            xaxt="n", yaxt="n")
  title(ylab = expression(paste("Selection Intensity /", "10"^"4")), 
        mgp = c(2, 0, 0))
  title(xlab = "Gene", mgp = c(2, 0, 0))
  
  axis(1, mgp = c(0, 0.2, 0), 
       at = seq(1.5 , 20 , 2), 
       labels = genes, 
       tick=FALSE , cex.axis=0.5)
  
  axis(2, at = yticks * 1e4, las=2,
       labels = yticks, cex.axis=0.6)
  
  # Add the grey vertical lines
  for(i in seq(0.5 , 21 , 2)){ 
    abline(v=i, lty=1, col="grey")
  }
  
  # Add a legend
  legend("topright", legend = group_names, 
         col=palette, 
         pch = 15, bty = "n", pt.cex = 2, cex = 1,  horiz = F)
}

pdf("stage_box0511.pdf", width = 8, height = 6)
si_boxplot(stage_plot_data, c("Stage I", "Stages II, III, IV"), 
           top_gene_stage, "Greens", yticks = seq(0, 40, 2.5))
dev.off()

pdf("size_med_box0515.pdf", width = 8, height = 6)
size_plt <- si_boxplot(size_plot_data, c("Size 0 ~ 5", "Size > 5"), 
           top_gene_size, "Blues", yticks = seq(0, 40, 2.5))
dev.off()

# ---- Extract outliers from the box plot ----
pick_outlier <- function(data, stats) {
  outlier <- which(data$selection_intensity < stats[1] | data$selection_intensity > stats[5])
  return(data[outlier, ])
}

pairs <- unique(size_plot_data[, c("gene", "group")]) %>%
  arrange(gene, group)

myplt <-
  boxplot(selection_intensity ~ group*gene, data = size_plot_data, plot = F)

df <- data.frame()
for (i in 1:20) {
  data <- filter(size_plot_data, gene == pairs[i, 1], 
         group == pairs[i, 2])
  stats <- myplt$stats[, i]
  outlier <- pick_outlier(data, stats)
  df <- rbind(df, outlier)
}
df <- df %>% arrange(gene, group)
write.csv(df, "size_med_outliers.csv", row.names = F, quote = F)

