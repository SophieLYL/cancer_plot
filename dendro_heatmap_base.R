# -----------------------------------------------
# Heatmap with Dendrogram for Tumor and Signature
# Last Modified: Apr.14, 2020
# -----------------------------------------------

library(dplyr)
library(RColorBrewer)

load("cesanal_bystage.RData")
load("cesanal_bysize.RData")

working_data <- kirc_02  # stage data
# working_data <- analysis_size  # size data

# ---- Extract signature data ----
tumor2sig <- 
  working_data@trinucleotide_mutation_weights$signatures_output_list

sig_df <- data.frame()
for (record in tumor2sig) {
  sig_df <- rbind(sig_df, record$signatures_output$weights)
}
weight.prevalence <- apply(sig_df, 2, function(x) {length(which(x > 0))})

# Exclude signatures with zeros
weight.prevalence <- weight.prevalence[weight.prevalence > 0]
sig_df <- sig_df %>%
  select(c(names(weight.prevalence), "SBS27"))

# ---- Generate dendrogram ----
dd.col <- as.dendrogram(hclust(dist(sig_df)))
dd.row <- as.dendrogram(hclust(dist(t(sig_df))))

# ---- Generate heatmap ----
heatmap(as.matrix(sig_df), labRow = "", 
        xlab = "Signatures", ylab = "Tumors",
        Colv = dd.row, Rowv = dd.col, scale = "none",
        col= colorRampPalette(brewer.pal(8, "Greens"))(25))  # Green for stage
        # col= colorRampPalette(brewer.pal(8, "Blues"))(25))  # Blue for size
