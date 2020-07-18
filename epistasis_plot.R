
# library(cancereffectsizeR)
library(reshape2)
library(ggplot2)
library(dplyr)

# ------- Generate plot for individual pair of genes -------
#' @param ces_obj CESAnalysis object with `gene_epistasis_results`
#' @param i row number for the pair of genes [default: 1]
#' @param colors color for gene1 and gene2

# Last update: July. 18, 2020

plot_pair <- function(ces_obj, i = 1, colors = c("red", "darkgreen")) {
  epistasis_results <- ces_obj@gene_epistasis_results
  genes <- epistasis_results[i, c("gene_1", "gene_2")]
  gene1_2 <- data.frame(substitutions = c(0:2),
                        gene_1=as.numeric(c(0, epistasis_results[i, c(3,6)])),
                        gene_2=as.numeric(c(0, epistasis_results[i, c(4,5)])))
  colnames(gene1_2) <- c("substitutions", genes)
  # Cumulative selection intensity
  gene1_2[3,2] <- gene1_2[2,2] + gene1_2[3,2]
  gene1_2[3,3] <- gene1_2[2,3] + gene1_2[3,3]
  
  gene1_2.m <- reshape2::melt(gene1_2, id=c("substitutions"))
  sub_01 <- gene1_2.m %>% filter(substitutions < 2)
  sub_12 <- gene1_2.m %>% filter(substitutions >= 1)
  sub_12$variable <- rev(sub_12$variable)
  
  plt <- ggplot() +
    geom_line(data = sub_01, aes(x = substitutions, y = value, color = variable)) +
    geom_point(data = sub_01, aes(x = substitutions, y = value, color = variable)) +
    geom_line(data = sub_12, aes(x = substitutions, y = value, color = variable)) +
    geom_point(data = sub_12, aes(x = substitutions, y = value, color = variable)) +
    theme_bw()+
    scale_colour_manual(values=colors)+
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(legend.title=element_blank(), legend.position="bottom") +
    xlab("Substitutions") + ylab("Cumulative Selection") +
    scale_x_continuous(breaks=c(0,1,2))
  
  return(plt)
}

# ---- Usage ----
# Load CES object (with epistasis results)
load("kirc_epistasis_stage_new.RData")
plot_pair(stage_epistasis)
ggsave("PVRM1_VHL.png", width = 5, dpi=300, height = 4)

