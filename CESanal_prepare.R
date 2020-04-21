# Preprocess

# set working directory

maf_file <- read.csv("TCGA.KIRC.MAF_hg19_with_Sizes.csv", stringsAsFactors = F)
library(cancereffectsizeR)
# Create CESAnalysis object and define the chronological tumor progression stages

ces_anal = CESAnalysis(genome = "hg19", 
                       progression_order = c("Size >3", "Size 0-3"))

analysis_size <- load_maf(ces_anal, maf = maf_file, 
                          progression_col = "tumor_size")
analysis_size = 
  calc_baseline_mutation_rates(analysis_size, covariate_file = "kirc_pca")
# Calculate selection intensities and produce human-readable results
# If you have multiple computing cores and the parallel library, you can parallelize the operation
analysis_size = ces_snv(analysis_size, include_genes_without_recurrent_mutations = T)



