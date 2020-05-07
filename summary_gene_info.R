# -----------------------------------------------------
# Summarise information of genes with multiple variants
# Output: table of #variant, mean_SI, std, top variant
#         for each gene (by stage/size)
# Last Modified: May.5, 2020
# -----------------------------------------------------

library(dplyr)

summary_multivar <- function(data, group_name) {
  data_clean <- data %>% 
    filter(group == group_name) %>%
    arrange(desc(selection_intensity)) %>%
    filter(selection_intensity > 1)
  
  # Summarise information of gene with multiple variants
  info1 <- data_clean %>% group_by(gene) %>%
    summarise(cum_si = sum(selection_intensity), # change sum to mean and sd
              mean_si = mean(selection_intensity),
              sd = sd(selection_intensity),
              max_si = max(selection_intensity),
              n_variant = n_distinct(variant)) %>%
    filter(n_variant > 1) 
  
  top_variant <- data_clean %>%
    group_by(gene) %>% filter(row_number() == 1)
  
  merge_info <- merge(info1, top_variant[, -3], by.x = "gene") %>%
    arrange(desc(cum_si), desc(n_variant))
  return(merge_info)
}

# load data
load("cesanal_bysize.RData")  # size data (`analysis_size`)
load("cesanal_bystage.RData") # stage data (`kirc_02`)

size_data <- 
  data.frame(variant = analysis_size@selection_results$variant,
             gene = analysis_size@selection_results$gene,
             selection_intensity = analysis_size@selection_results$selection_intensity,
             group = analysis_size@selection_results$progression)

stage_data <- 
  data.frame(variant = kirc_02@selection_results$variant,
             gene = kirc_02@selection_results$gene,
             selection_intensity = kirc_02@selection_results$selection_intensity,
             group = kirc_02@selection_results$progression)

large_info <- summary_multivar(size_data, group_name = "Size >3")
small_info <- summary_multivar(size_data, group_name = "Size 0-3")

stage1_info <- summary_multivar(stage_data, group_name = "Stage I")
stage234_info <- summary_multivar(stage_data, group_name = "Stages234")

fill_na <- function(x, fill = 0) {
  x = ifelse(is.na(x), fill, x)
  return(x)
}

# ---- size data ----
sd_info <- size_data %>% 
  filter(selection_intensity > 1) %>%
  group_by(gene) %>% summarise(sd = sd(selection_intensity))

size_merge <- merge(large_info, small_info, by = "gene", all = T, 
                    suffixes = c(".large", ".small")) %>%
  mutate_at(c("cum_si.small", "cum_si.large", 
              "mean_si.small", "mean_si.large", 
              "sd.small", "sd.large",
              "n_variant.large", "n_variant.small"), fill_na) %>%
  # mutate(cum_si_total = cum_si.small + cum_si.large) %>%
  mutate(n_variant_total = n_variant.small + n_variant.large,
         mean_si_total = (cum_si.small + cum_si.large) / n_variant_total) 

# size_merge <- merge(size_merge, sd_info, by.x = "gene") %>%
#   arrange(desc(n_variant_total))

size_merge <- size_merge %>%
  mutate(n_variant = paste(as.character(n_variant.small), 
                           as.character(n_variant.large), sep = "|"),
         topvar.large = ifelse(is.na(variant.large), "NA", 
                               paste0(variant.large, "(", 
                                      round(max_si.large/1e4, 1),"x10^4)")),
         topvar.small = ifelse(is.na(variant.small), "NA", 
                               paste0(variant.small, "(", 
                                      round(max_si.small/1e4, 1),"x10^4)")),
         topvar = paste(topvar.small, topvar.large, sep = "|"),
         mean_si = paste(round(mean_si.small/1e4, 2), round(mean_si.large/1e4, 2), sep = "|"),
         mean_si_total = round(mean_si_total/1e4, 2),
         sd = paste(round(sd.small/1e4, 2), round(sd.large/1e4, 2), sep = "|")
  ) %>%
  arrange(desc(n_variant_total))

size_merge_sub <- size_merge %>%
  select("gene", "n_variant", "mean_si", "sd", "topvar") %>%
  # arrange(desc(cum_si_total)) %>%
  mutate(format = "Small Size | Large Size")

names(size_merge_sub) <- 
  c("gene", "#variant", "mean SI (10^4)", "sd (10^4)", "TopVariant", "format")
write.csv(size_merge_sub, file = "gene_multivar_0505.csv", quote = F, row.names = F)

# ---- stage data ----
sd_info <- stage_data %>% 
  filter(selection_intensity > 1) %>%
  group_by(gene) %>% summarise(sd = sd(selection_intensity))

stage_merge <- merge(stage1_info, stage234_info, by = "gene", all = T, 
                     suffixes = c(".1", ".234")) %>%
  mutate_at(c("cum_si.1", "cum_si.234", 
              "mean_si.1", "mean_si.234", 
              "sd.1", "sd.234",
              "n_variant.1", "n_variant.234"), fill_na) %>%
  mutate(n_variant_total = n_variant.1 + n_variant.234, 
    mean_si_total = (cum_si.1 + cum_si.234) / n_variant_total) %>%
  arrange(desc(n_variant_total))


stage_merge <- stage_merge %>%
  mutate(n_variant = paste(as.character(n_variant.1), 
                           as.character(n_variant.234), sep = "|"),
         topvar.1 = ifelse(is.na(variant.1), "NA", 
                           paste0(variant.1, "(", 
                                  round(max_si.1/1e4, 1),"x10^4)")),
         topvar.234 = ifelse(is.na(variant.234), "NA", 
                             paste0(variant.234, "(", 
                                    round(max_si.234/1e4, 1),"x10^4)")),
         topvar = paste(topvar.1, topvar.234, sep = "|"),
         mean_si = paste(round(mean_si.1/1e4, 2), round(mean_si.234/1e4, 2), sep = "|"),
         mean_si_total = round(mean_si_total/1e4, 2),
         sd = paste(round(sd.1/1e4, 2), round(sd.234/1e4, 2), sep = "|")
  )

stage_merge_sub <- stage_merge %>%
  select("gene", "n_variant", "mean_si", "sd", "topvar") %>%
  # arrange(desc(cum_si_total)) %>%
  mutate(format = "Early Stage | Late Stage")

names(stage_merge_sub) <- 
  c("gene", "#variant", "mean SI (10^4)", "sd (10^4)", "TopVariant", "format")
write.csv(stage_merge_sub, file = "stage_gene_multivar_0505.csv", quote = F, row.names = F)
