# --------------------------------------------------
# Find common/similar variants between different 
# stages or sizes
# Input: table of genes with multiple variants (CSV)
# Output: table of common/similar variants (CSV)
# Last Modified: Apr.7, 2020
# --------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)

# load data (genes with multiple variants)
var_size <- read.csv("gene_multi_var_bysize.csv")
var_stage <- read.csv("gene_multi_var_bystage.csv")

gene_list <- c("VHL", "SPEN", "MAML1", "PBRM1", "CHIT1", "PTPRU", 
               "CARD11", "TRIOBP", "ELP2")
extract_size <- 
  var_size %>% filter(gene %in% gene_list, 
                      selection_intensity > 1000) %>%
  arrange(gene, varient)

var_tmp <- 
  as.data.frame(str_split_fixed(extract_size$varient, " ", 4)) %>%
  mutate(striped = paste(V1, V2))

extract_size$var_strip <- var_tmp$striped
common_var <- extract_size %>% group_by(var_strip) %>%
  summarise(n_grp = n_distinct(group)) %>% filter(n_grp > 1)
common_var_size <- extract_size %>% 
  filter(var_strip %in% common_var$var_strip)

extract_stage <- var_stage %>% filter(gene %in% gene_list,
                                      selection_intensity > 1000)

var_tmp <- 
  as.data.frame(str_split_fixed(extract_stage$varient, " ", 4)) %>%
  mutate(striped = paste(V1, V2))

extract_stage$var_strip <- var_tmp$striped
common_var <- extract_stage %>% group_by(var_strip) %>%
  summarise(n_grp = n_distinct(group)) %>% filter(n_grp > 1)
common_var_stage <- extract_stage %>% 
  filter(var_strip %in% common_var$var_strip)

extract_common <- rbind(common_var_size, common_var_stage) %>%
  arrange(gene, varient, group)%>%
  mutate(info = paste0(varient, "(", round(selection_intensity/1e4, 1),
                       "x10^4)"))

common_wide <- 
  pivot_wider(extract_common, names_from = group, values_from = info) %>%
  select(-varient, -selection_intensity)
write.csv(common_wide, "common_tab.csv", row.names = F, quote = F)


