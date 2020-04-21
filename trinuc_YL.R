# -----------------------------------------------
# Heatmap for Upstream & Downstream Mutation 
# (by tumor size/stage)
# Last Modified: Apr.21, 2020
# -----------------------------------------------

library(cancereffectsizeR)
library(cowplot)
library(dplyr)
library(ggplot2)


# Set pathname and load data
load("cesanal_bysize.RData")  # dataset named "analysis_size" 
load("cesanal_bystage.RData")  # dataset named "kirc_02"

# Extracting trinuc data
KIRC_trinuc <- analysis_size@trinucleotide_mutation_weights$trinuc_proportion_matrix

# Seperating trinuc data into Size_3_and_less and Size_3.1_and_up
tumor_names <- unique(analysis_size@annotated.snv.maf$Unique_Patient_Identifier)

tumor_Size_3_and_less <- c()
tumor_Size_3.1_and_up <- c()

for (x in tumor_names){
  if (analysis_size@progressions@progression_by_tumor[[x]] == 1){
    tumor_Size_3_and_less <- c(tumor_Size_3_and_less, x)
  }
  else if (analysis_size@progressions@progression_by_tumor[[x]] == 2){
    tumor_Size_3.1_and_up <- c(tumor_Size_3.1_and_up, x)
  }
}

tumor_trinuc_Size_3_and_less <- KIRC_trinuc[tumor_Size_3_and_less, ]
tumor_trinuc_Size_3.1_and_up <- KIRC_trinuc[tumor_Size_3.1_and_up, ]

#Averaging trinuc data
tumor_trinuc_Size_3_and_less_avg <- 
  apply(tumor_trinuc_Size_3_and_less, 2, mean)
tumor_trinuc_Size_3.1_and_up_avg <- 
  apply(tumor_trinuc_Size_3.1_and_up, 2, mean)

# Ordering trinuc average data
tumor_trinuc_Size_3_and_less_avg_ordered <- 
  data.frame(average = tumor_trinuc_Size_3_and_less_avg,
             mutation = names(tumor_trinuc_Size_3_and_less_avg)) %>%
  mutate(upstream = substr(mutation, 1, 1),
         downstream = substr(mutation, 7, 7),
         mut_from = substr(mutation, 3, 3),
         mut_to = substr(mutation, 5, 5),
         mutation_name = paste0(mut_from, "\u2192", mut_to)) %>%
  arrange(downstream, upstream)

tumor_trinuc_Size_3.1_and_up_avg_ordered <- 
  data.frame(average = tumor_trinuc_Size_3.1_and_up_avg,
             mutation = names(tumor_trinuc_Size_3_and_less_avg)) %>%
  mutate(upstream = substr(mutation, 1, 1),
         downstream = substr(mutation, 7, 7),
         mut_from = substr(mutation, 3, 3),
         mut_to = substr(mutation, 5, 5),
         mutation_name = paste0(mut_from, "\u2192", mut_to)) %>%
  arrange(downstream, upstream)

trinuc.mutation_data <- tumor_trinuc_Size_3_and_less_avg_ordered[, -1]

# ---- Genereate heat map ----
KIRC_trinuc_heatmap_data <- 
  data.frame(deconstrucSig = trinuc.mutation_data$mutation,
             Upstream = trinuc.mutation_data$upstream,
             Downstream = trinuc.mutation_data$downstream,
             mutated_from = trinuc.mutation_data$mut_from,
             mutated_to = trinuc.mutation_data$mut_to,
             trinuc_Size_3_and_less = tumor_trinuc_Size_3_and_less_avg_ordered$average,
             trinuc_Size_3.1_and_up = tumor_trinuc_Size_3.1_and_up_avg_ordered$average) %>%
  mutate(mutation = paste0(mutated_from, "\u2192", mutated_to))

levels(KIRC_trinuc_heatmap_data$mutation) <- 
  c("C\u2192A", "C\u2192G", "C\u2192T", "T\u2192A", "T\u2192C", "T\u2192G")
save(KIRC_trinuc_heatmap_data, file="KIRC_trinuc_heatmap_data_sizes.RData")

KIRC_trinuc_heatmap_Size_3_and_less <- 
  ggplot(data=KIRC_trinuc_heatmap_data, aes(x=Downstream, Upstream)) +
  geom_tile(aes(fill = trinuc_Size_3_and_less*100), color="white") + 
  scale_fill_gradient(low="white", high="steelblue", name="Percent") +
  facet_grid(.~mutation)+
  geom_text(aes(label = round(trinuc_Size_3_and_less, 4)*100), size=2) +
  labs(title="Size_3_and_less", x="Downstream", y="Upstream") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text = element_text(size=15),
                     axis.title.x = element_text(size=15),
                     axis.title.y = element_text(size=15),
                     axis.text.x = element_text(size=12),
                     axis.text.y = element_text(size=12))

KIRC_trinuc_heatmap_Size_3.1_and_up <- 
  ggplot(data=KIRC_trinuc_heatmap_data, aes(x=Downstream, Upstream)) +
  geom_tile(aes(fill = trinuc_Size_3.1_and_up*100), color="white") + 
  scale_fill_gradient(low="white", high="steelblue", name="Percent") +
  facet_grid(.~mutation)+
  geom_text(aes(label = round(trinuc_Size_3.1_and_up, 4)*100), size=2) +
  labs(title="Size_3.1_and_up", x="Downstream", y="Upstream") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text = element_text(size=15),
                     axis.title.x = element_text(size=15),
                     axis.title.y = element_text(size=15),
                     axis.text.x = element_text(size=12),
                     axis.text.y = element_text(size=12))

# Save to file
KIRC_trinuc_heatmap_Size_3_and_less
ggsave("KIRC_trinuc_heatmap_Size_3_and_less.png", width=8, height=2.25)

KIRC_trinuc_heatmap_Size_3.1_and_up
ggsave("KIRC_trinuc_heatmap_Size_3.1_and_up.png", width=8, height=2.25)

combined_trinuc_heatmap <- 
  plot_grid(KIRC_trinuc_heatmap_Size_3_and_less, 
            KIRC_trinuc_heatmap_Size_3.1_and_up, 
            labels = c("A", "B"), align="h", ncol=1)
combined_trinuc_heatmap
ggsave("combined_trinuc_heatmap_sizes.png", width=8, height=4.5)


# ---- heatmap by stage ----
# kirc_02: stage data generated by the `cancereffectsizeR` package
KIRC_trinuc_stage <- 
  kirc_02@trinucleotide_mutation_weights$trinuc_proportion_matrix

#Seperating trinuc data into Size_3_and_less and Size_3.1_and_up
tumor_names <- unique(kirc_02@annotated.snv.maf$Unique_Patient_Identifier)

stage1 <- c()
stage234 <- c()

for (x in tumor_names){
  if (kirc_02@progressions@progression_by_tumor[[x]] == 1){
    stage1 <- c(stage1, x)
  }
  else if (kirc_02@progressions@progression_by_tumor[[x]] == 2){
    stage234 <- c(stage234, x)
  }
}

tumor_stage1 <- KIRC_trinuc_stage[stage1, ]
tumor_stage234 <- KIRC_trinuc_stage[stage234, ]

#Averaging trinuc data
tumor_stage1_avg <- 
  apply(tumor_stage1, 2, mean)
tumor_stage234_avg <- 
  apply(tumor_stage234, 2, mean)

# Ordering trinuc average data
tumor_stage1_avg_ordered <- 
  data.frame(average = tumor_stage1_avg,
             mutation = names(tumor_stage1_avg)) %>%
  mutate(upstream = substr(mutation, 1, 1),
         downstream = substr(mutation, 7, 7),
         mut_from = substr(mutation, 3, 3),
         mut_to = substr(mutation, 5, 5),
         mutation_name = paste0(mut_from, "\u2192", mut_to)) %>%
  arrange(downstream, upstream)

tumor_stage234_avg_ordered <- 
  data.frame(average = tumor_stage234_avg,
             mutation = names(tumor_stage1_avg)) %>%
  mutate(upstream = substr(mutation, 1, 1),
         downstream = substr(mutation, 7, 7),
         mut_from = substr(mutation, 3, 3),
         mut_to = substr(mutation, 5, 5),
         mutation_name = paste0(mut_from, "\u2192", mut_to)) %>%
  arrange(downstream, upstream)

trinuc.mutation_data <- tumor_stage1_avg_ordered[, -1]

# ---- Genereate heatmap ----

KIRC_trinuc_stage_heatmap_data <- 
  data.frame(deconstrucSig=trinuc.mutation_data$mutation,
             Upstream=trinuc.mutation_data$upstream,
             Downstream=trinuc.mutation_data$downstream,
             mutated_from=trinuc.mutation_data$mut_from,
             mutated_to=trinuc.mutation_data$mut_to,
             trinuc_Size_3_and_less=tumor_stage1_avg_ordered$average,
             trinuc_Size_3.1_and_up=tumor_stage234_avg_ordered$average) %>%
  mutate(mutation = paste0(mutated_from, "\u2192", mutated_to))


levels(KIRC_trinuc_stage_heatmap_data$mutation) <- 
  c("C\u2192A", "C\u2192G", "C\u2192T", "T\u2192A", "T\u2192C", "T\u2192G")

# backup heatmap data
# save(KIRC_trinuc_stage_heatmap_data, 
#      file="KIRC_trinuc_heatmap_data_stage.RData")

# Stage 1
KIRC_trinuc_stage1_heatmap <- 
  ggplot(data=KIRC_trinuc_stage_heatmap_data, aes(x=Downstream, Upstream)) +
  geom_tile(aes(fill = trinuc_Size_3_and_less*100), color="white") + 
  scale_fill_gradient(low="white", high="dark green", name="Percent") +
  facet_grid(.~mutation)+
  geom_text(aes(label = round(trinuc_Size_3_and_less, 4)*100), size=2) +
  labs(title="Stage 1", x="Downstream", y="Upstream") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text = element_text(size=15),
                     axis.title.x = element_text(size=15),
                     axis.title.y = element_text(size=15),
                     axis.text.x = element_text(size=12),
                     axis.text.y = element_text(size=12))

# Stage 2,3,4
KIRC_trinuc_stage234_heatmap <- 
  ggplot(data=KIRC_trinuc_stage_heatmap_data, aes(x=Downstream, Upstream)) +
  geom_tile(aes(fill = trinuc_Size_3.1_and_up*100), color="white") + 
  scale_fill_gradient(low="white", high="dark green", name="Percent") +
  facet_grid(.~mutation)+
  geom_text(aes(label = round(trinuc_Size_3.1_and_up, 4)*100), size=2) +
  labs(title="Stage 2,3,4", x="Downstream", y="Upstream") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text = element_text(size=15),
                     axis.title.x = element_text(size=15),
                     axis.title.y = element_text(size=15),
                     axis.text.x = element_text(size=12),
                     axis.text.y = element_text(size=12))
# Save to file
KIRC_trinuc_stage1_heatmap
ggsave("KIRC_trinuc_heatmap_stage1.png", width=8, height=2.25)

KIRC_trinuc_stage234_heatmap
ggsave("KIRC_trinuc_heatmap_stage234.png", width=8, height=2.25)

combined_trinuc_heatmap_stage <- 
  plot_grid(KIRC_trinuc_stage1_heatmap, KIRC_trinuc_stage234_heatmap, 
            labels = c("A", "B"), align="h", ncol=1)
combined_trinuc_heatmap_stage
ggsave("combined_trinuc_heatmap_sizes.png", width=8, height=4.5)
