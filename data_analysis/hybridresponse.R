###### input: output deseq2, hybrid response
###### output: rose venn diagram, MA and scatter plots. Candidate lists, filtered on directionality and 2lfc.

library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(RColorBrewer)
library(ggrepel)
library(grid)
library(ggvenn)
library(tidyverse)
library(readxl)
library(reshape2)


rm(list = ls())

# 18 oC deseq2 output, hybrid response
atxaa_at <- read_excel("output/2025/hybridresponse/base_analysis/big-overview/atxaa_at_18.xlsx") %>% as.data.frame()
atxaa_aa <- read_excel("output/2025/hybridresponse/base_analysis/big-overview/atxaa_aa_18.xlsx") %>% as.data.frame()
atxal_at <- read_excel("output/2025/hybridresponse/base_analysis/big-overview/atxal_at_18.xlsx") %>% as.data.frame()
atxal_al <- read_excel("output/2025/hybridresponse/base_analysis/big-overview/atxal_al_18.xlsx") %>% as.data.frame()

atxal_atxaa <- read_excel("output/2025/hybridresponse/base_analysis/big-overview/atxal_atxaa_18.xlsx")
atxal_atxaa <- atxal_atxaa %>%
  dplyr::filter(padj < 0.05)
# Gene descriptions from araport
araport <- read_tsv("input/genedescriptions/all_genedescriptions.tsv", col_names = FALSE)
colnames(araport) <- c("Gene", "Symbol", "Genename", "Description")

# Ortholog lists
ortho_aa <- read_tsv("input/orthologs/At_Aa_orthologs.tsv")
ortho_al <- read_tsv("input/orthologs/At_Al_orthologs.tsv")

common_orthos <- intersect(ortho_aa$Ath, ortho_al$Ath)

# Extracting only interesting data, and keeping in mind merging of datasets will happen. Only analyzing 1:1:1 orthologs
atxaa_at <- atxaa_at %>%
  dplyr::select(Gene, baseMean, log2FoldChange, padj) %>%
  dplyr::filter(Gene %in% common_orthos) %>% 
  dplyr::filter(padj < 0.05)
  
atxaa_aa <- atxaa_aa %>%
  dplyr::select(Gene, baseMean, log2FoldChange, padj) %>%
  dplyr::filter(Gene %in% common_orthos) %>% 
  dplyr::filter(padj < 0.05)

atxal_at <- atxal_at %>%
  dplyr::select(Gene, baseMean, log2FoldChange, padj) %>%
  dplyr::filter(Gene %in% common_orthos) %>% 
  dplyr::filter(padj < 0.05)

atxal_al <- atxal_al %>%
  dplyr::select(Gene, baseMean, log2FoldChange, padj) %>%
  dplyr::filter(Gene %in% common_orthos) %>% 
  dplyr::filter(padj < 0.05)

# Venn diagram
ggvenn(
  list(atxaa_at = atxaa_at$Gene, atxaa_aa = atxaa_aa$Gene, atxal_al = atxal_al$Gene, atxal_at = atxal_at$Gene),
  fill_color = c("cyan", "magenta", "yellow", "green")
)

# Gene lists with genes from interesting intersections in last venn diagram
common_all <- Reduce(intersect, list(atxaa_at$Gene, atxaa_aa$Gene, atxal_at$Gene, atxal_al$Gene))
common_aa <- intersect(atxaa_at$Gene, atxaa_aa$Gene)
exclusive_aa <- setdiff(common_aa, union(atxal_at$Gene, atxal_al$Gene))
common_lyr <- intersect(atxal_at$Gene, atxal_al$Gene)
exclusive_lyr <- setdiff(common_lyr, union(atxaa_at$Gene, atxaa_aa$Gene))



#### ONE-SIDED HYBRRES
## Arenosa
# In merging, .x will always refer to comparisons against thaliana!
aahybr_oneside <- merge(atxaa_at, atxaa_aa, by = "Gene")
aahybr_oneside <- aahybr_oneside %>%
  dplyr::filter(Gene %in% exclusive_aa) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

ggvenn(
  list(aa_oneside = aahybr_oneside$Gene, atxal_v_atxaa = atxal_atxaa$Gene),
  fill_color = c("cyan", "magenta")
)



# Directionality filter  
aa_hybr <- aahybr_oneside %>%
  dplyr::filter((log2FoldChange.x > 0 & log2FoldChange.y > 0) | (log2FoldChange.x < -0 & log2FoldChange.y < -0)) 

aa_hybr <- merge(araport, aa_hybr, by = "Gene")
aa_hybr$Description <- trimws(gsub(";\\(.*$", "", aa_hybr$Description))


#### Making plots
aahybr_oneside <- aahybr_oneside %>%
  dplyr::mutate(log2_baseMean.x = log2(baseMean.x + 1)) %>%
  dplyr::mutate(log2_baseMean.y = log2(baseMean.y + 1))

# Is the gene a significant DEG in hybrid-hybrid comparison?
aahybr_oneside <- aahybr_oneside %>%
  mutate(hybridcomp = case_when(
    
    Gene %in% atxal_atxaa$Gene ~ "yes",
    
    TRUE ~ "no"
  ))

# Setting labels for genes regulated at least 2 lfc in the same direction in both comparisons
aahybr_oneside <- aahybr_oneside %>%
  mutate(regulated = case_when(
    log2FoldChange.x > 2 & log2FoldChange.y > 2 ~ "up",
    log2FoldChange.x < -2 & log2FoldChange.y < -2 ~ "down",
    
    TRUE ~ "no"
  ))


aahybr_oneside <- aahybr_oneside %>%
  mutate(candidate = case_when(
    hybridcomp == "yes" & regulated != "no" ~ "yes",
    TRUE ~ "no"
    
  ))
 
sum(aahybr_oneside$hybridcomp == "yes")
sum(aahybr_oneside$regulated != "no")
sum(aahybr_oneside$candidate == "yes")

# MA arenosa hybrid vs at
ggplot(aahybr_oneside, aes(x = log2_baseMean.x, y = log2FoldChange.x)) +
  geom_point(aes(color = regulated, shape = hybridcomp, size = regulated)) +
  scale_color_manual(values = c("up" = "#FF6F3C", "down" = "#1E615E", "no" = "darkgrey")) +
  scale_shape_manual(values = c("yes" = 17, "no" = 1)) +
  scale_size_manual(values = c("up" = 1.4, "down" = 1.4, "no" = 0.8), guide = "none") + 
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "log2BaseMean", 
       y = "log2FC", 
       title = "One-sided hybrid response: A.t x A.a vs A.t") +
  scale_x_continuous(limits = c(0, ceiling(max(aahybr_oneside$log2_baseMean.x))), 
                     breaks = seq(0, ceiling(max(aahybr_oneside$log2_baseMean.x)), by = 2)) +
  scale_y_continuous(limits = c(-6, 12),  
                     breaks = seq(-6, 12, by = 2)) +  
  theme_light() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



# MA arenosa hybrid vs aa
ggplot(aahybr_oneside, aes(x = log2_baseMean.y, y = log2FoldChange.y)) +
  geom_point(aes(color = regulated, shape = hybridcomp, size = regulated)) +
  scale_color_manual(values = c("up" = "#FF6F3C", "down" = "#1E615E", "no" = "darkgrey")) +
  scale_shape_manual(values = c("yes" = 17, "no" = 1)) +
  scale_size_manual(values = c("up" = 1.5, "down" = 1.5, "no" = 0.8), guide = "none") + 
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "log2BaseMean", 
       y = "log2FC", 
       title = "One-sided hybrid response: A.t x A.a vs A.a") +
  scale_x_continuous(limits = c(0, ceiling(max(aahybr_oneside$log2_baseMean.y))), 
                     breaks = seq(0, ceiling(max(aahybr_oneside$log2_baseMean.y)), by = 2)) +
  scale_y_continuous(limits = c(-12, 12),  
                     breaks = seq(-12, 12, by = 2)) +  
  theme_light() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


## lyrata
alhybr_oneside <- merge(atxal_at, atxal_al, by = "Gene")
alhybr_oneside <- alhybr_oneside %>%
  dplyr::filter(Gene %in% exclusive_lyr) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

# Directionality filter
al_hybr <- alhybr_oneside %>%
  dplyr::filter(Gene %in% exclusive_lyr) %>%
  dplyr::filter((log2FoldChange.x > 0 & log2FoldChange.y > 0) | (log2FoldChange.x < -0 & log2FoldChange.y < -0)) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

al_hybr <- merge(araport, al_hybr, by = "Gene")
al_hybr$Description <- trimws(gsub(";\\(.*$", "", al_hybr$Description))


#### Making plots
alhybr_oneside <- alhybr_oneside %>%
  dplyr::mutate(log2_baseMean.x = log2(baseMean.x + 1)) %>%
  dplyr::mutate(log2_baseMean.y = log2(baseMean.y + 1)) 

# Is the gene a significant DEG in hybrid-hybrid comparison?
alhybr_oneside <- alhybr_oneside %>%
  mutate(hybridcomp = case_when(
    
    Gene %in% atxal_atxaa$Gene ~ "yes",
    
    TRUE ~ "no"
  ))

# Setting labels for genes regulated at least 2 lfc in the same direction in both comparisons
alhybr_oneside <- alhybr_oneside %>%
  mutate(regulated = case_when(
    log2FoldChange.x > 2 & log2FoldChange.y > 2 ~ "up",
    log2FoldChange.x < -2 & log2FoldChange.y < -2 ~ "down",
    
    TRUE ~ "no"
  ))

alhybr_oneside <- alhybr_oneside %>%
  mutate(candidate = case_when(
    hybridcomp == "yes" & regulated != "no" ~ "yes",
    TRUE ~ "no"
    
  ))

sum(alhybr_oneside$hybridcomp == "yes")
sum(alhybr_oneside$regulated != "no")
sum(alhybr_oneside$candidate == "yes")

# MA lyrata hybrid
ggplot(alhybr_oneside, aes(x = log2_baseMean.x, y = log2FoldChange.x)) +
  geom_point(aes(color = regulated, shape = hybridcomp, size = regulated)) +
  scale_color_manual(values = c("up" = "#FF6F3C", "down" = "#1E615E", "no" = "darkgrey")) +
  scale_shape_manual(values = c("yes" = 17, "no" = 1)) +
  scale_size_manual(values = c("up" = 1.4, "down" = 1.4, "no" = 0.8), guide = "none") + 
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "log2BaseMean", 
       y = "log2FC", 
       title = "One-sided hybrid response: A.t x A.l vs A.t") +
  scale_x_continuous(limits = c(0, ceiling(max(alhybr_oneside$log2_baseMean.x))), 
                     breaks = seq(0, ceiling(max(alhybr_oneside$log2_baseMean.x)), by = 2)) +
  scale_y_continuous(limits = c(-6, 12),  
                     breaks = seq(-6, 12, by = 2)) +  
  theme_light() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# MA lyrata hybrid vs lyrata
ggplot(alhybr_oneside, aes(x = log2_baseMean.y, y = log2FoldChange.y)) +
  geom_point(aes(color = regulated, shape = hybridcomp, size = regulated)) +
  scale_color_manual(values = c("up" = "#FF6F3C", "down" = "#1E615E", "no" = "darkgrey")) +
  scale_shape_manual(values = c("yes" = 17, "no" = 1)) +
  scale_size_manual(values = c("up" = 1.2, "down" = 1.2, "no" = 0.5), guide = "none") + 
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "log2BaseMean", 
       y = "log2FC", 
       title = "One-sided hybrid response: A.t x A.l vs A.l") +
  scale_x_continuous(limits = c(0, ceiling(max(alhybr_oneside$log2_baseMean.y))), 
                     breaks = seq(0, ceiling(max(alhybr_oneside$log2_baseMean.y)), by = 2)) +
  scale_y_continuous(limits = c(-12, 12),  
                     breaks = seq(-12, 12, by = 2)) +  
  theme_light() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#### COMMON HYBRRES
aa <- merge(atxaa_at, atxaa_aa, by = "Gene")
aa <- aa %>%
  dplyr::rename(baseMean_atxaa_at = baseMean.x,
                baseMean_atxaa_aa = baseMean.y,
                log2FoldChange_atxaa_at = log2FoldChange.x,
                log2FoldChange_atxaa_aa = log2FoldChange.y,
                padj_atxaa_at = padj.x,
                padj_atxaa_aa = padj.y)

al <- merge(atxal_at, atxal_al, by = "Gene")
al <- al %>%
  dplyr::rename(baseMean_atxal_at = baseMean.x,
                baseMean_atxal_al = baseMean.y,
                log2FoldChange_atxal_at = log2FoldChange.x,
                log2FoldChange_atxal_al = log2FoldChange.y,
                padj_atxal_at = padj.x,
                padj_atxal_al = padj.y)

common <- merge(aa, al, by = "Gene")

common <- common %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

common <- merge(araport, common, by = "Gene")


# Directionality filter
common_filtered <- common %>%
  dplyr::filter(Gene %in% common_all) %>%
  dplyr::filter((log2FoldChange_atxaa_at > 0 & log2FoldChange_atxaa_aa > 0) | 
                  (log2FoldChange_atxaa_at < -0 & log2FoldChange_atxaa_aa < -0)) %>% # Directionality filter, atxaa: removing any rows where atxaa_at and atxaa_aa show regulation in antagonistic manner
  dplyr::filter((log2FoldChange_atxal_at > 0 & log2FoldChange_atxal_al > 0) | 
                  (log2FoldChange_atxal_at < -0 & log2FoldChange_atxal_al < -0)) # Directionality filter, atxal

common_filtered$regulated <- apply(common_filtered[, c("log2FoldChange_atxaa_at", 
                                     "log2FoldChange_atxaa_aa", 
                                     "log2FoldChange_atxal_at", 
                                     "log2FoldChange_atxal_al")], 
                          1, 
                          function(x) {
                            log2FoldChange_atxaa_at <- x[1]
                            log2FoldChange_atxal_at <- x[3]
                            difference <- abs(log2FoldChange_atxaa_at - log2FoldChange_atxal_at)
                            
                            if (all(x >= 2 | x <= -2)) {
                              return("yes")
                            } else if ((log2FoldChange_atxaa_at * log2FoldChange_atxal_at < 0) && 
                                       (difference >= 2)) {
                              return("antagonistic")
                            } else {
                              return("no")
                            }
                          })

filter <- common_filtered %>%
  dplyr::select(Gene, regulated)

common <- left_join(common, filter, by = "Gene")
common <- common %>%
  dplyr::mutate(regulated = ifelse(is.na(regulated), "no", regulated))

common$color <- with(common, ifelse(regulated == "no", "darkgrey", 
                                    ifelse(regulated == "antagonistic", "#FFCC00", 
                                           ifelse(log2FoldChange_atxaa_at > 2, "#FF6F3C", 
                                                  ifelse(log2FoldChange_atxaa_at < -2, "#1E615E", "darkgrey")))))


common <- common %>%
  mutate(hybridcomp = case_when(
    
    Gene %in% atxal_atxaa$Gene ~ "yes",
    
    TRUE ~ "no"
  ))

common_filtered <- common_filtered %>%
  mutate(hybridcomp = case_when(
    
    Gene %in% atxal_atxaa$Gene ~ "yes",
    
    TRUE ~ "no"
  ))

sum(common$hybridcomp == "yes")
sum(common$regulated == "antagonistic")

# plot
ggplot(common, aes(x = log2FoldChange_atxaa_at, y = log2FoldChange_atxal_at)) +
  geom_point(aes(color = color, shape = hybridcomp, size = regulated), alpha = 0.8) + 
  scale_size_manual(values = c("yes" = 1.2, "antagonistic" = 1.2, "no" = 0.5), guide = "none") + 
  scale_shape_manual(values = c("yes" = 17, "no" = 1)) +
  scale_color_identity() + 
  geom_hline(yintercept = 0, color = "black") + 
  geom_vline(xintercept = 0, color = "black") +
  labs(x = "A.t x A.a vs A.t (log2FC)", 
       y = "A.t x A.l vs A.t (log2FC)", 
       title = "Common hybrid response") +
  scale_x_continuous(breaks = seq(floor(min(common$log2FoldChange_atxaa_at)), ceiling(max(common$log2FoldChange_atxaa_at)), by = 2)) + 
  scale_y_continuous(breaks = seq(floor(min(common$log2FoldChange_atxal_at)), ceiling(max(common$log2FoldChange_atxal_at)), by = 2)) +
  theme_light() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# plot with aa and al comparisons instead of at
ggplot(common, aes(x = log2FoldChange_atxaa_aa, y = log2FoldChange_atxal_al)) +
  geom_point(aes(color = color, shape = hybridcomp, size = regulated), alpha = 0.8) + 
  scale_size_manual(values = c("yes" = 1.2, "antagonistic" = 1.2, "no" = 0.5), guide = "none") + 
  scale_shape_manual(values = c("yes" = 17, "no" = 1)) +
  scale_color_identity() + 
  geom_hline(yintercept = 0, color = "black") + 
  geom_vline(xintercept = 0, color = "black") +
  labs(x = "A.t x A.a vs A.a (log2FC)", 
       y = "A.t x A.l vs A.l (log2FC)", 
       title = "Common hybrid response") +
  scale_x_continuous(breaks = seq(floor(min(common$log2FoldChange_atxaa_aa)), ceiling(max(common$log2FoldChange_atxaa_aa)), by = 2)) + 
  scale_y_continuous(breaks = seq(floor(min(common$log2FoldChange_atxal_al)), ceiling(max(common$log2FoldChange_atxal_al)), by = 2)) +
  theme_light() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())




# Write lists
write.xlsx(aa_hybr, file = "output/2025/hybridresponse/filtered_results/2025-06-07_onesided_atxaa.xlsx")
write.xlsx(al_hybr, file = "output/2025/hybridresponse/filtered_results/2025-06-07_onesided_atxal.xlsx")
write.xlsx(common_filtered, file = "output/2025/hybridresponse/filtered_results/2025-06-07_common_18.xlsx")

# Filtering candidates on 2lfc
aa_hybr <- aa_hybr %>%
  dplyr::filter((log2FoldChange.x > 2 & log2FoldChange.y > 2) | (log2FoldChange.x < -2 & log2FoldChange.y < -2))

al_hybr <- al_hybr %>%
  dplyr::filter((log2FoldChange.x > 2 & log2FoldChange.y > 2) | (log2FoldChange.x < -2 & log2FoldChange.y < -2))

common_filtered <- common_filtered %>%
  dplyr::filter(str_detect(regulated, "antagonistic"))


# Write lists
write.xlsx(aa_hybr, file = "output/2025/hybridresponse/filtered_results/2025-06-07_onesided_atxaa_2lfc.xlsx")
write.xlsx(al_hybr, file = "output/2025/hybridresponse/filtered_results/2025-06-07_onesided_atxal_2lfc.xlsx")
write.xlsx(common_filtered, file = "output/2025/hybridresponse/filtered_results/2025-06-07_antagonists.xlsx")

