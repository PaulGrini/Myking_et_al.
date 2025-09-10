## Making temperature filters from homozygous samples
Differentially expressed genes between 26 and 18 oC in homozygous samples regulated at least 1 log fold are removed from temperature analysis

```R
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

rm(list = ls())

at <- read_excel("output/2025/temperatureeffect/controls/at_26v18.xlsx")
aa <- read_excel("output/2025/temperatureeffect/controls/aa_26v18.xlsx")
al <- read_excel("output/2025/temperatureeffect/controls/al_26v18.xlsx")

at_filter <- at %>%
  filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>%
  dplyr::select(Gene)

aa_filter <- aa %>%
  filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>%
  dplyr::select(Gene)

al_filter <- al %>%
  filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>%
  dplyr::select(Gene)

write.table(at_filter, "output/2025/temperatureeffect/controls/temperature_filter_at.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(aa_filter, "output/2025/temperatureeffect/controls/temperature_filter_aa.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(al_filter, "output/2025/temperatureeffect/controls/temperature_filter_al.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
```
## Filtering outputs from DESeq2 from hybrids, 26 vs 18 oC
Filtering based on homozygous temperature effects and significance, and making MA plots (Figure 6B and SFigure 13)

```R
# temperature effects
arenosahybr <- read_excel("output/2025/temperatureeffect/unfiltered_results/atxaa_26v18.xlsx") %>% as.data.frame()
lyratahybr <- read_excel("output/2025/temperatureeffect/unfiltered_results/atxal_26v18.xlsx") %>% as.data.frame()

# temperature controls: files with all genes regulated at least 1 logFC in homozygous comparisons between 26 and 18 degrees
at <- read_tsv("output/2025/temperatureeffect/controls/temperature_filter_at.tsv")
aa <- read_tsv("output/2025/temperatureeffect/controls/temperature_filter_aa.tsv")
al <- read_tsv("output/2025/temperatureeffect/controls/temperature_filter_al.tsv")

# orthologs
ortho_aa <- read_tsv("input/orthologs/At_Aa_orthologues.tsv")
ortho_al <- read_tsv("input/orthologs/At_Al_orthologues.tsv")

ortho <- intersect(ortho_aa$Ath, ortho_al$Ath)

filter <- bind_rows(at, aa, al)

# filter on common orthologs and on parental temperature effect
atxaa <- arenosahybr %>%
  distinct(Gene, .keep_all = TRUE) %>%
  dplyr::filter(Gene %in% ortho) %>%
  anti_join(at, by = "Gene") %>%
  anti_join(aa, by = "Gene")

atxal <- lyratahybr %>%
  distinct(Gene, .keep_all = TRUE)  %>%
  dplyr::filter(Gene %in% ortho) %>%
  anti_join(at, by = "Gene") %>%
  anti_join(al, by = "Gene")

# extracting one-sided temperature effects
os_atxaa <- atxaa %>%
  anti_join(atxal, by = "Gene")

os_atxal <- atxal %>%
  anti_join(atxaa, by = "Gene")



############################################### setting up hybridresponse candidates (endosperm specific) to highlight in plots
# Load candidate gene lists (18 degrees)
aa_hybrres <- read_excel("output/2025/hybridresponse/filtered_results/2025-06-07_2lfc_onesided_atxaa.xlsx")
al_hybrres <- read_excel("output/2025/hybridresponse/filtered_results/2025-06-07_2lfc_onesided_atxal.xlsx")
antag_hybrres <- read_excel("output/2025/hybridresponse/filtered_results/2025-06-07_common-antagonists.xlsx")

# Load deseq between the hybrids
atxal_atxaa <- read_excel("output/2025/hybridresponse/base_analysis/2025-04-22_atxal_vs_atxaa_18.xlsx")

# Load filter: endosperm
endosperm <- read_tsv("input/imprinting/combined-endosperm-filter.tsv")

aa_hybrres <- aa_hybrres %>%
  dplyr::rename(log2FoldChange_atxaa_at = log2FoldChange.x,
                log2FoldChange_atxaa_aa = log2FoldChange.y)

al_hybrres <- al_hybrres %>%
  dplyr::rename(log2FoldChange_atxal_at = log2FoldChange.x,
                log2FoldChange_atxal_al = log2FoldChange.y)

# Filtering on genes that are different between the two hybrids
aa_hybrres <- aa_hybrres %>%
  dplyr::filter(Gene %in% atxal_atxaa$Gene) %>%
  dplyr::filter(Gene %in% endosperm$Gene)

al_hybrres <- al_hybrres %>%
  dplyr::filter(Gene %in% atxal_atxaa$Gene) %>%
  dplyr::filter(Gene %in% endosperm$Gene)

antag_hybrres <- antag_hybrres %>%
  dplyr::filter(Gene %in% atxal_atxaa$Gene) %>%
  dplyr::filter(Gene %in% endosperm$Gene)

all_hybrres <- bind_rows(aa_hybrres, al_hybrres, antag_hybrres)
##########################################################################



# log-trransforming baseMean for MA-plot
os_atxaa <- os_atxaa %>%
  dplyr::mutate(log2_baseMean = log2(baseMean + 1))

os_atxaa <- os_atxaa %>%
  mutate(regulated = case_when(
    
    log2FoldChange < -2 ~ "down",
    log2FoldChange > 2 ~ "up",
    
    Gene %in% all_hybrres$Gene ~ "hybridresponse",
    
    TRUE ~ "no"
  ))

# MA - atxaa
ggplot(os_atxaa, aes(x = log2_baseMean, y = log2FoldChange)) +
  geom_point(aes(color = ifelse(log2FoldChange < -2, "#1E615E", 
                                ifelse(log2FoldChange > 2, "#FF6F3C", "darkgrey"))), 
             size = 0.5, alpha = 0.8) + 
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "log2BaseMean", 
       y = "log2FC", 
       title = "One-sided temperature effect: A.t x A.a") +
  scale_x_continuous(limits = c(2, 20),
                     breaks = seq(2, 20, by = 2)) + 
  scale_y_continuous(breaks = seq(floor(min(os_atxaa$log2FoldChange)), 
                                  ceiling(max(os_atxaa$log2FoldChange)), by = 2)) +
  scale_color_identity() +
  theme_light() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  theme(legend.position = "none")


os_atxal <- os_atxal %>%
  dplyr::mutate(log2_baseMean = log2(baseMean + 1))

# MA - atxal
ggplot(os_atxal, aes(x = log2_baseMean, y = log2FoldChange)) +
  geom_point(aes(color = ifelse(log2FoldChange < -2, "#1E615E", 
                                ifelse(log2FoldChange > 2, "#FF6F3C", "darkgrey"))), 
             size = 0.5, alpha = 0.8) + 
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "log2BaseMean", 
       y = "log2FC", 
       title = "One-sided temperature effect: A.t x A.l") +
  scale_x_continuous(limits = c(2, 20),
                     breaks = seq(2, 20, by = 2)) + 
  scale_y_continuous(breaks = seq(floor(min(os_atxal$log2FoldChange)), 
                                  ceiling(max(os_atxal$log2FoldChange)), by = 2)) +
  scale_color_identity() +
  theme_light() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  theme(legend.position = "none")


## common temperature effect genes
common_temp <- merge(atxaa, atxal, by = "Gene")
common_temp <- common_temp %>%
  dplyr::rename(log2FoldChange_atxaa = log2FoldChange.x,
                log2FoldChange_atxal = log2FoldChange.y)

# filtering check
sum((common_temp$log2FoldChange_atxaa > 2 | common_temp$log2FoldChange_atxaa < -2) & 
      (common_temp$log2FoldChange_atxal > 2 | common_temp$log2FoldChange_atxal < -2))

# setting colors for scatterplot
common_temp <- common_temp %>%
  mutate(regulated = case_when(
    
    (log2FoldChange_atxaa < -2 | log2FoldChange_atxaa > 2) & 
      (log2FoldChange_atxal < -2 | log2FoldChange_atxal > 2) ~ "yes",
    
    Gene %in% all_hybrres$Gene ~ "hybridresponse",
    
    TRUE ~ "no"
  ))


common_temp$color <- with(common_temp, ifelse(regulated == "no", "darkgrey", 
                                              ifelse(regulated == "yes" & log2FoldChange_atxaa > 0, "#FF6F3C", 
                                                     ifelse(regulated == "hybridresponse", "#F2A900",
                                                            ifelse(regulated == "yes" & log2FoldChange_atxaa < 0, "#1E615E", 
                                                                   ifelse(log2FoldChange_atxaa > 2, "#FF6F3C", 
                                                                          ifelse(log2FoldChange_atxal < -2, "#1E615E", "darkgrey")))))))



ggplot(common_temp, aes(x = log2FoldChange_atxaa, y = log2FoldChange_atxal)) +
  geom_point(aes(color = color, shape = regulated, size = regulated), alpha = 0.8) + 
  scale_shape_manual(values = c("yes" = 16, "no" = 16, "hybridresponse" = 17)) +
  scale_size_manual(values = c("yes" = 1.2, "hybridresponse" = 1.4, "no" = 0.5), guide = "none") +
  geom_hline(yintercept = 0, color = "black") + 
  geom_vline(xintercept = 0, color = "black") +
  scale_x_continuous(breaks = seq(floor(min(common_temp$log2FoldChange_atxaa)), ceiling(max(common_temp$log2FoldChange_atxaa)), by = 1)) +
  scale_y_continuous(breaks = seq(floor(min(common_temp$log2FoldChange_atxal)), ceiling(max(common_temp$log2FoldChange_atxal)), by = 1)) +
  labs(
    x = "Log2FC A.t x A.a 26vs18",
    y = "Log2FC A.t x A.l 26vs18"
  ) +
  theme_light() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_color_identity()

# Temperature candidates: common regulation of at least 2 logFC
common_2lfc <- common_temp %>%
  dplyr::filter(log2FoldChange_atxaa < -2 & log2FoldChange_atxal < -2 | log2FoldChange_atxaa > 2 & log2FoldChange_atxal > 2)

endosperm <- read_tsv("input/imprinting/combined-endosperm-filter.tsv")

candidates <- common_2lfc %>%
  dplyr::filter(Gene %in% endosperm$Gene)

write.xlsx(candidates, file = "output/2025/temperatureeffect/filtered_results/2025-06-12_2lfc_temperature_candidates.xlsx")
```
