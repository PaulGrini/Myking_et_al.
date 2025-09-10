### Feulgen seeds: barplot and statistical tests

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)

# Spring cleaning
rm(list = ls())

# Read input file (long format)
df <- read_excel("input/feulgen/feulgen_stats.xlsx", sheet = "Endosperm")

# total counts pr cross
df_counts <- as.data.frame(table(df$Cross, df$Endosperm))
colnames(df_counts) <- c("Cross", "Endosperm", "Count")
df_counts$total <- ave(df_counts$Count, df_counts$Cross, FUN = sum)

# Calculating the percentage of each stage
df_counts$percentage <- df_counts$Count / df_counts$total * 100
df_counts$Endosperm <- factor(df_counts$Endosperm, levels = c("Complete", "Partially complete", "Periferal", "Micropylar", "Syncytial"))

# Read input file 2 (long format)
df_2 <- read_excel("input/feulgen/feulgen_stats.xlsx", sheet = "Embryo")

# total counts pr cross
df2_counts <- as.data.frame(table(df_2$Cross, df_2$Embryo))
colnames(df2_counts) <- c("Cross", "Embryo", "Count")
df2_counts$total <- ave(df2_counts$Count, df2_counts$Cross, FUN = sum)

# Calculating the percentage of each stage
df2_counts$percentage <- df2_counts$Count / df2_counts$total * 100
df2_counts$Embryo <- factor(df2_counts$Embryo, levels = c("Trans", "Globular", "Octant"))

df_counts$Type <- "EN"
df2_counts$Type <- "EM"

# Renaming colums prior to merge
df_counts_long <- df_counts %>% dplyr::rename(Category = Endosperm)
df2_counts_long <- df2_counts %>% dplyr::rename(Category = Embryo)

# Combine the dfs
combined_df <- bind_rows(df_counts_long, df2_counts_long)

combined_df <- combined_df %>%
  mutate(Cross = recode(Cross,
                        'Al18' = 'Lyrata 18',
                        'Al26' = 'Lyrata 26',
                        'AtxAl18' = 'A.t x A.l 18',
                        'AtxAl26' = 'A.t x A.l 26',
                        'Col18' = 'Col-0 18',
                        'Col26' = 'Col-0 26'))

# Set the correct order
combined_df$Cross <- factor(combined_df$Cross, levels = c("A.t x A.l 18", "A.t x A.l 26", "Col-0 18", "Col-0 26", "Lyrata 18", "Lyrata 26"))
combined_df$Type <- factor(combined_df$Type, levels = c("EN", "EM"))

# Stick with the pretty colors
palette <- c("#495044", "#1E615E", "#76A5A3", "#B2DCD9", "#E3F2F1", "#FFE5CC", "#FFB27A", "#FF9655")

# Plot
ggplot(combined_df, aes(x = Type, y = percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "", x = "", fill = "") +
  theme_light() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  scale_fill_manual(values = palette) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  facet_wrap(~Cross, ncol = length(unique(combined_df$Cross)), scales = "free_x")


#### Statistics

# Make table for chi-test
EN <- table(df$Cross, df$Endosperm)

# Big Chi? Simulate p-value to deal with all the 0 in the table
cs_EN <- chisq.test(EN, simulate.p.value = TRUE)
print(cs_EN)

# Pairwise Chi-squared tests
pairwise_chisq_test <- function(contingency_table, p.adjust.method = "bonferroni") {
  groups <- rownames(contingency_table)
  results <- list()
  
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      # Create a temporary table for the specific pairs
      temp_table <- contingency_table[c(groups[i], groups[j]), , drop = FALSE]
      
      # Perform Chi-squared test
      test <- chisq.test(temp_table)
      
      # Store results in a df
      results[[paste(groups[i], "vs", groups[j])]] <- list(statistic = test$statistic,
                                                           p.value = test$p.value)
    }
  }
  result_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(statistic = x$statistic, p.value = x$p.value)
  }))
  
  # Adjust p-values
  result_df$adjusted_p.value <- p.adjust(result_df$p.value, method = p.adjust.method)
  
  return(as.data.frame(result_df))
}

# run pairwise chi function, endosperm
EN_pairs <- pairwise_chisq_test(EN)
print(EN_pairs)



# Same thing for embryo data
EM <- table(df_2$Cross, df_2$Embryo)

cs_EM <- chisq.test(EM, simulate.p.value = TRUE)
print(cs_EM)

# run pairwise chi function, embryo
EM_pairs <- pairwise_chisq_test(EM)
print(EM_pairs)


