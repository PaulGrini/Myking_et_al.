### Processing raw limma outputs
### Three filtering steps
### Stacked barplot as summary

library(dplyr)
library(readxl)
library(ggplot2)
library(grid)
library(ggvenn)
library(tidyr)
library(tidyverse)
library(ggrepel)

# Spring cleaning
rm(list = ls())

# loading gene descriptions
descriptions <- read_tsv("input/genedescriptions/all_genedescriptions.tsv")

# loading limma output
atxaa18 <- read.csv("output/2025/limma/atxaa/AtxAa18.counts.csv.filtered.final.csv", header = FALSE)
colnames(atxaa18) <- c("Gene", 
                       "Mat-BR1", "Mat-BR2", "Mat-BR3", "Mat-BR4", 
                       "Pat-BR1", "Pat-BR2", "Pat-BR3", "Pat-BR4",
                       "linenumber", "logFC_limma", "AveExpr", "t.stat", "p_limma", "padj_limma", "Blog-odds", "FC_limma")

# Removing mitochondria and chloroplast genes, as well as small touchups
atxaa18 <- atxaa18 %>%
  dplyr::select(-linenumber) %>%
  dplyr::filter(!str_detect(Gene, "ATM|ATC"))
atxaa18 <- merge(descriptions, atxaa18, by = "Gene")
atxaa18$Description <- trimws(gsub(";\\(.*$", "", atxaa18$Description))

# Do the same thing with all limma outputs
atxaa26 <- read.csv("output/2025/limma/atxaa/AtxAa26.counts.csv.filtered.final.csv", header = FALSE)
colnames(atxaa26) <- c("Gene", 
                       "Mat-BR1", "Mat-BR2", "Mat-BR3", "Mat-BR4", 
                       "Pat-BR1", "Pat-BR2", "Pat-BR3", "Pat-BR4",
                       "linenumber", "logFC_limma", "AveExpr", "t.stat", "p_limma", "padj_limma", "Blog-odds", "FC_limma")

atxaa26 <- atxaa26 %>%
  dplyr::select(-linenumber) %>%
  dplyr::filter(!str_detect(Gene, "ATM|ATC"))
atxaa26 <- merge(descriptions, atxaa26, by = "Gene")
atxaa26$Description <- trimws(gsub(";\\(.*$", "", atxaa26$Description))


atxal18 <- read.csv("output/2025/limma/atxal/AtxAl18.counts.csv.filtered.final.csv", header = FALSE) 
colnames(atxal18) <- c("Gene", 
                       "Mat-BR1", "Mat-BR2", "Mat-BR3", "Mat-BR4", 
                       "Pat-BR1", "Pat-BR2", "Pat-BR3", "Pat-BR4",
                       "linenumber", "logFC_limma", "AveExpr", "t.stat", "p_limma", "padj_limma", "Blog-odds", "FC_limma")

atxal18 <- atxal18 %>%
  dplyr::select(-linenumber) %>%
  dplyr::filter(!str_detect(Gene, "ATM|ATC"))
atxal18 <- merge(descriptions, atxal18, by = "Gene")
atxal18$Description <- trimws(gsub(";\\(.*$", "", atxal18$Description))


atxal26 <- read.csv("output/2025/limma/atxal/AtxAl26.counts.csv.filtered.final.csv", header = FALSE) 
colnames(atxal26) <- c("Gene", 
                       "Mat-BR1", "Mat-BR2", "Mat-BR3", "Mat-BR4", 
                       "Pat-BR1", "Pat-BR2", "Pat-BR3", "Pat-BR4",
                       "linenumber", "logFC_limma", "AveExpr", "t.stat", "p_limma", "padj_limma", "Blog-odds", "FC_limma")

atxal26 <- atxal26 %>%
  dplyr::select(-linenumber) %>%
  dplyr::filter(!str_detect(Gene, "ATM|ATC"))
atxal26 <- merge(descriptions, atxal26, by = "Gene")
atxal26$Description <- trimws(gsub(";\\(.*$", "", atxal26$Description))


##### Adding filtering labels

### STEP 1: thresholds
step_1 <- function(limma_data, threshold_total, threshold_replicate) {
  # Calculate total reads across all samples: columns 5 to 12
  limma_data$total_reads <- rowSums(limma_data[, 5:12], na.rm = TRUE)
  
  # Calculate reads per replicate
  limma_data$replicate_1_reads <- limma_data[, 5] + limma_data[, 9]
  limma_data$replicate_2_reads <- limma_data[, 6] + limma_data[, 10]
  limma_data$replicate_3_reads <- limma_data[, 7] + limma_data[, 11]
  limma_data$replicate_4_reads <- limma_data[, 8] + limma_data[, 12]
  
  # Initialize filter_status column
  limma_data$filter_status <- NA
  
  # Update filter_status based on thresholds
  limma_data$filter_status <- ifelse(
    limma_data$total_reads < threshold_total |
      limma_data$replicate_1_reads < threshold_replicate |
      limma_data$replicate_2_reads < threshold_replicate |
      limma_data$replicate_3_reads < threshold_replicate |
      limma_data$replicate_4_reads < threshold_replicate,
    "low_counts", limma_data$filter_status
  )
  
  # Return the modified data frame, keeping all rows
  # Remove the additional columns used for filtering
  limma_data <- limma_data[, -which(names(limma_data) %in% 
                                      c("total_reads", "replicate_1_reads", 
                                        "replicate_2_reads", "replicate_3_reads", 
                                        "replicate_4_reads"))]
  
  return(limma_data)
}


# Set thresholds: genes with lower read counts will be considered to have low expression
threshold_total <- 100       # This is the threshold for the number of total reads pr cross
threshold_replicate <- 20    # This is the threshold for the number of reads pr replica


atxaa18 <- step_1(atxaa18, threshold_total, threshold_replicate)
atxaa26 <- step_1(atxaa26, threshold_total, threshold_replicate)
atxal18 <- step_1(atxal18, threshold_total, threshold_replicate)
atxal26 <- step_1(atxal26, threshold_total, threshold_replicate)


#### STEP 2: significance. Status "beg" and "weak_bias" are set.
step_2 <- function(df) {

  df$filter_status[df$logFC_limma >= -0.5 & df$logFC_limma <= 0.5 & is.na(df$filter_status)] <- "beg"
  df$filter_status[df$padj_limma > 0.05 & df$logFC_limma > -0.5 & df$logFC_limma < 0.5 & is.na(df$filter_status)] <- "weak_bias"
  df$filter_status[df$padj_limma < 0.05 & is.na(df$filter_status) & 
                     ((df$logFC_limma <= -0.5 & df$logFC_limma >= -1) | 
                        (df$logFC_limma >= 0.5 & df$logFC_limma <= 1))] <- "weak_bias"
  df$filter_status[df$padj_limma > 0.05 & is.na(df$filter_status)] <- "weak_bias"
  
  return(df)
}

atxaa18 <- step_2(atxaa18)
atxaa26 <- step_2(atxaa26)
atxal18 <- step_2(atxal18)
atxal26 <- step_2(atxal26)

#### STEP 3: only endosperm specific genes can be megs
# Load list of endosperm specific genes
endosperm <- read_tsv("input/imprinting/combined-endosperm-filter.tsv")

step_3 <- function(df, filter) {
  # Check that the required columns exist in the dataframe
  if (!("logFC_limma" %in% colnames(df)) || !("Gene" %in% colnames(df))) {
    stop("The dataframe must contain 'logFC_limma' and 'Gene' columns.")
  }
  
  # Update filter_status based on logFC_limma and filter condition 
  # Only for rows that do not already have "beg" or "low_counts" in filter_status
  df$filter_status <- ifelse(
    !(df$filter_status %in% c("beg", "low_counts", "weak_bias")),  # Keep existing
    ifelse(df$logFC_limma <= -1 & df$padj_limma < 0.05,
           ifelse(is.na(df$filter_status), "peg", df$filter_status),  # Set "peg" labels
           ifelse(df$logFC_limma >= 1 & df$padj_limma < 0.05 & df$Gene %in% filter$Gene,
                  ifelse(is.na(df$filter_status), "meg", df$filter_status),  # Set "meg" labels
                  ifelse(is.na(df$filter_status), "not_specific", df$filter_status)  # Set "not_specific" labels for maternally biased genes that are not endosperm enriched
           )
    ),
    df$filter_status  # Do not overwrite any labels already set in previous steps!!
  )
  
  return(df)
}

atxaa18 <- step_3(atxaa18, endosperm)
atxaa26 <- step_3(atxaa26, endosperm)
atxal18 <- step_3(atxal18, endosperm)
atxal26 <- step_3(atxal26, endosperm)

# Check output!!!
sum(atxaa18$filter_status == "peg")
sum(atxaa18$filter_status == "meg")
sum(atxaa18$filter_status == "meg" | atxaa18$filter_status == "peg")
sum(is.na(atxaa18$filter_status))

sum(atxaa26$filter_status == "peg")
sum(atxaa26$filter_status == "meg")
sum(is.na(atxaa26$filter_status))

sum(atxal18$filter_status == "peg")
sum(atxal18$filter_status == "meg")
sum(atxal18$filter_status == "meg" | atxal18$filter_status == "peg")
sum(is.na(atxal18$filter_status))

sum(atxal26$filter_status == "peg")
sum(atxal26$filter_status == "meg")
sum(is.na(atxal26$filter_status))

##### barplot summary
filtered_dfs <- list(
  df1 = atxaa18,
  df2 = atxaa26,
  df3 = atxal18,
  df4 = atxal26
)

# Combine all data frames into one
combined_df <- do.call(rbind, lapply(names(filtered_dfs), function(name) {
  df <- filtered_dfs[[name]]
  df$source <- name
  return(df)
}))

# Count occurrences of each filter_status
status_counts <- combined_df %>%
  group_by(source, filter_status) %>%
  summarise(count = n(), .groups = 'drop')

status_counts$filter_status <- factor(status_counts$filter_status, 
                                      levels = c("low_counts", "beg", "weak_bias", "not_specific", "meg", "peg"))

# Create stacked bar plot
ggplot(status_counts, aes(x = source, y = count, fill = filter_status)) +
  geom_bar(stat = "identity") + 
  labs(title = "", 
       x = "Cross+temperature", 
       y = "Number of genes", 
       fill = "") +
  scale_fill_manual(values = c("lightgrey", "#B2DCD9", "#76A5A3", "#1E615E", "#FFB27A", "#E76F51")) +
  theme_minimal()



