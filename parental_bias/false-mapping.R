#############################################################
## Mapping-check: are genes mapping to the correct parent? ##
## Example script for A. thaliana mapped to hybrid genomes ##
#############################################################

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(DESeq2)
library(pheatmap)
library(glmpca)
library(apeglm)
library(genefilter)
library(ggrepel)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(openxlsx)
library(scales)

# Spring cleaning
rm(list = ls())

##################################
# Preparation of raw mapping data
##################################

# Define the path to the directory containing your files
path_to_files <- "input/readsprgene_controltest/"

# List all the files matching your pattern
files <- list.files(path = path_to_files, pattern = "*.tab", full.names = TRUE)

# Function to read a file and extract the sample name
read_sample_file <- function(file) {
  # Extract samples based on common file naming (default, but not necessary if all files from a directory are included anyway)
  sample_name <- str_replace(basename(file), "\\_ReadsPerGene\\.out\\.tab$", "")
  # Read the file with appropriate column names
  df <- read_tsv(file, col_names = c("Gene", sample_name), col_types = cols(.default = "c"))
  df <- df[1:2]
  df[[sample_name]] <- as.numeric(df[[sample_name]])
  return(df)
}

# Read all files and store them in a list
list_of_dfs <- lapply(files, read_sample_file)

# Reduce the list of data frames into a single data frame by merging on the "Gene" column
merged_df <- purrr::reduce(list_of_dfs, full_join, by = "Gene")

# Convert the result to a tibble
result_tibble <- as_tibble(merged_df)

# View the resulting tibble
print(result_tibble)
colnames(result_tibble)

# Find Genes with row sum equal to 0
result_tibble %>%
  filter(rowSums(select(., -Gene)) == 0) # %>% write_csv("Genes_0_count.csv")

result_tibble[is.na(result_tibble)] <- 0

write.xlsx(result_tibble, file = "2024-11-21_doumapping_controls.xlsx", rowNames = FALSE)


#### A. thaliana mapped to A.t x A.a genome ######
At_aa <- result_tibble %>%
  dplyr::filter(grepl("^AT|^aa", Gene)) %>%
  dplyr::select(Gene, dplyr::matches("^2.*Aa-")) %>%
  dplyr::filter(rowSums(dplyr::select(., -Gene)) > 0)

tail(At_aa)



ATgenes <- At_aa %>%
  filter(str_detect(Gene, "^AT"))

aagenes <- At_aa %>%
  filter(str_detect(Gene, "^aa"))

ortho_aa <- read_tsv("input/orthologues/At_Aa_orthologues.tsv")
ortho_al <- read_tsv("input/orthologues/At_Al_orthologues.tsv")

# Extracting ATG codes from the ortho-list
ortholist_AT <- unique(ortho_aa$Ath)
ortholist_Aa <- unique(ortho_aa$Aaa)

# Filter on known orthos
ATgenes <- ATgenes[ATgenes$Gene %in% ortholist_AT, ]
aagenes <- aagenes[aagenes$Gene %in% ortholist_Aa, ]

# Translate aa_jg into their corresponding ATGs
aagenes$TranslatedGene <- NA

# Use match() to find the indices of the translation in "ortho_aa" for each gene id in "aagenes$Gene"
translation_indices <- match(aagenes$Gene, ortho_aa$Aaa)

# Using the translation indices, update the "TranslatedGene" column with the corresponding Ath values
aagenes$TranslatedGene <- ortho_aa$Ath[translation_indices]

# Move "TranslatedGene" column to the very left
aagenes <- aagenes[, c("TranslatedGene", names(aagenes)[names(aagenes) != "TranslatedGene"])]
aagenes <- subset(aagenes, select = -Gene)
colnames(aagenes)[colnames(aagenes) == "TranslatedGene"] <- "Gene"

# Now kiss
merged_genes <- merge(ATgenes, aagenes, by = "Gene", all = TRUE, check.names = FALSE)
merged_genes[is.na(merged_genes)] <- 0

colnames(merged_genes)

results_At_aa <- merged_genes %>%
  mutate(
    x_sum = rowSums(dplyr::select(., ends_with(".x"))),
    y_sum = rowSums(dplyr::select(., ends_with(".y"))),
    total = x_sum + y_sum,
    x_percentage = (x_sum / total) * 100,
    y_percentage = (y_sum / total) * 100
  ) %>%
  dplyr::select(Gene, x_percentage, y_percentage, x_sum, y_sum)

colnames(results_At_aa)[2:5] <- c("percent_AT", "percent_aa", "mapped_AT", "mapped_aa")

### for lowest number of counts from one replicate: merging parental reads 
# making new, merged columns based on .x and .y
columns_x <- grep("\\.x$", colnames(merged_genes), value = TRUE)
columns_y <- grep("\\.y$", colnames(merged_genes), value = TRUE)
base_names_x <- sub("\\.x$", "", columns_x)
base_names_y <- sub("\\.y$", "", columns_y)
base_names <- intersect(base_names_x, base_names_y)
for (base_name in base_names) {
  merged_genes[[base_name]] <- merged_genes[[paste0(base_name, '.x')]] + merged_genes[[paste0(base_name, '.y')]]
}

# Get the column names that end with ".x" or ".y"
remove_cols <- grep("\\.x$|\\.y$", names(merged_genes), value = TRUE)

# Remove the columns, so only merged columns are left. 
merged_genes <- merged_genes[, !(names(merged_genes) %in% remove_cols)]

colnames(merged_genes)

merged_genes$lowest18 <- apply(merged_genes[, 2:4], 1, min)
merged_genes$lowest26 <- apply(merged_genes[, 5:7], 1, min)
lowest_counts <- merged_genes %>%
  dplyr::select(Gene, lowest18, lowest26)

results_At_aa <- merge(results_At_aa, lowest_counts, by ="Gene")

print(sum(results_At_aa$mapped_AT > 20, na.rm = TRUE))


##### A. thaliana mapped to A.t x A.l genome ########

At_al <- result_tibble %>%
  dplyr::filter(grepl("^AT|^al", Gene)) %>%
  dplyr::select(Gene, dplyr::matches("^3.*Al-")) %>%
  dplyr::filter(rowSums(dplyr::select(., -Gene)) > 0)

tail(At_al)



ATgenes <- At_al %>%
  filter(str_detect(Gene, "^AT"))

algenes <- At_al %>%
  filter(str_detect(Gene, "^al"))

# Extracting ATG codes from the ortho-list
ortholist_AT <- unique(ortho_al$Ath)
ortholist_Al <- unique(ortho_al$Ape)

# Filter on known orthos
ATgenes <- ATgenes[ATgenes$Gene %in% ortholist_AT, ]
algenes <- algenes[algenes$Gene %in% ortholist_Al, ]

# Translate aa_jg into their corresponding ATGs
algenes$TranslatedGene <- NA

# Use match() to find the indices of the translation in "ortho_aa" for each gene id in "aagenes$Gene"
translation_indices <- match(algenes$Gene, ortho_al$Ape)

# Using the translation indices, update the "TranslatedGene" column with the corresponding Ath values
algenes$TranslatedGene <- ortho_al$Ath[translation_indices]

# Move "TranslatedGene" column to the very left
algenes <- algenes[, c("TranslatedGene", names(algenes)[names(algenes) != "TranslatedGene"])]
algenes <- subset(algenes, select = -Gene)
colnames(algenes)[colnames(algenes) == "TranslatedGene"] <- "Gene"

# Now kiss
merged_genes <- merge(ATgenes, algenes, by = "Gene", all = TRUE, check.names = FALSE)
merged_genes[is.na(merged_genes)] <- 0

colnames(merged_genes)


results_At_al <- merged_genes %>%
  mutate(
    x_sum = rowSums(dplyr::select(., ends_with(".x"))),
    y_sum = rowSums(dplyr::select(., ends_with(".y"))),
    total = x_sum + y_sum,
    x_percentage = (x_sum / total) * 100,
    y_percentage = (y_sum / total) * 100
  ) %>%
  dplyr::select(Gene, x_percentage, y_percentage, x_sum, y_sum)

colnames(results_At_al)[2:5] <- c("percent_AT", "percent_al", "mapped_AT", "mapped_al")

print(sum(results_At_al$mapped_AT > 10, na.rm = TRUE))

### for lowest number of counts from one replicate: merging parental reads 
# making new, merged columns based on .x and .y
columns_x <- grep("\\.x$", colnames(merged_genes), value = TRUE)
columns_y <- grep("\\.y$", colnames(merged_genes), value = TRUE)
base_names_x <- sub("\\.x$", "", columns_x)
base_names_y <- sub("\\.y$", "", columns_y)
base_names <- intersect(base_names_x, base_names_y)
for (base_name in base_names) {
  merged_genes[[base_name]] <- merged_genes[[paste0(base_name, '.x')]] + merged_genes[[paste0(base_name, '.y')]]
}

# Get the column names that end with ".x" or ".y"
remove_cols <- grep("\\.x$|\\.y$", names(merged_genes), value = TRUE)

# Remove the columns, so only merged columns are left. 
merged_genes <- merged_genes[, !(names(merged_genes) %in% remove_cols)]

colnames(merged_genes)

merged_genes$lowest18 <- apply(merged_genes[, 2:4], 1, min)
merged_genes$lowest26 <- apply(merged_genes[, 5:7], 1, min)
lowest_counts <- merged_genes %>%
  dplyr::select(Gene, lowest18, lowest26)

results_At_al <- merge(results_At_al, lowest_counts, by ="Gene")




###################################################
# Thaliana_arenosa

At_aa <- result_tibble %>%
  dplyr::filter(grepl("^AT|^aa", Gene)) %>%
  dplyr::select(Gene, dplyr::ends_with("At_aa")) %>%
  dplyr::select(Gene, dplyr::matches("Col")) %>%
  dplyr::filter(rowSums(dplyr::select(., -Gene)) > 0)

head(At_aa)
tail(At_aa)



ATgenes <- At_aa %>%
  filter(str_detect(Gene, "^AT"))

aagenes <- At_aa %>%
  filter(str_detect(Gene, "^aa"))

# Extracting ATG codes from the ortho-list
ortholist_AT <- unique(ortho_aa$Ath)
ortholist_Aa <- unique(ortho_aa$Aaa)

# Filter on known orthos
ATgenes <- ATgenes[ATgenes$Gene %in% ortholist_AT, ]
aagenes <- aagenes[aagenes$Gene %in% ortholist_Aa, ]

# Translate aa_jg into their corresponding ATGs
aagenes$TranslatedGene <- NA

# Use match() to find the indices of the translation in "ortho_aa" for each gene id in "aagenes$Gene"
translation_indices <- match(aagenes$Gene, ortho_aa$Aaa)

# Using the translation indices, update the "TranslatedGene" column with the corresponding Ath values
aagenes$TranslatedGene <- ortho_aa$Ath[translation_indices]

# Move "TranslatedGene" column to the very left
aagenes <- aagenes[, c("TranslatedGene", names(aagenes)[names(aagenes) != "TranslatedGene"])]
aagenes <- subset(aagenes, select = -Gene)
colnames(aagenes)[colnames(aagenes) == "TranslatedGene"] <- "Gene"

# Now kiss
merged_genes <- merge(ATgenes, aagenes, by = "Gene", all = TRUE, check.names = FALSE)
merged_genes[is.na(merged_genes)] <- 0

colnames(merged_genes)


results_At_aa <- merged_genes %>%
  mutate(
    x_sum = rowSums(dplyr::select(., ends_with(".x"))),
    y_sum = rowSums(dplyr::select(., ends_with(".y"))),
    total = x_sum + y_sum,
    x_percentage = (x_sum / total) * 100,
    y_percentage = (y_sum / total) * 100
  ) %>%
  dplyr::select(Gene, x_percentage, y_percentage, x_sum, y_sum)

colnames(results_At_aa)[2:5] <- c("percent_AT", "percent_aa", "mapped_AT", "mapped_aa")

print(sum(results_At_aa$mapped_AT < 15, na.rm = TRUE))

### for lowest number of counts from one replicate: merging parental reads 
# making new, merged columns based on .x and .y
columns_x <- grep("\\.x$", colnames(merged_genes), value = TRUE)
columns_y <- grep("\\.y$", colnames(merged_genes), value = TRUE)
base_names_x <- sub("\\.x$", "", columns_x)
base_names_y <- sub("\\.y$", "", columns_y)
base_names <- intersect(base_names_x, base_names_y)
for (base_name in base_names) {
  merged_genes[[base_name]] <- merged_genes[[paste0(base_name, '.x')]] + merged_genes[[paste0(base_name, '.y')]]
}

# Get the column names that end with ".x" or ".y"
remove_cols <- grep("\\.x$|\\.y$", names(merged_genes), value = TRUE)

# Remove the columns, so only merged columns are left. 
merged_genes <- merged_genes[, !(names(merged_genes) %in% remove_cols)]

colnames(merged_genes)

merged_genes$lowest18 <- apply(merged_genes[, 2:5], 1, min)
merged_genes$lowest26 <- apply(merged_genes[, 6:8], 1, min)
lowest_counts <- merged_genes %>%
  dplyr::select(Gene, lowest18, lowest26)

results_At_aa <- merge(results_At_aa, lowest_counts, by ="Gene")



write.xlsx(results_At_aa, file = "output/2025/controlmapping/2025-05-29_controlmapping_Col_at_aa.xlsx", rowNames = FALSE)

## Thaliana_lyrata
At_al <- result_tibble %>%
  dplyr::filter(grepl("^AT|^al", Gene)) %>%
  dplyr::select(Gene, dplyr::ends_with("At_al")) %>%
  dplyr::select(Gene, dplyr::matches("Col")) %>%
  dplyr::filter(rowSums(dplyr::select(., -Gene)) > 0)

head(At_al)
tail(At_al)



ATgenes <- At_al %>%
  filter(str_detect(Gene, "^AT"))

algenes <- At_al %>%
  filter(str_detect(Gene, "^al"))

# Extracting ATG codes from the ortho-list
ortholist_AT <- unique(ortho_al$Ath)
ortholist_Al <- unique(ortho_al$Ape)

# Filter on known orthos
ATgenes <- ATgenes[ATgenes$Gene %in% ortholist_AT, ]
algenes <- algenes[algenes$Gene %in% ortholist_Al, ]

# Translate aa_jg into their corresponding ATGs
algenes$TranslatedGene <- NA

# Use match() to find the indices of the translation in "ortho_aa" for each gene id in "aagenes$Gene"
translation_indices <- match(algenes$Gene, ortho_al$Ape)

# Using the translation indices, update the "TranslatedGene" column with the corresponding Ath values
algenes$TranslatedGene <- ortho_al$Ath[translation_indices]

# Move "TranslatedGene" column to the very left
algenes <- algenes[, c("TranslatedGene", names(algenes)[names(algenes) != "TranslatedGene"])]
algenes <- subset(algenes, select = -Gene)
colnames(algenes)[colnames(algenes) == "TranslatedGene"] <- "Gene"

# Now kiss
merged_genes <- merge(ATgenes, algenes, by = "Gene", all = TRUE, check.names = FALSE)
merged_genes[is.na(merged_genes)] <- 0

colnames(merged_genes)


results_At_al <- merged_genes %>%
  mutate(
    x_sum = rowSums(dplyr::select(., ends_with(".x"))),
    y_sum = rowSums(dplyr::select(., ends_with(".y"))),
    total = x_sum + y_sum,
    x_percentage = (x_sum / total) * 100,
    y_percentage = (y_sum / total) * 100
  ) %>%
  dplyr::select(Gene, x_percentage, y_percentage, x_sum, y_sum)

colnames(results_At_al)[2:5] <- c("percent_AT", "percent_al", "mapped_AT", "mapped_al")

print(sum(results_At_al$mapped_AT < 15, na.rm = TRUE))

### for lowest number of counts from one replicate: merging parental reads 
# making new, merged columns based on .x and .y
columns_x <- grep("\\.x$", colnames(merged_genes), value = TRUE)
columns_y <- grep("\\.y$", colnames(merged_genes), value = TRUE)
base_names_x <- sub("\\.x$", "", columns_x)
base_names_y <- sub("\\.y$", "", columns_y)
base_names <- intersect(base_names_x, base_names_y)
for (base_name in base_names) {
  merged_genes[[base_name]] <- merged_genes[[paste0(base_name, '.x')]] + merged_genes[[paste0(base_name, '.y')]]
}

# Get the column names that end with ".x" or ".y"
remove_cols <- grep("\\.x$|\\.y$", names(merged_genes), value = TRUE)

# Remove the columns, so only merged columns are left. 
merged_genes <- merged_genes[, !(names(merged_genes) %in% remove_cols)]

colnames(merged_genes)

merged_genes$lowest18 <- apply(merged_genes[, 2:5], 1, min)
merged_genes$lowest26 <- apply(merged_genes[, 6:8], 1, min)
lowest_counts <- merged_genes %>%
  dplyr::select(Gene, lowest18, lowest26)

results_At_al <- merge(results_At_al, lowest_counts, by ="Gene")


write.xlsx(results_At_al, file = "output/2025/controlmapping/2025-05-29_controlmapping_Col_at_al.xlsx", rowNames = FALSE)
