
# Prepping hybrid samples for limma
Input files: ReadsPerGene.out.tab files from mapping with STAR  

## Gather all samples and calculate normalization factors  

```R
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

# Spring cleaning
rm(list = ls())

##################################
# Preparation of raw mapping data
##################################

# Define the path to the directory containing your files
path_to_files <- "input/readsprgene/"

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

# NA to 0
result_tibble[is.na(result_tibble)] <- 0

# Extract the right samples and genes.
hybrids <- result_tibble %>% 
  filter(str_detect(Gene, "^(?i)a")) %>%
  dplyr::select(Gene, dplyr::matches("At_aa|At_al")) 

colnames(hybrids)
head(hybrids)
tail(hybrids)


### Normalization factor for limma
# Calculate total reads for each sample
total_reads <- colSums(hybrids[, -1], na.rm = TRUE)

# Calculate the average reads across every four samples
average_cross <- sapply(seq(1, length(total_reads), by=4), function(i) {
  mean(total_reads[i:min(i+3, length(total_reads))], na.rm = TRUE)
})
average_cross <- rep(average_cross, each = 4)


# Create new rows for the DataFrame
df <- rbind(
  total_reads,
  average_cross
)

df <- as.data.frame(df)

norm_factor <- df[2, ] / df[1, ]

df_normhybrids <- rbind(df, norm_factor)
rownames(df_normhybrids)[3] <- 'norm_factor'
```

## Make .tsv input files fit for limma

```R
# Extracting reads 
ATgenesaa <- hybrids %>%
  filter(str_detect(Gene, "^AT")) %>%
  dplyr::select(-contains("At-x-Al"))

aagenes <- hybrids %>%
  filter(str_detect(Gene, "^aa")) %>%
  dplyr::select(-contains("At-x-Al"))

ATgenesal <- hybrids %>%
  filter(str_detect(Gene, "^AT")) %>%
  dplyr::select(-contains("At-x-Aa"))

algenes <- hybrids %>%
  filter(str_detect(Gene, "^al")) %>%
  dplyr::select(-contains("At-x-Aa"))


# Import ortolog overview
ortho_aa <- read_tsv("input/orthologues/At_Aa_orthologues.tsv")
ortho_al <- read_tsv("input/orthologues/At_Al_orthologues.tsv")



## ARENOSA HYBRIDS
# Extracting ATG codes from the ortho-list
ortholist_ATaa <- unique(ortho_aa$Ath)
ortholist_Aa <- unique(ortho_aa$Aaa)

# Filter on known orthos
ATgenesaa <- ATgenesaa[ATgenesaa$Gene %in% ortholist_ATaa, ]
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
aagenes <- aagenes[aagenes$Gene %in% ATgenesaa$Gene, ]


######################################################
###### Making tsvs for each sample - arenosa #########
######################################################

create_sample_dfs_legacy_with_saving <- function(ATgenesaa, aagenes, output_dir = "input/limma") {
  # Get the sample names (excluding the 'Gene' column)
  sample_names <- colnames(ATgenesaa)[-1]
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (sample in sample_names) {
    # Extract counts from ATgenesaa
    at_counts <- data.frame(
      gene = ATgenesaa$Gene,
      allele = "At",
      count = ATgenesaa[[sample]]
    )
    
    # Extract counts from aagenes
    aa_counts <- data.frame(
      gene = aagenes$Gene,
      allele = "Aa",
      count = aagenes[[sample]]
    )
    
    # Combine both counts into one data frame
    combined_counts <- rbind(at_counts, aa_counts)
    
    # Order the combined data frame alphabetically by gene
    combined_counts <- combined_counts[order(combined_counts$gene), ]
    
    # Assign the new data frame to a variable with the name of the sample
    assign(sample, combined_counts, envir = .GlobalEnv)
    
    # Save the combined data frame as a .tsv file
    write.table(combined_counts, 
                file = file.path(output_dir, paste0(sample, ".tsv")), 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
  }
}

# Run the function to create new data frames and save them as .tsv files for the ATgenesaa and aagenes
create_sample_dfs_legacy_with_saving(ATgenesaa, aagenes)



##### LYRATA HYBRIDS
# Extracting ATG codes from the ortho-list
ortholist_ATal <- unique(ortho_al$Ath)
ortholist_Al <- unique(ortho_al$Ape)

# Filter on known orthos
ATgenesal <- ATgenesal[ATgenesal$Gene %in% ortholist_ATal, ]
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

algenes <- algenes[algenes$Gene %in% ATgenesal$Gene, ]



######################################################
###### Making tsvs for each sample - lyrata ##########
######################################################

create_sample_dfs_with_saving <- function(ATgenesal, algenes, output_dir = "input/limma") {
  # Get the sample names (excluding the 'Gene' column)
  sample_names <- colnames(ATgenesal)[-1]
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (sample in sample_names) {
    # Extract counts from ATgenesal
    at_counts <- data.frame(
      gene = ATgenesal$Gene,
      allele = "At",
      count = ATgenesal[[sample]]
    )
    
    # Extract counts from algenes
    al_counts <- data.frame(
      gene = algenes$Gene,
      allele = "Al",
      count = algenes[[sample]]
    )
    
    # Combine both counts into one data frame
    combined_counts <- rbind(at_counts, al_counts)
    
    # Order the combined data frame alphabetically by gene
    combined_counts <- combined_counts[order(combined_counts$gene), ]
    
    # Assign the new data frame to a variable with the name of the sample
    assign(sample, combined_counts, envir = .GlobalEnv)
    
    # Save the combined data frame as a .tsv file
    write.table(combined_counts, 
                file = file.path(output_dir, paste0(sample, ".tsv")), 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
  }
}

# Run the function to create new data frames and save them as .tsv files
create_sample_dfs_with_saving(ATgenesal, algenes)
```

  
## Make model.tsv for each hybrid using normalization factors calculated in first step

```R
norm_factor_aa <- as.numeric(df_normhybrids[3, 1:8])

model_aa <- data.frame(
  cross = c(rep("AtxAa18", 4), rep("AtxAa26", 4)),
  mat = rep("At", 8),
  pat = rep("Aa", 8),
  matmult = rep(0.5, 8),
  replicate = rep(colnames(df_normhybrids)[-length(colnames(df_normhybrids))], length.out = 8)
)

model_aa$norm <- norm_factor_aa

print(model_aa)

write.table(model_aa, file = "input/limma/atxaa/model.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


norm_factor_al <- as.numeric(df_normhybrids[3, 9:16])

model_al <- data.frame(
  Cross = c(rep("AtxAl26", 4), rep("AtxAl18", 4)),
  Mat = rep("At", 8),
  Pat = rep("Al", 8),
  MatMult = rep(0.5, 8),
  Replicate = rep(colnames(df_normhybrids)[9:16], length.out = 8)
)

model_al$norm <- norm_factor_al

print(model_al)
write.table(model_al, file = "input/limma/atxal/model.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
```

## Genes pass filter: only include genes able to map to the correct parent 
Input: mapping rates from homozygous samples mapped to hybrid genomes

```R
col_atxaa <- read_excel("output/2025/controlmapping/2025-04-05_controlmapping_Col_at_aa.xlsx") %>% as.data.frame()
aa_atxaa <- read_excel("output/2025/controlmapping/2025-03-11_controlmapping_aa.xlsx") %>% as.data.frame()
col_atxal <- read_excel("output/2025/controlmapping/2025-04-05_controlmapping_Col_at_al.xlsx") %>% as.data.frame()
al_atxal <- read_excel("output/2025/controlmapping/2025-03-11_controlmapping_al.xlsx") %>% as.data.frame()

col_atxaa <- col_atxaa %>%
  filter(mapped_AT >= 83.33333)

aa_atxaa <- aa_atxaa %>%
  filter(mapped_aa >= 83.33333)

col_atxal <- col_atxal %>%
  filter(mapped_AT >= 83.33333)

al_atxal <- al_atxal %>%
  filter(mapped_al >= 83.33333)


genes_pass_filter_aa <- intersect(col_atxaa$Gene, aa_atxaa$Gene)
genes_pass_filter_al <- intersect(col_atxal$Gene, al_atxal$Gene)

write.table(genes_pass_filter_aa, file = "input/limma/atxaa/genes_pass_filter.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(genes_pass_filter_al, file = "input/limma/atxal/genes_pass_filter.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

