#############################
## PCA plots, all samples ##
############################

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

##################################
# Preparation of raw mapping data
##################################

# Spring cleaning
rm(list = ls())

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

# Swap out NA with 0
result_tibble[is.na(result_tibble)] <- 0

# Keep all samples, only remove rows with no reads
hybrids <- result_tibble %>%
  dplyr::filter(rowSums(dplyr::select(., -Gene)) > 0)

colnames(hybrids)

head(hybrids)
tail(hybrids)

########################################
## Merge reads from orthologous genes ##
########################################

# Extract the thaliana genes
ATgenes <- hybrids_18 %>%
  filter(str_detect(Gene, "^AT"))

# Extract the arenosa genes
aagenes <- hybrids_18 %>%
  filter(str_detect(Gene, "^aa"))

# Extract the lyrata genes
algenes <- hybrids_18 %>%
  filter(str_detect(Gene, "^al"))

# Read ortholog keys
ortho_aa <- read_tsv("input/orthologs/At_Aa_orthologs.tsv")
ortho_al <- read_tsv("input/orthologs/At_Al_orthologs.tsv")

# Merge the keys to one df
ortho <- merge(ortho_aa, ortho_al, by = "Ath")
ortho <- ortho %>%
  dplyr::select(Ath, Aaa, Ape)

# Extracting ATG codes from key
ortholist_AT <- unique(ortho$Ath)
ortholist_Aa <- unique(ortho$Aaa)
ortholist_Al <- unique(ortho$Ape)

# Filter on known orthos
ATgenes <- ATgenes[ATgenes$Gene %in% ortholist_AT, ]
aagenes <- aagenes[aagenes$Gene %in% ortholist_Aa, ]
algenes <- algenes[algenes$Gene %in% ortholist_Al, ]

# Translate arenosa and lyrata gene ids into their thaliana ortholog gene id
aagenes$TranslatedGene <- NA
algenes$TranslatedGene <- NA

translation_indices_aa <- match(aagenes$Gene, ortho$Aaa)
translation_indices_al <- match(algenes$Gene, ortho$Ape)

aagenes$TranslatedGene <- ortho$Ath[translation_indices_aa]
algenes$TranslatedGene <- ortho$Ath[translation_indices_al]

# Move "TranslatedGene" column to the left
aagenes <- aagenes[, c("TranslatedGene", names(aagenes)[names(aagenes) != "TranslatedGene"])]
aagenes <- subset(aagenes, select = -Gene)
colnames(aagenes)[colnames(aagenes) == "TranslatedGene"] <- "Gene"

algenes <- algenes[, c("TranslatedGene", names(algenes)[names(algenes) != "TranslatedGene"])]
algenes <- subset(algenes, select = -Gene)
colnames(algenes)[colnames(algenes) == "TranslatedGene"] <- "Gene"
# removing lyrata samples from the arenosa genes and arenosa samples from the lyrata genes, as they are empty
aagenes <- dplyr::select(aagenes, -contains("x-Al"))
algenes <- dplyr::select(algenes, -contains("x-Aa"))

# Now kiss, once. Thaliana columns will have suffix .x, arenosa columns will have suffix .y
merged_genes <- merge(ATgenes, aagenes, by = "Gene", all = TRUE, check.names = FALSE)
colnames(merged_genes)

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

# Remove the columns, so only the merged columns are left
merged_genes <- merged_genes[, !(names(merged_genes) %in% remove_cols)]

colnames(merged_genes)

head(merged_genes)
tail(merged_genes)

# Now kiss, twice. Lyrata columns have the suffix .y
merged_genes <- merge(merged_genes, algenes, by = "Gene", all = TRUE, check.names = FALSE)

colnames(merged_genes)

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

# Remove the columns, left are just merged columns
merged_genes <- merged_genes[, !(names(merged_genes) %in% remove_cols)]

colnames(merged_genes)

head(merged_genes)
tail(merged_genes)

# Swap out NA with 0 just in case
merged_genes[is.na(merged_genes)] <- 0

# Rearrange 
merged_genes <- merged_genes %>%
  dplyr::select(Gene, 18:24, 2:17, 25:36)


###############
## Metadata ##
##############

# Read full metadata file and make sure R knows how to read each column.
metadata <- read_excel("input/metadata/temperaturehybrids-metadata.xlsx") %>% as.data.frame()

row.names(metadata) <- metadata$Sample_name
metadata$Temperature <- as.factor(metadata$Temperature)
metadata$Replicates <- as.factor(metadata$Replicates)
metadata$Sample_type1 <- as.factor(metadata$Sample_type1)
metadata$Sample_type2 <- as.factor(metadata$Sample_type2)

# Extract count data, or just load of already obtained!
count_data <- merged_genes %>%
  column_to_rownames(var = "Gene") %>%
  as.matrix()

# Check that the metadata corresponds with the samples
all(rownames(metadata) == colnames(count_data))

# does it give 'FALSE'? Run this, which tells you which rows and columns are not matching.
#error_indices <- which(!rownames(metadata) == colnames(count_data))

# Create a DESeq2 dataset. In design, the last argument is what will be compared in differential expression. Defined data is count data and metadata. 
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, 
                              design=~Replicates+Temperature)

###############################
### Perform DESeq2 analysis ###
###############################

# Calculate the row sums (total counts per gene) and filter out rows with low counts
min_counts <- 1  # Define the minimum count threshold
keep_genes <- rowSums(counts(dds) >= min_counts) > 1
dds <- dds[keep_genes, ]

# Perform DESeq2 analysis
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

# Variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)


### MAKE PCAs

# Extract PCA data
pcaData <- plotPCA(vsd, intgroup = "Temperature", returnData = TRUE)

# PCA
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = Temperature, shape = Sample_type2)) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_color_manual(values = c("#76A5A3", "#FFB27A")) +
  scale_shape_manual(values = c(16, 1, 15, 0, 17)) +
  scale_x_continuous(breaks = c(seq(-30, 40, by = 10), 0))

# only 18
pcaData_filtered <- pcaData[pcaData$Temperature != 26, ]

percentVar <- round(100 * attr(pcaData_filtered, "percentVar"))
ggplot(pcaData_filtered, aes(x = PC1, y = PC2, color = Sample_type2)) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_color_manual(values = c("#1E615E", "#76A5A3", "#FFB27A", "#E76F51", "gold")) +
  scale_x_continuous(breaks = c(seq(-30, 40, by = 10), 0))

