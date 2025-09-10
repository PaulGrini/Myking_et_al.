### Plot to set thresholds for low counts in limma outputs

library(dplyr)
library(readxl)
library(ggplot2)
library(grid)
library(ggvenn)
library(tidyr)
library(scales)
library(tidyverse)

# Spring cleaning
rm(list = ls())

# Load mapping numbers (made in false-mapping.R)
col_aa <- read_excel("output/2025/controlmapping/2025-05-29_controlmapping_Col_at_aa.xlsx")
col_al <- read_excel("output/2025/controlmapping/2025-05-29_controlmapping_Col_at_al.xlsx")
aa <- read_excel("output/2025/controlmapping/2025-05-29_controlmapping_aa.xlsx")
al <- read_excel("output/2025/controlmapping/2025-05-29_controlmapping_al.xlsx")

# arrange by most reads
col_aa <- col_aa %>%
  mutate(sum_mapped = .[[4]] + .[[5]]) %>%
  arrange(sum_mapped) %>%    
  mutate(row_mapped = row_number())

# For col_al
col_al <- col_al %>%
  mutate(sum_mapped = .[[4]] + .[[5]]) %>%
  arrange(sum_mapped) %>%    
  mutate(row_mapped = row_number())

# For aa
aa <- aa %>%
  mutate(sum_mapped = .[[4]] + .[[5]]) %>%
  arrange(sum_mapped) %>%    
  mutate(row_mapped = row_number())

# For al
al <- al %>%
  mutate(sum_mapped = .[[4]] + .[[5]]) %>%
  arrange(sum_mapped) %>%    
  mutate(row_mapped = row_number())

# arrange by lowest readcount in a single replicate
col_aa <- col_aa %>%
  arrange(lowest18) %>%    
  mutate(row_low18 = row_number()) %>%
  arrange(lowest26) %>%
  mutate(row_low26 = row_number())

col_al <- col_al %>%
  arrange(lowest18) %>%    
  mutate(row_low18 = row_number()) %>%
  arrange(lowest26) %>%
  mutate(row_low26 = row_number())

aa <- aa %>%
  arrange(lowest18) %>%    
  mutate(row_low18 = row_number()) %>%
  arrange(lowest26) %>%
  mutate(row_low26 = row_number())

al <- al %>%
  arrange(lowest18) %>%    
  mutate(row_low18 = row_number()) %>%
  arrange(lowest26) %>%
  mutate(row_low26 = row_number())


### adding fold difference. Putting on some .00 pseudocounts to combat 0s
col_aa <- col_aa %>%
  mutate(fd = (mapped_AT + 0.1) / (mapped_aa + 0.1)) %>%
  arrange(fd) %>%
  mutate(row_fd = row_number())

col_al <- col_al %>%
  mutate(fd = (mapped_AT + 0.001) / (mapped_al + 0.001)) %>%
  arrange(fd) %>%
  mutate(row_fd = row_number())

aa <- aa %>%
  mutate(fd = (mapped_aa + 0.001) / (mapped_AT + 0.001)) %>%
  arrange(fd) %>%
  mutate(row_fd = row_number())

al <- al %>%
  mutate(fd = (mapped_al + 0.001) / (mapped_AT + 0.001)) %>%
  arrange(fd) %>%
  mutate(row_fd = row_number())

ggplot() +
  # Plot for sum_mapped values
  geom_line(data = col_aa, aes(x = row_mapped, y = ifelse(sum_mapped < 1, NA, sum_mapped)), color = "#A0D6D2", linetype = "solid") +  
  geom_line(data = col_al, aes(x = row_mapped, y = ifelse(sum_mapped < 1, NA, sum_mapped)), color = "#A0D6D2", linetype = "dashed") + 
  geom_line(data = aa, aes(x = row_mapped, y = ifelse(sum_mapped < 1, NA, sum_mapped)), color = "#008B77", linetype = "solid") +  
  geom_line(data = al, aes(x = row_mapped, y = ifelse(sum_mapped < 1, NA, sum_mapped)), color = "#008B77", linetype = "dashed") + 
  
  # Plot for lowest18 values
  geom_line(data = col_aa, aes(x = row_low18, y = ifelse(lowest18 < 1, NA, lowest18)), color = "#E76F51", linetype = "solid") +  
  geom_line(data = col_al, aes(x = row_low18, y = ifelse(lowest18 < 1, NA, lowest18)), color = "#E76F51", linetype = "dashed") + 
  geom_line(data = aa, aes(x = row_low18, y = ifelse(lowest18 < 1, NA, lowest18)), color = "#FFA500", linetype = "solid") +  
  geom_line(data = al, aes(x = row_low18, y = ifelse(lowest18 < 1, NA, lowest18)), color = "#FFA500", linetype = "dashed") + 
  
  geom_hline(yintercept = 100, linetype = "dashed", color = "#008B77") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "#FF4500") +
  
  labs(title = "",
       x = "Placement",
       y = "total reads (log scale)") +
  
  scale_y_continuous(trans = 'log10',
                     limits = c(1, NA),        # Limit the y-axis starting at 1
                     breaks = trans_breaks("log10", function(x) 10^x),  # Specify breaks
                     labels = trans_format("log10", math_format(10^.x))) +  # Format labels
  
  coord_cartesian(ylim = c(1, NA)) +  # Set y-axis limits to start at 1
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())
