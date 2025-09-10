# Data analysis

## Hybrid response genes
[hybridresponse.R](hybridresponse.R) loads all total output files from hybrid samples compared to their parents in DESeq2 (18 oC), filters the output based on significance and transgressive expression criteria, and creates MA-plots (figure 1B and 1C, SFigure 4).  

## Temperature effect
[feulgen_statistics.R](feulgen_statistics.R) makes barplots (Figure 5) and performs statistical analysis on endosperm and embryo phenotypes.  
[temperatureeffect.md](temperatureeffect.md) loads all total output files from hybrid samples 26 oC compared to 18 oC, creates and uses filters based on parental temperature effect and significance, and creates MA-plots (Figure 6B, SFigure 13).  

## Parental bias  
[setting-limma-thresholds.R](setting-limma-thresholds.R) determines the minimum count thresholds using output files from [false-mapping.R](../parental_bias/false-mapping.R). [filtering_limma_outputs.R](filtering_limma_outputs.R) uses these thresholds and defines bias status (BEGs, MEGs and PEGs) for all genes in output files from parental bias analysis.
