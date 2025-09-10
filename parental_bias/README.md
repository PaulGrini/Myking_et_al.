# Parental bias analysis using Limma  
[false-mapping.R](false-mapping.R) calculates the ability of each gene in the homozygous samples mapped to hybrid reference genomes, to map to their "correct" parent.  

[prepare-for-limma.md](prepare-for-limma.md) calculates normalization factors for all hybrid samples, formats .tsv files for input to limma, creates model.tsv (like [this](model.tsv) file) and genes-pass-filter.tsv, a list of genes mapping at least 5-fold to the correct parent, based on the files created with false-mapping.R.  

  
Running [run_prepare_run_stats_V2.sh](run_prepare_run_stats_V2.sh) calls on all other scripts (apply_filter_to_counts.py, prepare_heterozygous_counts.py, statistics_from_counts_py and limma.foldchange.R) to perform differential expression analysis  
between the two orthologs. Input files are .tsv files from prepare-for-limma.md: formatted readcounts from hybrid samples, model.tsv and genes-pass-filter.tsv. Output files will be in .csv format. 
