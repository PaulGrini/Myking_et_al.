#!/bin/sh

#SBATCH --account=nn9525k
#SBATCH --job-name=Prep&Stats
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G  # 16 GB total
#SBATCH --cpus-per-task=4  # 4 cpu is optimal for 4 threads

ACCOUNT=nn9525k

echo MODULES
module --force purge
module load StdEnv
module load GCC/11.3.0
module load Python/3.10.4-GCCcore-11.3.0
# module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2 # incompatible with GCC 11.3.0
module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1
module list

p=python3
s='.'  # for now, read scripts from this directory

# This script will generate csv files.
# There shouldn't be any csv files in this directory.
# If there are, they are assumed to be input and lead to redundant outputs.
# In case we rerun this script, we need to get rid of the previously generated csv files.
echo
echo clean up
rm -v *.log
rm -v *.csv

echo
echo prepare
# Collate all replicates into one row per gene like MAT,MAT,MAT,PAT,PAT,PAT.
# Generate several files like MATxPAT.counts.csv WITH
# AT4G38740.1,26,36,39,33,37,30
# AT1G10760.1,87,343,286,134,291,274
# ...
$p ${s}/prepare_heterozygous_counts.py --debug model.tsv
echo -n $?
echo " exit status"


echo
echo filter
function filter() {
    echo "Filter $1 with gene list $2"
    # For Yuri, here we filtered out plastid genes by ATM and ATC IDs.
    # Otherwise, just filter by counts.
    ${p} ${s}/apply_filter_to_counts.py ${2} ${1} > ${1}.filtered
    echo "Statistics $1"
    ${p} ${s}/statistics_from_counts.py --debug ${s}/limma.foldchange.r ${1}.filtered
}
filter MxS.counts.csv genes_pass_filter.tsv
filter SxM.counts.csv genes_pass_filter.tsv

echo
echo clean up
#rm -v *.json
#rm -v *.noplastid
#rm -v *.counts.csv
#rm -v *.filtered

# statistics_from_counts.py leaves behind timestamps in *json files
echo Uncomment these lines to remove unnecessary files.
#echo clean up
#rm -v *.json
#rm -v *.noplastid
#rm -v *.counts.csv
#rm -v *.filtered


echo
echo "Look here for P-value statistics."
ls -l *.final.csv

echo
echo done
