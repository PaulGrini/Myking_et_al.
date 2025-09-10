library(methods)

# Expect 9 columns per line: gene M M M M P P P P
# Four replicates of two conditions.

args<-commandArgs(trailingOnly=TRUE)
gene.counts<-args[1]
gene.sigs<-paste(args[1], "de", sep=".")

dax=read.table(gene.counts, header=FALSE, sep=",")
#da=data.matrix(dax)
#da=dax
pseudocount=0.01 # 0.01 gave astronomical fold changes
pseudocount=1.0 # this was used through 2019
pseudocount=0.0 # as of 2020, inputs are already normalized including pseudocounts
da<-log2(dax[,2:9]+pseudocount)

library(limma)
Group <- factor(c("M","M","M","M", "P","P","P","P"))
design <- model.matrix(~0 + Group)
colnames(design) <- c("M","P")

fit <- lmFit(da[,1:8],design)

contrast.matrix<-makeContrasts(M-P, levels=design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
result1<-topTable(fit2,number=25000)   # use large number so all genes re listed

library(gtools)
x<-logratio2foldchange(result1$logFC, base=2)
res1<-cbind(result1,x)
write.csv(res1,gene.sigs)

# The output is missing the gene name but the first column is the ordinal of the original row.
# cat outfile | tr -d '"' | sort -k1,1n
# Then paste to the input file.

# Output columns:
# gene
# Col,Col,Col,other,other,other - six counts from three biological replicates
# line number - starts at 1 and increases by one
# logFC - estimate of the log2-fold-change corresponding to the effect or contrast
# AveExpr - average log2-expression for the probe over all arrays and channels
# t moderated t-statistic
# P.Value raw p-value
# adj.P.Value adjusted p-value
# B log-odds that the gene is differentially expressed 
# FC - fold change