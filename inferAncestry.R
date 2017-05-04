#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)


calculateAncestryPCA <- function(filteredVcf,outPrefix){
library(SNPRelate) 


###This fxn is based on the SNPRelate tutorial (www.codearray.sourceforge.net/tutorials/SNPRelate)
##It accepts as input a .vcf file and outputs 1.) an RDS file and 2.) text file of the genotype PCA results


##Make list of output file names
extensions <- c(".gds", "_pca.rds", '_pca.pdf', "_pca.txt")
outputNames <- sapply(extensions, function(x){return(paste0(outPrefix, x))})
reformatted <- snpgdsVCF2GDS(filteredVcf, outputNames[[1]])


##"open" the output file
genoFile <- snpgdsOpen(reformatted)

###Prune the snps - this picks so as to minimize the effects linkage disequilibirum (ie if 3 snps are ###always inherited together they get disproportional influence relative to how muc hthey should ###inform relatedness if you count them 3x )
snpset <- snpgdsLDpruning(genoFile, ld.threshold=0.2)
snpset.id <- unlist(snpset)

###PCA on pruned set of SNP
pca <- snpgdsPCA(genoFile, snp.id = snpset.id, num.thread=2)

##Write out .rds file with PCA results
saveRDS(pca, outputNames[[2]])

##Pull out PCs 1 and 2 to plot
tab <- data.frame(sample.id = pca$sample.id,
EV1 <- pca$eigenvect[,1],
EV2 <- pca$eigenvect[,2],
  stringsAsFactors = FALSE)

eigenVect <- pca$eigenvect
table <- data.frame(sample.id = pca$sample.id,pca$eigenvect, stringsAsFactors = FALSE)

##Plot pc1 and pc2, labelling point with sample name, write out as a pdf
pdf(outputNames[[3]])
plot(tab$EV1, tab$EV2, xlab = 'PC1', ylab = 'PC2')
text(tab$EV1, tab$EV2, labels = tab$sample.id, cex = .5, pos = 1, offset = .2)
dev.off()

##Write out .txt
write.table(table, outputNames[[4]], quote =FALSE)
return(pca)
}

calculateAncestryPCA(args[1],args[2])
