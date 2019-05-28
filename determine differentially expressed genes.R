# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
# 
# install.packages("openxlsx")
# install.packages("ggplot2")

library(openxlsx)
library(edgeR)
library(ggplot2)


differentially_expressed_genes <- function(){
  # open file
  countMatrix <- read.csv("countMatrix.csv", header = TRUE)
  rownames(countMatrix) <- countMatrix[,"Gene"]
  
  # create DGE list
  exp   <- c("E.coli", "E.coli", "B.subtillis", "B.subtillis")
  group <- factor(exp)
  y     <- DGEList(counts = countMatrix[,3:6], group = group)
  
  # select all genes with at least 50 counts per million (cpm) in two samples
  keep.genes <- rowSums(cpm(y) > 50) >= 2
  y <- y[keep.genes,]
  
  # recalculate the library size
  y$samples$lib.size <- colSums(y$counts)

  # normalization
  y <- calcNormFactors(y, method = "TMM")
  
  # create design matrix (samples grouped by conditions)
  design <- model.matrix(~0+group, data = y$samples)
  colnames(design) <- levels(y$samples$group)
  print(design)
  
  # estimate dispersion
  y <- estimateGLMCommonDisp( y, design)
  y <- estimateGLMTrendedDisp(y, design, method = "power")
  y <- estimateGLMTagwiseDisp(y, design)
  
  # plot normalized data
  pdf("Normalization_Results.pdf")
  plotMDS(y)
  plotBCV(y)
  dev.off()
  
  # Determine differentially expressed genes
  fit <- glmFit(y, design)
  mc  <- makeContrasts(exp.r=E.coli-B.subtillis, levels = design)
  fit <- glmLRT(fit, contrast = mc)
  res <- topTags(fit, n = 100000, p.value = 0.05)
  print(res)
  
  result_df <- as.data.frame(topTags(fit, n = 100000, p.value = 0.05))
  
  write.csv(res ,"diff_expressed_genes.csv", row.names = TRUE)
  
}


main <- function(){
  differentially_expressed_genes()
 
  
}

main()

