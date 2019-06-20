# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
# 
# install.packages("openxlsx")
# install.packages("ggplot2")

library(openxlsx)
library(edgeR)
library(ggplot2)
library(data.table)
library(grid)
library(gridExtra)

#Function to create countmatrix from seperate text files.
create_countmatrix <- function(){
  
  #Open all count files
  count85 <- (read.csv("counts_SRR5832185.txt", header = FALSE, sep = "\t"))
  count86 <- (read.csv("counts_SRR5832186.txt", header = FALSE, sep = "\t"))
  count87 <- (read.csv("counts_SRR5832187.txt", header = FALSE, sep = "\t"))
  
  count94 <- (read.csv("counts_SRR5832194.txt", header = FALSE, sep = "\t"))
  count95 <- (read.csv("counts_SRR5832195.txt", header = FALSE, sep = "\t"))
  count96 <- (read.csv("counts_SRR5832196.txt", header = FALSE, sep = "\t"))
  
  
  #split the chromosome and gene information
  split1 <- do.call(rbind ,strsplit(as.character(count85[,1]), ":", fixed = TRUE))
  split2 <- do.call(rbind ,strsplit(as.character(split1[,2]) , ".", fixed = TRUE))
  
  #Add chromosome and gene information  to count85 dataframe and removes old column
  count85 <- cbind(split2, count85)
  count85[,"V1"] <- NULL
  
  
  #creates count matrix 
  countMatrix <- do.call(cbind, list(count85, count86$V2, count87$V2, count94$V2, count95$V2, count96$V2))
  colnames(countMatrix) <- c("Chromosome", "Gene", "E.coli_85", "E.coli_86", "E.coli_87", "B.subtillis_94", "B.subtillis_95", "B.subtillis_96")
  
  write.csv(countMatrix ,"countMatrix.csv", row.names = FALSE)
  return(countMatrix)
}


filter_and_normalize <- function(countMatrix){
  # Create rownames
  rownames(countMatrix) <- countMatrix[,"Gene"]
  
  # create DGE list
  exp   <- c("E.coli", "E.coli", "E.coli", "B.subtillis", "B.subtillis", "B.subtillis")
  group <- factor(exp)
  y     <- DGEList(counts = countMatrix[,3:8], group = group)
  
  # select all genes with at least 50 counts per million (cpm) in two samples
  keep.genes <- rowSums(cpm(y) > 50) >= 2
  y <- y[keep.genes,]

  # recalculate the library size
  y$samples$lib.size <- colSums(y$counts)
  
  # normalization
  y <- calcNormFactors(y, method = "TMM")
  return(y)
}


dispersion <- function(y){
  # create design matrix (samples grouped by conditions)
  design <- model.matrix(~0+group, data = y$samples)
  colnames(design) <- levels(y$samples$group)
  
  # estimate dispersion
  y <- estimateGLMCommonDisp( y, design)
  y <- estimateGLMTrendedDisp(y, design, method = "power")
  y <- estimateGLMTagwiseDisp(y, design)
  
  #create PCA plot
  df <- as.data.frame(y$counts)
  df_pca <- prcomp(df)
  
  #create new dataframe to revise pca plot
  df_out <- as.data.frame(df_pca$x)
  exp    <- c("E.coli", "E.coli", "E.coli", "B.subtillis", "B.subtillis", "B.subtillis")
  df_out$exp <- sapply(strsplit(row.names(df), "_"), "[[", 1 )
  
  p <- ggplot(df_out,aes(x=PC1,y=PC2, color=exp ))
  p <- p+geom_point()
  
  # plot normalized data
  pdf("Normalization_Results.pdf") 
  plotMDS(y)
  plotBCV(y)
  plot(df_pca$x[,1], df_pca$x[,2], main="PCA plot",
       xlab="PCA1",ylab="PCA2")

  dev.off()
  
  return(design)
}

# Determine differentially expressed genes
determine_sig_genes <- function(design){
  fit <- glmFit(y, design)
  mc  <- makeContrasts(exp.r=E.coli-B.subtillis, levels = design)
  fit <- glmLRT(fit, contrast = mc)
  res <- topTags(fit, n = 100000, p.value = 0.05)
  result_df <- as.data.frame(topTags(fit, n = 100000, p.value = 0.05))
  
  #write csv file 
  write.csv(res ,"diff_expressed_genes.csv", row.names = TRUE)

  return(res)
}

#creates a text file out of a csv file. 
create_gsea_input <- function(res){ 
  subSet <<- as.data.frame(res[,4:5])
  setDT(subSet, keep.rownames = TRUE)[]
  
  write.table(subSet, "gsea_input.txt", sep = "\t", row.names = FALSE)
}


main <- function(){
  countMatrix <- create_countmatrix()
  y <- filter_and_normalize(countMatrix)
  design <- dispersion(y)
  res <- determine_sig_genes(design)
  create_gsea_input(res)
}

main()
