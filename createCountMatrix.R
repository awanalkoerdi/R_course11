#Open all count files
count85 <- (read.csv("counts_SRR5832185.txt", header = FALSE, sep = "\t"))
count86 <- (read.csv("counts_SRR5832186.txt", header = FALSE, sep = "\t"))

count94 <- (read.csv("counts_SRR5832194.txt", header = FALSE, sep = "\t"))
count95 <- (read.csv("counts_SRR5832195.txt", header = FALSE, sep = "\t"))

#split the chromosome and gene information
split1 <- do.call(rbind ,strsplit(as.character(count85[,1]), ":", fixed = TRUE))
split2 <- do.call(rbind ,strsplit(as.character(split1[,2]),  ".", fixed = TRUE))

#Add chromosome and gene information to count85 dataframe and removes old column
count85 <- cbind(split2, count85)
count85[,"V1"] <- NULL

#creates count matrix 
countMatrix <- cbind(cbind(count85, count86$V2),cbind(count94$V2, count95$V2))
colnames(countMatrix) <- c("Chromosome", "Gene", "E.coli_85", "E.coli_86", "B.subtillis_94", "B.subtillis_95")


write.csv(countMatrix ,"countMatrix.csv", row.names = FALSE)


