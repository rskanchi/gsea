# Use this file to input the data files, and to run the functions in gsea.R that implement the GSEA algorithm

################### DATA INPUT ###############################

setwd("") # set the working directory here containing the input files and the R functions file gsea.R

# two input data files: pathways.txt and leukemia.txt
# Input File 1: leukemia.txt consisting of gene expression data
D <- read.delim("leukemia.txt",header=FALSE,row.names=1,stringsAsFactors = FALSE)

# Input File 2: pathways.txt consisting of gene set information
nCol <- max(count.fields("pathways.txt", sep = "\t"),na.rm=TRUE) # 449
pathways <- read.delim("pathways.txt",header=FALSE,fill=TRUE,col.names=1:nCol)
#dim(pathways) # 522 449
pathwayInfo <- pathways[,1:2]; colnames(pathwayInfo) <- c("pathway","Annotation") # saving the annotations separately
rownames(pathways) <- pathways[,1]
pathways <- pathways[,-(1:2)] # removing annotations 
# dim(pathways) # 522 447

################# END OF DATA INPUT #########################

source("gsea.R") # source the functions file

# get the expression data and phenotypic labels separately for ease of work
tempD <- format.D(D)
exprData <- tempD$exprData; phenLabels <- tempD$phenLabels

# get pathways with at least "minGenes" genes common with the expression data gene list 
minGenes <- 50
minGenePathways <- get.minPathways(pathways,rownames(exprData),minGenes=minGenes)
minPathways <- minGenePathways$minPathways
# 39 pathways/ gene sets with at least 50 genes from the expr data

# normalized scores
nperm <- 100
testNES <- compute.NES(exprData, phenLabels, pathways=minPathways, minGenes=minGenes, rankMetric="t-statistic",
                       p=1, nperm=nperm, fixPerm=TRUE, pi=NULL,computeMinpathways=FALSE)
# p-values based on random permutations
testPerm <- perm.ES(exprData, phenLabels, pathways=minPathways, minGenes=minGenes, rankMetric="t-statistic",
                    p=1, nperm=nperm, fixPerm=FALSE,pi=NULL,computeMinpathways=FALSE)

# put the NES+FWER+FDR and permutations-based p value together
output <- merge(testNES$NormalizedES,testPerm$ES,by="row.names")
output <- output[order(output[,"NES"],decreasing = TRUE),]
output[,c(2,5)] <- round(output[,c(2,5)],3)
write.csv(output,"Output.csv",quote=FALSE)
