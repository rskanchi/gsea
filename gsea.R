# File contains functions to implement Gene Set Enrichment Analysis (GSEA)
# Mon/Year: Apr/2017
# Author: R Kanchi
# Assignment for the Genialis Bioinformatics position

# compute.NES computes the computes the normalized enrichment scores, FWER p-value, FDR q value
# inputs are the expression data (exprData), phenotypic labels (phenLabels for ex. AML, ALL..), pathways file 
# with pathway names and the gene members, minimum number of genes in a pathway that should be available in the 
# expression data (minGenes), rank metric (t-statistic or correlation etc) to rank the genes 
# based on the association of the expr data with the phenotype, weight of the step (p default 1) to control the 
# step size of enrichment score when gene in a pathway is "hit", number of permutations (nperm default 1000) in generating 
# the null distribution of enrichment scores, the fixed permutations matrix (pi). If pi is not supplied to 
# the function as an input a matrix is generated within the function which will be returned as part of 
# the function output for reproducibility of results

compute.NES <- function(exprData, phenLabels, pathways, minGenes=25, rankMetric="t-statistic",p=1,
                        nperm=1000, fixPerm=TRUE, pi=NULL,computeMinpathways=FALSE){
  permRes <- perm.ES(exprData, phenLabels, pathways, minGenes, rankMetric, p, nperm, fixPerm,pi,computeMinpathways)
  minPathways <- permRes$minPathways; NgeneSets <- permRes$NgeneSets 
  L <- permRes$L; observed.ES <- permRes$observedES
  pi <- permRes$pi; piES <- permRes$piES
  
  # compute normalized ES based on the +ve or -ve signs of the ES values..divide by the associated mean
  message("Normalizing the observed and null distribution of ES using the positive and negative parts separately")
  posMeans <- sapply(1:NgeneSets,FUN=function(x) mean(piES[x,][piES[x,]>=0]))
  negMeans <- sapply(1:NgeneSets,FUN=function(x) abs(mean(piES[x,][piES[x,]<0])))
  #observed normalized ES: NES values - function output!
  norm.ES <- sapply(1:NgeneSets,function(x) 
    ifelse(observed.ES[x]>=0,observed.ES[x]/posMeans[x],observed.ES[x]/negMeans[x]))

  # the normalized ES values for the fixed-permutatations based piES
  norm.piES <- sapply(1:nperm,FUN=function(x){
                  sapply(1:NgeneSets,FUN=function(y) 
                    ifelse(piES[y,x]>=0,piES[y,x]/posMeans[y],piES[y,x]/negMeans[y])) })
  # compute FWER p-value - function output!
  message("computing the FWER p-values")
  posExtreme <- sapply(1:nperm,FUN=function(x) {y <- norm.piES[,x]; ifelse(length(y[y >= 0])>0,max(y[y >= 0]),NA) })
  posExtreme <- posExtreme[!is.na(posExtreme)]
  negExtreme <- sapply(1:nperm,FUN=function(x) {y <- norm.piES[,x]; ifelse(length(y[y < 0])>0,min(y[y < 0]),NA) })
  negExtreme <- negExtreme[!is.na(negExtreme)]
  FWERp <- sapply(1:NgeneSets,function(x) ifelse(norm.ES[x]>=0, sum(posExtreme >= norm.ES[x])/length(posExtreme),
                                                 sum(negExtreme <= norm.ES[x])/length(negExtreme)))
  # compute FDR q-value - function output!
  message("computing the FDR q-values")
  null.numer <- sapply(1:NgeneSets,FUN=function(x){ y <- norm.piES[x,]
                    ifelse(norm.ES[x] >= 0, sum(y[y>=0] >= norm.ES[x])/length(y[y>=0]),
                          sum(y[y<0] <= norm.ES[x])/length(y[y<0]))})

  obs.denom <- sapply(1:NgeneSets,FUN=function(x) ifelse(norm.ES[x]>=0,
                                    sum(norm.ES[norm.ES>=0]>=norm.ES[x])/length(norm.ES[norm.ES>=0]),
                                    sum(norm.ES[norm.ES<0]<=norm.ES[x])/length(norm.ES[norm.ES<0]) ))
  FDRq <- null.numer/obs.denom
  resMat <- matrix(c(norm.ES,FWERp,FDRq),nr=NgeneSets,nc=3)
  rownames(resMat) <- rownames(minPathways)
  colnames(resMat) <- c("NES","FWER-pval","FDR-qval")
  resMat <- resMat[order(resMat[,"NES"],decreasing=TRUE),]
  return(list("NormalizedES"=resMat,"observedES"= observed.ES,"rankedL"=L,"Nperm"=nperm,"minGenes"=minGenes,
              "rankMetric"=rankMetric,"weightP"=p,"fixedPermutations"=pi))  # returning fixed permutations for reproducibility of results
} # end of compute.NES

###############################################################################
# function to get the pathways with at least a specified number (minGenes) of genes in each pathway that
# are in the gene list (genes for which expression data are available for two phenotype groups
get.minPathways <- function(pathways,geneList,minGenes=25){ 
  # compute the number of genes in each gene set common with the gene list 
  Ncommon <-  sapply(1:nrow(pathways),FUN=function(x) { 
    length(intersect(pathways[x,][pathways[x,] != ""],geneList)) })
  # determine the gene sets with at least minGenes genes common with the gene list
  minPathways <- pathways[which(Ncommon >= minGenes),]
  NgeneSets <- nrow(minPathways)
  message("Number of pathways with atleast ", minGenes," genes common with the gene list is ",NgeneSets)
  return(list("minPathways"=minPathways,"NgeneSets"=NgeneSets))
} #end of function get.minPathways

################################################################################
# function to estimate significance of the observed ES using permutations
# inputs are expression data (exprData), phenotypic labels (phenLabels for ex. AML, ALL..),
# gene sets or pathways, minimum number of genes in a pathway that should be available in the 
# expression data (minGenes),rankMetric, weight (p) to control the step size of the "hit" genes while computing ES
# number of permutations (nperm), logical (TRUE/FALSE) on whether to fix the permutations across the pathways
# fixperm = TRUE for computing NES,FWER and FDR; and FALSE for permutations based p value significance of ES
# matrix of fixed permutaions (pi) can be supplied which will be computed if not
perm.ES <- function(exprData, phenLabels, pathways, minGenes=25, rankMetric="t-statistic", p=1, 
                    nperm=1000, fixPerm=FALSE, pi=NULL,computeMinpathways=FALSE){
  if (computeMinpathways){
        minGenePathways <- get.minPathways(pathways,rownames(exprData),minGenes)
        minPathways <- minGenePathways$minPathways
        NgeneSets <- minGenePathways$NgeneSets
        rm(minGenePathways)
  } else {minPathways <- pathways; NgeneSets <- nrow(minPathways)}
  # observed ES for the pathways
  message("Computing the observed ranked list of genes (L)")
  L <- compute.L(exprData, phenLabels,rankMetric)
  message("Computing the observed enrichment scores ES")
  observed.ES <- sapply(1:NgeneSets,FUN=function(y) {Sy <- minPathways[y,][minPathways[y,] != ""]
  compute.ES(Sy,L,p)})
  names(observed.ES) <- rownames(minPathways)

  # a matrix of nperm permutations of pheno labels
  k <- length(phenLabels)
  # the fixed phenotype permutations for all S if fixPerm is TRUE (required for computing NES, FWER, FDR)
  # the fixPerm would need to be FALSE if computing permutations based p-value for the observed ES
  if (fixPerm & is.null(pi)) pi <- replicate(nperm,sample(k,k)) 
  
  # compute ES for all genesets for each permutation (fixed or otherwise)
  message("Starting permutations..buiding a null distribution of enrichment scores")
  piES <- sapply(1:nperm,FUN=function(x){
    if (fixPerm) pheno.x <- phenLabels[pi[,x]] else pheno.x <- phenLabels[sample(k,k)]
    Lx <- compute.L(exprData,pheno.x,rankMetric)
    sapply(1:NgeneSets,FUN=function(y) {Sy <- minPathways[y,][minPathways[y,] != ""]
    compute.ES(Sy,Lx,p)})
  }) # end of computing piES (matrix of nrow = NgeneSets and ncol = nperm)

  # p-value using the positive or negative Es values depending on the sign of the observed ES
  if (!(fixPerm)){  pval <- sapply(1:NgeneSets, FUN=function(x){
      if (observed.ES[x] >= 0) { signed.nullES <- piES[x,][which(piES[x,] >= 0)]
          pval <- sum(signed.nullES >= observed.ES[x])/length(signed.nullES)
          } else {
                    signed.nullES <- piES[x,][which(piES[x,] <= 0)]
                    pval <- sum(signed.nullES <= observed.ES[x])/length(signed.nullES)
                  } # end of if-else    
    }) # end of computing pval
    resMat <- matrix(c(observed.ES,pval),nc=2)
    rownames(resMat) <- rownames(minPathways); colnames(resMat) <- c("ES","Perm-pval")
  } # end of if in fixPerm

  # if fixPerm is TRUE, returning two matrices
  # pi (for reproducibility of results - the permutations were fixed for all S) and piES (for further computations) 
  if (fixPerm)  return(list("observedES"= observed.ES,"rankedL"=L,"minPathways"=minPathways,"NgeneSets"=NgeneSets,
                            "Nperm"=nperm,"rankMetric"=rankMetric,"weightP"=p,"pi"=pi,"piES"=piES))
  else return(list("ES"=resMat,"rankedL"=L,"minPathways"=minPathways,
                   "NgeneSets"=NgeneSets,"Nperm"=nperm,"rankMetric"=rankMetric,"weightP"=p))  
} # end of function perm.ES

################################################################################
# compute.ES computes the enrichment score ES for a given gene set S using
# the ranked list of genes L consisting of the ranking metric (t-statistic or corr etc) and 
# the weight p to control the step size when a gene in S is hit during the walk (scan) from top to bottom of L 
compute.ES <- function(S,L,p=1){
  N <- length(L) # number of genes in the gene list
  NH <- length(S)  # number of genes in the gene set S
  contrib <- rep(-1/(N-NH), N)  # vector of contributions from each gene not in the set S
  # now replacing the values in the positions of genes in set S with their contributions 
  S.indices <- which(names(L) %in% S) # positions of S-genes in L
  NR <-  sum(abs(L[S.indices])^p)
  contrib[S.indices] <- abs(L[S.indices])^p/NR
  Esi <- cumsum(contrib)
  #plot(Esi,type = "l")
  computedES <- Esi[which.max(abs(Esi))] # max deviation of Phit-Pmiss from 0; i.e. max abs value with its' sign
  return(computedES)
} # end of function compute.ES

################################################################################
# compute.L computes a ranked list of genes using three input arguments
# 1. the expression data matrix (N genes x k samples and no pheno labels) 
# to rank order the N genes in it to form L = {g1,g2,...gN} according to a
# 2. ranking metric (t-statistic, correlation etc) where the gene expression is the response variable and 
# 3. the phenotype (here AML/ALL) is the explanatory variable
# using switch so that ranking metric(s) can be added in future 

compute.L <- function(exprData,phenLabels,rankMetric="t-statistic"){
    L <- switch(rankMetric,
           "t-statistic" = apply(exprData,1,function(x) t.test(x~phenLabels)$statistic),
           "corr" = apply(exprData,1,function(x) { summ <- summary(lm(x~phenLabels)) 
                                                   sign(summ$coefficients[2,1])*sqrt(summ$r.squared) })
             ) # end of switch
  L <- sort(L,decreasing=TRUE) # sorted in decreasing order to get the ranked list of genes
  return(L)
} # end of compute.L

################################################################################
# function to extract the expression data, phenotype labels and gene labels separately from D
format.D <- function(data){
  colnames(data) <- paste("p",1:ncol(data),sep="") 
  phenLabels <- unlist(data[1,])
  data <- as.matrix(data[-1,])
  class(data) <- "numeric"
  geneLabels <- rownames(data) # apply(phenLabels,1,table) ALL 47; AML 25
  return(list("exprData"=data,"phenLabels"=phenLabels,"geneLabels"=geneLabels))
} # end of function format.D







