## Synopsis
This repository includes R functions to implement Gene Set Enrichment Analysis (GSEA) following the algorithm described in *Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles*, Subramanian et al. 2005 PNAS 102(43): 15545-15550.
## Overview of GSEA
Genomewide expression data are increasing in volume and variety. The problem today is in developing powerful computational tools and techniques to analyze these data and in extracting biological insight from such information. GSEA overcomes many analytical challenges (refer article) by evaluating the data at the level of gene sets, that is, groups of genes that share common biological function, chromosomal location or regulation.
### Big picture: GSEA Method...
* Considers genomewide expression profiles from samples belonging to TWO phenotypic classes
* Uses an *a priori* defined set of genes S based on biological knowledge (published information about biochemical pathways or coexpression in previous experiments)
* Computes a ranked list L of genes according to their differential expression between the phenotypic classes using any suitable metric (e.g. correlation) 
* Determines whether the members of a gene set S are randomly distributed throughout L or primarily found at the top or bottom of L in which case the gene set is associated with the phenotypic distinction.
### More details: Three key elements of GSEA method…
#### 1. Calculation of Enrichment Score (ES)
ES of a gene set S represents the degree to which genes in S are overrepresented at the extremes (top or bottom) of the ranked list L. First, contributions of each gene in L are computed with a weighted bump for genes in set S and a constant dip for genes not in set S. The weights for genes in S depend on the magnitude of association of the gene with the phenotype. The ES(S) is the maximum deviation (from zero) of the running sum of these contributions. 
#### 2. Estimation of Significance level of ES(S)
Nominal p-value for each ES(S) is computed using permutation test procedure. The phenotype labels are permuted (thus preserving the complex correlation structure of the gene expression data or the gene-gene correlations) and the L and ES recomputed for each permutation. The empirical, nominal p-value of the observed ES is computed using this null distribution of ES, for each S.
#### 3. Adjustment for Multiple Hypothesis Testing
When multiple gene sets are evaluated, the ES for each gene set S is normalized to account for the size of the set to give normalized ES (NES). For each NES, the false discovery rate (FDR) is computed by comparing the observed and null distributions of NES. The null  distributions of NES are generated using a fixed set of permutations for all S.
## Input Files
Two example input files are part of this repo:
1. leukemia.txt: The first row of this file contains the phenotype labels of the two phenotypic states for *k* samples. The gene expression data (*N* genes x *k* samples) are available from the second row onwards with the gene name in the first column followed by *k* expression values for the samples.
2. pathways .txt: This is a tab delimited file with each row corresponding to a pathway/gene set starting with the pathway/gene set name followed by a tab space, a description of the pathway/gene set followed by a tab space and then the gene names in that pathway/gene set.

These two input files would need to be read into R and converted to a more general form to be passed on as input arguments to the most important functions in this repository.

```R
setwd("...") 	# set path to your working directory/project
source("gsea.R") # make the functions in gsea.R available for R
# Input File 1: expression data file with phenotypic labels in the first row
data <- read.delim("leukemia.txt",header=FALSE,row.names=1,stringsAsFactors = FALSE)
tempData <- get.ExpressionData(data) # get.ExpressionData() is a function in gsea.R
exprData <- tempData$exprData; phenLabels <- tempData$phenLabels
```

The exprData matrix (*N* genes x *k* samples) with gene names as row names and the vector phenLabels (*k* samples) containing the two phenotypic states are required as inputs.

```R
# Input File 2: pathways.txt consisting of gene set information
nCol <- max(count.fields("pathways.txt", sep = "\t"),na.rm=TRUE)
pathways <- read.delim("pathways.txt",header=FALSE,fill=TRUE,col.names=1:nCol)
rownames(pathways) <- pathways[,1]
pathways <- pathways[,-(1:2)] # removing annotations 
```

Pathways matrix with pathway/gene set names as row names, and each row containing names of genes in the corresponding pathway is required as an input.

## The compute.NES() function
The ```compute.NES(...)``` function can be used to get the output statistics described in the article (Subramanian et al. 2005). Expression data, phenotype labels and gene sets are required as input arguments in the format discussed in the previous section on input files. Except these three, all other input arguments have default values assigned. More details on the inputs can be found in gsea.R in the function definition.
The following list of named data objects are returned by ```compute.NES() ```:

***```NES```:** Matrix with rows corresponding to pathway/gene sets given by its’ row names. The columns are ES, permutation based nominal p-value for ES, normalized enrichment score (NES), FWER p-value, FDR q value.

**```rankedL```:** Vector with values of the rank metric (correlation, t-statistic etc) used to measure the association of each gene with the phenotype. The *N* genes in the expression data matrix form the Names of the vector. The vector is sorted in decreasing order of the values of the rank metric.

**```Nperm```:** Scalar giving the number of permutations used in building the null distributions in the analysis. Number of permutations is an input argument with a default of 1000.

**```minGenes```:** Scalar giving the minimum number of genes common to the pathway/gene set and the expression data set. ```minGenes``` is an input argument with a default of 15. Pathways/gene sets with at least ```minGenes``` genes can be selected for robust analyses.

**```rankMetric```:** Scalar with the character value of the rank metric (ex. “Corr” or “t-statistic”). This is an input argument.

**```weightP```:** Scalar giving the value of the weight used in computing the ES. This is an input argument.

**```fixedPermutations```:** Matrix giving the indices for fixed permutations of the phenotype labels in generating null distribution of ES. This matrix has ```Nperm``` columns each with *k* (number of samples) rows giving the specific permutation of pheno labels.

For each pathway of interest, the observed enrichment scores (ES) and permutation based p-value significance are computed. Further, using a fixed set of permutations, a null distribution of enrichment scores is generated to compute the normalized enrichment scores for each gene set. The normalized null distributions are also used to compute the FWER p-values and FDR q values for each gene set.

#### Usage example
```R
myRes <- compute.NES(exprData, phenLabels, pathways, minGenes=15, rankMetric="t-statistic",p=1, nperm=1000, pi=NULL,computeMinpathways=TRUE)
```
The above output objects can be accessed as ```myRes$NES```, ```myRes$rankedL``` etc., and printed to file formats of choice.
