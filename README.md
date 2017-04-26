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
1. leukemia.txt: The first row of this file contains the phenotype labels for the two phenotypic states for *k* samples. The gene expression data (*N* genes and *k* samples) are available from the second row onwards with the gene name in the first column followed by *k* expression values for the samples.
2. pathways .txt: This is a tab delimited file with each row corresponding to a pathway/gene set starting with the pathway/gene set name followed by a tab space, a description of the pathway/gene set followed by a tab space and then the gene names in that pathway/gene set. 
## R functions in this repository
Coming Soon..
## Output
For each pathway of interest, the observed enrichment scores (ES) and permutation based p-value significance are computed. Further, using a fixed set of permutations, a null distribution of enrichment scores is generated to compute the normalized enrichment scores for each gene set. The normalized null distributions are also used to compute the FWER p-values and FDR q values for each gene set.
