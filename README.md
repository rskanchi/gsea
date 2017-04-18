# gsea
This project includes R functions to implement the gene set enrichment analysis (GSEA) algorithm based on Subramanian et al. 2005 PNAS.

For each pathway of interest, the observed enrichment scores (ES) and permutation based p-value significance are computed. Further, using a fixed set of permutations, a null distribution of enrichment scores is generated to compute the normalized enrichment scores for each gene set. The normalized null distributions are also used to compute the FWER p-values and FDR q values for each gene set.

**** INPUT FILES ****
Two input files are required: gene expression data file and the pathways/gene sets file. Format of the two files can be discussed here.
