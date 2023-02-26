# Differential-gene-expression-analysis

## 🐥 Chicken Analysis Reports
Organism: Gallus gallus

Tissues: Intestine, Meat, Reproductive Tissues, Macrophage and Monocyte

Total Samples: 2 replicates for each tissue

## Method

For each pairwise analysis, we filtered out non-expressed and low-expressed genes.

The surrogate variables (batch effects) were estimated by the SVA package (Leek JT et al., 2022)
and adopted in the DE analysis model. The DE analysis and normalization data were performed
using limma package (Ritchie ME et al., 2015). We considered a DE significant gene showed FDR < 0.05.

Note: 
-	Surrogate variables are covariates constructed directly from high-dimensional data (like gene expression/RNA sequencing/methylation/brain imaging data) that can be used in subsequent analyses to adjust for unknown, unmodeled, or latent sources of noise. 
-	n.sv is the number of latent factors that need to be estimated (number of surrogate variables) 

## Visualization

- PCA
- Volcano Plot
