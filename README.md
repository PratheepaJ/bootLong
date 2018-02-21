# bootLong
The Block Bootstrap Method for Longitudinal  Microbiome Data


`SamplingIndices()` Creates moving block bootstrap samples within each subject: given a block length, this function creates blocks of repeated samples from each subject. Then, this function resampling blocks within each subject to create bootstrap samples.

`Bootphyloseq()` Creates moving block bootstrap phyloseq: given the bootstrap samples from `SamplingIndices()`, this method creates bootstrap phyloseq object.

`MBB_method()` Computes the adjusted p-values and confidence sets: given a block length, this fucntion calls `Bootphyloseq()` to create bootstrap samples of phyloseq objects. Then, call an appropriate function on each bootstrap samples to compute corresponding statistic. Finally, this function computes the adjusted p-values and confidence sets.

`MBB_blocklength()` Computes the optimal block length for the given second level parameter (bias or variance or one/two-sided probability) of interest. Given an initial block length and percentage of repeated observations to make subsamples, this function computes the total mean squared error (MSE) in estimating the second level parameter with different block lengths. Finally, we choose the optimal block length based on the computed MSE.

# Installation

```{r}
# Install the the development version from GitHub:
devtools::install_github("PratheepaJ/bootLong")
```

# Usage

```{r}
#   simulate a longitudinal microbiome data 

#   make a sample data

#   make a phyloseq

#   Choose factors for differential analysis
#   factors <- "FACTOR_NAME"

#   Choose taxonomy level for the inference
#   taxlevel <- "Species"

#   Choose FDR
#   FDR <- .1

#   Assign R and RR (for double bootstrap)
#   R <- 200
#   RR <- 50

#   Compute the optimal block length
#   opt.length <- MBB_blocklength(ps,R,FDR,factors,taxlevel="Species",m=1,RR=RR)

#   Optimal block length
#   l.O <- opt.length[2]

#   MBB inference: compute the adjusted p-values and confidence limits
#   mbb.res <- MBB_method(ps,b=l.O,R,FDR,factors,taxlevel="Species",m=1,RR=RR)
```
