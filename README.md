# bootLong
The Block Bootstrap Method Inference for Longitudinal Microbiome Data

`bootLongIndices()` Creates sample indices to construct a moving block bootstrap realization. 

`bootLongPhyloseq()` Creates moving block bootstrap phyloseq realizations: given the block size, this function creates blocks of temporally contiguous observations for each subject by calling `bootLongIndices()`.

`bootLongMethod()` Computes the adjusted p-values and confidence sets: given a block length, this fucntion calls `bootLongPhyloseq()` to create moving block bootstrap phyloseq realizations. Then, it calls ``computeStat`` on each bootstrap realization to compute the statistic and pivotal statistic. Finally, this function computes the adjusted p-values and confidence sets.

`bootLongSubsampling()` Computes the MSE of two-sided probability using subsampling procedure. Given an initial block size and percentage of repeated observations to make subsamples, this function computes the total mean squared error (MSE) in estimating the two-sided probability with different block sizes 1:$l_{I}$. Finally, we choose the optimal block length based on the computed MSE.

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
