# bootLong
The Block Bootstrap Method Inference for Longitudinal Microbiome Data

`bootLongIndices()` Creates sample indices to construct a moving block bootstrap realization. 

`bootLongPhyloseq()` Creates moving block bootstrap phyloseq realizations: given the block size, this function creates blocks of temporally contiguous observations for each subject by calling `bootLongIndices()`.

`bootLongMethod()` Computes the adjusted p-values and confidence sets: given a block length, this fucntion calls `bootLongPhyloseq()` to create moving block bootstrap phyloseq realizations. Then, it calls ``computeStat`` on each bootstrap realization to compute the statistic and pivotal statistic. Finally, this function computes the adjusted p-values and confidence sets.

`bootLongSubsampling()` Computes the MSE of two-sided probability using subsampling procedure. Given an initial block size and percentage of repeated observations to make subsamples, this function computes the total mean squared error (MSE) in estimating the two-sided probability with different block sizes 1:(lI-1). 

# Installation

```{r}
# Install the the development version from GitHub:
devtools::install_github("PratheepaJ/bootLong")
```

# Workflow
See vignettes for detailed workflow.
