# bootLong

bootLong implements the block bootstrap method inference for longitudinal microbiome data. In addition, bootLong provides  tools to explore the sampling schedule and within-subject dependence.

bootLong strats with a `phyloseq` object and returns the confidence intervals and adjusted p-values for the differential abundance analysis. 

The input `phyloseq` must have the count data. 


##  Installation

Install the development version from GitHub with:
```{r}
# install.packages("devtools")
devtools::install_github("PratheepaJ/bootLong")
```

## Workflow
See [vignettes](https://github.com/PratheepaJ/bootLong/blob/master/vignettes/bootLong.Rmd) for the detailed workflow.
