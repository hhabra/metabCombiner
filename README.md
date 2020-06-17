# metabCombiner

This is an R package for aligning a pair of disparately-acquired untargeted LC-MS metabolomics . metabCombiner takes peak-picked and conventionally aligned untargeted LC-MS datasets and determines the overlapping <mass-to-charge (m/z), retention time (rt)>  features, concatenating their measurements to form a combined table of sample mass spectral measurements.

## installation

```r
#if devtools is not already installed:
install.packages("devtools")

devtools::install_github("hhabra/metabCombiner")
library(metabCombiner)
```

## How to Use metabCombiner

In the [/inst directory](/inst/) directory, we have an example metabCombiner script showing a demonstration of the package's utilities. See the vignette, e.g. type browseVignettes("metabCombiner"), for a more detailed description of how to use metabCombiner.







