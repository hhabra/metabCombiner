# metabCombiner

This is an R package for aligning a pair of disparately-acquired untargeted LC-MS metabolomics. *metabCombiner* takes peak-picked and conventionally aligned untargeted LC-MS datasets and determines the overlapping <mass-to-charge (m/z), retention time (rt)>  features, concatenating their measurements to form a combined table of sample mass spectral measurements.

## installation

```r
#if devtools is not already installed:
install.packages("devtools")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

#for R versions 4.0 and later
devtools::install_github("hhabra/metabCombiner", build_vignettes = TRUE)
library(metabCombiner)

#for R versions between 3.5 and 4
devtools::install_github("hhabra/metabCombiner", build_vignettes = TRUE, ref = "R3.5")
```

## How to Use metabCombiner

In the [/demo](/demo/) directory, we have an example R script showing a demonstration of the package's utilities. See the vignette by entering browseVignettes("metabCombiner") for a more detailed description of how to use *metabCombiner*.







