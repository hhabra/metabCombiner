# metabCombiner: LC-MS metabolomics Feature-Matching 

This is an R package for cross-dataset alignment of disparately-acquired untargeted LC-MS metabolomics datasets. metabCombiner takes a pair of peak-picked and aligned untargeted LC-MS datasets and determines overlapping <mass-to-charge (m/z), retention time (rt)>  features, concatenating their measurements to form a combined table of sample measurements.

##installation

```r
#if devtools is not already installed:
install.packages("devtools")

devtools::install_github("hhabra/metabCombiner")
library(metabCombiner)
```
##How to Use R Package

In the [/inst directory](/inst/) directory, we have an example script showing an example of how to use this package. A longer description can be obtained by viewing the vignette, e.g. using browseVignettes("metabCombiner").







