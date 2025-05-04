
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ShennongObject

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/zerostwo/shennong-object/graph/badge.svg)](https://app.codecov.io/gh/zerostwo/shennong-object)
[![R-CMD-check](https://github.com/zerostwo/shennong-object/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zerostwo/shennong-object/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of `ShennongObject` is to provide an S4 class, `Shennong`,
designed for storing and analyzing multi-modal single-cell and bulk
omics data. It builds upon the `MultiAssayExperiment` framework,
offering a structured way to manage multiple assays (e.g., RNA, ATAC),
data layers within assays (e.g., counts, normalized data, scaled data),
dimensional reductions (e.g., PCA, UMAP), sample metadata, and analysis
provenance.

## Installation

To install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("zerostwo/shennong-object")
```

## Quick Start Example

This is a basic example showing how to create a `Shennong` object:

``` r
library(ShennongObject)

# 1. Create dummy RNA counts matrix (genes x samples)
genes <- paste0("Gene", 1:100)
samples <- paste0("Sample", 1:10)
rna_counts <- matrix(
  rpois(1000, lambda = 10),
  nrow = 100,
  dimnames = list(genes, samples)
)
# Make it sparse for demonstration
rna_counts[sample(1000, 300)] <- 0
rna_counts <- as(rna_counts, "dgCMatrix")

# 2. Create dummy sample metadata
sample_metadata <- data.frame(
  condition = sample(c("Control", "Treated"), 10, replace = TRUE),
  batch = factor(sample(c("A", "B"), 10, replace = TRUE)),
  row.names = samples
)

# 3. Create the Shennong object
sn_obj <- sn_initialize_shennong_object(
  counts = rna_counts,
  metadata = sample_metadata,
  project = "QuickStart",
  assay = "RNA" # Name for the primary assay
)

# Print basic summary
print(sn_obj)
#> A Shennong object:  QuickStart 
#> Organism: NA 
#> A Shennong object of 1 listed
#>  experiment with a user-defined name and respective class.
#>  Containing an ExperimentList class object of length 1:
#>  [1] RNA: SummarizedExperiment with 100 rows and 10 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
#> 
#> --- Shennong Specific Slots ---
#> Active Assay: RNA 
#> Active Ident: 1 levels (first few: QuickStart ) 
#> 0 dimensional reduction(s) present: None 
#> Misc data: None 
#> Command history stored: No 
#> Tools data stored: None 
#> Object created with ShennongObject version: 0.1.0
```

For more detailed usage, please see the package vignettes, starting with
the “Quick Start Guide”:

## Code of Conduct

Please note that the ShennongObject project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
