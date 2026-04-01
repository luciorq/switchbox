# switchbox

<!-- badges: start -->
[![R-CMD-check](https://github.com/luciorq/switchbox/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/luciorq/switchbox/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

K-Top-Scoring-Pair (KTSP) binary classifiers based on pairwise feature
comparisons. The KTSP algorithm uses the relative ordering of feature values
(e.g. gene expression levels) to make predictions, producing robust and
transparent decision rules.

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("luciorq/switchbox")
```

## Quick Start

``` r
library(switchbox)

# Load example breast cancer data
data(trainingData)
data(testingData)

# Train a KTSP classifier
classifier <- swap_train_ktsp(matTraining, trainingGroup, k_range = 3:8)
classifier$TSPs

# Classify test samples
predictions <- swap_ktsp_classify(matTesting, classifier)
table(predictions, testingGroup)

# Evaluate performance
swap_prediction_stats(predictions, testingGroup)
```

## Features

- Train KTSP classifiers with automatic K selection (`swap_train_ktsp()`)
- Wilcoxon-based feature filtering (`swap_filter_wilcoxon()`)
- Restricted pair selection for incorporating prior knowledge
- Tie-aware scoring
- K-fold and leave-one-out cross-validation (`swap_ktsp_cv()`, `swap_ktsp_loo()`)
- ggplot2-based diagnostic plots (ROC curves, vote heatmaps, gene pair plots)
- C-backed signed score computation for performance

## Citation

If you use switchbox in published research, please cite:

> Afsari B, Fertig E, Geman D, Marchionni L. switchBox: an R package for
> k-Top Scoring Pairs classifier development. Bioinformatics. 2015 Jan
> 15;31(2):273-4. doi:10.1093/bioinformatics/btu622

## License

GPL-2 | GPL-3
