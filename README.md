
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CellPhe

<!-- badges: start -->
<!-- badges: end -->

CellPhe provides functions to accompany the paper (TODO Add reference
when available) to phenotype cells from time-lapse videos.

## Installation

You can install the latest version of CellPhe from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("uoy-research/CellPhe")
```

## Example

Included with the package is an example dataset to demonstrate CellPhe’s
capabilities, this data is available in `data` and comprises 3 parts:

-   The time-lapse stills as TIFF images
    (`data/05062019_B3_3_imagedata`)
-   Existing pre-extracted features
    (`data/05062019_B3_3_Phase-FullFeatureTable.csv`)
-   Region-of-interest boundaries already demarked
    (`data/05062019_B3_3_Phase.zip`)

``` r
library(CellPhe)
trial_name = "05062019_B3_3"
basedir <- "data"
```

The first step is to read in the raw images and normalise the pixel
values to \[0, 255\].

``` r
imagedata =  sprintf("%s/%s_imagedata", basedir, trial_name)
frames = readTiffs(imagedata)
normalised_frames = lapply(frames, normaliseImage, lower = 0, upper = 255)
```

Optionally, if there are already features available for each cell
tracked they can be loaded in and copied. This step isn’t obligatory.
(TODO What format should this CSV be in? 1 row per frame, with at least
FrameId and CellId columns?)

``` r
input_feature_table =  sprintf("%s/%s_Phase-FullFeatureTable.csv", basedir, trial_name)
min_frames = 50
feature_table = copyFeatures(input_feature_table, min_frames)
numcells_over_thresh = feature_table[[1]]
original_IDs = feature_table[[2]]
missing_frames = feature_table[[3]]
features = feature_table[[4]]
```

Given a directory containing ROI files (TODO what format should these be
in?), the `extractFeatures` function can be used to extract time-series
features for each cell.

``` r
roi_files =  sprintf("%s/%s_Phase", basedir, trial_name)
new_features = extractFeatures(roi_files, original_IDs, missing_frames, normalised_frames, min_frames)
#> err epscell, contact the developer please,liujianfei@pku.edu.cn
```

The original features and new features can be combined together in the
`varsFromTimeSeries` functions to extract a final list of variables.
TODO What exactly is this doing differently to extractFeatures?

``` r
tsvariables = varsFromTimeSeries(features, new_features, trial_name, original_IDs)
```

TODO classification/clustering