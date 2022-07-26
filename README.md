
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
-   Region-of-interest (ROI) boundaries already demarked in ImageJ
    format (`data/05062019_B3_3_Phase.zip`)

``` r
library(CellPhe)
trial_name <- "05062019_B3_3"
basedir <- "data"
```

The first step is to read in the raw images.

``` r
imagedata <- sprintf("%s/%s_imagedata", basedir, trial_name)
frames <- readTiffs(imagedata)
```

Optionally, if there are already features available for each cell
tracked they can be loaded in and copied. For example, the features
requiring phase information can be copied in from the feature table
output by PhaseFocus Livecyte software using the function
`copyPhaseFeatures()`. Only the values for the features volume and
sphericity, which require phase information, are copied as all other
features can be calculated by the function `extractFeatures()`. Only
cells that are tracked for a minimum of min_frames are copied into the
new feature table. The output is a data frame with each row
corresponding to a cell tracked in a given frame. Only cells that are
present in the specified minimum number of frames (`min_frames`). are
included. The first 2 columns provide numerical identifiers for the
frame and cell, while any remaining columns are features extracted from
the PhaseFocus feature table (currently just volume and sphericity).

``` r
# COPY PHASE FEATURES FROM PHASEFOCUS FEATURE TABLE FOR EACH CELL TRACKED FOR MORE THAN MIN_FRAMES AND FIND MISSING FRAMES:
min_frames <- 50
input_feature_table <- sprintf("%s/%s_Phase-FullFeatureTable.csv", basedir, trial_name)
feature_table <- copyPhaseFeatures(input_feature_table, min_frames)
```

Alternatively, pre-calculated features from any other tracking software
can be provided as long as it is organised in the same manner, with a
dataframe with `FrameID` and `CellID` as the first 2 columns and any
features present in the remaining columns.

If no pre-existing features are available, a data frame with just
`FrameID` and `CellID` columns will need to be made available.

In addition to any pre-calculated features, the `extractFeatures()`
function reads in information on cell boundaries (`roi_folder`) and
together with the images for every frame (`frames`). The images and
boundary information are used to calculate a total of 72 descriptive
features for each cell on every frame based on size, shape, texture and
movement including the local cell density The output is a dataframe
comprising the `FrameID` and `CellID` columns, the 72 features as
columns, and any additional features that may be present (such as from
`copyPhaseFeatures()`) in further columns.

``` r
roi_files <- sprintf("%s/%s_Phase", basedir, trial_name)
new_features <- extractFeatures(feature_table, frames, roi_files, min_frames, framerate=0.0028)
```

Variables are calculated from the time series for any pre-existing
features as well as the output of `extractFeatures()`, providing both
summary statistics and indicators of time-series behaviour at different
levels of detail obtained via wavelet analysis. 15 summary scores are
calculated for each feature, in addition to the cell trajectory, thereby
resulting in a default output of 1081 features (15x72 + 1). These are
output in the form of a dataframe with the first column being the
`CellID` used previously.

``` r
tsvariables = varsFromTimeSeries(new_features)
```

TODO classification/clustering
