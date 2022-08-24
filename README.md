
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
```

The first step is to prepare a dataframe containing metadata and any
pre-existing attributes. If PhaseFocus Livecyte or Trackmate software
has been used to generate the region-of-interest (ROI) files, then a
helper function is available to create the required metadata format:
`copyFeatures`. The dataframe format comprises each row corresponding to
a cell tracked in a given frame, indexed by columns `FrameID` and
`CellID` which contain numerical identifiers (NB: `FrameID` must be in
ascending chronological order). The only other required field is
`ROI_filename`, which specifies the filename of the ROI file
corresponding to the frame-cell combination. Any features can be
provided in additional columns, `copyFeatures` returns volume and
sphericity from PhaseFocus software.

The example below creates the metadata dataframe from a PhaseFocus
experimental setup, only including cells that were tracked for at least
50 frames.

``` r
trial_name <- "05062019_B3_3"
basedir <- "data"
min_frames <- 50
input_feature_table <- sprintf("%s/%s_Phase-FullFeatureTable.csv", basedir, trial_name)
feature_table <- copyFeatures(input_feature_table, min_frames, source="Phase")
```

In addition to any pre-calculated features, the `extractFeatures()`
function generates 74 descriptive features for each cell on every frame
using the frame images and pre-generated cell boundaries, based on size,
shape, texture, and the local cell density. The output is a dataframe
comprising the `FrameID`, `CellID`, and `ROI_filename` identifying
columns, the 74 features as columns, and any additional features that
may be present (such as from `copyFeatures()`) in further columns. The
program expects frames to be named according to the scheme
`<experiment name>-<frameid>.tif`, where `<frameid>` is a 4 digit
zero-padded integer corresponding to the `FrameID` column, and located
in the `frame_folder` directory, while ROI files are named according to
the `ROI_filename` column and located in the `roi_folder` directory.

``` r
roi_folder <- sprintf("%s/%s_Phase", basedir, trial_name)
image_folder <- sprintf("%s/%s_imagedata", basedir, trial_name)
new_features <- extractFeatures(feature_table, roi_folder, image_folder, framerate=0.0028)
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
tsvariables <- varsFromTimeSeries(new_features)
#> j=1
#> j=2
#> j=3
#> j=4
#> j=5
#> j=6
#> j=7
#> j=8
#> j=9
#> j=10
#> j=11
#> j=12
#> j=13
#> j=14
#> j=15
#> j=16
#> j=17
#> j=18
#> j=19
#> j=20
#> j=21
#> j=22
#> j=23
#> j=24
#> j=25
#> j=26
#> j=27
#> j=28
#> j=29
#> j=30
#> j=31
#> j=32
#> j=33
#> j=34
#> j=35
#> j=36
#> j=37
#> j=38
#> j=39
#> j=40
#> j=41
#> j=42
#> j=43
#> j=44
#> j=45
#> j=46
#> j=47
#> j=48
#> j=49
#> j=50
#> j=51
#> j=52
#> j=53
#> j=54
#> j=55
#> j=56
#> j=57
#> j=58
#> j=59
#> j=60
#> j=61
#> j=62
#> j=63
#> j=64
#> j=65
#> j=66
#> j=67
#> j=68
#> j=69
#> j=70
#> j=71
#> j=72
#> j=73
#> j=74
#> j=75
#> j=76
#> j=77
#> j=78
#> j=79
#> j=80
#> j=81
#> j=82
#> j=83
#> j=84
#> j=85
#> j=86
#> j=87
#> j=88
#> j=89
#> j=90
#> j=91
#> j=92
#> j=93
#> j=94
#> j=95
#> j=96
```

TODO classification/clustering
