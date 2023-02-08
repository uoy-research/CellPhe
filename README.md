
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CellPhe

<!-- badges: start -->
<a href="https://zenodo.org/badge/latestdoi/449769672"><img src="https://zenodo.org/badge/449769672.svg" alt="DOI"></a>
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
capabilities, this data is available in `example_data.zip` and comprises
3 parts:

-   The time-lapse stills as TIFF images (`05062019_B3_3_imagedata`)
-   Existing pre-extracted features
    (`05062019_B3_3_Phase-FullFeatureTable.csv`)
-   Region-of-interest (ROI) boundaries already demarked in ImageJ
    format (`05062019_B3_3_Phase`)

These should be extracted into a suitable location before proceeding
with the rest of the tutorial.

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
min_frames <- 50
input_feature_table <- "05062019_B3_3_Phase-FullFeatureTable.csv"
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
roi_folder <- "05062019_B3_3_Phase"
image_folder <- "05062019_B3_3_imagedata"
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
```
