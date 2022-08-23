devtools::load_all()
library(CellPhe)
library(tidyverse)
trial_name <- "05062019_B3_3"
basedir <- "data"

min_frames <- 200
input_feature_table <- sprintf("%s/%s_Phase-FullFeatureTable.csv", basedir, trial_name)
feature_table <- copyFeatures(input_feature_table, min_frames, source="Phase")

roi_folder <- sprintf("%s/%s_Phase", basedir, trial_name)
frame_folder <- sprintf("%s/%s_imagedata", basedir, trial_name)
new_features <- extractFeatures(feature_table, roi_folder, frame_folder, framerate=0.0028)
tsvariables <- varsFromTimeSeries(new_features)


############## Benchmark extractFeatures
old_feat <- readRDS("tests/expected_output2/extractFeatures_output.rds")
old_ids <- unlist(readRDS("tests/expected_output2/originalIDs.rds"))
# Need to convert to dataframe for comparison
old_feat_df <- bind_rows(setNames(lapply(old_feat[[2]], as.data.frame), old_ids), .id='CellID')

# NB: It's possible to have more cells in old dataframe than new
# The old code counts the number of frames as the difference
# between the max and min frame id, rather than how many rows there
# are per cell.
# This shows how many cells would be found with the old code,
# which matches up with the expected output.
# NB: the > rather than >= which I think is a bug
read_csv("data/05062019_B3_3_Phase-FullFeatureTable.csv", skip=1) |>
    group_by(`Tracking ID`) |> 
    summarise(frame_diff = max(Frame) - min(Frame) + 1) |>
    filter(frame_diff > 200) |>
    nrow()

# This test will only compare the cells present in both datasets
old_feat_df <- old_feat_df |> filter(CellID %in% unique(feature_table$CellID))

# Both functions return slightly different numbers of columns
# The old function doesn't return:
#   - FrameID
#   - the Phase features of Volume and Sphericity
#   - x-y centres
# Density also has a different name
old_feat_df <- old_feat_df |> rename(dens=Den)

# So now can investigate the number of frames issue
# Which I think is to do with this missing frames business
missing_frames <- readRDS("tests/expected_output2/missing_frames.rds")
# Subset this to the CellIDs that have the required number of frames available
missing_frames <- missing_frames[old_ids %in% unique(feature_table$CellID)]
old_ids <- old_ids[old_ids %in% unique(feature_table$CellID)]

# NB: The old results have fully NA rows for frames with missing cells
# This code removes those
missing_df <- map_dfr(setNames(missing_frames, old_ids), 
                      function(x) tibble(framenum=1:length(x), frame_missing=x), 
                      .id="CellID") |> 
              filter(frame_missing == 0)
  
old_feat_df_nomissing <- old_feat_df |>
  group_by(CellID) |>
  mutate(framenum = row_number()) |>
  inner_join(missing_df, by=c("CellID", "framenum")) |>
  select(-frame_missing) |>
  mutate(CellID = as.numeric(CellID)) |>
  ungroup()

# Reorder to put frame number second for ease of inspection
cols <- colnames(old_feat_df_nomissing)
cols <- c("CellID", "framenum", setdiff(cols, c("CellID", "framenum")))
old_feat_df_nomissing <- old_feat_df_nomissing[, cols]

# New dataset has those 4 more features: ROI_filename, volume, sphericity, x, y
assertthat::are_equal(nrow(old_feat_df_nomissing), nrow(new_features))
assertthat::are_equal(ncol(old_feat_df_nomissing), ncol(new_features)-5)

# Now can compare on the feature values after ordering them the same
# Assuming the order of frames in the old dataframe is the same 
# as for the new one (that has explicit frameIDs)
old_feat_df_nomissing <- old_feat_df_nomissing |> arrange(CellID, framenum)
new_features <- new_features |> arrange(CellID, FrameID)
cols_to_compare <- setdiff(colnames(old_feat_df_nomissing), c("CellID", "framenum"))

# This shouldn't raise any exception
# Using are_equal to allow for a tolerance in floats
comparison <- sapply(cols_to_compare, function(col) {
  assertthat::are_equal(new_features[[col]], old_feat_df_nomissing[[col]])
})

# The only failure is density, but this is because the algorithm changed
table(comparison)
comparison[!comparison]

############## Benchmark varsFromTimeSeries
tsvariables <- tsvariables |> dplyr::arrange(CellID)
old_ts <- readRDS("tests/expected_output2/tsvariables_output.rds")
# Rename columns to be consistent with old output
colnames(tsvariables) <- gsub("CellID", "ID", colnames(tsvariables))
colnames(tsvariables) <- gsub("Volume_", "Vol_", colnames(tsvariables))
colnames(tsvariables) <- gsub("Sphericity_", "Sph_", colnames(tsvariables))
colnames(tsvariables) <- gsub("dens_", "Den_", colnames(tsvariables))

# Restrict analysis to CellIDs that appear in both
# The old CellIDs are stored separately
ids_to_use <- old_ids %in% unique(new_features$CellID)

# Raises error if not equal
# Again using are_equal to allow for tolerance in floats
comparison <- sapply(colnames(tsvariables), function(col) {
  assertthat::are_equal(tsvariables[[col]], old_ts[, col])
})

# The only failure is density, but again this is because the algorithm changed
table(comparison)
comparison[!comparison]

# Finally compare on trajArea, which is now calculated inside varsFromTimeSeries()
# Test passed
assertthat::are_equal(tsvariables$trajArea, old_ts[, "trajArea"])
