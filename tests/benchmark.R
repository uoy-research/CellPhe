# Expects data unzipped and placed in data/
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
results <- c()
old_feat <- readRDS("tests/expected_output/extractFeatures_output.rds")
old_ids <- unlist(readRDS("tests/expected_output/originalIDs.rds"))
# Need to convert from list to dataframe for comparison
old_feat_df <- bind_rows(setNames(lapply(old_feat[[2]], as.data.frame), old_ids), .id='CellID')

# NB: It's possible to have more cells in old dataframe than new
# The old code counts the number of frames per cell as the difference
# between the max and min frame id, rather than how many rows there
# are per cell.
# This test will only compare the cells present in both datasets
old_feat_df <- old_feat_df |> filter(CellID %in% unique(feature_table$CellID))

# Density used to be called Den, it's now dens
old_feat_df <- old_feat_df |> rename(dens=Den)

# There's also a discrepancy between the number of frames in both outputs
# This is caused by the way in which 'missing_frames' are handled.
# This was previously saved as a binary vector
missing_frames <- readRDS("tests/expected_output/missing_frames.rds")
# Subset this to the CellIDs that have the required number of frames available
# using the new logic
missing_frames <- missing_frames[old_ids %in% unique(feature_table$CellID)]
old_ids <- old_ids[old_ids %in% unique(feature_table$CellID)]

# NB: The old features have fully NA rows for frames with missing cells
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

# Now can compare on the feature values after ordering them the same
# Assuming the order of frames in the old dataframe is the same 
# as for the new one (that has explicit frameIDs)
old_feat_df_nomissing <- old_feat_df_nomissing |> arrange(CellID, framenum)
new_features <- new_features |> arrange(CellID, FrameID)
cols_to_compare <- setdiff(colnames(old_feat_df_nomissing), c("CellID", "framenum"))

cat("\nextractFeatures\n~~~~~~~~~~~~~~~\n")
cat(sprintf("Same number of rows (%d vs %d)\n", nrow(old_feat_df_nomissing), nrow(new_features)))
res <- assertthat::are_equal(nrow(old_feat_df_nomissing), nrow(new_features))
results <- c(results, res)
res

cat(sprintf("Same number of columns after accounting for ROI_filename, volume, sphericity, xpos, and ypos (%d vs %d)\n", ncol(old_feat_df_nomissing), ncol(new_features)))
res <- assertthat::are_equal(ncol(old_feat_df_nomissing), ncol(new_features)-5)
results <- c(results, res)
res

# Compare the feature values themselves
# Using are_equal to allow for a tolerance in floats
comparison <- sapply(cols_to_compare, function(col) {
  assertthat::are_equal(new_features[[col]], old_feat_df_nomissing[[col]])
})

# The only failure is density and Cooc02 features, as we changed Density calculation
# and CoocO2 was incorrectly returning wrong value before
cat(sprintf("Columns have the same feature values after accounting for fixes to Cooc02 (14 features), Density, minBox (2), minBox:Area, rectangularity (%d/%d)\n", sum(comparison), length(comparison)))
res <- assertthat::are_equal(sum(comparison), 53)
results <- c(results, res)
res

############## Benchmark varsFromTimeSeries
tsvariables <- tsvariables |> dplyr::arrange(CellID)
old_ts <- readRDS("tests/expected_output/tsvariables_output.rds")
# Rename columns to be consistent with old output
colnames(tsvariables) <- gsub("CellID", "ID", colnames(tsvariables))
colnames(tsvariables) <- gsub("Volume_", "Vol_", colnames(tsvariables))
colnames(tsvariables) <- gsub("Sphericity_", "Sph_", colnames(tsvariables))
colnames(tsvariables) <- gsub("dens_", "Den_", colnames(tsvariables))

# Restrict analysis to CellIDs that appear in both
# The old CellIDs are stored separately
ids_to_use <- old_ids %in% unique(new_features$CellID)

# Check feature values
comparison <- sapply(colnames(tsvariables), function(col) {
  assertthat::are_equal(tsvariables[[col]], old_ts[, col])
})

cat("\nvarsFromTimeSeries\n~~~~~~~~~~~~~~~~~~\n")
cat(sprintf("Columns have the same feature values after accounting for the 19 features (x15 time-series variables) that have had bug-fixes (%d/%d)\n", sum(comparison), length(comparison)))
res <- assertthat::are_equal(sum(comparison), 828)
results <- c(results, res)
res

cat("\nResults\n~~~~~~~\n")
cat(sprintf("%d/%d tests passed\n", sum(results), length(results)))
status <- as.integer(!all(results))
quit(save="no", status=status)
