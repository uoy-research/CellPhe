devtools::load_all()
library(CellPhe)
trial_name <- "05062019_B3_3"
basedir <- "data"

imagedata <- sprintf("%s/%s_imagedata", basedir, trial_name)
frames <- readTiffs(imagedata)
min_frames <- 200
input_feature_table <- sprintf("%s/%s_Phase-FullFeatureTable.csv", basedir, trial_name)
feature_table <- copyPhaseFeatures(input_feature_table, min_frames)

roi_files <- sprintf("%s/%s_Phase", basedir, trial_name)
new_features <- extractFeatures(feature_table, frames, roi_files, min_frames, framerate=0.0028)

tsvariables = varsFromTimeSeries(new_features)

############## Benchmark extractFeatures
# TODO Need to rerun the original code (initial-release branch) and 
# also save the OriginalIDs list
old_feat <- readRDS("tests/expected-output/extractFeatures_output.rds")
# Need to convert to dataframe for comparison
old_feat_df <- do.call('rbind', old_feat[[2]])

############## Benchmark varsFromTimeSeries
tsvariables <- tsvariables |> dplyr::arrange(CellID)
old <- readRDS("tests/expected-output/tsvariables_output.rds")
colnames(tsvariables) <- gsub("CellID", "ID", colnames(tsvariables))
colnames(tsvariables) <- gsub("Volume_", "Vol_", colnames(tsvariables))
colnames(tsvariables) <- gsub("Sphericity_", "Sph_", colnames(tsvariables))
colnames(tsvariables) <- gsub("dens_", "Den_", colnames(tsvariables))
comp <- sapply(colnames(tsvariables), function(x) all(tsvariables[[x]] == old_ts[, x]))
# Nope this has changed a fair bit!
table(comp)
