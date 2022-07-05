library(CellPhe)

# TODO I don't have access to this file
#name = '18112020_A3_4'
name = '05062019_B3_3'
basedir <- sprintf("DATA/%s", name)
imagedata =  sprintf("%s/%s_imagedata", basedir, name)
input_feature_table =  sprintf("%s/%s_Phase-FullFeatureTable.csv", basedir, name)
roi_files =  sprintf("%s/%s_Phase", basedir, name)
outputdata =  sprintf("%s/%s_outputdata.csv", basedir, name)

# READ IN IMAGE FILES, NORMALISING TO [0, 255]:
frames = CellPhe::readTiffs(imagedata)
normalised_frames = lapply(frames, CellPhe::normaliseImage, lower = 0, upper = 255)

# COPY EXISTING FEATURES FOR EACH CELL TRACKED FOR MORE THAN MIN_FRAMES AND FIND MISSING FRAMES:
# TODO USER MIGHT WANT TO DO THEIR OWN THING HERE FOR THEIR OWN DATASET
min_frames = 50
feature_table = CellPhe::copyFeatures(input_feature_table, min_frames)
numcells_over_thresh = feature_table[[1]]
original_IDs = feature_table[[2]]
missing_frames = feature_table[[3]]
features = feature_table[[4]]

# CALCULATE NEW FEATURES FOR EACH CELL TRACKED FOR MORE THAN MIN_FRAMES AND FIND MISSING FRAMES:
new_features = CellPhe::extractFeatures(roi_files, original_IDs, missing_frames, normalised_frames, min_frames)

# CALCULATE VARIABLES FROM EACH FEATURE'S TIME-SERIES
tsvariables = CellPhe::varsFromTimeSeries(features, new_features, name, original_IDs)

# OUTPUT TIME SERIES VARIABLES TO FILE:
write.csv(tsvariables, outputdata)

