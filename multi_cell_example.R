library('Rcpp')
library('CellPhe')

# READ IN IMAGE FILES, NORMALISING TO [0, 255]:
frames = CellPhe::readTiffs('25112020_D4_5_imagedata')
normalised_frames = lapply(frames, CellPhe::normaliseImage, lower = 0, upper = 255)

# COPY EXISTING FEATURES FOR EACH CELL TRACKED FOR MORE THAN MIN_FRAMES AND FIND MISSING FRAMES:
min_frames = 50
feature_table = copyFeatures('25112020_D4_5_Phase-FullFeatureTable.csv', min_frames)
numcells_over_thresh = feature_table[[1]]
original_IDs = feature_table[[2]]
missing_frames = feature_table[[3]]
features = feature_table[[4]]

# CALCULATE NEW FEATURES FOR EACH CELL TRACKED FOR MORE THAN MIN_FRAMES AND FIND MISSING FRAMES:
new_features = extractFeatures('25112020_D4_5', original_IDs, missing_frames, normalised_frames, min_frames)

# CALCULATE VARIABLES FROM EACH FEATURE'S TIME-SERIES
tsvariables = varsFromTimeSeries(features, new_features, '25112020_D4_5', original_IDs)

# OUTPUT TIME SERIES VARIABLES TO FILE:
write.csv(tsvariables, "25112020_D4_5_outputdata.csv")

