prepareFeatureTable = function(featureTableFile, cellID) {
  featureTableData = read.csv(featureTableFile, header = FALSE, skip = 2)
  return as.matrix(subset(featureTableData, featureTableData[ , 3] == cellID))
}

# TODO: add number of frames calculation
