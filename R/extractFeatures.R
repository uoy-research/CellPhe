prepareFeatureTable = function(featureTableFile, cellID) {
  featureTableData = read.csv(featureTableFile, check.names=FALSE, skip = 1)
  colnames(featureTableData) = unlist(lapply(names(featureTableData), function(x) gsub(' \\(.+$', '', x)))
  return (as.matrix(subset(featureTableData, featureTableData$"Tracking ID" == cellID)))
}

# TODO: add number of frames calculation
