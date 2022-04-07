# TODO: add number of frames calculation

prepareFeatureTable = function(featureTableFile, cellID) {
  featureTableData = read.csv(featureTableFile, check.names = FALSE, skip = 1)
  colnames(featureTableData) = unlist(lapply(names(featureTableData), function(x) gsub(' \\(.+$', '', x)))
  return (as.matrix(subset(featureTableData, featureTableData$"Tracking ID" == cellID)))
}

getRois = function(roiArchive, cellID) {
  roiFullFileList = unzip(roiArchive, list = TRUE)

  grepString = paste('*-', cellID, '.roi', sep = '')
  roiSubsetFileList = grep(grepString, roiFullFileList$"Name", value = TRUE)
  unzip(roiArchive, files = roiSubsetFileList)
  
  orderedRoiSubsetFileList = roiSubsetFileList[order(as.numeric(sub("([0-9]+)\\-.*\\.roi", "\\1", roiSubsetFileList)))]
  rois = lapply(orderedRoiSubsetFileList, RImageJROI::read.ijroi)
  lapply(orderedRoiSubsetFileList, unlink)
  return (rois)
}

transformRoi = function(roiData) {
  frameID = sub("-.*", "", roiData[["name"]])
  cellID = sub(".*-", "", roiData[["name"]])
  numberOfCoordinates = roiData[["n"]]
  coordinates = c(t(roiData[["coords"]]))
  
  return (as.integer(c(frameID, cellID, numberOfCoordinates, coordinates)))
}

prepareBoundaryCoordinates = function(roiArchive, cellID) {
  return (lapply(getRois(roiArchive, cellID), transformRoi))
}

normalise = function(values, lower, upper) {
  scale = (upper - lower) / (max(values) - min(values))
  offset = upper - (scale * max(values))

  return (scale * values + offset)
}

readTiffs = function(directory) {
  tiffFileList = list.files(directory, pattern = "*.tif$", full.name = TRUE)
  return (lapply(tiffFileList, tiff::readTIFF))
}

applyRoiMask = function(roi, frameNumber, frames) {
  # FOR min(X) TO max(X) IN roi
    # FOR min(Y) TO max(Y) IN roi
      # IF [x, y] IS in OR on surface OF polygon described by ROI
        # ADD [x, y] value TO list of intensities to return
      # IF [x, y] IS NOT in polygon described by ROI
        # ADD -1 TO list of intensities to return
  # DROP NaN VALUES from intensities to return
}
