prepareFeatureTable = function(featureTableFile, cellID) {
  featureTableData = read.csv(featureTableFile, check.names = FALSE, skip = 1)
  colnames(featureTableData) = unlist(lapply(names(featureTableData), function(x) gsub(' \\(.+$', '', x)))
  return (as.matrix(subset(featureTableData, featureTableData$"Tracking ID" == cellID)))
}

prepareBoundaryCoordinates = function(roiArchive, cellID) {
  roiFullFileList = unzip(roiArchive, list = TRUE)

  grepString = paste('*-', cellID, '.roi', sep = '')
  roiSubsetFileList = grep(grepString, roiFullFileList$"Name", value = TRUE)
  unzip(roiArchive, files = roiSubsetFileList)
  
  orderedRoiSubsetFileList = roiSubsetFileList[order(as.numeric(sub("([0-9]+)\\-.*\\.roi", "\\1", roiSubsetFileList)))]

  prepareRoiOutput = function(roiFile) {
    roi = RImageJROI::read.ijroi(roiFile)
    
    frameID = sub("([0-9]+)\\-.*", "\\1", roi$"name")
    numberOfCoordinates = roi$"n"
    coordinates = c(t(roi$"coords"))
    
    #return (sapply(c(frameID, cellID, numberOfCoordinates, coordinates), as.integer, USE.NAMES = FALSE))
    return (as.integer(c(frameID, cellID, numberOfCoordinates, coordinates)))
  }

  boundaryCoordinates = lapply(orderedRoiSubsetFileList, prepareRoiOutput)
  lapply(orderedRoiSubsetFileList, unlink)

  return (boundaryCoordinates)
}
# TODO: add number of frames calculation
