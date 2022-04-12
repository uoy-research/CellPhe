# TODO: add number of frames calculation

prepareFeatureTable = function(featureTableFile, cellId) {
  featureTableData = read.csv(featureTableFile, check.names = FALSE, skip = 1)
  colnames(featureTableData) = unlist(lapply(names(featureTableData), function(x) gsub(' \\(.+$', '', x)))
  return (as.matrix(subset(featureTableData, featureTableData[["Tracking ID"]] == cellId)))
}

readRois = function(roiArchive, cellId) {
  roiFullFileList = unzip(roiArchive, list = TRUE)

  grepString = paste('*-', cellId, '.roi', sep = '')
  roiSubsetFileList = grep(grepString, roiFullFileList[["Name"]], value = TRUE)
  unzip(roiArchive, files = roiSubsetFileList)
  
  orderedRoiSubsetFileList = roiSubsetFileList[order(as.numeric(sub("([0-9]+)\\-.*\\.roi", "\\1", roiSubsetFileList)))]
  rois = lapply(orderedRoiSubsetFileList, RImageJROI::read.ijroi)
  rois = lapply(rois, function(roi) {roi[["cellId"]] = as.numeric(sub(".*-", "", roi[["name"]])); roi[["frameId"]] = as.numeric(sub("-.*", "", roi[["name"]])); return (roi)})
  lapply(orderedRoiSubsetFileList, unlink)

  return (rois)
}

simplifyRoi = function(roi) {
  return (as.integer(c(roi[["frameId"]], roi[["cellId"]], roi[["n"]], c(t(roi[["coords"]])))))
}

prepareBoundaryCoordinates = function(roiArchive, cellId) {
  return (lapply(readRois(roiArchive, cellId), simplifyRoi))
}

normalise = function(values, lower, upper) {
  scale = (upper - lower) / (max(values) - min(values))
  offset = upper - (scale * max(values))

  return (scale * values + offset)
}

readTiffs = function(directory) {
  tiffFileList = list.files(directory, pattern = "*.tif$", full.name = TRUE)

  return (lapply(tiffFileList, function(filename) {
         return (list("frameId" = as.numeric(sub(".*-([0-9]{4})\\.tif", "\\1", filename)), "intensities" = tiff::readTIFF(filename)))
  }))
}

applyRoiMask = function(roi, frame) {
  intensities = matrix(ncol = 3)

  for (i in roi[["left"]] : roi[["right"]]) {
    for (j in roi[["top"]] : roi[["bottom"]]) {
      if (ptinpoly::pip2d(roi[["coords"]], cbind(i, j)) %in% c(0, 1)) {
        intensities = rbind(intensities, cbind(i, j, frame[["intensities"]][j + 1, i + 1]))
      } else if (ptinpoly::pip2d(roi[["coords"]], cbind(i, j)) == -1) {
        intensities = rbind(intensities, cbind(i, j, -1))
      }
    }
  }

  intensities = na.omit(intensities)
  attributes(intensities)[["na.action"]] = NULL
  
  intensities[ , 1] = intensities[ , 1] - roi[["left"]]
  intensities[ , 2] = intensities[ , 2] - roi[["top"]]

  intensities = cbind(intensities, apply(intensities, 1, function(row) return (row[[1]] + (roi[["width"]] * row[[2]]))))
  intensities = intensities[order(intensities[ , 4]), ]
  
  return (as.integer(c(frame[["frameId"]], roi[["width"]] + 1, roi[["height"]] + 1, intensities[ , 3])))
}

prepareMiniImage = function(rois, frames) {
  lapply(rois, function(roi) {CellPhe::applyRoiMask(roi, frames[which(vapply(frames, "[[", 1, "frameId") == roi[["frameId"]])][[1]])})
}
