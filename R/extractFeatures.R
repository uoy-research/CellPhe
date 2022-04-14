#' Prepare Feature Table
#'
#' Reads a full feature table into a matrix, suitable for hand-off to feature extraction.
#'
#' @param featureTableFile File name of a full feature table CSV, as a character vector.
#' @param cellId The tracking ID of the cell you want to evaluate, as an integer.
#' 
#' @return A matrix containing full feature table data for the passed cell ID.
#'
#' @export
prepareFeatureTable = function(featureTableFile, cellId) {
  featureTableData = read.csv(featureTableFile, check.names = FALSE, skip = 1)
  colnames(featureTableData) = unlist(lapply(names(featureTableData), function(x) gsub(' \\(.+$', '', x)))
  return (as.matrix(subset(featureTableData, featureTableData[["Tracking ID"]] == cellId)))
}

#' Read ROI Files
#' 
#' Reads a collection of ImageJ ROI files, selecting those relevant to a specific cell ID.
#'
#' @param directory File name of a directory containing ROI files, as a character vector.
#' @param cellId The tracking ID of the cell you want to evaluate, as an integer.
#'
#' @return A list of ROI objects, including additional parameters `cellId` and `frameId`
#'
#' @export
readRois = function(directory, cellId) {
  roiFileList = list.files(directory, pattern = paste('.*-', cellId, '.roi', sep = ''), full.name = TRUE)
  orderedRoiFileList = roiFileList[order(as.numeric(sub("([0-9]+)\\-.*\\.roi", "\\1", basename(roiFileList))))]
  rois = lapply(orderedRoiFileList, RImageJROI::read.ijroi)
  rois = lapply(rois, function(roi) {roi[["cellId"]] = as.numeric(sub(".*-", "", roi[["name"]])); roi[["frameId"]] = as.numeric(sub("-.*", "", roi[["name"]])); return (roi)})
  return (rois)
}

#' Simplify an ROI Object
#'
#' Simplifies an ROI object to its frame ID, cell ID, and region of interest coordinates.
#'
#' @param roi ImageJ ROI data, as an ijroi object.
#'
#' @return ROI frame ID, cell ID, and coordinates, as an integer vector.
simplifyRoi = function(roi) {
  return (as.integer(c(roi[["frameId"]], roi[["cellId"]], roi[["n"]], c(t(roi[["coords"]])))))
}

#' Prepare ROI Boundary Coordinates
#'
#' Simplifies a collection of ImageJ ROI data into frame ID, cell ID, and coordinates.
#'
#' @param rois A collection of ImageJ ROI data, as a list.
#'
#' @return A collection of simplified ROI data, as a list of integer vectors.
#'
#' @export
prepareBoundaryCoordinates = function(rois) {
  return (lapply(rois, simplifyRoi))
}

#' Normalise an Image
#'
#' Normalises an image to a specified range.
#'
#' @param values The image to normalise, as a matrix.
#' @param lower The lower bound of the target normalisation range, as an integer.
#' @param upper The upper boundo of the target normalisation range, as an integer.
#'
#' @return The normalised image, as a matrix.
#'
#' @export
normaliseImage = function(values, lower, upper) {
  scale = (upper - lower) / (max(values) - min(values))
  offset = upper - (scale * max(values))

  return (scale * values + offset)
}

#' Read a Collection of TIFF images
#'
#' Reads a directory of TIFF images, sorting by sequence number in the file name.
#'
#' @param directory File name of a directory containing TIFF files, as a character vector.
#'
#' @return A collection of TIFF images, as a list of matrices.
#'
#' @export
readTiffs = function(directory) {
  tiffFileList = list.files(directory, pattern = "*.tif$", full.name = TRUE)
  orderedTiffFileList = tiffFileList[order(as.numeric(sub(".*-([0-9]{4})\\.tif", "\\1", basename(tiffFileList))))]

  return (lapply(tiffFileList, tiff::readTIFF))
}

#' Apply an ImageJ ROI Mask to a TIFF Image
#'
#' Masks a section of a TIFF image corresponding to an ImageJ ROI, including those pixels on the boundary of the ROI.
#'
#' @param roi An ImageJ ROI, as an ijroi object.
#' @param frame A TIFF image, as a matrix.
#'
#' @return The frame ID, width, height, and pixel values of the masked region, as an integer vector.
applyRoiMask = function(roi, frame) {
  intensities = matrix(ncol = 3)

  for (i in roi[["left"]] : roi[["right"]]) {
    for (j in roi[["top"]] : roi[["bottom"]]) {
      if (ptinpoly::pip2d(roi[["coords"]], cbind(i, j)) %in% c(0, 1)) {
        intensities = rbind(intensities, cbind(i, j, frame[j + 1, i + 1]))
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
  
  return (as.integer(c(roi[["frameId"]], roi[["width"]] + 1, roi[["height"]] + 1, intensities[ , 3])))
}

#' Prepare Mini Image
#'
#' Masks a collection of TIFF images to a region of interest for a specific cell.
#'
#' @param rois A collection of ImageJ ROIs pertaining to a specific cell, as a list of ijroi objects.
#' @param frames A collection of TIFF images, as a list of matrices.
#'
#' @return A collection of frame IDs, widths, heights, and pixel values for each TIFF image, as a list of integer vectors.
#'
#' @export
prepareMiniImage = function(rois, frames) {
  lapply(rois, function(roi) {applyRoiMask(roi, frames[[roi[["frameId"]]]])})
}
