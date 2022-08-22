#' Prepare Feature Table
#'
#' Reads a full feature table into a matrix, suitable for hand-off to feature extraction.
#'
#' @param featureTableFile File name of a full feature table CSV, as a character vector.
#' @param cellId The tracking ID of the cell you want to evaluate, as an integer.
#'
#' @return A matrix containing full feature table data for the passed cell ID.
prepareFeatureTable = function(featureTableFile, cellId) {
  featureTableData = utils::read.csv(featureTableFile, check.names = FALSE, skip = 1)
  colnames(featureTableData) = sapply(names(featureTableData), function (name) {
    gsub(' \\(.+$', '', name)
  }, USE.NAMES = FALSE)
  return (as.matrix(subset(
    featureTableData, featureTableData[["Tracking ID"]] == cellId
  )))
}

#' Read ROI Files
#'
#' Reads a collection of ImageJ ROI files, selecting those relevant to a specific cell ID.
#'
#' @param directory File name of a directory containing ROI files, as a character vector.
#' @param cellId The tracking ID of the cell you want to evaluate, as an integer.
#'
#' @return A list of ROI objects, including additional parameters `cellId` and `frameId`
readRois = function(directory, cellId) {
  roiFileList = list.files(
    directory,
    pattern = paste('.*-', cellId, '.roi', sep = ''),
    full.names = TRUE
  )
  orderedRoiFileList = roiFileList[order(as.numeric(sub(
    "([0-9]+)\\-.*\\.roi", "\\1", basename(roiFileList)
  )))]
  rois = lapply(orderedRoiFileList, RImageJROI::read.ijroi)
  rois = lapply(rois, function(roi) {
    roi[["cellId"]] = as.numeric(sub(".*-", "", roi[["name"]]))
    roi[["frameId"]] = as.numeric(sub("-.*", "", roi[["name"]]))
    return (roi)
  })
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
  return (as.integer(c(roi[["frameId"]], roi[["cellId"]], roi[["n"]], c(t(
    roi[["coords"]]
  )))))
}

#' Prepare ROI Boundary Coordinates
#'
#' Simplifies a collection of ImageJ ROI data into frame ID, cell ID, and coordinates.
#'
#' @param rois A collection of ImageJ ROI data, as a list.
#'
#' @return A collection of simplified ROI data, as a list of integer vectors.
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
  tiffFileList = list.files(directory, pattern = "*.tif$", full.names = TRUE)
  orderedTiffFileList = tiffFileList[order(as.numeric(sub(
    ".*-([0-9]{4})\\.tif", "\\1", basename(tiffFileList)
  )))]
  
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
  
  for (i in roi[["left"]]:roi[["right"]]) {
    for (j in roi[["top"]]:roi[["bottom"]]) {
      if (ptinpoly::pip2d(roi[["coords"]], cbind(i, j)) %in% c(0, 1)) {
        intensities = rbind(intensities, cbind(i, j, frame[j + 1, i + 1]))
      } else if (ptinpoly::pip2d(roi[["coords"]], cbind(i, j)) == -1) {
        intensities = rbind(intensities, cbind(i, j,-1))
      }
    }
  }
  
  intensities = stats::na.omit(intensities)
  attributes(intensities)[["na.action"]] = NULL
  
  intensities[, 1] = intensities[, 1] - roi[["left"]]
  intensities[, 2] = intensities[, 2] - roi[["top"]]
  
  intensities = cbind(intensities, apply(intensities, 1, function(row)
    return (row[[1]] + (roi[["width"]] * row[[2]]))))
  intensities = intensities[order(intensities[, 4]),]
  
  return (as.integer(c(roi[["frameId"]], roi[["width"]] + 1, roi[["height"]] + 1, intensities[, 3])))
}

#' Prepare Mini Image
#'
#' Masks a collection of TIFF images to a region of interest for a specific cell.
#'
#' @param rois A collection of ImageJ ROIs pertaining to a specific cell, as a list of ijroi objects.
#' @param frames A collection of TIFF images, as a list of matrices.
#'
#' @return A collection of frame IDs, widths, heights, and pixel values for each TIFF image, as a list of integer vectors.
prepareMiniImage = function(rois, frames) {
  lapply(rois, function(roi) {
    applyRoiMask(roi, frames[[roi[["frameId"]]]])
  })
}

#' Calculates cell features from timelapse videos
#'
#' Calculates 1109 features related to size, shape, texture and movement for each
#' cell on every non-missing frame, as well as the cell density around each
#' cell on each frame and a measure describing the trajectory of the cell
#' over all frames.
#'
#' @details
#' Reads in information on cell boundaries from ROI files, the cell identifiers
#' from the original tracking (i.e. before cells tracked for less than
#' min_frames were removed), original_IDs, information on missing frames for
#' each cell in a list of vectors with missing frames indicated by 1 and
#' non-missing frames by 0, missing_frames, and the normalised images for every frame, normalised_frames.
#' @param df DataFrame where every row corresponds to a combination of a cell
#' tracked in a frame. It must have at least columns \code{CellID} and \code{FrameID},
#' along with any additional features.
#' @param roi_folder A path to a directory containing multiple Report Object Instance
#' (ROI) files named in the format \code{cellid}-\code{frameid}.roi
#' @param frame_folder A path to a directory containing multiple frames in TIFF format.
#' It is assumed these are named under the pattern \code{<experiment name>-<frameid>.tif}, where 
#' \code{<frameid>} is a 4 digit zero-padded integer.
#' @return A dataframe with 76 columns and 1 row per cell per frame it's present in:
#' \itemize{
#'   \item{\code{frameID}: the numeric frameID}
#'   \item{\code{cellID}: the numeric cellID}
#'   \item{\code{...}: 72 frame specific features
#'   \item{\code{xcentres}: the x-coordinate of the cell in that frame}
#'   \item{\code{ycentres}: the y-coordinate of the cell in that frame}
#'   \item{\code{...}: Any other data columns that were present in \code{df}}
#' }
#' @export
extractFeatures = function(df,
                           roi_folder,
                           frame_folder,
                           framerate = 1) {
  n_cells <- length(unique(df$CellID))
  all_features <- vector(mode = "list", length = n_cells)
  centroids <- vector(mode = "list", length = n_cells)
  RandA <- vector(mode = "list", length = n_cells)
  meanr = rep(NA, n_cells)
  
  mfeature_cols <- c("Dis", "Trac", "D2T", "Vel")
  bfeature_cols <- c(
    "Rad",
    "VfC",
    "Curv",
    "Len",
    "Wid",
    "Area",
    "A2B",
    "Box",
    "Rect",
    "poly1",
    "poly2",
    "poly3",
    "poly4"
  )
  
  tfeature_cols <- c(
    "FOmean",
    "FOsd",
    "FOskew",
    "Cooc01ASM",
    "Cooc01Con",
    "Cooc01IDM",
    "Cooc01Ent",
    "Cooc01Cor",
    "Cooc01Var",
    "Cooc01Sav",
    "Cooc01Sen",
    "Cooc01Den",
    "Cooc01Dva",
    "Cooc01Sva",
    "Cooc01f13",
    "Cooc01Sha",
    "Cooc01Pro",
    "Cooc12ASM",
    "Cooc12Con",
    "Cooc12IDM",
    "Cooc12Ent",
    "Cooc12Cor",
    "Cooc12Var",
    "Cooc12Sav",
    "Cooc12Sen",
    "Cooc12Den",
    "Cooc12Dva",
    "Cooc12Sva",
    "Cooc12f13",
    "Cooc12Sha",
    "Cooc12Pro",
    "Cooc02ASM",
    "Cooc02Con",
    "Cooc02IDM",
    "Cooc02Ent",
    "Cooc02Cor",
    "Cooc02Var",
    "Cooc02Sav",
    "Cooc02Sen",
    "Cooc02Den",
    "Cooc02Dva",
    "Cooc02Sva",
    "Cooc02f13",
    "Cooc02Sha",
    "Cooc02Pro",
    "IQ1",
    "IQ2",
    "IQ3",
    "IQ4",
    "IQ5",
    "IQ6",
    "IQ7",
    "IQ8",
    "IQ9"
  )
  
  # SHAPE FEATURES FROM BOUNDARIES:
  bfeatures = matrix(NA, nrow = nrow(df), ncol = length(bfeature_cols))
  colnames(bfeatures) = bfeature_cols
  
  # TEXTURE FEATURES FROM INTERIOR PIXELS:
  tfeatures = matrix(NA, nrow = nrow(df), ncol = length(tfeature_cols))
  colnames(tfeatures) <- tfeature_cols
  
  xcentres <- rep(NA, nrow(df))
  ycentres <- rep(NA, nrow(df))
  ids <- matrix(NA, nrow=nrow(df), ncol=2)
  colnames(ids) <- c("CellID", "FrameID")
  
  all_frame_ids <- unique(df$FrameID)
  row_num <- 1
  for (frame_id in all_frame_ids) {
    # Find all Cells that feature in this FrameID
    cell_ids <-
      df |> dplyr::filter(FrameID == frame_id) |> dplyr::distinct(CellID) |> dplyr::pull(CellID)
    n_cells <- length(cell_ids)
    
    # Load frame into memory
    tiff_fn <- list.files(frame_folder, pattern = sprintf(".*-%04d.tif$", frame_id), full.names = TRUE)
    if (length(tiff_fn) != 1) {
      stop("Cannot find tif for frame id %d in %s. Check that the filename convention is as expected (?extractFeatures).", frame_id, frame_folder)
    }
    frame <- normaliseImage(tiff::readTIFF(tiff_fn), lower=0, upper=255)
  
    for (cell_id in cell_ids) {
      ids[row_num, ] <- c(cell_id, frame_id)
      roi_fn <-
        df |> dplyr::filter(CellID == cell_id, FrameID == frame_id) |> dplyr::distinct(ROI_filename) |> dplyr::pull(ROI_filename)
      if (length(roi_fn) > 1) {
        stop("Error: found more than one ROI filename for Cell %d and Frame %d",
             cell_id,
             frame_id)
      }
      roi = RImageJROI::read.ijroi(sprintf("%s/%s.roi", roi_folder, roi_fn))
      # It is possible to have negative coordinates
      roi$coords[which(roi$coords < 0)] = 0
      # check whether cell too small in either direction
      w = max(roi$coords[, 1]) - min(roi$coords[, 1]) + 1
      h = max(roi$coords[, 2]) - min(roi$coords[, 2]) + 1
      if ((w < 8) | (h < 8)) {
        next
      }
    
      # EXTRACT SUB-IMAGE, MASK, CELL PIXELS AND BOUNDARY COORDINATES FOR SPECIFIC CELL:
      sub_image_info = subImageInfo(roi, frame)
      cell_pixels = sub_image_info[[5]]
      interior_pixels = sub_image_info[[6]]
      boundary_coordinates = sub_image_info[[7]]
      xcentres[row_num] = sub_image_info[[8]][1]
      ycentres[row_num] = sub_image_info[[8]][2]
      
      # SHAPE FEATURES
      vfc = varFromCentre(boundary_coordinates)
      # AVERAGE RADIUS
      bfeatures[row_num, 1] = vfc[1]
      # VARIANCE ON BOUNDARY PIXEL DISTANCES FROM CELL CENTRE:
      bfeatures[row_num, 2] = vfc[2]
      # MEASURE OF BOUNDARY CURVATURE:
      bfeatures[row_num, 3] = curvature(boundary_coordinates, 4)
      # WIDTH AND HEIGHT OF MINIMAL BOX:
      box = minBox(boundary_coordinates)
      bfeatures[row_num, 4] = box[1]
      bfeatures[row_num, 5] = box[2]
      # AREA, CALCULATED FROM THE MINI-IMAGE:
      bfeatures[row_num, 6] = nrow(cell_pixels)
      # AREA TO BOUNDARY RATIO:
      bl = boundary_coordinates$length
      bfeatures[row_num, 7] = bfeatures[row_num, 6] / (bl * bl)
      # MINIMAL BOX TO AREA RATIO:
      bfeatures[row_num, 8] = (box[1] * box[2]) / bfeatures[row_num, 3]
      # RECTANGULARITY:
      m = max(box[1], box[2])
      bfeatures[row_num, 9] = m / (box[1] + box[2])
      # FITTED POLYGON FEATURES (MAX_SIDE, MIN_ANGLE, ANGLE_VARIANCE, SIDE_LENGTH_VARIANCE)
      for (k in 1:4) {
        bfeatures[row_num, (k + 9)] = polyClass(boundary_coordinates)[k]
      }
      
      # TEXTURE FEATURES
      # FIRST ORDER FEATURES
      tfeatures[row_num, 1] = mean(cell_pixels[, 3])
      tfeatures[row_num, 2] = sqrt(stats::var(cell_pixels[, 3]))
      tfeatures[row_num, 3] = e1071::skewness(cell_pixels[, 3], type = 2)
      # CALCULATE COOCCURRENCES MATRICES
      cooccurrence_levels = 10
      cooc = cooccur(sub_image_info, cooccurrence_levels)
      # HARALICK FEATURES
      haralickfeatures01 <-
        calculateHaralickFeatures(cooc$cooc01)
      haralickfeatures12 <-
        calculateHaralickFeatures(cooc$cooc12)
      haralickfeatures02 <-
        calculateHaralickFeatures(cooc$cooc12)
      for (k in 1:14) {
        tfeatures[row_num, (k + 3)] = haralickfeatures01[k]
        tfeatures[row_num, (k + 17)] = haralickfeatures12[k]
        tfeatures[row_num, (k + 31)] = haralickfeatures02[k]
      }
      # INTENSITY QUANTILE FEATURES
      quantileVars = intensityQuantiles(boundary_coordinates, interior_pixels)
      for (k in 1:9) {
        tfeatures[row_num, (k + 45)] = quantileVars[[k]]
      }
    
      row_num <- row_num + 1
    }
  }
  # Combine and add movement features
  res <- cbind(ids, bfeatures, tfeatures, xpos=xcentres, ypos=ycentres) |> as.data.frame()
  res <- res |> 
    dplyr::arrange(CellID, FrameID) |>
    dplyr::group_by(CellID) |> 
    dplyr::mutate(startx = xpos[which.min(FrameID)],
           starty = ypos[which.min(FrameID)],
           Dis = sqrt((xpos - startx)**2 + (ypos - starty)**2),
           dist_timestamp = sqrt((xpos - dplyr::lag(xpos))**2 + (ypos - dplyr::lag(ypos))**2),
           dist_timestamp = ifelse(FrameID == min(FrameID), 0, dist_timestamp),
           Trac = cumsum(dist_timestamp),
           D2T = Dis / Trac,
           D2T = ifelse(is.infinite(D2T) | is.nan(D2T), 0, D2T),
           Vel = (framerate * dist_timestamp) / (FrameID - lag(FrameID, default=0))) |>
    dplyr::ungroup() |>
    dplyr::select(-startx, -starty, -dist_timestamp)
    
  meanrad <- res |> 
    dplyr::group_by(CellID) |> 
    summarise(meanr = mean(Rad, na.rm=T)) |>
    ungroup() |>
    summarise(meanrad = mean(meanr, na.rm=T)) |>
    pull(meanrad)
  
  ## CALCULATE DENSITY FOR EACH CELL
  # TODO can I refactor this to just take columns?
  #res <- do.call('rbind', lapply(1:n_cells, function(j) {
  #  feat_df <- as.data.frame(all_features[[j]])
  #  feat_df$dens <- densityCalc(j, centroids, RandA, meanrad)
  #  feat_df <- cbind(feat_df, centroids[[j]])
  #  feat_df
  #}))
  # TODO add the movement features
  # TODO do the density stuff
  ############################################################
  
  df |> dplyr::inner_join(res, by = c("CellID", "FrameID"))
}

subImageInfo = function(roi, frame) {
  sub_image = frame[(min(roi$coords[, 2]) + 1):(max(roi$coords[, 2]) + 1), (min(roi$coords[, 1]) +
                                                                              1):(max(roi$coords[, 1]) + 1)]
  
  width = max(roi$coords[, 1]) - min(roi$coords[, 1]) + 1
  height = max(roi$coords[, 2]) - min(roi$coords[, 2]) + 1
  pixel_type = as.list(rep(1, width * height))
  
  # EXTRACT BOUNDARY PIXELS
  length = roi$n
  x = roi$coords[, 1] - min(roi$coords[, 1]) + 1
  y = roi$coords[, 2] - min(roi$coords[, 2]) + 1
  newbcs = cbind(x, y)
  xycoords = list(length, x, y)
  names(xycoords) <- c("length", "x", "y")
  
  # CELL CENTROID
  cx = mean(roi$coords[, 1])
  cy = mean(roi$coords[, 2])
  
  centroid = cbind(cx, cy)
  
  # ADD BOUNDARY TO PIXEL_TYPE
  pixel_type = as.list(rep(1, width * height))
  inds <- (xycoords[[2]] - 1) * height + xycoords[[3]]
  pixel_type[inds] = 0
  
  # FILL IN MASK FROM EDGES OF IMAGE
  matpix_type <-
    matrix(unlist(pixel_type), nrow = height, ncol = width)
  matpix_type = apply(matpix_type, 2, fill_mask)
  matpix_type = t(apply(matpix_type, 1, fill_mask))
  
  # CHECK NEIGHBOURS OF MASK PIXELS
  n = 1
  while (n > 0) {
    n = 0
    for (j in 2:(height - 1)) {
      for (i in 2:(width - 1)) {
        if (matpix_type[j, i] == -1) {
          if (matpix_type[j - 1, i] == 1) {
            matpix_type[j - 1, i] = -1
            n = 1
          }
          if (matpix_type[j, i - 1] == 1) {
            matpix_type[j, i - 1] = -1
            n = 1
          }
          if (matpix_type[j, i + 1] == 1) {
            matpix_type[j, i + 1] = -1
            n = 1
          }
          if (matpix_type[j + 1, i] == 1) {
            matpix_type[j + 1, i] = -1
            n = 1
          }
        }
      }
    }
  }
  
  n = 1
  intensities = matrix(nrow = width * height, ncol = 4)
  for (i in 1:width) {
    for (j in 1:height) {
      intensities[n, ] = c(i, j, sub_image[j, i], matpix_type[j, i])
      n = n + 1
    }
  }
  
  # IDENTIFY INTERIOR PIXELS (WITHIN, BUT NOT ON THE CELL BOUNDARY)
  cellpixels = intensities[which(intensities[, 4] != -1), 1:3]
  interiorpixels = intensities[which(intensities[, 4] == 1), 1:3]
  
  # MASK NON-CELL PIXELS
  mask = rep(0, width * height)
  mask[which(intensities[, 4] >= 0)] = 1
  mask = matrix(mask, nrow = height)
  
  return (
    list(
      width,
      height,
      sub_image,
      mask,
      cellpixels,
      interiorpixels,
      xycoords,
      centroid
    )
  )
}

varFromCentre = function(bc) {
  dist = c(1:bc$length)
  for (j in 1:bc$length) {
    cx = mean(bc$x)
    cy = mean(bc$y)
    x = bc$x[j]
    y = bc$y[j]
    dist[j] = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy))
  }
  vdist = stats::var(dist)
  rad = mean(dist)
  return(c(rad, stats::var(dist)))
}

curvature = function(bc, gap) {
  curv = 0.0
  for (j in 1:bc$length) {
    x = bc$x[j]
    y = bc$y[j]
    xmgap = bc$x[(j - gap + bc$length) %% bc$length]
    
    ymgap = bc$y[(j - gap + bc$length) %% bc$length]
    
    if (j == gap) {
      xmgap = bc$x[bc$length]
      ymgap = bc$y[bc$length]
    }
    xpgap = bc$x[(j + gap + bc$length) %% bc$length]
    
    ypgap = bc$y[(j + gap + bc$length) %% bc$length]
    
    if (j == (bc$length - gap)) {
      xpgap = bc$x[bc$length]
      ypgap = bc$y[bc$length]
    }
    i1 = (x - xmgap) * (x - xmgap) + (y - ymgap) * (y - ymgap)
    i2 = (x - xpgap) * (x - xpgap) + (y - ypgap) * (y - ypgap)
    i3 = (xmgap - xpgap) * (xmgap - xpgap) + (ymgap - ypgap) * (ymgap - ypgap)
    newval = sqrt(i1) + sqrt(i2) - sqrt(i3)
    curv = curv + newval
  }
  return(curv / bc$length)
}

minBox = function(bc) {
  # FIND MAXIMUM DISTANCE BETWEEN ANY TWO BOUNDARY POINTS
  d = cbind(bc$x, bc$y)
  dd = as.matrix(stats::dist(d))
  k1 = which(dd == max(dd), arr.ind = T)[1]
  k2 = which(dd == max(dd), arr.ind = T)[2]
  keepx1 = bc$x[k1]
  
  keepy1 = bc$y[k1]
  
  keepx2 = bc$x[k2]
  
  keepy2 = bc$y[k2]
  
  alpha = atan2((keepy1 - keepy2), (keepx1 - keepx2))
  # rotating points by -alpha makes keepx1-keepx2 lie along x-axis
  roty = keepy1 - sin(alpha) * (bc$x - keepx1) + cos(alpha) * (bc$y - keepy1)
  rotx = keepx1 + cos(alpha) * (bc$x - keepx1) + sin(alpha) * (bc$y - keepy1)
  minx = min(rotx)
  maxx = max(rotx)
  miny = min(roty)
  maxy = max(roty)
  xlength = maxx - minx
  ylength = maxy - miny
  return(c(xlength, ylength))
}

polyClass = function(bc) {
  points = polygon(bc)
  points1 = c(points[2:length(points)], points[1])
  points2 = c(points1[2:length(points1)], points1[1])
  mat01 = cbind(bc$x[points], bc$y[points], bc$x[points1], bc$y[points1])
  mat02 = cbind(bc$x[points], bc$y[points], bc$x[points2], bc$y[points2])
  mat12 = cbind(bc$x[points1], bc$y[points1], bc$x[points2], bc$y[points2])
  Asq = apply(mat01, 1, sqreucdist)
  Bsq = apply(mat12, 1, sqreucdist)
  Csq = apply(mat02, 1, sqreucdist)
  squaredLengths = rbind(Asq, Bsq, Csq)
  maxLength = max(sqrt(Asq))
  varLength = stats::var(sqrt(Asq))
  angles = apply(squaredLengths, 2, polyAngle)
  minAngle = min(angles)
  varAngle = stats::var(angles)
  output = c(maxLength, minAngle, varAngle, varLength)
  return(output)
}

polygon = function(bc) {
  thresh = 2.5
  # FIND MAXIMUM DISTANCE FROM FIRST BOUNDARY POINT
  d = cbind(bc$x, bc$y)
  dd = as.matrix(stats::dist(d))
  indkeep = which(dd[, 1] == max(dd[, 1]))
  numpoints = 2
  pointArray = as.vector(c(1, indkeep))
  alldone = 0
  while (alldone == 0) {
    tempArray = 1
    n = 0
    alldone = 1
    for (k in 2:numpoints) {
      x1 = bc$x[pointArray[k - 1]]
      y1 = bc$y[pointArray[k - 1]]
      x2 = bc$x[pointArray[k]]
      y2 = bc$y[pointArray[k]]
      # COEFFICIENTS IN EQUATION POF THE LINE
      a = y2 - y1
      b = x2 - x1
      c = -a * x1 + b * y1
      denom = sqrt(a * a + b * b)
      lcoeffs = c(a, b, c, denom)
      # BOUNDARY POINTS TO CHECK DISTANCE FROM
      if ((pointArray[k] - pointArray[k - 1]) > 4) {
        x0 = bc$x[(pointArray[k - 1] + 1):(pointArray[k] - 1)]
        y0 = bc$y[(pointArray[k - 1] + 1):(pointArray[k] - 1)]
        v = cbind(x0, y0)
        v = stats::na.omit(v)
        dist = apply(v, 1, pointttolinedist, lcoeffs)
        indkeep = which(dist == max(dist))
        if (max(dist) > thresh) {
          tempArray = c(tempArray, (pointArray[k - 1] + indkeep[1]))
          alldone = 0
          n = n + 1
        }
      }
      tempArray = c(tempArray, pointArray[k])
    }
    # NOW DO LINE BACK TO START
    x1 = bc$x[pointArray[numpoints]]
    y1 = bc$y[pointArray[numpoints]]
    x2 = bc$x[pointArray[1]]
    y2 = bc$y[pointArray[1]]
    # COEFFICIENTS IN EQUATION POF THE LINE
    a = y2 - y1
    b = x2 - x1
    c = -a * x1 + b * y1
    denom = sqrt(a * a + b * b)
    lcoeffs = c(a, b, c, denom)
    # BOUNDARY POINTS TO CHECK DISTANCE FROM
    if ((bc$length - pointArray[numpoints]) > 4) {
      x0 = bc$x[(pointArray[numpoints] + 1):bc$length]
      y0 = bc$y[(pointArray[numpoints] + 1):bc$length]
      v = cbind(x0, y0)
      v = stats::na.omit(v)
      dist = apply(v, 1, pointttolinedist, lcoeffs)
      indkeep = which(dist == max(dist))
      if (max(dist) > thresh) {
        tempArray = c(tempArray, (pointArray[numpoints] + indkeep[1]))
        alldone = 0
        n = n + 1
      }
    }
    pointArray = tempArray
    numpoints = numpoints + n
  }
  return(pointArray)
}

pointttolinedist = function(v, lc) {
  x = v[1]
  y = v[2]
  numer = abs(lc[1] * x - lc[2] * y + lc[3])
  denom = lc[4]
  dist = numer / denom
  return (dist)
}

sqreucdist = function(v) {
  dist = (v[1] - v[3]) * (v[1] - v[3]) + (v[2] - v[4]) * (v[2] - v[4])
  return (dist)
}

polyAngle = function(v) {
  angle = 2.0 * pi
  a = sqrt(v[1])
  b = sqrt(v[2])
  if (abs((v[1] + v[2] - v[3]) / (2.0 * a * b) - 1.0) > 0.001) {
    angle = acos((v[1] + v[2] - v[3]) / (2.0 * a * b))
  }
  return(angle)
}

cooccur = function(mini_image_info, nc) {
  image = mini_image_info[[3]]
  mask = mini_image_info[[4]]
  lev1image = waveTran2D(image)
  lev2image = waveTran2D(lev1image)
  lev1image = doubleImage(lev1image)
  lev1image = lev1image[1:dim(image)[1], 1:dim(image)[2]]
  lev2image = doubleImage(doubleImage(lev2image))
  lev2image = lev2image[1:dim(image)[1], 1:dim(image)[2]]
  # calculate co-occurrence matrix between image and first level wavelet approximation
  cooc01 = getCoocMatrix(image, lev1image, mask, nc)
  # calculate co-occurrence matrix between first and second level wavelet approximation
  cooc12 = getCoocMatrix(lev1image, lev2image, mask, nc)
  # calculate co-occurrence matrix between image and second level wavelet approximation
  cooc02 = getCoocMatrix(image, lev2image, mask, nc)
  cooc = list(cooc01, cooc12, cooc02)
  names(cooc) <- c("cooc01", "cooc12", "cooc02")
  return(cooc)
}

waveTran2D = function(image) {
  temp = image
  width = ncol(image)
  height = nrow(image)
  ww = 0
  if (height %% 2 != 0)
    ww = 1
  if (ww == 1)
    temp = rbind(temp, temp[height, ])
  # do one level y transform
  w = apply(temp, 2, daub2, (height + ww), 1)
  # do one level x transform
  temp = w[1:((height + ww) / 2), ] / sqrt(2)
  ww = 0
  if (width %% 2 != 0)
    ww = 1
  if (ww == 1)
    temp = cbind(temp, temp[, width])
  w = apply(temp, 1, daub2, (width + ww), 1)
  output = w[1:((width + ww) / 2), ] / sqrt(2)
  return(t(output))
}

daub2 = function(a, n, isign) {
  D0 = 0.70710678
  D1 = 0.70710678
  wa = c(1:n)
  if (n >= 4) {
    nh = n / 2
    if (isign == 1) {
      i = 1
      for (j in seq(1, n, 2)) {
        wa[i] = D0 * a[j] + D1 * a[j + 1]
        wa[i + nh] = D1 * a[j] - D0 * a[j + 1]
        i = i + 1
      }
    }
    else if (isign == -1) {
      j = 1
      for (i in 1:nh) {
        wa[j] = D0 * a[i] + D1 * a[i + nh]
        wa[j + 1] = D1 * a[i] - D0 * a[i + nh]
        j = j + 2
      }
    }
    return(wa)
  }
}

doubleImage = function(image) {
  width = ncol(image)
  height = nrow(image)
  s = apply(image, 2, doublevector)
  output = apply(s, 1, doublevector)
  return(t(output))
}

doublevector = function(v) {
  n = length(v)
  s = rep(0, (2 * n))
  for (i in 1:n) {
    s[(2 * i) - 1] = v[i]
    s[(2 * i)] = v[i]
  }
  return(s)
}

getCoocMatrix = function(image1, image2, mask, nc) {
  cooc = matrix(0, nrow = nc, ncol = nc)
  image1 = rescale(image1, nc)
  image2 = rescale(image2, nc)
  for (i in 1:nrow(image1)) {
    for (j in 1:ncol(image1)) {
      if (mask[i, j] != 0) {
        cooc[image1[i, j], image2[i, j]] = cooc[image1[i, j], image2[i, j]] + 1
      }
    }
  }
  return(cooc)
}

rescale = function(image, nc) {
  scale = nc / (max(image) - min(image))
  image = (image - min(image)) * scale
}

calculateHaralickFeatures = function(glcm)
{
  o.hara <- matrix(0, 14)
  pglcm <- glcm / sum(glcm)
  nx <- ncol(pglcm)
  ny <- nrow(pglcm)
  px <- colSums(pglcm)
  py <- rowSums(pglcm)
  pxpy <- matrix(px, nx, ny) * t(matrix(py, nx, ny))
  px_y <- matrix(0, nx + ny)
  pxmy <- matrix(0, (nx + ny) / 2)
  vx <- 1:nx
  vy <- 1:ny
  mx <- sum(px * vx)
  my <- sum(py * vx)
  stdevx <- sum(px * (vx - mx) ^ 2)
  stdevy <- sum(py * (vy - my) ^ 2)
  hxy1_0 <- matrix(0, nx, ny)
  hxy2_0 <- matrix(0, nx, ny)
  hxy1_0 <- pglcm * log10(pxpy)
  hxy2_0 <- (pxpy) * log10(pxpy)
  hx <- -sum(px * log10(px), na.rm = TRUE)
  hy <- -sum(py * log10(py), na.rm = TRUE)
  hxy1 <- -sum(hxy1_0, na.rm = TRUE)
  hxy2 <- -sum(hxy2_0, na.rm = TRUE)
  op <- matrix(1:nx, nx, ny)
  oq <- t(op)
  spq <- matrix(1:nx, nx, ny) + t(matrix(1:ny, nx, ny))
  dpq <- abs(matrix(1:nx, nx, ny) - t(matrix(1:ny, nx,
                                             ny)))
  o.hara[1, 1] <- sum(pglcm ^ 2)
  o.hara[2, 1] <- sum(dpq ^ 2 * pglcm)
  o.hara[3, 1] <- sum(pglcm / (1 + dpq ^ 2))
  o.hara[4, 1] <- -sum(pglcm * log10(pglcm), na.rm = TRUE)
  o.hara[5, 1] <-
    sum((op - mx) * (oq - my) * pglcm / (sqrt(stdevx * stdevy)))
  if (stdevx * stdevy == 0)
    o.hara[5, 1] <- 0
  o.hara[6, 1] <- sum((op - ((mx + my) / 2)) ^ 2 * pglcm)
  o.hara[7, 1] <- sum(spq * pglcm)
  sen <- array(0, (2 * nx))
  den.1 <- array(0, nx)
  den.2 <- array(0, nx)
  pglcm2 <- cbind(pglcm[, nx:1])
  for (i in 2:nx) {
    sen[i] <- sum(diag(pglcm2[1:i, (nx - i + 1):nx]))
    den.1[i] <- sum(diag(pglcm[1:i, (nx - i + 1):nx]))
    sen[1] <- pglcm2[1, nx]
    den.1[1] <- pglcm[1, nx]
  }
  for (i in 1:(nx - 2)) {
    sen[i + nx] <- sum(diag(pglcm2[(i + 1):nx, 1:(nx -
                                                    i)]))
    den.2[nx - i] <- sum(diag(pglcm[(i + 1):nx, 1:(nx -
                                                     i)]))
  }
  sen[nx + nx - 1] <- pglcm2[nx, 1]
  den.2[1] <- pglcm[nx, 1]
  o.hara[8, 1] <- -sum(sen * log10(sen), na.rm = TRUE)
  den <- den.1 + den.2
  o.hara[9, 1] <- -sum(den * log10(den), na.rm = TRUE)
  o.hara[10, 1] <- sum(((dpq - o.hara[9]) ^ 2) * pglcm)
  o.hara[11, 1] <- sum(((spq - o.hara[8]) ^ 2) * pglcm)
  o.hara[12, 1] <- sqrt(1 - exp(-2 * abs(hxy2 - o.hara[4])))
  o.hara[13, 1] <- sum((spq - mx - my) ^ 3 * pglcm)
  o.hara[14, 1] <- sum((spq - mx - my) ^ 4 * pglcm)
  
  rownames(o.hara) <- c(
    "asm",
    "con",
    "idm",
    "ent",
    "cor",
    "var",
    "sav",
    "sen",
    "den",
    "dva",
    "sva",
    "f13",
    "sha",
    "pro"
  )
  return(o.hara)
}

intensityQuantiles = function(bc, intensities)
{
  qs = stats::quantile(intensities[, 3], probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  quantileVars <- vector(mode = "list", length = length(qs))
  for (i in 1:length(qs)) {
    s = subset(intensities, intensities[, 3] >= qs[i])
    d = as.vector(stats::dist(s[, 1:2]))
    quantileVars[[i]] = stats::var(d) / mean(d)
  }
  return(quantileVars)
}

densityCalc = function(jj, centroids, RandA, avrad) {
  ncells = length(centroids)
  nframesjj = length(centroids[[jj]][, 1])
  diam = 2 * avrad
  density = rep(0, nframesjj)
  for (i in 1:nframesjj) {
    f = RandA[[jj]][i, 1]
    if (is.na(f) == T) {
      density[i] = NA
    }
    else{
      x = centroids[[jj]][i, 1]
      y = centroids[[jj]][i, 2]
      r = RandA[[jj]][i, 2]
      #  DISTANCE TO OTHER CELLS
      for (k in 1:ncells) {
        if (k != jj) {
          nframesk = length(centroids[[k]][, 1])
          for (j in 1:nframesk) {
            f1 = RandA[[k]][j, 1]
            if ((f1 == f) & (is.na(f1) == F)) {
              x1 = centroids[[k]][j, 1]
              y1 = centroids[[k]][j, 2]
              r1 = RandA[[k]][j, 2]
              dist = sqrt(((x - x1) ^ 2) + ((y - y1) ^ 2))
              standard = r + diam + r1
              # d IS THE HEIGHT OF THE CIRCULAR SEGMENT OF THE NEARBY CELL TO BE INCLUDED IN THE DENSITY CALCULATION
              d = standard - dist
              if (dist < standard) {
                swap = 0
                if (d > r1) {
                  d = d - r1
                  swap = 1
                }
                # ANGLE SUBTENDED AT CENTRE
                alpha = abs(2 * cos((r1 - d) / r1))
                # PROPORTION OF CELL WITHIN CLOSE REGION
                prop = 1.0 / (2 * pi * (alpha - sin(alpha)))
                # SEE IF MORE THAN HALF THE CELL IS IN CLOSE REGION
                if (swap == 1) {
                  prop = 1 - prop
                }
                # IF ENTIRE CELL INSIDE CLOSE REGION
                if (standard > (dist + r1))
                  prop = 1.0
                A1 = RandA[[k]][j, 3]
                # AREA AROUND FIRST CELL
                A = 2 * pi * diam * (diam + 2 * r)
                density[i] = density[i] + 100.0 * prop * A1 / A
              }
            }
          }
        }
      }
    }
  }
  return(density)
}

calculateTrajArea <- function(x, y)
{
  xCentres <- stats::na.omit(x)
  yCentres <- stats::na.omit(y)
  numframes = length(xCentres)
  trajArea = ((max(xCentres) - min(xCentres)) * (max(yCentres) - min(yCentres))) /
    numframes
  return(trajArea)
}

fill_mask = function(v) {
  inds = which(v == 0)
  firstzero = inds[1]
  lastzero = inds[length(inds)]
  if (v[1] != 0)
    v[1:(firstzero - 1)] = -1
  if (v[length(v)] != 0)
    v[(lastzero + 1):length(v)] = -1
  return(v)
}
