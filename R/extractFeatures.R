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

#' Calculates cell features from timelapse videos
#'
#' Calculates 74 features related to size, shape, texture and movement for each
#' cell on every non-missing frame, as well as the cell density around each
#' cell on each frame.
#' NB: while the ROI filenames are expected to be provided in \code{df} and found in 
#' \code{roi_folder}, the frame filenames are just expected to follow the naming convention
#' \code{<some text>-<FrameID>.tiff}, where FrameID is a 4 digit leading zero-padded number, corresponding 
#' to the \code{FrameID} column in \code{df}.
#'
#' @param df DataFrame where every row corresponds to a combination of a cell
#' tracked in a frame. It must have at least columns \code{CellID}, \code{FrameID} and
#' \code{ROI_filename} along with any additional features.
#' @param roi_folder A path to a directory containing multiple Report Object Instance
#' (ROI) files named in the format \code{cellid}-\code{frameid}.roi
#' @param frame_folder A path to a directory containing multiple frames in TIFF format.
#' It is assumed these are named under the pattern \code{<experiment name>-<frameid>.tif}, where 
#' \code{<frameid>} is a 4 digit zero-padded integer.
#' @param framerate The frame-rate, used to provide a meaningful measurement unit for velocity,
#'    otherwise a scaleless unit is implied with \code{framerate=1}.
#' @return A dataframe with 77+N columns (where N is the number of imported features)
#' and 1 row per cell per frame it's present in:
#' \itemize{
#'   \item{\code{FrameID}: the numeric frameID}
#'   \item{\code{CellID}: the numeric cellID}
#'   \item{\code{ROI_filename}: the ROI filename}
#'   \item{\code{...}: 74 frame specific features}
#'   \item{\code{...}: Any other data columns that were present in \code{df}}
#' }
#' @export
extractFeatures = function(df,
                           roi_folder,
                           frame_folder,
                           framerate = 1) {

  required_cols <- c("FrameID", "CellID", "ROI_filename")
  if (! all(required_cols %in% colnames(df))) {
    stop("Columns FrameID, CellID, and ROI_filename must be present in df.")
  }

  n_cells <- length(unique(df$CellID))
  
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
      stop(sprintf("Cannot find tif for frame id %d in %s. Check that the filename convention is as expected (?extractFeatures).", frame_id, frame_folder))
    }
    frame <- normaliseImage(tiff::readTIFF(tiff_fn), lower=0, upper=255)
  
    for (cell_id in cell_ids) {
      ids[row_num, ] <- c(cell_id, frame_id)
      roi_fn <-
        df |> dplyr::filter(CellID == cell_id, FrameID == frame_id) |> dplyr::distinct(ROI_filename) |> dplyr::pull(ROI_filename)
      if (length(roi_fn) > 1) {
        stop(sprintf("Error: found more than one ROI filename for Cell %d and Frame %d",
             cell_id,
             frame_id))
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
      # TODO vectorisable
      bl = boundary_coordinates$length
      bfeatures[row_num, 7] = bfeatures[row_num, 6] / (bl * bl)
      # MINIMAL BOX TO AREA RATIO:
      # TODO vectorisable
      bfeatures[row_num, 8] = (box[1] * box[2]) / bfeatures[row_num, 3]
      # RECTANGULARITY:
      # TODO vectorisable
      m = max(box[1], box[2])
      bfeatures[row_num, 9] = m / (box[1] + box[2])
      # FITTED POLYGON FEATURES (MAX_SIDE, MIN_ANGLE, ANGLE_VARIANCE, SIDE_LENGTH_VARIANCE)
      bfeatures[row_num, 10:13] <- polyClass(boundary_coordinates)
      
      # TEXTURE FEATURES
      # FIRST ORDER FEATURES
      tfeatures[row_num, 1] = mean(cell_pixels[, 3])
      tfeatures[row_num, 2] = sqrt(stats::var(cell_pixels[, 3]))
      tfeatures[row_num, 3] = e1071::skewness(cell_pixels[, 3], type = 2)
      # CALCULATE COOCCURRENCES MATRICES
      cooccurrence_levels = 10
      # TODO refactor to separate mask and image rather than passing in whole sub_image_info
      cooc = cooccur(sub_image_info, cooccurrence_levels)
      # HARALICK FEATURES
      haralickfeatures01 <-
        calculateHaralickFeatures(cooc$cooc01)
      haralickfeatures12 <-
        calculateHaralickFeatures(cooc$cooc12)
      haralickfeatures02 <-
        calculateHaralickFeatures(cooc$cooc02)
      for (k in 1:14) {
        tfeatures[row_num, (k + 3)] = haralickfeatures01[k]
        tfeatures[row_num, (k + 17)] = haralickfeatures12[k]
        tfeatures[row_num, (k + 31)] = haralickfeatures02[k]
      }
      # INTENSITY QUANTILE FEATURES
      tfeatures[row_num, 46:54] <- intensityQuantiles(boundary_coordinates, interior_pixels)
    
      row_num <- row_num + 1
    }
  }
  # Combine and add movement features
  res <- cbind(ids, bfeatures, tfeatures, xpos=xcentres, ypos=ycentres) |> as.data.frame()
  res <- res |> 
    dplyr::filter(!is.na(CellID), !is.na(FrameID)) |>  # Removes any frames with missing ROIs
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
           Vel = (framerate * dist_timestamp) / (FrameID - dplyr::lag(FrameID, default=0))) |>
    dplyr::ungroup() |>
    dplyr::select(-startx, -starty, -dist_timestamp) |>
    dplyr::select(Dis, Trac, D2T, Vel, dplyr::everything())

  # Calculate density 
  dens <- densityCalc(res)
  # Merge back in, cells-frames with NA mean they have 0 density
  res <- res |> 
    dplyr::left_join(dens, by=c("CellID"="cell1", "FrameID")) |>
    dplyr::mutate(dens = ifelse(is.na(dens), 0, dens))
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
  while (TRUE) {
    # find all values with -1
    neg <- which(matpix_type == -1)
    # generate coordinates of neighbours as x+1, x-1, x+nrow, x-nrow
    neighbours <- c(neg, neg+1, neg-1, neg+height, neg-height)
    neighbours <- unique(neighbours[neighbours > 0 & neighbours <= length(matpix_type)])  # As might find neighbours in first row/col
    positive_neighbours <- neighbours[matpix_type[neighbours] == 1]
    # set positive neighbours to be negative if any
    if (length(positive_neighbours) == 0) {
      break
    }
    matpix_type[positive_neighbours] <- -1
  }
  
  intensities <- unname(as.matrix(
    # expand.grid iterates over first column first, we want opposite
    cbind(expand.grid(1:height, 1:width)[, c(2, 1)],
    as.numeric(sub_image),
    as.numeric(matpix_type))
  ))
  
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
  
  mat_all <- array(c(mat01, mat12, mat02), dim=c(dim(mat01), 3))
  squaredLengths <- (mat_all[, 1, ] - mat_all[, 3, ])**2 + (mat_all[, 2, ] - mat_all[, 4, ])**2
  maxLength = max(sqrt(squaredLengths[, 1]))
  varLength = stats::var(sqrt(squaredLengths[, 1]))
  angles <- polyAngle(squaredLengths)
  minAngle = min(angles)
  varAngle = stats::var(angles)
  output = c(maxLength, minAngle, varAngle, varLength)
  return(output)
}

polygon = function(bc) {
  bc_df <- cbind(x=bc$x, y=bc$y)
  # FIND MAXIMUM DISTANCE FROM FIRST BOUNDARY POINT
  dd = as.matrix(stats::dist(bc_df))
  indkeep = which.max(dd[, 1])
  pointArray = as.vector(c(1, indkeep))
  previousArray = 1
  while (length(pointArray) > length(previousArray)) {
    tempArray = 1
    numpoints <- length(pointArray)
    
    for (k in 2:numpoints) {
      out <- poly_distance(bc_df, pointArray, k-1, k, FALSE)
      if (!is.na(out)) {
        tempArray <- c(tempArray, out)
      }
      # TODO better way to do this than vector concatenation?
      tempArray = c(tempArray, pointArray[k])
    }
    # NOW DO LINE BACK TO START
    out <- poly_distance(bc_df, pointArray, numpoints, 1, TRUE)
    if (!is.na(out)) {
      tempArray <- c(tempArray, out)
    }
    previousArray <- pointArray
    pointArray <- tempArray
  }
  return(pointArray)
}

poly_distance <- function(df, points, i1, i2, final, thresh=2.5) {
    if (final) {
      end <- nrow(df)
      end_v <- end
    } else {
      end <- points[i2]
      end_v <- end - 1
    }
  
    if ((end - points[i1]) <= 4) {
      return(NA)
    }
    
    c1 <- df[points[i1], ]
    c2 <- df[points[i2], ]
    diff_rows <- rev(c2 - c1)
    diff_rows[1] <- -diff_rows[1]
    c_vals <- sum(diff_rows * c1)
    denom2 <- sqrt(sum(diff_rows**2))
    v_df <- df[(points[i1]+1) : end_v, ]
    v_df <- v_df[complete.cases(v_df), ]  # TODO is this needed or can it be replaced?
    v_diff <- -v_df %*% diag(diff_rows)
    dist2 <- abs(v_diff[, 1] + v_diff[, 2] + c_vals) / denom2
    
    if (max(dist2) > thresh) {
      points[i1] + which.max(dist2)
    } else {
      NA
    }
}

polyAngle = function(v) {
  a <- sqrt(v[, 1])
  b <- sqrt(v[, 2])
  calc <- (v[, 1] + v[, 2] - v[, 3]) / (2.0 * a * b)
  ifelse(abs(calc - 1.0) > 0.001, acos(calc), 2.0*pi)
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
  w = apply(temp, 2, daub2, 1)
  # do one level x transform
  temp = w[1:((height + ww) / 2), ] / sqrt(2)
  ww = 0
  if (width %% 2 != 0)
    ww = 1
  if (ww == 1)
    temp = cbind(temp, temp[, width])
  w = apply(temp, 1, daub2, 1)
  output = w[1:((width + ww) / 2), ] / sqrt(2)
  return(t(output))
}

daub2 = function(a, isign) {
  n <- length(a)
  if (n < 3) return()
  D0 = 0.70710678
  D1 = 0.70710678
  wa = numeric(n)
  nh = n / 2
  D0_a <- D0 * a
  D1_a <- D1 * a
  odd <- seq(1, n, 2)
  even <- odd + 1
  if (isign == 1) {
    wa[1:nh] <- D0_a[even] + D1_a[odd]
    wa[(nh+1):n] <- D1_a[odd] - D0_a[even]
  }
  else if (isign == -1) {
    wa[even] <- D0_a[1:nh] + D1_a[(nh+1):n]
    wa[odd] <- D1_a[1:nh] - D0_a[(nh+1):n]
  }
  wa
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
  # TODO vectorize
  for (i in 1:n) {
    s[(2 * i) - 1] = v[i]
    s[(2 * i)] = v[i]
  }
  return(s)
}

getCoocMatrix = function(image1, image2, mask, nc) {
  image1 = rescale(image1, nc)
  image2 = rescale(image2, nc)
  cooc = matrix(0, nrow = nc, ncol = nc)
  masked <- mask == 1
  xs <- floor(image1[masked])
  ys <- floor(image2[masked])
  positive_coordinates <- xs > 0 & ys > 0
  flattened <- xs[positive_coordinates] + nc * (ys[positive_coordinates] - 1)
  counts <- tabulate(flattened)
  cooc[seq_along(counts)] <- counts
  cooc
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
  # I've spent ages trying to optimise this by calculating the distance matrix
  # up front to save calculating it 10 times, but the indexing or required matrix conversion
  # is so awkward that it doesn't result in a speed up
  qs = stats::quantile(intensities[, 3], probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  vapply(qs, function(thresh) {
    d = stats::dist(intensities[intensities[, 3] >= thresh, 1:2])
    stats::var(d) / mean(d)
  }, numeric(1))
}

densityCalc = function(df, radius_threshold=6) {
  # Calculate distance matrices between each cell for each Frame
  dists <- do.call('rbind', lapply(unique(df$FrameID), function(id) {
    sub_df <- df |> dplyr::filter(FrameID == id)
    dists <- as.matrix(stats::dist(sub_df[, c('xpos', 'ypos')], diag=TRUE, upper=TRUE))
    dists <- as.data.frame(dists)
    colnames(dists) <- sub_df$CellID
    dists$cell1 <- sub_df$CellID
    dists |> 
      tidyr::pivot_longer(-cell1, names_to="cell2", values_to="dist") |> 
      dplyr::mutate(FrameID = id)
  }))
  # Use the Rad column from the original dataframe to subset to cells that are close to each other
  dists |>
    dplyr::inner_join(df |> dplyr::select(FrameID, CellID, Rad), by=c("cell1"="CellID", "FrameID")) |>
    dplyr::filter(dist < radius_threshold * Rad, dist > 0) |>
    dplyr::group_by(FrameID, cell1) |>
    dplyr::summarise(dens = sum(1/dist)) |>
    dplyr::ungroup()
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

fill_mask <- function(v){
  inds = which(v == 0)
  if (length(inds) > 0){
    firstzero = inds[1]
    lastzero = inds[length(inds)]
    if (v[1] != 0) v[1:(firstzero-1)] = -1
    if (v[length(v)] != 0) v[(lastzero+1):length(v)] = -1
  }
  else{
    v = -v
  }
  return(v)
}
