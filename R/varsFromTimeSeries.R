#' Calculate variables from the time series of each extracted cell feature.
#'
#' Reads in the time series for any pre-existing features as well as the 
#' output of \code{extractFeatures}. Variables are calculated from the time series, 
#' providing both summary statistics and indicators of time-series behaviour 
#' at different levels of detail obtained via wavelet analysis. 
#' @param df A data frame, the same as the output from \code{extractFeatures}.
#' @return A data frame with cells in rows and variables describing their behaviour over time in columns.
#' @export
varsFromTimeSeries = function(df) {
  cell_ids <- unique(df$CellID)
  num_cells = length(cell_ids)
  
  n_ts_vars <- 15
  n_old_features <- 72
  # Dataframe at minimum is 72 features, 3 id cols, 2 coord cols
  n_new_features <- ncol(df) - (n_old_features + 3 + 2)
  # Have a time-series summary for each feature + CellID + trajArea
  numcols <- n_ts_vars * (n_old_features + n_new_features) + 1 + 1
  output = matrix(NA, nrow = num_cells, ncol = numcols)
  
  for (j in 1:num_cells) {
    cell_id <- cell_ids[j]
    timeseries = df[ df$CellID == cell_id, setdiff(colnames(df), c("CellID", "FrameID", "ROI_filename", "xcentres", "ycentres"))]
    frame_ids <- df$FrameID[df$CellID == cell_id]
    numvars = dim(timeseries)[2]
    
    stats <- matrix(NA, nrow = 3, ncol = numvars)
    eleVars <- matrix(NA, nrow = 3, ncol = numvars)
    wVars <- matrix(NA, nrow = 9, ncol = numvars)
    
    # Check to see if need to interpolate any missing frames
    if (any(diff(frame_ids) > 1)) {
      # CALCULATE SUMMARY STATISTICS
      stats = apply(timeseries, 2, summaryStats)
      
      # Add NAs for missing frames
      missing_frame_ids <- setdiff(seq(min(frame_ids), max(frame_ids)), frame_ids)
      missing_frames <- data.frame(FrameID=missing_frame_ids)
      for (col in colnames(timeseries)) {
        missing_frames[[col]] <- NA
      }
      timeseries <- cbind(FrameID=frame_ids, timeseries)
      # merge missing frames into main df and drop FrameID
      timeseries <- rbind(timeseries, missing_frames)
      timeseries <- timeseries[order(timeseries$FrameID), ]
      timeseries <- timeseries[, -1]
      
      # INTERPOLATE MISSING VALUES
      timeseries = apply(timeseries, 2, interpolate)
    } else {
      # CALCULATE SUMMARY STATISTICS
      stats = apply(timeseries, 2, summaryStats)
    }
    
    # CALCULATE ASCENT, DESCENT AND MAX
    eleVars = apply(timeseries, 2, elevationVars)
    # CALCULATE VARIABLES FROM 3 LEVELS OF WAVELET DETAIL COEFFICIENTS
    wVars = apply(timeseries, 2, waveVars)
    
    vars = rbind(stats, eleVars, wVars)
    output[j, 1] <- cell_id
    output[j, 2:(numcols - 1)] = as.vector(vars)
    output[j, numcols] = calculateTrajArea(df$xcentres[df$CellID == cell_id], df$ycentres[df$CellID == cell_id])
  }
  cn = colnames(vars)
  rn = c(
    "mean",
    "std",
    "skew",
    "asc",
    "des",
    "max",
    "l1_asc",
    "l1_des",
    "l1_max",
    "l2_asc",
    "l2_des",
    "l2_max",
    "l3_asc",
    "l3_des",
    "l3_max"
  )
  name_perms <- expand.grid(rn, cn)
  names <- paste(name_perms$Var2, name_perms$Var1, sep = "_")
  names <- c(names, "trajArea")
  colnames(output) = c("CellID", names)
  res_df <- as.data.frame(output)
  res_df
}

summaryStats = function(v){
    m = mean(v)
    std = sqrt(stats::var(v))
    skew = 0
    if (max(v) > 0) skew = e1071::skewness(v, type = 2)                

	return(c(m, std, skew))
}

elevationVars = function(v){
    asc = 0
    des = 0
    max = v[1]
	for (k in 2:length(v)){
		diff = v[k] - v[k-1]
        if (diff > 0) {asc = asc + diff}
        else {des = des + diff}
        if (v[k] > max) max = v[k]
	}
    asc = asc/length(v)
    des = des/length(v)

	return(c(asc, des, max))
}

waveVars <- function(v)
{
 	wl = 0
 	# CALCULATE WAVELET DETAIL COEFFICIENTS 
	w = waveTran(v)
    for (i in 1:3){
 		 wl =  c(wl, detVars(w[[i]]))
    }
    wl = wl[-1]
    return(wl)
}

waveTran = function(v){
  n = length(v)
  ww = 0
  w <- vector(mode = "list", length = 3)
  if (n %% 2 != 0) ww = 1
  n = n - ww
  for (k in 1:3){
     # DO TRANSFORM
     newv = daub2(v, n, 1)
     w[[k]] = newv[((n/2) + 1):n]
     if (k < 3){
	     v = newv[1:(n/2)]/sqrt(2)
    	 n = n/2
         ww = 0
     	 if (n %% 2 != 0) ww = 1
         n = n - ww
      }
  } 
  return(w) 
}

detVars = function(v){
    asc = 0
    des = 0
    max = v[1]
	for (k in 1:length(v)){
        if (v[k] > 0) {asc = asc + v[k]}
        else {des = des + v[k]}
        if (v[k] > max) max = v[k]
	}
    asc = asc/length(v)
    des = des/length(v)

	return(c(asc, des, max))
}

interpolate <- function(v)
{
    nframes <- length(v)
    while (is.na(v[nframes]) == T){
      v = v[-nframes]
      nframes = nframes - 1
    }
    while (is.na(v[1]) == T){
      v = v[-1]
      nframes = nframes - 1
    }
   
    for (k in 1:nframes){
      r = 0
      if (is.na(v[k]) == T){
        r = 1
        while (((k+r) < nframes) & (is.na(v[k+r])) == T) r = r+1
      }
      value = (v[k+r]-v[k-1])/(r+1)
      if (r > 0){
        for (m in 1:r){
          v[k-1+m] = v[k-1] + m*value
        }
      }
      k = k+r
    }
    return(v[1:nframes])
}
