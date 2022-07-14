#' Calculate variables from the time series of each extracted cell feature.
#'
#' Reads in the time series for any pre-existing features as well as the 
#' output of \code{extractFeatures}. Variables are calculated from the time series, 
#' providing both summary statistics and indicators of time-series behaviour 
#' at different levels of detail obtained via wavelet analysis. Name is an 
#' experiment identifier, used together with the cell identifiers in 
#' original_IDs to provide rownames in the output matrix.
#' @param features A list of matrices, the same as the output from \code{copyPhaseFeatures}
#' @param new_features A list of two lists, the same as the output from \code{copyFeatures}
#' @param expname A experiment identifier as a string
#' @param originalIDs A list of integers
#' @return A matrix with cells in rows and variables describing their behaviour over time in columns.
#' @export
varsFromTimeSeries = function(features, new_features, expname, originalIDs){
	trajArea = new_features[[1]]

	num = length(new_features[[2]])
	timeseries <- vector(mode = "list", length = num)
	
	  numcols <- 1109 + ncol(features[[1]])
    output = matrix(NA, nrow = num, ncol = numcols)
	
	for (j in 1:num){
     	timeseries[[j]] = cbind(features[[j]], new_features[[2]][[j]])
        numvars = dim(timeseries[[j]])[2]
 
		stats <- matrix(NA, nrow = 3, ncol = numvars)
		eleVars <- matrix(NA, nrow = 3, ncol = numvars)
		wVars <- matrix(NA, nrow = 9, ncol = numvars)

        # Use first column of newfeatures which may have additional frames considered "missing"
        firstnewfeature = ncol(features[[j]])+1
    	missingInd = which(is.na(timeseries[[j]][,firstnewfeature]) == T)
        
        if (length(missingInd) > 0){
		    # add NAs for original features
		    timeseries[[j]][missingInd,1:ncol(features[[j]])] = NA
			# INTERPOLATE MISSING VALUES
			timeseries[[j]] = apply(timeseries[[j]], 2, interpolate)
  			# CALCULATE SUMMARY STATISTICS 
			stats = apply(timeseries[[j]][-missingInd,], 2, summaryStats)
  	    }
  	    else {
  			# CALCULATE SUMMARY STATISTICS 
			stats = apply(timeseries[[j]], 2, summaryStats)  	    
  	    }
 
  		# CALCULATE ASCENT, DESCENT AND MAX
  		eleVars = apply(timeseries[[j]], 2, elevationVars)
  		# CALCULATE VARIABLES FROM 3 LEVELS OF WAVELET DETAIL COEFFICIENTS
   		wVars = apply(timeseries[[j]], 2, waveVars)  

	    vars = rbind(stats, eleVars, wVars)
	    cn = colnames(vars)
  		rn = c("mean", "std", "skew", "asc", "des", "max", "l1_asc", "l1_des", "l1_max","l2_asc", "l2_des", "l2_max", "l3_asc", "l3_des", "l3_max")
        names = rep(NA, length = (length(cn)*length(rn)+1))
        m = 1
        for (i in 1:length(cn)){
        	for (k in 1:length(rn)){
               names[m] = paste(cn[i], rn[k], sep = "_")
               m = m+1
            }
        }
        names[m] = "trajArea"
	    output[j,1:(numcols-1)] = as.vector(vars)
        output[j,numcols] = trajArea[[j]]
        colnames(output) = names
        row.names(output) = paste(expname, originalIDs, sep = "_")
	}
    return(output)
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
	nframes = length(v)
  	while (is.na(v[nframes]) == T){
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
