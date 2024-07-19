#' Determine optimal separation threshold
#'
#'
#' @description Determines the optimal separation threshold using the method described in the CellPhe paper
#'
#'
#' @param group1data Feature table for group 1
#' @param group2data Feature table for group 2
#' @param group1name A name for group 1 cells
#' @param group2name A name for group 2 cells
#' 
#'
#' @return an optimal separation threshold
#' 
#' @export
optimalSepThreshold<-function(group1data, group2data, seps){
  varinds = seps[,1]
  scores = seps[,3]
  ordered = scores[order(-scores)]
  mins = min(ordered)
  maxs = max(ordered)
  n = c(0:(length(ordered)-1))
  d = rep(0,length(ordered))
  for (i in 1:length(ordered)){
    d[i] = ((mins-maxs)*n[i])/(length(ordered)-1) + maxs - ordered[i]
  }
  ind = which(d == max(d))
  thresh = ordered[ind]
  chosen = seps[, 2][which(scores > thresh)]
  return(chosen)
}
