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
#' @return A character vector containing the names of the optimal set of features
optimalSepThreshold<-function(group1data, group2data, seps){
  scores = seps[,3]
  ordered = scores[order(-scores)]
  mins = min(ordered)
  maxs = max(ordered)
  indices = 0:(length(ordered)-1)
  dists = ((mins-maxs)*indices)/(length(ordered)-1) + maxs - ordered
  ind = which.max(dists)
  thresh = ordered[ind]
  seps[, 2][which(scores > thresh)]
}
