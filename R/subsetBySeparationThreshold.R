#' Only include features above a set separation threshold
#'
#'
#' @description Subsets a data file to only include features above a set separation threshold
#'
#'
#' @param outputfile Feature table to be subset
#' @param separationscores list of separation scores, each element with a different separation threshold
#' @param t which index from separation scores to use
#'
#' @export
subsetBySeparationThreshold<-function(outputfile,separationscores,t)
{
  subbedOutputFile<-outputfile[ ,separationscores[[t]][,2]]
  return(subbedOutputFile)
}
