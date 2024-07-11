#' Identify heterogeneous cell clusters
#'
#' @description This function performs hierarchical clustering on a given data set in order to identify \code{k} heterogeneous cell clusters, where \code{k} is pre-defined.
#' @param dataset feature table for clustering
#' @param dataID a list of cell identifiers (e.g. cell IDs)
#' @param k number of desired clusters
#'
#' @return This function plots a dendogram of the obtained hierarchical clustering results, and outputs a matrix where each column \code{n} lists the cells that were assigned cluster \code{n}.
#' @export
identifyHeterogeneousClusters<-function(dataset, dataID, k)
{
  d<-stats::dist(scale(dataset), method = "euclidean")
  hierclust<-factoextra::hcut(d, hc_func = "agnes", hc_method = "ward.D", hc_metric = "euclidean", k = k)
  p <- factoextra::fviz_dend(hierclust, show_labels = FALSE)+ggplot2::theme_classic(base_size = 20)
  print(p)
  
  cellclusters<-matrix(nrow = dim(dataset)[1], ncol = k)
  for(j in 1:k)
  {
    l = 1
    for (i in 1:dim(dataset)[1]) 
    {
      if(hierclust$cluster[i] == j)
      {
        cellclusters[l,j] = as.character(dataID[i])
        l=l+1
        
      }
    }
    
  }
  return(cellclusters)
}
