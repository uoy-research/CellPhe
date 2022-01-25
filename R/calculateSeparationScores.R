#' Calculate a table containing separation scores for each variable
#' 
#' @description This function can be used to calculate the separation score for each variable, where a feature selection threshold can be defined by the user to identify discriminatory variables
#'
#' @param group1data Output file of variables for group 1
#' @param group2data Output file of variables for group 2
#' @param threshold Separation threshold, only features achieving separation greater than or equal to this threshold will be output
#'
#' @return A data frame of separation values where column 1 is a list of variable indices, column 2 a list of variable names, and column 3 the corresponding separation scores
#' @export
#'
#' @examples separationscores<-calculateSeparationScores(Untreated, Treated, 0.1)
calculateSeparationScores<-function(group1data, group2data, threshold)
{
  group1data<-as.data.frame(group1data)
  group2data<-as.data.frame(group2data)
  separationscores<-data.frame(nrow = length(group1data), ncol = 3)
  for(i in c(1:length(group1data)))
  {
    separationscores[i,1]=i
    separationscores[i,2]=colnames(group1data)[i]
    Vw = (((dim(group1data)[1]-1)*var(group1data[,i]))+((dim(group2data)[1]-1)*var(group2data[,i])))/(dim(group1data)[1]+dim(group2data)[1]-2)
    overmean = (dim(group1data)[1]*mean(group1data[,i]) + (dim(group2data)[1]*mean(group2data[,i])))/(dim(group1data)[1]+dim(group2data)[1])
    Vb = (((dim(group1data)[1]*(mean(group1data[,i])-overmean)^2))+((dim(group2data)[1]*(mean(group2data[,i])-overmean)^2)))/(dim(group1data)[1]+dim(group2data)[1]-2)
    separation = Vb/Vw
    separationscores[i,3]=separation
  }
  separationscores[,3]=separationscores[,3]
  separationscores<-subset(separationscores, separationscores[,3] >= threshold)
  return(separationscores)
}