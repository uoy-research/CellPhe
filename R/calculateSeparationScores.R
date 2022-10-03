#' Calculate a table containing separation scores for each variable
#' 
#' @description This function can be used to calculate the separation score for each variable, where a feature selection threshold can be defined by the user to identify discriminatory variables.
#' Note that \code{group1data} and \code{group2data} should include only columns of variables (i.e remove cell identifier columns prior to function use)
#'
#' @param group1data Feature table for group 1
#' @param group2data Feature table for group 2
#' @param threshold Separation threshold, only features achieving separation greater than or equal to this threshold will be output
#'
#' @return A data frame of separation values where column 1 is a list of variable indices, column 2 a list of variable names, and column 3 the corresponding separation scores
#' @export
calculateSeparationScores<-function(group1data, group2data, threshold = 0, calculateOptimalThresh = FALSE)
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
  
  if(calculateOptimalThresh == TRUE)
  {
    optThresh<-optimalSepThreshold(group1data, group2data)
    separationscores<-subset(separationscores, separationscores[,3] >= optThresh)
  }
  
  else
  {
    separationscores<-subset(separationscores, separationscores[,3] >= threshold)
  }
  colnames(separationscores)<-c("VarIndex", "VarName", "SepScore")
  return(separationscores)
}