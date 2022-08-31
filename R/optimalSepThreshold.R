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

optimalSepThreshold<-function(group1data, group2data){
  thresholds = c(0,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5)
  
  separationscores<-lapply(thresholds, CellPhe::calculateSeparationScores, group1data = group1data, group2data = group2data)
  
  ErrRate = matrix(nrow = length(thresholds), ncol = 1)
  group1names<-rep("group1", dim(group1data)[1])
  group2names<-rep("group2", dim(group2data)[1])
  group1<-cbind(group1names,group1data)
  colnames(group1)[colnames(group1) == 'group1names'] <- 'Group'
  group2<-cbind(group2names, group2data)
  colnames(group2)[colnames(group2) == 'group2names'] <- 'Group'
  size1<-dim(group1)[1]
  size2<-dim(group2)[1]
  
  if(size1 > size2)
  {
    sample<-sample(1:size1, size2, replace = FALSE)
    group1<-group1[sample,]
  }
  
  if(size2 > size1)
  {
    sample<-sample(1:size2, size1, replace = FALSE)
    group2<-group2[sample,]
  }
  
  Training<-rbind(group1,group2)
  
  for (i in (1:length(thresholds)))
  {
    correct = 0
    if(length(separationscores[[i]][[1]]) > 5)
    {
      subtrain<-subsetBySeparationThreshold(Training, separationscores, i)
      classificationresults = cellPopulationClassification(subtrain, subtrain, as.factor(Training[,1]))
      for(j in c(1:dim(Training)[1]))
      {
        if(classificationresults[j,4] == Training[j,1])
        {
          correct = correct+1
        }
      }
    }
    ErrRate[i] = 1-(correct/(dim(Training)[1]))
  }
  
  
  
  # choosing optimal separation threshold
  
  for (i in c(1:(length(ErrRate)-1)))
  {
    increase = (ErrRate[i+1] - ErrRate[i])*100
    
    if((increase > 1) == TRUE)
    {
      break
    }
  }
  
  optimal = thresholds[i-1]
  return(optimal)
  
}