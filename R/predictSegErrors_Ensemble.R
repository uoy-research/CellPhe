#' Train and test an ensemble of decision trees for segmentation error prediction
#'
#'@description An alternative to \code{predictSegErrors()}. This function trains \code{K} sets of \code{num} decision trees and uses them in ensemble to predict whether or not new cells experience segmentation errors. Predictions from each of the \code{K} sets are obtained as described in \code{predictSegErrors()}. However this function adds further stringency by repeating \code{predictSegErrors()} a number of times and a cell is given a final classification of segmentation error if it receives a vote for this class in at least half of the repeated runs. This code outputs a list
#'of cells that were classified as segmentation error.
#'
#' @param segerrors A feature table of ground truth segmentation errors
#' @param correctsegs A feature table of ground truth correctly segmented cells
#' @param num Number of decision trees to be trained
#' @param K The number of repeated runs of predictSegErrors() to be performed
#' @param testset Test set for segmentation error predictions to be made
#' @param dataID List of test set identifiers (e.g. cell IDs)
#' @param proportion Proportion of votes needed for a final classification of segmentation error to be made (e.g. 0.7 if 70% of the votes are needed for segmentation error classification to be made)
#' @param dup_size Passed to smotefamily::SMOTE. The number of times to
#' duplicate the minority class.
#'
#' @return This function returns the list of identifiers that were predicted as segmentation errors, note these cells will have been predicted a segmentation error in at least \code{K/2} of the repeated classification runs
#' @export
predictSegErrors_Ensemble<-function(segerrors, correctsegs, num, K, testset, dataID, proportion, dup_size)
{ 
  votes<-list()
  for(i in c(1:K))
  {
    segtest<-predictSegErrors(segerrors, correctsegs, num, testset, dataID, proportion, dup_size)
    votes[[i]]<-segtest
  }
  
  list<-NULL
  for(i in c(1:K)) 
  { 
    list<-c(list, as.character(votes[[i]])) 
  } 

  segtest<-vector()
  j=1 
  freqdata<-as.data.frame(table(list)) 
  for(i in c(1:dim(freqdata)[1]))
  { 
    if(freqdata$Freq[i] >= K/2) 
    { 
      segtest[j] = as.character(freqdata$list[i]) 
      j=j+1 
    }
  }
  return(segtest)
}