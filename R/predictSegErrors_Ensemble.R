#' Train and test an ensemble of decision trees for segmentation error prediction
#'
#'@description An alternative to \code{CellPhe::predictSegErrors()}. This function trains \code{K} sets of \code{num} decision trees and uses them in ensemble to predict whether or not new cells experience segmentation errors. Predictions from each of the \code{K} sets are obtained as described in \code{predictSegErrors()}. However this function adds further stringency by repeating \code{predictSegErrors()} a number of times and a cell is given a final classification of segmentation error if it receives a vote for this class in at least half of the repeated runs. This code outputs a list
#'of cells that were classified as segmentation error.
#'
#' @param smalldata Output from \code{prepareSegmentationTrainingSet()} - feature table of either segmentation errors or correctly segmented cells, whichever has the \bold{lowest} sample size
#' @param bigdata Output from \code{prepareSegmentationTrainingSet()} - feature table of either segmentation errors or correctly segmented cells, whichever has the \bold{greatest} sample size
#' @param smallclass Output from \code{prepareSegmentationTrainingSet()} - list of class labels for the class with the smallest sample size
#' @param bigclass Output from \code{prepareSegmentationTrainingSet()} - list of class labels for the class with the largest sample size
#' @param num Number of decision trees to be trained
#' @param K The number of repeated runs of CellPhe::predictSegErrors() to be performed
#' @param testset Test set for segmentation error predictions to be made
#' @param dataID List of test set identifiers (e.g. cell IDs)
#' @param proportion Proportion of votes needed for a final classification of segmentation error to be made (e.g. 0.7 if 70% of the votes are needed for segmentation error classification to be made)
#'
#' @return This function returns the list of identifiers that were predicted as segmentation errors, note these cells will have been predicted a segmentation error in at least \code{K/2} of the repeated classification runs
#'
#' @examples segmentation_errors <- predictSegErrors(segerrordata, correctsegdata, segerrorlabels, correctseglabels, 50, 10, testset, testset$cellnames, 0.7)
#' @export
predictSegErrors_Ensemble<-function(smalldata, bigdata, smallclass, bigclass, num, K, testset, dataID, proportion) 
{ 
  votes<-list()
  for(i in c(1:K))
  {
    segtest<-predictSegErrors(smalldata, bigdata, smallclass, bigclass, num, testset, dataID, proportion)
    votes[[i]]<-segtest
  }
  
  View(votes)
  
  list<-NULL
  for(i in c(1:K)) 
  { 
    list<-c(list, as.character(votes[[i]])) 
  } 
  View(list)
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