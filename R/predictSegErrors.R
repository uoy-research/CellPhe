#' Train and test decision trees for segmentation error prediction
#'
#'
#'@description This function trains \code{num} decision trees and uses them to predict whether or not new cells experience segmentation errors. Final classifications are
#'made via a voting system, where a cell is classified as segmentation error if more than a defined \code{proportion} of decision trees predict it as such. This code outputs a list
#'of cells that were classified as segmentation error.
#'
#' @param smalldata Output from \code{prepareSegmentationTrainingSet()} - feature table of either segmentation errors or correctly segmented cells, whichever has the \bold{lowest} sample size
#' @param bigdata Output from \code{prepareSegmentationTrainingSet()} - feature table of either segmentation errors or correctly segmented cells, whichever has the \bold{greatest} sample size
#' @param smallclass Output from \code{prepareSegmentationTrainingSet()} - list of class labels for the class with the smallest sample size
#' @param bigclass Output from \code{prepareSegmentationTrainingSet()} - list of class labels for the class with the largest sample size
#' @param num Number of decision trees to be trained
#' @param testset Test set for segmentation error predictions to be made
#' @param dataID List of test set identifiers (e.g. cell IDs)
#' @param proportion Proportion of votes needed for a final classification of segmentation error to be made (e.g. 0.7 if 70% of the votes are needed for segmentation error classification to be made)
#'
#' @return This function returns the list of identifiers that were predicted as segmentation errors
#'
#' @examples segmentation_errors <- predictSegErrors(segerrordata, correctsegdata, segerrorlabels, correctseglabels, 50, testset, testset$cellnames, 0.7)
#' @export
predictSegErrors<-function(smalldata, bigdata, smallclass, bigclass,
                           num, testset, dataID, proportion) 
{ 
  n1 = length(smallclass)
  n2 = length(bigclass)
  treelist = list()
  for (i in 1:num)
  { 
    inds = sample.int(n2, n1)
    data = rbind(bigdata[inds,],smalldata)
    class = c(bigclass[1:n1], smallclass)
    data = data.frame(class, data)
    mytree = tree::tree(as.factor(class)~., data=data)
    treelist[[i]] <- mytree 
  } 
  
  trees<-treelist
  
  pred = matrix(" ", nrow = num, ncol = nrow(testset)) 
  for (i in 1:num)
  { 
    x = predict(trees[[i]], testset, type = "class")
    pred[i,] <- x 
  } 
  predictions<-pred
  
  testnumseg = vector(mode = "integer", length = ncol(predictions)) 
  for (i in 1:ncol(predictions))
  { 
    testnumseg[i] = length(which(predictions[,i] == 2)) 
  } 
  testvote = rep("nonseg", length = ncol(predictions))
  for (i in 1:ncol(predictions))
  { 
    if (testnumseg[i] >  proportion * nrow(predictions)) testvote[i] = "seg" 
  } 
  ind = which(testvote == "seg")
  segtest = dataID[ind]
  return(segtest)
}