#' Train and test decision trees for segmentation error prediction
#'
#'
#'@description This function trains \code{num} decision trees and uses them to predict whether or not new cells experience segmentation errors. Final classifications are
#'made via a voting system, where a cell is classified as segmentation error if more than a defined \code{proportion} of decision trees predict it as such. This code outputs a list
#'of cells that were classified as segmentation error.
#'
#' @param segerrors A feature table of ground truth segmentation errors
#' @param correctsegs A feature table of ground truth correctly segmented cells
#' @param num Number of decision trees to be trained
#' @param testset Test set for segmentation error predictions to be made
#' @param dataID List of test set identifiers (e.g. cell IDs)
#' @param proportion Proportion of votes needed for a final classification of segmentation error to be made (e.g. 0.7 if 70% of the votes are needed for segmentation error classification to be made)
#'
#' @return This function returns the list of identifiers that were predicted as segmentation errors
#'
#' @export
predictSegErrors<-function(segerrors, correctsegs,
                           num, testset, dataID, proportion) 
{ 
  seginfo<-prepareSegmentationTrainingSet(segerrors, correctsegs)
  smalldata = seginfo[[1]]
  bigdata = seginfo[[2]]
  smallclass = seginfo[[3]]
  bigclass = seginfo[[4]]
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