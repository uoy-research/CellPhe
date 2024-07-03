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
#' @param dup_size Passed to smotefamily::SMOTE. The number of times to
#' duplicate the minority class.
#'
#' @return This function returns the list of identifiers that were predicted as segmentation errors

#' @export
predictSegErrors<-function(segerrors, correctsegs,
                           num, testset, dataID, proportion, dup_size=1)
{ 
  seginfo<-prepareSegmentationTrainingSet(segerrors, correctsegs, dup_size)
  smalldata = seginfo$minority_data
  bigdata = seginfo$majority_data
  smallclass = seginfo$minority_class
  bigclass = seginfo$majority_class
  n1 = nrow(smalldata)
  n2 = nrow(bigdata)
  treelist = list()
  for (i in 1:num)
  { 
    inds = sample.int(n2, n1)
    data = rbind(bigdata[inds,],smalldata)
    class = c(rep(bigclass, n1), rep(smallclass, n1))
    data = data.frame(class, data)
    mytree = tree::tree(as.factor(class)~., data=data)
    treelist[[i]] <- mytree 
  } 
  
  trees<-treelist
  
  pred = matrix(" ", nrow = num, ncol = nrow(testset)) 
  for (i in 1:num)
  { 
    x = stats::predict(trees[[i]], testset, type = "class")
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