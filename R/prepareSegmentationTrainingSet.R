#' Prepare a training set for classification of cell segmentation errors
#'
#'
#' @description This function prepares a training set that can be used for classification of cell segmentation errors. This includes use of the
#'\code{smotefamily::SMOTE()} function to over-sample the minority class to avoid biased classifier training. This function outputs feature tables
#'for ground truth segmentation errors and correctly segmented cells (with smote data included for the previously under-sampled class), as 
#'well as a list of ground truth labels for each class.
#'
#'
#' @param segerrors A feature table of ground truth segmentation errors
#' @param correctsegs A feature table of ground truth correctly segmented cells
#'
#' 

prepareSegmentationTrainingSet<-function(segerrors,correctsegs)
{
  Segerrortraining <- segerrors
  Correctsegtraining <- correctsegs
  
  ## collate correct segmentation and segmentation error data into one data frame
  data = rbind(Segerrortraining, Correctsegtraining) 
  segerrordata<-Segerrortraining[,-1] 
  correctsegdata<-Correctsegtraining[,-1]
  alldata = rbind(segerrordata, correctsegdata)
  
  ## first column lists ground truth data labels
  class1 = rep("segerror", dim(segerrordata)[1])
  class2 = rep("correct", dim(correctsegdata)[1]) 
  class = c(class1, class2)
  all = data.frame(class, alldata)
  
  ## SMOTE() to over-sample segmentation errors 
  smoteseg<-smotefamily::SMOTE(all[,-1], all$class, K = 3, dup_size = 1)
  
  ## add smote segmentation errors to the training sets and update segmentation error data table
  smotesegsyn_data<-smoteseg$syn_data[,c((dim(segerrors)[2]), 1:(dim(segerrors)[2]-1))] 
  all<-rbind(all, smoteseg$syn_data) 
  segdata<-rbind(segerrordata, smotesegsyn_data[,-1]) 
  newclass<-c(class, smotesegsyn_data[,1]) 
  class1 = rep("segerror", dim(segerrordata)[1])
  seginfo<-list(segerrordata, correctsegdata, class1, class2)
  return(seginfo)
}