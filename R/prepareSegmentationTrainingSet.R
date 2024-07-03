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
#' @param dup_size Passed to smotefamily::SMOTE. The number of times to
#' duplicate the minority class.
#'
#' @export
prepareSegmentationTrainingSet<-function(segerrors,correctsegs, dup_size=1)
{
  alldata = rbind(segerrors, correctsegs)[, -1]
  
  ## first column lists ground truth data labels
  class1 = rep("segerror", nrow(segerrors))
  class2 = rep("correct", nrow(correctsegs))
  class = c(class1, class2)
  
  ## SMOTE() to over-sample segmentation errors 
  smoteseg<-smotefamily::SMOTE(alldata, class, K = 3, dup_size = dup_size)
  
  ## add smote segmentation errors to the training sets and update segmentation error data table
  smotesegsyn_data<-smoteseg$syn_data[,c((dim(segerrors)[2]), 1:(dim(segerrors)[2]-1))] 
  list(
    minority_data = rbind(
      smoteseg$orig_P,
      smoteseg$syn_data
    )[, -ncol(smoteseg$orig_P)],
    majority_data=smoteseg$orig_N[, -ncol(smoteseg$orig_N)],
    minority_class=smoteseg$orig_P$class[1],
    majority_class=smoteseg$orig_N$class[1]
  )
}
