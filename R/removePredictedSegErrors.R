## keep only correctly segmented cells
#' Remove predicted segmentation errors from test set
#'
#'@description This function can be used to automate removal of predicted segmentation errors from the test set
#'
#' @param testset Test set for segmentation error predictions to be made
#' @param k Column index for the column of cell identifiers within the test set, e.g. if column 2 of the test set is a column of cell identifiers then here k = 2
#' @param predictedSegErrors Output from either CellPhe::predictSegErrors() or CellPhe::predictSegErrors_Ensemble(), a list of cell identifiers for cells classified as segmentation error
#'
#' @examples removePredictedSegErrors(testset, 2, predicted_seg_errors)
removePredictedSegErrors<-function(testset, k, predictedSegErrors)
{
  testset<-subset(testset, testset[,k] %in% predictedSegErrors == FALSE)
  return(testset)
}