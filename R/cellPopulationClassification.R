#' Train and test an ensemble of classifiers for cell population classification
#'
#'@description This function trains three classifiers (Linear Disriminant Analysis, Random Forest and Support Vector Machine) and uses these in ensemble to obtain final predictions for cell type classification of a test set. The function returns a matrix of classification results, including predictions from each of the classifiers individually and when used in ensemble. Note that if training and test sets are to be feature selected ahead of classification, their subsetted form should be used as input here.
#'
#' @param TrainingSet feature table to be used as a training set  
#' @param TestSet feature table to be used as a test set
#' @param TrainingLabels list of training set ground truth labels
#'
#' @return This function returns a matrix of classification results. Columns 1, 2, 3 and 4 contain test set predictions from LDA, RF, SVM and the ensemble respectively.
#' @export
cellPopulationClassification<-function(TrainingSet, TestSet, TrainingLabels)
{
  dataforscaling<-rbind(TrainingSet, TestSet)
  dataforscaling<-scale(dataforscaling)
  TrainingSet<-dataforscaling[c(1:nrow(TrainingSet)),]
  TestSet<-dataforscaling[-c(1:nrow(TrainingSet)),]
  
  ## classifier training
  ldamodel<-MASS::lda(TrainingSet, TrainingLabels)
  rfmodel <- randomForest::randomForest(TrainingLabels~., data = TrainingSet, ntree=200, mtry=5, importance=TRUE, norm.votes = TRUE)
  svmmodel<-e1071::svm(TrainingSet, TrainingLabels, kernel = 'radial', probability = TRUE)
  
  ## classifier testing
  ldapred = stats::predict(ldamodel, TestSet)
  rfpred = stats::predict(rfmodel, TestSet)
  svmpred = stats::predict(svmmodel, TestSet)
  
  ## ensemble classification, final predicted label based on majority vote
  classificationvotes<-cbind(as.character(ldapred$class), as.character(rfpred), as.character(svmpred))
  
  classificationvotes<-as.data.frame(classificationvotes)
  
  for(i in c(1:dim(TestSet)[1])) 
  {
    if((classificationvotes[i,1]==unique(TrainingLabels)[1] && classificationvotes[i,2]==unique(TrainingLabels)[1]) ||(classificationvotes[i,1]==unique(TrainingLabels)[1] && classificationvotes[i,3]==unique(TrainingLabels)[1]) 
       ||(classificationvotes[i,2]==unique(TrainingLabels)[1] && classificationvotes[i,3]==unique(TrainingLabels)[1])) 
    {
      classificationvotes[i,4] = unique(TrainingLabels)[1] 
    }
    else
    {
      classificationvotes[i,4] = unique(TrainingLabels)[2] 
    }
  }
  
  colnames(classificationvotes)<-c("LDA", "RF", "SVM", "Ensemble")
  return(classificationvotes)
}