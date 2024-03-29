% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellPopulationClassification.R
\name{cellPopulationClassification}
\alias{cellPopulationClassification}
\title{Train and test an ensemble of classifiers for cell population classification}
\usage{
cellPopulationClassification(TrainingSet, TestSet, TrainingLabels)
}
\arguments{
\item{TrainingSet}{feature table to be used as a training set}

\item{TestSet}{feature table to be used as a test set}

\item{TrainingLabels}{list of training set ground truth labels}
}
\value{
This function returns a matrix of classification results. Columns 1, 2, 3 and 4 contain test set predictions from LDA, RF, SVM and the ensemble respectively.
}
\description{
This function trains three classifiers (Linear Discriminant Analysis, Random Forest and Support Vector Machine) and uses these in ensemble to obtain final predictions for cell type classification of a test set. The function returns a matrix of classification results, including predictions from each of the classifiers individually and when used in ensemble. Note that if training and test sets are to be feature selected ahead of classification, their subsetted form should be used as input here.
}
