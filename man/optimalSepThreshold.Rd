% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimalSepThreshold.R
\name{optimalSepThreshold}
\alias{optimalSepThreshold}
\title{Determine optimal separation threshold}
\usage{
optimalSepThreshold(group1data, group2data, seps)
}
\arguments{
\item{group1data}{Feature table for group 1}

\item{group2data}{Feature table for group 2}

\item{group1name}{A name for group 1 cells}

\item{group2name}{A name for group 2 cells}
}
\value{
A character vector containing the names of the optimal set of features
}
\description{
Determines the optimal separation threshold using the method described in the CellPhe paper
}
