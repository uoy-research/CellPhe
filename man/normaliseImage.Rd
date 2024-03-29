% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractFeatures.R
\name{normaliseImage}
\alias{normaliseImage}
\title{Normalise an Image}
\usage{
normaliseImage(values, lower, upper)
}
\arguments{
\item{values}{The image to normalise, as a matrix.}

\item{lower}{The lower bound of the target normalisation range, as an integer.}

\item{upper}{The upper boundo of the target normalisation range, as an integer.}
}
\value{
The normalised image, as a matrix.
}
\description{
Normalises an image to a specified range.
}
