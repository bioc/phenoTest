\name{export2CSV}
\alias{export2CSV}
\docType{data}
\title{
Export object to comma-separated text file.
}
\description{
Saves object as comma-separated text file (CSV), using \code{write.csv}.
}
\usage{
export2CSV(x, file, row.names=FALSE, ...)
}
\arguments{
\item{x}{ object to be exported. Currently methods for objects of class epheno (produced with
  ExpressionPhenoTest function) are implemented.}
\item{file}{ Name of the file where the results are to be saved}
\item{row.names}{Passed on to \code{write.csv}}
\item{...}{Other arguments to be passed on to \code{write.csv}}
}
\keyword{datasets}
