\name{summary.gseaSignificance}
\alias{summary.gseaSignificanceSign}
\alias{summary.gseaSignificanceVar}
\docType{data}
\title{
Obtain a data.frame with the pvalues and fdr for all signatures and
variables of a gseaSignificanceSign or gseaSignificanceVar object.
}
\description{
Builds a data.frame object that can easily be written to a csv file
containing the ES, NES, pval.ES, pval.NES and FDR.
}
\usage{
summary.gseaSignificanceSign(object,...)
}
\arguments{
\item{object}{ object of class \code{gseaSignificanceSign} or \code{gseaSignificanceVar}.}
\item{...}{ Arguments to be passed to \code{summary}.}
}
\references{
Aravind Subramanian, (October 25, 2005) \emph{Gene Set Enrichment Analysis}. 
\url{www.pnas.org/cgi/doi/10.1073/pnas.0506580102}
}
\examples{
#for examples see the help file of gseaSigntaures: ?gseaSignatures
}
\author{
Evarist Planet
}
\keyword{datasets}
