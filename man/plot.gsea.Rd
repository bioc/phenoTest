\name{plot.gseaData}
\alias{plot.gseaData}
\docType{data}
\title{
GSEA-like Plot.
}
\description{
Builds a GSEA plot using a \code{gseaData} object.
\code{gseaData} object can be obtained with the \code{gsea} function.
}
\usage{
plot.gseaData(x,selGsets,selVars,...)
}
\arguments{
\item{x}{ this has to be of class \code{gseaData}}
\item{selGsets}{ object of class \code{character} containing the names
  of the gene sets that we want to plot.}
\item{selVars}{ object of class \code{character} containing the names of
  the variables that we want to plot.}
\item{...}{ Arguments to be passed to \code{plot}.}
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
