\name{plot.gseaSignatures}
\alias{plot.gseaSignaturesSign}
\alias{plot.gseaSignaturesVar}
\docType{data}
\title{
GSEA-like Plot.
}
\description{
Builds a GSEA plot using a gseaSignature object (one of
\code{gseaSignaturesSign} or \code{gseaSignaturesVar} obtained with the
\code{gseaSignatures} function) and a gseaSignificance object (one of
\code{gseaSignificanceSign} or \code{gseaSignificanceVar} obtained with
the \code{gseaSignificance} function).
}
\usage{
plot.gseaSignaturesSign(x,gseaSignificance,es.ylim,nes.ylim,es.nes="both",...)
}
\arguments{
\item{x}{ object of class \code{gseaSignaturesSign} or
  \code{gseaSignaturesVar}.}
\item{gseaSignificance}{ object of class \code{gseaSignificanceSign} or \code{gseaSignificanceVar}.}
\item{es.ylim}{ ylim values for the ES plot.}
\item{nes.ylim}{ ylim values for the NES plot.}
\item{es.nes}{ can be "es" if we want to plot enrichment score, "nes" if
  we want to plot normalised enrichment scores or "both"if we want to
  plot them both.}
\item{...}{ Arguments to be passed to \code{plot}.}
}
\references{
Aravind Subramanian, (October 25, 2005) \emph{Gene Set Enrichment Analysis}. 
\url{www.pnas.org/cgi/doi/10.1073/pnas.0506580102}
}
\examples{
#for examples see the help file of gseaSigntaures: ?gseaSignatures
}
\seealso{
plot.gseaSignaturesSign, plot.gseaSignaturesVar
}
\author{
Evarist Planet
}
\keyword{datasets}
