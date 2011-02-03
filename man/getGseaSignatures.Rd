\name{get gseaSignatures' elements}
\alias{getEs}
\alias{getEsSim}
\alias{getFcHr}
\docType{data}
\title{
Substract element's of a gseaSignaturesSign or gseaSignaturesVar object (obtained using the
gseaSignatures function). 
}
\description{
getEs returns ES (enrichment scores) getEsNes returns simulated ES
(needed to compute pvals) and getFcHr returns the fold changes or hazard
ratios used to compute the ES and simulated ES.
}
\usage{
getEs(x)
getEsSim(x)
getFcHr(x)
}
\arguments{
\item{x}{ an \code{gseaSignaturesSign} or \code{gseaSignaturesVar} object. Those
  objects are obtained using the \code{gseaSignatures} function.}
}
\keyword{datasets}
\author{
Evarist Planet
}