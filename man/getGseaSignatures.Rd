\name{get gseaSignatures' elements}
\alias{getEs}
\alias{getEsSim}
\alias{getNes}
\alias{getFcHr}
\docType{data}
\title{
Substract element's of a gseaSignaturesSign or gseaSignaturesVar object (obtained using the
gseaSignatures function). 
}
\description{
getEs returns ES (enrichment scores) getEsSim returns simulated ES
(needed to compute pvals), getNes returns NES (normalized enrichment
scores) and getFcHr returns the fold changes or hazard used to compute
the ES, simulated ES and NES.
}
\usage{
getEs(x)
getEsSim(x)
getNes(x)
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