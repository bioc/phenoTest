\name{gseaData-class}
\Rdversion{1.1}
\docType{class}
\alias{gseaData-class}
\alias{show,gseaData-method}
\alias{getEs,gseaData-method}
\alias{getEsSim,gseaData-method}
\alias{getFcHr,gseaData-method}

\title{Class "gseaData"}
\description{
This class is an ES (enrichment score) and ES.sim (simulated enrichment
score) container that will be used in the GSEA (Gene Set Enrichment
Analysis) process.
There is one container for every gene signature.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gseaData", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} .}
    \item{\code{gseaSignaturesSign}:}{Object of class
      \code{"gseaSignaturesSign"} or \code{"gseaSignaturesVar"}.}
    \item{\code{gseaSignificanceSign}:}{Object of class
      \code{"gseaSignificanceSign"} or \code{"gseaSignificanceVar"}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{getEs}{\code{signature(x = "gseaData")}: Returns the
      enrichment scores.}
    \item{getEsSim}{\code{signature(x = "gseaData")}: Returns
      the simulated enrichemt scores (the ones obtained after
      permutations).} 
    \item{getFcHr}{\code{signature(x = "gseaData")}: Returns
      the fold change and/or the hazard ratio that were used to compute
      the enrichment scores. }
	 }
}
\author{Evarist Planet}
\examples{
showClass("gseaSignaturesSign")
}
\keyword{classes}
