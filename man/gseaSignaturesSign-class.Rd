\name{gseaSignaturesSign-class}
\Rdversion{1.1}
\docType{class}
\alias{gseaSignaturesSign-class}
\alias{getEs,gseaSignaturesSign-method}
\alias{getEsSim,gseaSignaturesSign-method}
\alias{getFcHr,gseaSignaturesSign-method}
\alias{gseaSignificance,gseaSignaturesSign-method}

\title{Class "gseaSignaturesSign"}
\description{
This class is an ES (enrichment score) and ES.sim (simulated enrichment
score) container that will be used in the GSEA (Gene Set Enrichment
Analysis) process.
There is one container for every gene signature.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gseaSignaturesSign", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} .}
    \item{\code{gseaSignatures}:}{Object of class
      \code{"gseaSignatures"} This is the object that will contain the
      ES and ES.sim. }
    \item{\code{es.sim.gam}:}{Object of class \code{"matrix"} enrichment
      scores computed with the gam method. }
    \item{\code{fc.hr}:}{Object of class \code{"character"} fold change
      or hazard ratio used to compute the enrichment scores. }
    \item{\code{s}:}{Object of class \code{"logical"} The subset of
      genes we are interested in.} 
    \item{\code{test}:}{Object of class \code{"character"} The
      statistical test that will be used. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{getEs}{\code{signature(x = "gseaSignaturesSign")}: Returns the
      enrichment scores.}
    \item{getEsSim}{\code{signature(x = "gseaSignaturesSign")}: Returns
      the simulated enrichemt scores (the ones obtained after
      permutations).} 
    \item{getFcHr}{\code{signature(x = "gseaSignaturesSign")}: Returns
      the fold change and/or the hazard ratio that were used to compute
      the enrichment scores. }
    \item{gseaSignificance}{\code{signature(x = "gseaSignaturesSign")}:
      This is the next step in the process of performing GSEA. This
      funtion will test if the gene sets are enriched. }
	 }
}
\author{Evarist Planet}
\examples{
showClass("gseaSignaturesSign")
}
\keyword{classes}
