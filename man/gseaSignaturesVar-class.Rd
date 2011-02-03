\name{gseaSignaturesVar-class}
\Rdversion{1.1}
\docType{class}
\alias{gseaSignaturesVar-class}
\alias{getEs,gseaSignaturesVar-method}
\alias{getEsSim,gseaSignaturesVar-method}
\alias{getFcHr,gseaSignaturesVar-method}
\alias{gseaSignificance,gseaSignaturesVar-method}

\title{Class "gseaSignaturesVar"}
\description{This class is an ES (enrichment score) and ES.sim (simulated enrichment
score) container that will be used in the GSEA (Gene Set Enrichment
Analysis) process.
There is one container for every phenotype. Every one of this containers
(of class \code{gseaSignaturesSign}) is a container itsef and has the
enrichment scores of all signatures.
\code{GseaSignaturesVar} contains one element per phenotype (phenotypic
variable). Every one of this elements is of class
\code{gseaSignaturesSign} and contains one element per signature.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gseaSignaturesVar", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}. }
    \item{\code{gseaSignatures}:}{Object of class
      \code{"gseaSignaturesSign"} This object contains the enrichment
      scores and other elements that will be used in the GSEA process. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{getEs}{\code{signature(x = "gseaSignaturesVar")}: Returns the
      enrichment scores.}
    \item{getEsSim}{\code{signature(x = "gseaSignaturesVar")}: Returns
      the simulated enrichemt scores (the ones obtained after
      permutations).} 
    \item{getFcHr}{\code{signature(x = "gseaSignaturesVar")}: Returns
      the fold change and/or the hazard ratio that were used to compute
      the enrichment scores.}
    \item{gseaSignificance}{\code{signature(x = "gseaSignaturesVar")}:
      This is the next step in the process of performing GSEA. This
      funtion will test if the gene sets are enriched. } 
	 }
}
\author{Evarist Planet}
\examples{
showClass("gseaSignaturesVar")
}
\keyword{classes}
