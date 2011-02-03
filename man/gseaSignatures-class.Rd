\name{gseaSignatures-class}
\Rdversion{1.1}
\docType{class}
\alias{gseaSignatures-class}
\title{Class "gseaSignatures" ES and EsSim container.}
\description{This object contains de ES (enrichment scores) and
  simulated ES that will be used in the GSEA (Gene Set Enrichment
  Analysis) process.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gseaSignatures", ...)}. 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} .}
    \item{\code{es}:}{Object of class \code{"numeric"} Contains the
      observed enrichment scores. The ones that were compted from the
      data without permuting anything.  }
    \item{\code{es.sim}:}{Object of class \code{"numeric"} Contains the
      enrichment score that were obtained after permutations. }
    \item{\code{signature}:}{Object of class \code{"numeric"} The subset
      of genes we are interested in. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
No methods defined with class "gseaSignatures" in the signature.
}
\author{Evarist Planet}
\examples{
showClass("gseaSignatures")
}
\keyword{classes}
