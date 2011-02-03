\name{gseaSignificanceVar-class}
\Rdversion{1.1}
\docType{class}
\alias{gseaSignificanceVar-class}

\title{Class "gseaSignificanceVar"}
\description{This object contains the results of the test of enrichment
that was performed on each gene set and phenotype.
There is one container for every phenotype. Every one of this containers
(of class \code{gseaSignificanceSign}) is a container itsef and has the
results of the tests for all signatures.
\code{GseaSignificanceVar} contains one element per phenotype (phenotypic
variable). Every one of this elements is of class
\code{gseaSignificanceSign} and contains one element per signature.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gseaSignificanceVar", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} .}
    \item{\code{gseaSignificance}:}{Object of class
      \code{"gseaSignificanceSign"} This object contains the results of
      the tests. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
No methods defined with class "gseaSignificanceVar" in the signature.
}
\author{Evarist Planet}
\examples{
showClass("gseaSignificanceVar")
}
\keyword{classes}
