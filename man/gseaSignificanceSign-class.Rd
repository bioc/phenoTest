\name{gseaSignificanceSign-class}
\Rdversion{1.1}
\docType{class}
\alias{gseaSignificanceSign-class}

\title{Class "gseaSignificanceSign"}
\description{This object contains the results of the test of enrichment
that was performed on each gene set.
There is one container for every gene signature.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gseaSignificanceSign", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} .}
    \item{\code{gseaSignificance}:}{Object of class \code{"matrix"}
      Contains the statistics. Use the \code{summary} method to access
      this information.}
    \item{\code{p.adjust.method}:}{Object of class \code{"character"}
      The p-value adjustment method that was used.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
No methods defined with class "gseaSignificanceSign" in the signature.
}
\author{Evarist Planet}
\examples{
showClass("gseaSignificanceSign")
}
\keyword{classes}
