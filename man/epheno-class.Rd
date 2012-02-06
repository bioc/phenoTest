\name{epheno-class}
\Rdversion{1.1}
\docType{class}
\alias{epheno-class}
\alias{[,epheno,ANY,ANY-method}
\alias{dim,epheno-method}
\alias{getFc,epheno-method}
\alias{getHr,epheno-method}
\alias{getMeans,epheno-method}
\alias{getPvals,epheno-method}
\alias{getPostProbs,epheno-method}
\alias{getSignif,epheno-method}
\alias{getSummaryDif,epheno-method}
\alias{logFcHr,epheno-method}
\alias{p.adjust.method,epheno-method}
\alias{phenoClass,epheno-method}
\alias{phenoNames,epheno-method}
\alias{approach,epheno-method}
\alias{show,epheno-method}

\title{Class "epheno"}
\description{
Object obtained with the ExpressionPhenoTest function. Contains FC, HR
and pvals from testing expression values of each gene against phenotypic
variables.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("epheno",
  assayData, phenoData, featureData, exprs, ...)}. 
}
\section{Slots}{
  \describe{
    \item{\code{p.adjust.method}:}{Object of class \code{"character"}
      containing the multiple testing adjustment method used (if one was
      used). }
    \item{\code{approach}:}{Object of class \code{"character"}
      containing 'frequentist' or 'bayesian' depending on the user's
      selecion.} 
    \item{\code{assayData}:}{Object of class \code{"AssayData"} that is
      inherited from the \code{ExpressionSet} object used to create the
      \code{epheno} object. }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}
      that contains information about the variables stored in the
      experimentData slot such as their class (continuous, categorical,
      etc) or type (mean, summaryDif, pval, etc).}
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} that is
      inherited from the \code{ExpressionSet} object used to create the
      \code{epheno} object. }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} that is
      inherited from the \code{ExpressionSet} object used to create the
      \code{epheno} object. }
    \item{\code{annotation}:}{Object of class \code{"character"} that is
      inherited from the \code{ExpressionSet} object used to create the
      \code{epheno} object. }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} that is
      inherited from the \code{ExpressionSet} object used to create the
      \code{epheno} object. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} that is
      inherited from the \code{ExpressionSet} object used to create the
      \code{epheno} object. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{ExpressionSet}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "ExpressionSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "ExpressionSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "ExpressionSet", distance 4.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "epheno", i = "ANY", j = "ANY")}:
      inherited from the \code{ExpressionSet} class. }
    \item{dim}{\code{signature(x = "epheno")}: inherited from the
      \code{ExpressionSet} class.}
    \item{export2CSV}{\code{signature(x = "epheno")}: ... }
    \item{getFc}{\code{signature(x = "epheno")}: getter for the fold
      changes.}
    \item{getHr}{\code{signature(x = "epheno")}: getter for the hazard
      ratios.}
    \item{getMeans}{\code{signature(x = "epheno")}: getter for the means.}
    \item{getSignif}{\code{signature(x = "epheno")}: getter for the
      pvalues or posterior probabilities.}
    \item{getPvals}{\code{signature(x = "epheno")}: getter for the
      pvalues.}
    \item{getPostProbs}{\code{signature(x = "epheno")}: getter for the
      posterior probabilities.}
    \item{getSummaryDif}{\code{signature(x = "epheno")}: getter that
      returns hazard ratios, fold changes and pvalues. }
    \item{gseaSignatures}{\code{signature(x = "epheno", signatures =
	"list")}: Used to compute GSEA. Please read the gseaSignatures
      manual.}
    \item{logFcHr}{\code{signature(x = "epheno")}: getter for the log of
      fold changes and hazard ratios.}
    \item{p.adjust.method}{\code{signature(x = "epheno")}: getter for
      the p value adjustment method that has been used. }
    \item{phenoClass}{\code{signature(x = "epheno")}: Returns the class
      off all variables. }
    \item{phenoNames}{\code{signature(x = "epheno")}: Returns the names
      of the tested phenotypes.}
    \item{show}{\code{signature(object = "epheno")}: Shows a brief
      overview of the object.}
	 }
}
\author{
Evarist Planet
}
\examples{
showClass("epheno")
}
\keyword{classes}
