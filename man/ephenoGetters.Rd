\name{getters for the epheno object}
\alias{getFc}
\alias{getHr}
\alias{getMeans}
\alias{getPvals}
\alias{getSummaryDif}
\alias{logFcHr}
\alias{p.adjust.method}
\alias{phenoClass}
\alias{phenoNames}
\docType{data}
\title{
Getters for the epheno object:
}
\description{
\code{getFc} gets the fold changes. \code{getHr} gets the hazard
ratios. \code{getMeans}  gets the means. \code{getPvals} gets the p
values. \code{getSummaryDif} gets fold changes and hazard
ratios. \code{lofFcHr} gets the fold changes and hazard ratios after log
scaling. \code{p.adjust.method} gets the p value adjustment method that
was used when creating the object. \code{phenoClass} returns a
data.frame telling the class (ordinal, continuous, categorical or
survival) of each phenotype. \code{phenoNames} gets the phenotype names.
}
\usage{
getFc(x)
getHr(x)
getMeans(x)
getPvals(x)
getSummaryDif(x)
logFcHr(x)
p.adjust.method(x)
phenoClass(x)
phenoNames(x)
}
\arguments{
\item{x}{ \code{epheno} object}
}
\keyword{datasets}
\author{
Evarist Planet
}