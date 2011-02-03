\name{pAdjust}
\alias{pAdjust}
\title{
Adjust p values of an epheno object.
}
\description{
Adjusts the p values of an epheno object. The \code{p.adjust} function will be
used. For more information read the \code{p.adjust} function's help.
}
\usage{
pAdjust(x, method = "BH")
}
\arguments{
  \item{x}{an \code{epheno} object. }
  \item{method}{the correction method that will be used. Se the
    \code{p.adjust} help for more info about the methods.}
}
\examples{
#load epheno object
data(epheno)
epheno

#Adjust pvalue
p.adjust.method(epheno)
epheno <- pAdjust(epheno,method='BH')
p.adjust.method(epheno)
}
\author{
Evarist Planet
}
