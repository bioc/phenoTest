\name{smoothCoxph}
\alias{smoothCoxph}
\title{
Plots the Cox proportional hazard smoothed by gene expression level.
}
\description{
Builds a plot showing how hazard behaves over different levels of
expression of a given gene. Confidence intervals are also provided. 
}
\usage{
smoothCoxph(time, event, x, xlim, ylim, xlab, ylab, logrisk=TRUE, ...)
}
\arguments{
  \item{time}{variable where time to survival is stored.}
  \item{event}{variable where survival event is stored.}
  \item{x}{\code{numeric} containing the expression levels of a given gene.}
  \item{xlim}{\code{xlim} for the plot.}
  \item{ylim}{\code{ylim} for the plot.}
  \item{xlab}{\code{xlab} for the plot.}
  \item{ylab}{\code{ylab} for the plot.}  
  \item{logrisk}{\code{logrisk} if we want to compute risk or logrisk
    estimates. By default this is TRUE, which has a better behaviour
    under small sample sizes.}  
  \item{\dots}{other arguments that will be passed to plot.}
}
\author{
David Rossell.
}
\examples{
#load eset
data(eset)

#make plot
smoothCoxph(pData(eset)$Months2Relapse,pData(eset)$Relapse,exprs(eset)[25,])
}
\keyword{datasets}
