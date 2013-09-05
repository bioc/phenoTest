\name{pca}
\alias{pca}
\title{
Principal components plot.
}
\description{
Creates a Principal Components plot where we can show paired samples,
and confidence intervals for the mean of every group of interest.
We can also choose the component or components we want to plot.
}
\usage{
pca(x, group, group2, pair, names, ellipse = FALSE, main = "", components = c(1, 2))
}
\arguments{
  \item{x}{
    An object of class ExpressionSet.
}
  \item{group}{
    Variable in pData(x) that contains the groups of interest. Samples
    of the same group will be plotted with the same color.
}
  \item{group2}{
    Variable in pData(x) that contains secondary groups of
    interest. Sample of the same secondary group of interest will be
    plotted with the same symbol.
}
  \item{pair}{
    Variable in pData(x) that contains the information about the pairs
    of data. Those pairs will be joined by a line.
}
  \item{names}{
    Variable in pData(x) that contains the information about the names
    of the samples.
}
  \item{ellipse}{
    If we want to plot ellipses with the 95 percent confidence intervals
    for every group.
}
  \item{main}{
    A title for the plot.
}
  \item{components}{
    Which components we want to plot. By default the first principal
    component will be plotted on the x axis and the second principal
    component will be plotted on the y axis.
    More than two components may be specified. If so multiple plots will
    be produced.
}
}
\author{
  Evarist Planet
}
\seealso{
  prcomp.
}
\examples{
data(eset)
pca(x=eset, group='Relapse', names='GEOaccession')
pca(x=eset, group='Relapse', names='GEOaccession', components=1:3)
}
\keyword{datasets}
    