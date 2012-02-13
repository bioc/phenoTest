\name{eset2genelevel}
\alias{eset2genelevel}
\title{
Filter ExpressionSet to keep one probeset per gene.
}
\description{
Only one probeset per gene will be kept and entrezid will be used as
gene identifier.
\code{nsFilter} from package \code{genefilter} is used to select the
probeset.
The selected probeset is the one with higher interquartilic range.
}
\usage{
eset2genelevel(x)
}
\arguments{
  \item{x}{an object of class \code{ExpressionSet}.}
  }
\author{
Evarist Planet
}
\seealso{
genefilter::nsFilter
}
\examples{
#data(eset)
#library(hgu133a.db)
#x <- eset2genelevel(eset)
#x
#head(featureNames(x))
}
\keyword{ ~kwd1 }

