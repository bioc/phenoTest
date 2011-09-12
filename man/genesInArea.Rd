\name{genesInArea}
\alias{genesInArea}
\title{
Find genes that are in given areas.
}
\description{
Combine the output of getEsPositions and findCopyNumber to see which
genes are in the enriched areas.

Given areas of enrichment (obtained with \code{findCopyNumber}) and a
set of genes or probes and their positions in the genome (obtained with
\code{getEsPositions}) the function tells which genes fall in each
area.
}
\usage{
genesInArea(x, regions)
}
\arguments{
  \item{x}{
An object of class \code{data.frame} with gene or probe identifiers as
row names and the following columns: es (the enrichment score), chr (the
chromosome where the gene or probe belong to) and pos (position in the
chromosome in megabases).
It can be obtained with the function \code{getEsPositions}.
}
  \item{regions}{
This is usually the output of \code{findCopyNumber} function.
}
}
\author{
Evarist Planet
}
\seealso{
getEsPositions, findCopyNumber
}
\examples{
data(epheno)
phenoNames(epheno)
mypos <- getEsPositions(epheno,'Relapse')
head(mypos)
#regions <- findCopyNumber(mypos)
#head(regions)
#genes <- genesInArea(mypos,regions)
#head(genes)
}
\keyword{ ~kwd1 }
