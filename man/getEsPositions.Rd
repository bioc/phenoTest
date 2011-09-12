\name{getEsPositions}
\alias{getEsPositions}
\title{
Obtain chromosome positions for each gene.
}
\description{
Given an object of class epheno obtain the gene positions on the genome.
}
\usage{
getEsPositions(epheno, phenoName, organism = "human", logEs = T, center = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epheno}{
An object of class \code{epheno} usually obtained with \code{ExpressionPhenoTest}
}
  \item{phenoName}{
The phenotype that we want to use. Has to be in \code{phenoNames(epheno)}
}
  \item{organism}{
Has to be 'human' or 'mouse'. The default is 'human'.
}
  \item{logEs}{
If the values have to be log scaled.
}
  \item{center}{
If the values have to be genome centered. If TRUE the genome average
will be substracted to every value.
}
}
\details{
The output will usually be passed to findCopyNumber.
}
\value{
An object of class \code{data.frame} will be returned containing 3
variables: es (enichment score for fold change or hazard ratio), chr
(chromosome), pos (position in Mb).
epheno's featureNames will be used as row names. 
}
\author{
Evarist Planet
}
\examples{
data(epheno)
phenoNames(epheno)
mypos <- getEsPositions(epheno,'Relapse')
head(mypos)
}
\keyword{ ~kwd1 }
