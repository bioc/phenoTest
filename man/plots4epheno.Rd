\name{barplotSignatures}
\alias{barplotSignatures}
\alias{barplotSignifSignatures}
\title{Summary plots for gene signature vs phenotype association}
\description{
Summarizes the univariate relationships between genes in one or more
signatures and several phenotype variables, as summarized in
\code{epheno} objects (which can be created with the
\code{ExpressionPhenoTest} function).

By default \code{barplotSignifSignatures} performs a binomial test
(\code{binom.test} from package \code{stats}) for each signature
to see if the number of up up regulated and down regulated genes is
different enough to be statistically different.
When a reference gene set is provided we test if the proportions of up
and down regulated genes of each gene set is different from the
proportions in the reference gene set. This has been done with a
chi-square test.
When a reference gene set is provided and parameter
\code{testUpDown} is \code{TRUE} (by default its \code{FALSE})
the number of genes corresponding to up and down regulated are
compared with those of the reference gene set separately.
}
\usage{
barplotSignatures(x, signatures, referenceSignature, alpha=.05,
p.adjust.method='none', ylab, cex.text=1, ...)
barplotSignifSignatures(x, signatures, referenceSignature, testUpDown=FALSE,
simulate.p.value = FALSE, B = 10^4, p.adjust.method='none', alpha=.05,
ylab, ylim=ylim, cex.text=1, ...)
}
\arguments{
\item{x}{ \code{epheno} object, as returned by
  \code{ExpressionPhenoTest}.} 
\item{signatures}{List with each element corresponding to a
  signature. The gene names in each signature must match those in
  \code{epheno}.}
\item{referenceSignature}{If specified, the average fold change in each
  signature is compared to the average fold change in the signature
  \code{referenceSignature}.}
\item{testUpDown}{If set to \code{TRUE}, bars corresponding to
  up and down-regulated genes are compared with those of
  \code{referenceSignature} separately. This argument is ignored if
  \code{referenceSignature} is not specified.}
\item{cex.text}{Character expansion for the text indicating the
  P-values. Ignored if referenceSignature is missing.}
\item{alpha}{Confidence levels for barplot error bars.}
\item{p.adjust.method}{P-value adjustment method, passed on to
  \code{p.adjust}.}
\item{simulate.p.value}{A logical indicating whether chi-square p-values
  should be computed by Monte Carlo simulation (passed on to \code{chisq.test}).}
\item{B}{Integer specifying the number of replicates in the Monte Carlo
  simulation (passed on to \code{chisq.test}).}
\item{ylab}{y-axis labels}
\item{ylim}{y-axis limits}
\item{...}{Other arguments to be passed on to boxplot.}
}
\value{
  When a single signature is provided as input, a single plot assessing the
  association of that signature with all phenotype variables is created.
  If several signatures are provided, one separate plot is created for
  each phenotype variable.
}
\examples{
#create epheno
data(epheno)

#construct two signatures
sign1 <- sample(featureNames(epheno))[1:20]
sign2 <- sample(featureNames(epheno))[1:15]
mySignature <- list(sign1,sign2)
names(mySignature) <- c('My first signature','My preferred signature')

#plot
barplotSignifSignatures(epheno[,'Relapse'],mySignature,alpha=0.05)
}
\author{
Evarist Planet
}
\keyword{graph}
