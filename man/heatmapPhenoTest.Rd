\name{heatmapPhenoTest}
\alias{heatmapPhenoTest}
\title{
Produce heatmap from phenotype data.
}
\description{
Show the associations between clusters that each sample belongs to and
each phenotype in a heatmap and/or a Kaplan-Meier plot.
}
\usage{
heatmapPhenoTest(x, signatures, vars2test, probes2genes = FALSE,
filterVar, filteralpha = 0.05, distCol = "pearson", nClust = 2, distRow
= "cor", p.adjust.method = "none", simulate.p.value = FALSE, B = 10^5,
linkage = "average", equalize = FALSE, center = TRUE, col, survCol,
heat.kaplan="both", ...)
}
\arguments{
  \item{x}{ExpressionSet with phenotype information stored in \code{pData(x)}.}
  \item{signatures}{Either character vector or list of character vectors
  with gene sets to be used to draw heatmaps (gene names should match
  those in \code{featureNames(x)}). A separate heatmap will be produced
  for each element in the list.}
  \item{vars2test}{ list with components 'continuous', 'categorical',
  'ordinal' and 'survival' indicating which phenotype variables should
  be tested. 'continuous', 'categorical' and 'ordinal' must be character
  vectors, 'survival' a matrix with columns named 'time' and
  'event'. The names must match names in \code{names(pData(x))}.}
  \item{probes2genes}{If set to \code{TRUE} a single probe is selected for
  each gene. \code{nsFilter} is used to select the probe with highest
  inter-quartile range.}
  \item{filterVar}{If specified, only genes with significant differences
  in the variable \code{filterVar} will be displayed in the
  heatmap. Note that this option will not affec the sample clustering,
  as this is obtained using both significant and non-significant genes.}
  \item{filteralpha}{Significance level for the filtering based on \code{filterVar}.}
  \item{distCol}{Distance metric used to cluster columns
  (e.g. patients/samples). Can take any value accepted by
  \code{dist}. Pearson and Spearman correlations are also 
  allowed. Write 'spearman' or 'pearson' to use them.}
  \item{nClust}{Number of desired clusters.}
  \item{distRow}{Distance metric used to cluster rows (e.g. genes). Can
  take any value accepted by \code{distancematrix}.}
  \item{p.adjust.method}{Method for P-value adjustment, passed on to
  \code{p.adjust}.}
  \item{simulate.p.value}{If set to FALSE the chi-square test p-values are
  computed using asymptotics, otherwise a simulation is used (see
  \code{chisq.test} for details).}
  \item{B}{An integer specifying the number of replicates used in the
  chi-square Monte Carlo test (passed on to \code{chisq.test}).}
  \item{linkage}{Linkage used for clustering. Must be either 'complete',
  'average' or 'minimum'.}
  \item{equalize}{Should color codes be equalized between genes, i.e. all
  genes present the same range of colors. Passed on to
  \code{heatmap_plus}.}
  \item{center}{centering is done by subtracting the column means
  (omitting NAs).}
  \item{col}{Color scheme to be used for heatmap. Defaults to a green/red
  scheme designed to look nice for microarray data.}
  \item{survCol}{Colors for the Kaplan-Meier survival curves.}
  \item{heat.kaplan}{can be "heat" if we want to plot a heatmap,
    "kaplan" if we want to plot a kaplan-meier or "both" if we want both
    of them.}
  \item{\dots}{Other arguments for the survival plot, e.g. lty etc.}
}
\details{
Makes two clusters of samples based on the expression levels of the
genes from the given signature and plots a heatmap and/or a Kaplan-Meier
showing the association  between belonging to one cluster or the other
and each phenotype.
  
For variables in vars2test\$continuous and vars2test\$ordinal a Kruskal-Wallis Rank Sum test is
used; for vars2test\$categorical a chi-square test (with exact p-value
if \code{simulate.p.value} is set to TRUE); for var2test\$survival a Cox proportional hazards likelihood-ratio test.
}
\examples{
#load data
data(eset)  
eset

#construct vars2test
survival <- matrix(c("Relapse","Months2Relapse"),ncol=2,byrow=TRUE)
colnames(survival) <- c('event','time')
vars2test <- list(survival=survival)
vars2test

#construct a signature
sign <- sample(featureNames(eset))[1:20]

#make plot
heatmapPhenoTest(eset,sign,vars2test=vars2test,heat.kaplan='heat')
heatmapPhenoTest(eset,sign,vars2test=vars2test,heat.kaplan='kaplan')
}
\author{
David Rossell
}
