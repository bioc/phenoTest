\name{ClusterPhenoTest}
\alias{ClusterPhenoTest}
\title{
Test association of clusters with phenotype.
}
\description{
Test the associations between clusters that each sample belongs to
(based on gene expression) and each phenotype.
}
\usage{
ClusterPhenoTest(x,cluster,vars2test,B=10^4,p.adjust.method='none')
}
\arguments{
  \item{x}{ExpressionSet with phenotype information stored in \code{pData(x)}.}
  \item{cluster}{variable of class \code{character} or \code{factor}
    telling at which cluster each sample belongs to.}
  \item{vars2test}{ list with components 'continuous', 'categorical',
  'ordinal' and 'survival' indicating which phenotype variables should
  be tested. 'continuous', 'categorical' and 'ordinal' must be character
  vectors, 'survival' a matrix with columns named 'time' and
  'event'. The names must match names in \code{names(pData(x))}.}
  \item{B}{An integer specifying the number of replicates used in the
  chi-square Monte Carlo test (passed on to \code{chisq.test}).}
  \item{p.adjust.method}{Method for P-value adjustment, passed on to
  \code{p.adjust}.}
}
\details{
Test association between the provided clusters and each phenotype.
  
For variables in vars2test\$continuous and vars2test\$ordinal a
Kruskal-Wallis Rank Sum test is used; for vars2test\$categorical a
chi-square test (with exact p-value if \code{simulate.p.value} is set to
TRUE); for var2test\$survival a Cox proportional hazards
likelihood-ratio test. 
}
\examples{
#load data
data(eset)  
eset

#construct vars2test
survival <- matrix(c("Relapse","Months2Relapse"),ncol=2,byrow=TRUE)
colnames(survival) <- c('event','time')
#add positive to have more than one category
pData(eset)[1:20,'lymph.node.status'] <- 'positive' 
vars2test <- list(survival=survival,categorical='lymph.node.status')
vars2test

#first half of the samples will be one cluster and the rest the other cluster
cluster <- c(rep('Cluster1',floor(ncol(eset)/2)),rep('Cluster2',ncol(eset)-floor(ncol(eset)/2)))

#test association
ClusterPhenoTest(eset,cluster,vars2test=vars2test)
}
\author{
David Rossell
}
