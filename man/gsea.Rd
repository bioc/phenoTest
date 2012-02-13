\name{gsea}
\alias{gsea}
\docType{data}
\title{
GSEA (Gene Set Enrichment Analysis).}
\description{
Computes the enrichment scores and simulated enrichment scores for
each variable and signature.
An important parameter of the function is \code{logScale}. Its default
value is TRUE which means that by default the provided scores (i.e. fold
changes, hazard ratios) will be log scaled. Remember to change this
parameter to FALSE if your scores are already log scaled.
The \code{getEs}, \code{getEsSim}, \code{getFc}, \code{getHr} and
\code{getFcHr} methods can be used to acces each subobject. For more
information please visit the man pages of each method.

It also computes the NES (normalized enrichment score), p values and fdr
(false discovery rate) for all variables and signatures.
For an overview of the output use the \code{summary} method.
}
\usage{
gsea(x,gsets,logScale=TRUE, absVals=FALSE, averageRepeats=FALSE, B=1000,
     mc.cores=1, test="perm",p.adjust.method="none",
     pval.comp.method="original",pval.smooth.tail=TRUE,minGenes=10,maxGenes=500) 
}
\arguments{
\item{x}{ \code{ePhenoTest}, \code{numeric} or \code{matrix} object
  containing hazard ratios or fold changes.}
\item{gsets}{ character or list object containing the names of the
  genes that belong to each signature.}
\item{logScale}{ if values should be log scaled.}
\item{absVals}{ if TRUE fold changes and hazard ratios that are negative
  will be turned into positive before starting the process. This is
  useful when genes can go in both directions.}
\item{averageRepeats}{ if x is of class numeric and has repeated names
  (several measures for some indivdual names) we can average the
  measures of the same names.}
\item{B}{ number of simulations to perform.}
\item{mc.cores}{ number of processors to use.}
\item{test}{ the test that will be used. 'perm' stands for the permutation
  based method, 'wilcox' stands for the wilcoxon test (this is the fastest
  one) and 'ttperm' stands for permutation t test.}
\item{p.adjust.method}{ p adjustment method to be used. Common options
  are 'BH', 'BY', 'bonferroni' or 'none'. All available options and
  their explanations can be found on the \code{p.adjust} function manual.}
\item{pval.comp.method}{ the p value computation method. Has to be one
  of 'signed' or 'original'. The default one is 'original'. See details for
  more information.}
\item{pval.smooth.tail}{ if we want to estimate the tail of the
  ditribution where the pvalues will be generated.}
\item{minGenes}{ gene sets with less than minGenes genes will be removed
  from the analysis.}
\item{maxGenes}{ gene sets with more than maxGenes genes will be removed
  from the analysis.}
}
\details{
The following preprocessing was done on the provided scores (i.e. fold
changes, hazard ratios) to avoid errors during the enrichment score
computation:
-When having two scores with the same name its average was used.
-Zeros were removed.
-Scores without names (which can not be in any signature) removed. 
-Non complete cases (i.e. NAs, NaNs) were removed.
ES score was calculated for each signature and variable (see
references). If parameter \code{test} is 'perm' the signature was 
permutted and the ES score was recalculated (this happened B times for
each variable, 1000 by default).
If \code{test} is 'wilcox' a wilcoxon test in which we test the fact
that the average value of the genes that do belong to our signtaure is
different from the average value of the genes that do not belong to our
signature will be performed.
If \code{test} is 'ttperm' a permutation t-test will be used.
Take into account that the final plot will be different when 'wilcox' is used.

The simulated enrichment scores and the calculated one
are used to find the p value.
 
P value calculation depends on the parameter
\code{pval.comp.method}. The default value is 'original'. In 'original' 
we are simply computing the proportion of anbolute simulated ES which
are larger than the observed absolute ES. In 'signed' we are computing
the proportion of simulated ES which are larger than the observed ES (in
case of having positive enrichment score) and the proportion of
simulated ES which are smaller than the observed ES (in case of having
negative enrichment score). 
}
\references{
Aravind Subramanian, (October 25, 2005) \emph{Gene Set Enrichment
  Analysis}.
\url{www.pnas.org/cgi/doi/10.1073/pnas.0506580102}

C.A. Tsai and J.J. Chen. \emph{Kernel estimation for adjusted p-values
  in multiple testing. Computational Statistics & Data Analysis}
\url{http://econpapers.repec.org/article/eeecsdana/v_3a51_3ay_3a2007_3ai_3a8_3ap_3a3885-3897.htm}
}
\keyword{datasets}
\seealso{
gsea.go, gsea.kegg
}
\examples{
#load epheno object
data(epheno)
epheno

#we construct two signatures
sign1 <- sample(featureNames(epheno))[1:20]
sign2 <- sample(featureNames(epheno))[50:75]
mySignature <- list(sign1,sign2)
names(mySignature) <- c('My first signature','My preferred signature')

#run gsea functions
gseaData <- gsea(x=epheno,gsets=mySignature,B=100,mc.cores=1)
my.summary <- summary(gseaData)
my.summary 
#plot(gseaData)
}
\author{
Evarist Planet
}