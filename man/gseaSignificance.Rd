\name{gseaSignificance}
\alias{gseaSignificance}
\docType{data}
\title{
ES' (enrichment scores) sifgnificance.
}
\description{
This function performs the second step in the process of obtaining a
GSEA-like plot.
It computes the NES (normalized enrichment score), p values and fdr
(false discovery rate) for all variables and signatures. A
\code{gseaSignaturesSign} or \code{gseaSignaturesVar} object will 
be needed as input (these objects can be obtained with the
\code{gseaSignatures} function).
For an overview of the output use the \code{summary} method.
The next step after using the \code{gseaSignificance} function would be using
the \code{plot} method.
}
\usage{
gseaSignificance(x,p.adjust.method='none',pval.comp.method='original',pval.smooth.tail=TRUE)
}
\arguments{
\item{x}{ \code{gseaSignaturesSign} or \code{gseaSignaturesVar} object obtained
  with the \code{gseaSignatures} method. This object contains the enrichment scores
  ,the simulated enrichment scores and the fold changes or hazard ratios.}
\item{p.adjust.method}{ p adjustment method to be used. Common options
  are 'BH', 'BY', 'bonferroni' or 'none'. All available options and
  their explanations can be found on the \code{p.adjust} function manual.}
\item{pval.comp.method}{ the p value computation method. Has to be one
  of 'signed' or 'original'. The default one is 'original'. See details for
  more information.}
\item{pval.smooth.tail}{ if we want to estimate the tail of the
  ditribution where the pvalues will be generated.}
}
\details{
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
\examples{
#for examples see the help file of gseaSigntaures: ?gseaSignatures
}
\author{
Evarist Planet
}
