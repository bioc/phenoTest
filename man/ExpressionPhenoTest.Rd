\name{ExpressionPhenoTest}
\alias{ExpressionPhenoTest}
\docType{data}
\title{Tests univariate association between a list of phenotype
  variables and gene expression.} 
\description{
Tests univariate association between a list of phenotype variables and
gene expression. 
}
\usage{
ExpressionPhenoTest(x, vars2test, adjustVars,
p.adjust.method='BH',continuousCategories=3,mc.cores,approach='frequentist') 
}
\arguments{
\item{x}{ ExpressionSet containing expression levels in \code{exprs(x)}
  and phenotype information in \code{pData(x)}.} 
\item{vars2test}{ list with components 'continuous', 'categorical',
  'ordinal' and 'survival' indicating which phenotype variables should 
  be tested. 'continuous', 'categorical' and 'ordinal' must be character
  vectors, 'survival' a matrix with columns named 'time' and
  'event'. The names must match names in names(pData(x)).} 
\item{adjustVars}{ variables that will be used as adjustment variables
  when fitting linear models and/or cox models. This variables have to
  exist in \code{colnames(pData(x))}.} 
\item{p.adjust.method}{ method for p-value adjustment, passed on to
  \code{p.adjust}. Valid values are c("holm", "hochberg", "hommel",
  "bonferroni", "BH", "BY", "fdr", "none").} 
\item{continuousCategories}{ number of categories used for continuous
  variables.} 
\item{mc.cores}{ the number of cores to use, i.e. how many processes
  will be spawned (at most).}
\item{approach}{ this can be either 'frequentist' or 'bayesian'. With
  frequentist pvalues will be computed. With 'bayesian' posterior
  probabilities will be computed.}
}
\details{
If approach is 'frequentist':
-The effect of both continuous, categorical and ordinal phenotype
variables on gene expression levels are tested via lmFit.
-For ordinal variables a single coefficient is used to test its effect
on gene expression (trend test), which is then used to obtain a P-value
(means for each category are reported in the output).
-Gene expression effects on survival are tested via Cox proportional
hazards model, as implemented in function 'coxph'.

If approach is bayesian posterior probabilities are computed comparing
the BIC of a model with the variable of interest as explanatory variable
against the BIC of the same model without the variable of interest as
explanatory variable.
}
\value{
The output is an \code{epheno} object, which basically extends an
\code{ExpressionSet} object. The means, fold changes, standarized hazard
ratios and pvalues are stored in the \code{experimentData} slot which is
accessible with the \code{exprs} method. Information about the kind of
information of each variable can be found in the \code{phenoData} slot
which is accessible with the \code{pData} method.

There are several methods that can be used to access the information
stored in an \code{epheno} object. For more information please type one
of the following: \code{getFc(x), getHr(x), getMeans(x), getSignif,
getPvals(x), getPostProbs, getSummaryDif(x), logFcHr(x),
p.adjust.method(x), phenoClass(x), phenoNames(x)}.
}
\references{
Kass R.E. and Wasserman L. A Reference Bayesian Test for Nested
Hypotheses and its Relationship to the Schwarz Criterion. Journal of the
American Statistical Association, 90, pp. 928-934.
}
\keyword{datasets}
\examples{
#load eset
data(eset)
eset

#prepare vars2test
survival <- matrix(c("Relapse","Months2Relapse"),ncol=2,byrow=TRUE)
colnames(survival) <- c('event','time')
vars2test <- list(survival=survival)

#run ExpressionPhenoTest
epheno <- ExpressionPhenoTest(eset,vars2test,p.adjust.method='none')
epheno
}
\author{
David Rossell
}
