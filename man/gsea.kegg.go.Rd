\name{gsea.kegg}
\alias{gsea.kegg}
\alias{gsea.go}
\docType{data}
\title{
Perform Gene Set Enrichment Analysis (GSEA) of Gene Ontologies (GO) and
Kegg gene sets.
}
\description{
The function obtains the GO or Kegg gene sets and perfomrs GSEA analysis
as implemented in the \code{gsea} function.
}
\usage{
gsea.go(x,species='Hs', ontologies='MF', logScale=TRUE, absVals=FALSE,
        averageRepeats=FALSE, B=1000, mc.cores=1, test="perm", 
        p.adjust.method="none", pval.comp.method="original",
        pval.smooth.tail=TRUE,minGenes=10,maxGenes=500) 
gsea.kegg(x,species='Hs', logScale=TRUE, absVals=FALSE,
        averageRepeats=FALSE, B=1000, mc.cores=1, test="perm", 
        p.adjust.method="none", pval.comp.method="original",
        pval.smooth.tail=TRUE,minGenes=10,maxGenes=500) 

}
\arguments{
\item{x}{ \code{ePhenoTest}, \code{numeric} or \code{matrix} object
  containing hazard ratios or fold changes.}
\item{species}{ a single character value specifying the species: "Dm"
  ("Drosophila_ melanogaster"), "Hs" ("Homo_sapiens"), "Rn"
  ("Rattus_norvegicus"), "Mm" ("Mus_musculus") or "Ce"
  ("Caenorhabditis_elegans")).}
\item{ontologies}{ a single character value or a character vector
  specifying an ontology or multiple ontologies. The current version
  provides the following choices: "BP", "CC" and "MF"}
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
This function relies on the following packages: GSEABase, GO.db.

For more information about how the gene sets are obtained see the man
page of the functions \code{getGo} and/or \code{getKegg}.
For more information about the implemented GSEA see the man page of th
function \code{gsea}.
}
\value{
a list of gene sets, with names as GO pathway IDs. Each gene set is a
character vector of Entrez gene identifiers.
}
\author{
Evarist Planet.
}
\seealso{
getGo
}
\examples{
##load libs
#library(KEGG.db)
#library(org.Hs.eg.db)

##get data
#data(eset.genelevel)
#eset.genelevel

##prepare vars2test
#survival <- matrix(c("Relapse","Months2Relapse"),ncol=2,byrow=TRUE)
#colnames(survival) <- c('event','time')
#vars2test <- list(survival=survival,categorical='ER.Status')

##run ExpressionPhenoTest
#epheno <- ExpressionPhenoTest(eset.genelevel,vars2test,p.adjust.method='none')
#epheno

##run gsea with kegg gene sets.
#gseaData <- gsea.kegg(epheno[,1],'Hs')
#summary(gseaData)
#plot(gseaData[[1]],gseaData[[2]],selGsets='hsa04062')
}
\keyword{ ~kwd1 }
