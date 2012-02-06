\name{epheno2html}
\alias{epheno2html}
\title{
Create html files and plots from an epheno object.
}
\description{
Creates html files and plots using an epheno object, which stores the
association between a list of variables and gene expression.
}
\usage{
epheno2html(x, epheno, outputdir, prefix = "", genelimit = 50, categories = 3, withPlots = TRUE, mc.cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{ExpressionSet} (used to generate
    the epheno object) containing expression levels in exprs(x),
    phenotype information in pData(x) and annotation in annotation(x).}
  \item{epheno}{an object produced by ExpressionPhenoTest. this object
    will contain univariate association between a list of phenotype
    variables and gene expression as weel as p-values.}
  \item{outputdir}{where to place files.}
  \item{prefix}{will be used to add a text to the beginning of the files
    that will be created.}
  \item{genelimit}{maximum number of genes on the list.}
  \item{categories}{Number of categories used for continuous
    variables. It has to be the same as the one used for
    ExpressionPhenoTest.}
  \item{withPlots}{when FALSE no plots will be produced. Makes the
    process faster.}
  \item{mc.cores}{number of cores that will be used to run the process.}
}
\author{
Evarist Planet
}
\examples{
#Example on building homology tables for human.
#mart <- useMart("ensembl", "hsapiens_gene_ensembl")
#homol.symbol <- getLDS(attributes = c("entrezgene"), 
#    mart = mart, attributesL = c("external_gene_id"), 
#    martL = mart, filters = "entrezgene", values = entrezid)
#mart <- useMart("ensembl", "hsapiens_gene_ensembl")
#homol.genename <- getLDS(attributes = c("entrezgene"), 
#    mart = mart, attributesL = c("description"), martL = mart, 
#    filters = "entrezgene", values = entrezid)
}