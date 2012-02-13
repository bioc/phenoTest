\name{gsea2html}
\alias{gsea2html}
\title{
Export an object of class gseaData to an html file.
}
\description{
Exports \code{gseaData} objects to html files with plots and links to online databases.
}
\usage{
gsea2html(gseaData, epheno, variable, title = "", path, file, digits = 3, plotEs = FALSE, limit=100)
}
\arguments{
  \item{gseaData}{ an object of class \code{gseaData}.}
  \item{epheno}{ the object of class \code{epheno} that was used to
    create \code{gseaData}.}
  \item{variable}{ varible that we are interested in.}
  \item{title}{ title that will be shown on top of the table.}
  \item{path}{ directory where we want to store the html files.}
  \item{file}{ filename.}
  \item{digits}{ Number of decimal digits that will be shown on the
    table.} 
  \item{plotEs}{ if this is TRUE enrichment score plots will be plotted
    instead of normalized enrichment score plots.}
  \item{limit}{ maximum number of gene sets that will be exported.}
}
\details{
This function produces a browseable version of the table that we can
obtain with \code{summary(gseaData)}.
We will obtain one plot per NES (or ES) and we will be able to see which
genes belong to each gene set and the values they have in the
\code{epheno} object.
}
\author{
Evarist Planet
}
\examples{
#WITH PROBESET AS IDENTIFIER
data(eset)
data(epheno)

set.seed(777)
sign1 <- sample(featureNames(eset))[1:20]
sign2 <- sample(featureNames(eset))[1:50]
mySignature <- list(sign1,sign2)
names(mySignature) <- c('My first signature','Another signature')
mySignature

mygsea <- gsea(x=epheno[,1],gsets=mySignature,B=100,p.adjust='BH')
summary(mygsea)

#following line has been commented to prevent the creation of files
#gsea2html(gseaData=mygsea,epheno=epheno,variable=phenoNames(epheno)[1],title='My test',path='~/Desktop',file='myGSEA.html')

#WITH ENTREZID AS IDENTIFIER
data(eset.genelevel)
eset.genelevel

set.seed(777)
sign1 <- sample(featureNames(eset.genelevel))[1:20]
sign2 <- sample(featureNames(eset.genelevel))[1:50]
mySignature.genelevel <- list(sign1,sign2)
names(mySignature.genelevel) <- c('My first signature','Another signature')
mySignature.genelevel

epheno.genelevel <- ExpressionPhenoTest(eset.genelevel,vars2test=list(categorical='lymph.node.status'))
mygsea.genelevel <- gsea(x=epheno.genelevel,gsets=mySignature.genelevel,B=100,p.adjust='BH')
summary(mygsea.genelevel)

#following line has been commented to prevent the creation of files
#gsea2html(gseaData=mygsea.genelevel,epheno=epheno.genelevel,variable=phenoNames(epheno.genelevel),title='My test (at genelevel)',path='~/Desktop',file='myGSEA_genelevel.html')
}
\keyword{ ~kwd1 }

