\name{getGo}
\alias{getGo}
\title{
Create a list of gene sets based on GO pathways terms.
}
\description{
This function creates a list of gene sets based on GO pathways
terms. It is species-specific, and returns a list of gene sets, each of
which is a character vector of Entrez gene identifiers.

This function is a wraper to the function \code{GOGeneSets} from
}
\usage{
getGo(species = "Dm", ontologies = "MF")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species}{ a single character value specifying the species: "Dm"
    ("Drosophila_ melanogaster"), "Hs" ("Homo_sapiens"), "Rn"
    ("Rattus_norvegicus"), "Mm" ("Mus_musculus") or "Ce"
    ("Caenorhabditis_elegans")).}
  \item{ontologies}{ a single character value or a character vector
    specifying an ontology or multiple ontologies. The current version
    provides the following choices: "BP", "CC" and "MF"}
}
\details{
This function relies on the following packages: GSEABase, GO.db.
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
#library(GO.db)
#library(org.Hs.eg.db)
#go.Hs <- getGo('Hs')
#str(go.Hs)
#go.Hs[1:2]
}
\keyword{ ~kwd1 }
