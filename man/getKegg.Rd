\name{getKegg}
\alias{getKegg}
\title{
Create a list of gene sets based on KEGG pathways terms.
}
\description{
This function creates a list of gene sets based on KEGG pathways
terms. It is species-specific, and returns a list of gene sets, each of
which is a character vector of Entrez gene identifiers.

This function is a wraper to the function \code{KeggGeneSets} from
package \code{HTSanalyzeR}.
}
\usage{
getKegg(species = "Dm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species}{ a single character value specifying the species: "Dm"
    ("Drosophila_ melanogaster"), "Hs" ("Homo_sapiens"), "Rn"
    ("Rattus_norvegicus"), "Mm" ("Mus_musculus") or "Ce"
    ("Caenorhabditis_elegans")).}
}
\details{
This function relies on the following packages: GSEABase, KEGG.db.
}
\value{
a list of gene sets, with names as KEGG pathway IDs. Each gene set is a
character vector of Entrez gene identifiers.
}
\author{
Evarist Planet.
}
\seealso{
getGo
}
\examples{
#library(KEGG.db)
#library(org.Hs.eg.db)
#kegg.Hs <- getKegg('Hs')
#str(kegg.Hs)
#kegg.Hs[1:2]
}
\keyword{ ~kwd1 }
