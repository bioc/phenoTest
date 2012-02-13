libLoaded <- function(lib) {
  loaded <- eval(parse(text=paste('"',lib,'" %in% loadedNamespaces()',sep='')))
  if (!loaded) stop(paste(lib,'has not been loaded.'))
}

getGo <- function(species = "Dm", ontologies = "MF") {
  libLoaded('GO.db')
  lib <- paste('org.',species,'.eg.db',sep='')
  libLoaded(lib)
  HTSanalyzeR::GOGeneSets(species, ontologies)
}

getKegg <- function(species = "Dm") {
  libLoaded('KEGG.db')
  lib <- paste('org.',species,'.eg.db',sep='')
  libLoaded(lib)
  HTSanalyzeR::KeggGeneSets(species)
}
