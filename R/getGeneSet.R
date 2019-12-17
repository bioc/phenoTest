libLoaded <- function(lib) {
  loaded <- eval(parse(text=paste('"',lib,'" %in% loadedNamespaces()',sep='')))
  if (!loaded) stop(paste(lib,'has not been loaded.'))
}

