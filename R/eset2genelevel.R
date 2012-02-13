eset2genelevel <- function(x) {
  stopifnot(class(x)=='ExpressionSet')
  annolib <- paste(annotation(x),'.db',sep='')
  if (!(annolib %in% loadedNamespaces())) stop(paste(annolib,'has not been loaded!'))
  #
  ans <- nsFilter(x,require.entrez=FALSE,remove.dupEntrez=TRUE,var.filter=FALSE)$eset
  text <- "as.character(unlist(AnnotationDbi::mget(featureNames(ans),eval(parse(text=paste(annotation(ans),'ENTREZID',sep=''))))))"
  featureNames(ans) <- eval(parse(text=text))
  text <- paste(annotation(ans),'ORGPKG',sep='')
  annotation(ans) <- eval(parse(text=text))
  ans
}
