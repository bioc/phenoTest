gsea2html <- function(gseaData,epheno,variable,title='',path,file,digits=3,plotEs=FALSE,limit=100) {
  #control errors
  stopifnot(class(gseaData)=='gseaData')
  stopifnot(class(epheno)=='epheno')
  stopifnot(variable %in% phenoNames(epheno))
  stopifnot(!missing(file))
  stopifnot(file.exists(path))

  #my label
  mylabel <- gsub('.html','',file)
  
  #preprocess
  id.entrezid <- 'org.' %in% substr(annotation(epheno),0,4)
  epheno <- epheno[,variable]
  sigs <- gseaData[[1]][[1]]$signatures

  gsets <- summary(gseaData)
  gsets <- data.frame(gsets[1:2],plot='',gsets[3:ncol(gsets)])

  genes <- exprs(epheno)
  sel <- pData(epheno)$phenoType=='signif'
  colnames(genes)[sel] <- paste(colnames(genes)[sel],ifelse(approach(epheno)=='frequentist','pvalue','postProb'),sep='.')

  fnames <- featureNames(epheno)
  fnames <- unique(fnames[fnames %in% unlist(sigs)])

  #gene annotation
  if (id.entrezid) {
    entrezid <- featureNames(epheno)
    orglib <- paste(annotation(epheno),'.db',sep='')
    if (!(orglib %in% loadedNamespaces())) stop('multicore library has not been loaded!')
    envir <- eval(parse(text=paste(annotation(epheno),'SYMBOL',sep='')))
    symbol <- unlist(AnnotationDbi::mget(featureNames(epheno),envir,ifnotfound=NA))
    envir <- eval(parse(text=paste(annotation(epheno),'GENENAME',sep='')))
    genename <- unlist(AnnotationDbi::mget(featureNames(epheno),envir,ifnotfound=NA))
  } else {
    entrezid <- unlist(AnnotationDbi::mget(fnames,AnnotationDbi::get(paste(annotation(epheno),"ENTREZID", sep = "")),ifnotfound=NA))
    symbol <- unlist(AnnotationDbi::mget(fnames,AnnotationDbi::get(paste(annotation(epheno),"SYMBOL", sep = "")),ifnotfound=NA))
    genename <- unlist(AnnotationDbi::mget(fnames,AnnotationDbi::get(paste(annotation(epheno),"GENENAME", sep = "")),ifnotfound=NA))
  }
  anno <- data.frame(entrezid=as.character(entrezid),symbol=as.character(symbol),genename=as.character(genename),stringsAsFactors=F)
  
  #filter gsets
  gsets <- gsets[order(gsets$fdr),]
  if (limit<nrow(gsets)) gsets <- gsets[1:limit,]
  
  #gene set annotation
  if (gseaData$gsetOrigin=='KEGG') {
    envir <- eval(parse(text='KEGGPATHID2NAME'))
    geneSetName <- unlist(AnnotationDbi::mget(gsub('mmu','',gsets[,'geneSet']),envir))
    gsets <- data.frame(gsets[,1:2],geneSetName=geneSetName,gsets[,3:ncol(gsets)])
  } else if (gseaData$gsetOrigin=='GO') {
    geneSetName <- Term(as.character(gsets[,'geneSet']))
    gsets <- data.frame(gsets[,1:2],geneSetName=geneSetName,gsets[,3:ncol(gsets)])
  }

  #set links
  linksout <- tiny.pic <- vector('list',ncol(gsets))
  linksout[[2]] <- paste('gsets_',mylabel,'/gset_',1:nrow(gsets),'.html',sep='')
  link2plot <- paste('plots_',mylabel,'/gset_',1:nrow(gsets),'.png',sep='')
  if (gseaData$gsetOrigin=='KEGG') {
    linksout[[3]] <- paste('http://www.genome.jp/dbget-bin/www_bget?map',gsub('mmu','',as.character(gsets$geneSet)),sep='')
    linksout[[4]] <- tiny.pic[[4]] <- link2plot
  } else if (gseaData$gsetOrigin=='GO') {
    linksout[[3]] <- paste('http://amigo.geneontology.org/cgi-bin/amigo/term_details?term=',as.character(gsets[,'geneSet']),sep='')
    linksout[[4]] <- tiny.pic[[4]] <- link2plot
  }

  #write html
  plotsdir <- file.path(path,paste('plots',mylabel,sep='_'))
  gsetsdir <- file.path(path,paste('gsets',mylabel,sep='_'))
  if (file.exists(plotsdir)) unlink(plotsdir, recursive = TRUE)
  if (file.exists(gsetsdir)) unlink(gsetsdir, recursive = TRUE)
  dir.create(plotsdir, showWarnings=F)
  dir.create(gsetsdir, showWarnings=F)
  write.html(gsets,file=file.path(path,file),links=linksout,tiny.pic=tiny.pic,title=title,digits=digits)
  for (i in 1:nrow(gsets)) {
    gsetName <- as.character(gsets[i,'geneSet'])
    dat <- data.frame(anno[match(sigs[[gsetName]],anno$entrezid),,drop=FALSE],genes[sigs[[gsetName]],,drop=FALSE])
    links <- vector('list',ncol(dat))
    links[[1]] <- paste('http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids=',anno[sigs[[gsetName]],'entrezid'],sep='')
    links[[2]] <- paste('http://www.genecards.org/index.php?path=/Search/keyword/',anno[sigs[[gsetName]],'symbol'],sep='')
    write.html(dat,links=links,file=paste(path,'/gsets_',mylabel,'/gset_',i,'.html',sep=''),title=title,digits=digits)
    png(file.path(path,paste('plots',mylabel,sep='_'),paste('gset_',i,'.png',sep='')))
    plot(gseaData,es.nes=ifelse(plotEs,'es','nes'),selGsets=as.character(gsetName),selVars=as.character(gsets[i,1]))
    dev.off()
  }
}
