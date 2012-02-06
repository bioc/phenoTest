ClusterPhenoTest <- function(x,cluster,vars2test,B=10^4,p.adjust.method='none') {
  if (length(cluster) != ncol(x)) stop('Dimensions of cluster and x do not match')
  if (class(cluster)=='character') cluster <- factor(cluster)
  if (!is.null(vars2test$continuous)) {
    if (sum(vars2test$continuous %in% names(pData(x))) != length(vars2test$continuous)) stop('Invalid variable names in vars2test$continuous')
  }
  if (!is.null(vars2test$categorical)) {
    if (sum(vars2test$categorical %in% names(pData(x))) != length(vars2test$categorical)) stop('Invalid variable names in vars2test$categorical')
  }
  if (!is.null(vars2test$ordinal)) {
    if (sum(vars2test$ordinal %in% names(pData(x))) != length(vars2test$ordinal)) stop('Invalid variable names in vars2test$ordinal')
  }  
  if (!is.null(vars2test$survival)) {
    if (any(!(as.character(vars2test$survival) %in% names(pData(x))))) stop('Invalid variable names in vars2test$survival')
  }

  summaryDif <- pval <- vector("list",4); names(summaryDif) <- names(pval) <- c('cont','categ','ordinal','surv')
  if (!is.null(vars2test$continuous)) {
    mysummary <- sapply(vars2test$continuous,function(y) by(pData(x)[,y,drop=FALSE],cluster,summary),simplify=FALSE)
    pval$cont <- sapply(vars2test$continuous,function(y) kruskal.test(pData(x)[,y],g=cluster)$p.value)
  }

  if (!is.null(vars2test$categorical)) {
    summaryDif$categ <- sapply(vars2test$categorical,function(y) table(factor(pData(x)[,y]),cluster),simplify=FALSE)
    pval$categ <- sapply(summaryDif$categ,function(y) chisq.test(y,simulate.p.value=TRUE,B=B)$p.value)
  }

  if (!is.null(vars2test$ordinal)) {
    summaryDif$ordinal <- sapply(vars2test$ordinal,function(y) table(factor(pData(x)[,y]),cluster),simplify=FALSE)
    pval$ordinal <- sapply(vars2test$ordinal,function(y) kruskal.test(pData(x)[,y],g=cluster)$p.value)
  }

  if (!is.null(vars2test$survival)) {
    summaryDif$surv <- apply(vars2test$survival,1,function(y) coxph(Surv(pData(x)[,y['time']],pData(x)[,y['event']]) ~ cluster))
    pval$surv <- sapply(summaryDif$surv,function(y) summary(y)$logtest['pvalue'])
  }

  allpvals <- p.adjust(c(pval$cont,pval$categ,pval$ordinal,pval$surv),method=p.adjust.method)
  names(allpvals) <- c(vars2test$continuous,vars2test$categorical,vars2test$ordinal,as.character(vars2test$survival[,'time']))
  return(list(summaryDif=summaryDif,p.value=allpvals))
}

setGeneric("heatmapPhenoTest",function(x, signatures, vars2test, probes2genes=FALSE, filterVar, filteralpha=.05, distCol='pearson', nClust=2, distRow='cor', p.adjust.method='none', simulate.p.value=FALSE, B=10^5, linkage='average',equalize=FALSE,center=TRUE,col,survCol,heat.kaplan='both',...) standardGeneric("heatmapPhenoTest"))

setMethod("heatmapPhenoTest",signature(x="ExpressionSet",signatures="character"), function(x, signatures, vars2test, probes2genes=FALSE, filterVar, filteralpha=.05, distCol='pearson', nClust=2, distRow='cor', p.adjust.method='none', simulate.p.value=FALSE, B=10^5, linkage='average',equalize=FALSE,center=TRUE,col,survCol,heat.kaplan='both',...) {
  if (!missing(filterVar)) {  if (is.null(pData(x)[,filterVar])) stop(paste(filterVar,'not found in pData(x)')) }
  idnotfound <- sum(signatures %in% featureNames(x))
  if ((idnotfound>0) & (idnotfound<length(signatures))) {
     warning('Some of the identifiers in signatures are not present in featureNames(x)')
     signatures <- signatures[signatures %in% featureNames(x)]
   } else if (idnotfound==0) {
     stop('None of the specified identifiers could be found in featureNames(x)')
   }
  if (!is.null(vars2test$continuous)) {
    if (sum(vars2test$continuous %in% names(pData(x))) != length(vars2test$continuous)) stop('Invalid variable names in vars2test$continuous')
  }
  if (!is.null(vars2test$categorical)) {
    if (sum(vars2test$categorical %in% names(pData(x))) != length(vars2test$categorical)) stop('Invalid variable names in vars2test$categorical')
  }
  if (!is.null(vars2test$survival)) {
    if (any(!(as.character(vars2test$survival) %in% names(pData(x))))) stop('Invalid variable names in vars2test$survival')    
  }

#Select single probe per gene
x <- x[signatures,]
if (probes2genes==TRUE) {
  if (length(annotation(x))==0) {
    warning('annotation(x) is empty. Cannot select single probe per gene. All rows in exprs(x) were used')
  } else {
    x <- nsFilter(x,require.entrez=FALSE,remove.dupEntrez=TRUE,var.filter=FALSE)$eset
  }
}
  
#Clustering
if (center) { exprs(x) <- exprs(x) - rowMeans(exprs(x)) }
if (distCol=='pearson') {
  hc1 <- hclust(as.dist((1 - cor(exprs(x)))/2), method=linkage)
} else if (distCol=='spearman') {
  hc1 <- hclust(as.dist((1 - cor(exprs(x),method='spearman'))/2), method=linkage)
} else {
  hc1 <- hclust(dist(t(exprs(x)), method=distCol), method=linkage)
}
clus1 <- cutree(hc1,k=nClust)

#Test cluster vs phenotype association
pvals <- ClusterPhenoTest(x,cluster=clus1,vars2test=vars2test,p.adjust.method=p.adjust.method)$p.value
pvals <- ifelse(pvals<.0001,'(P<0.0001)',paste('(P=',round(pvals,4),')',sep=''))
  
#Heatmap
varlabels <- c(vars2test$continuous, vars2test$categorical, vars2test$ordinal, vars2test$survival[,'event'])
addvar <- pData(x)[,varlabels,drop=FALSE]
addvar[is.na(addvar)] <- 0
if (!is.data.frame(addvar)) { addvar <- matrix(addvar,ncol=1); colnames(addvar) <- varlabels }
colnames(addvar) <- paste(colnames(addvar),pvals)
sel <- sapply(addvar,function(x) any(x==TRUE)) #select only logic, or 0/1 variables
if (any(!sel)) {
  warning(paste('The following variables in vars2test will not be ploted because they are not of type logic: ',paste(varlabels[!sel],collapse=', '),sep=''))
  addvar <- addvar[,sel & !is.na(sel),drop=FALSE]
}
eset4plot <- exprs(x); colnames(eset4plot) <- NULL
if (nrow(eset4plot)>100) rownames(eset4plot) <- NULL
  
if (!missing(filterVar)) {
  sel <- !is.na(pData(x)[,filterVar])
  lm1 <- lmFit(x[signatures,sel], model.matrix(~ pData(x)[,filterVar]))
  lm2 <- contrasts.fit(lm1,coefficients=colnames(coef(lm1))[grep('filterVar',colnames(coef(lm1)))])
  eb <- eBayes(lm2)
  if (sum(eb$F.p.value<filteralpha)<2) stop('Could not perform heatmap. When filtering by significant genes less than two were kept.')
  eset4plot <- eset4plot[eb$F.p.value<filteralpha,]
}

cluscol <- gray(1-0.5*((1:max(clus1)) %% 2))
if (distRow=='pearson') distRow <- 'cor'
dx <- as.dist(hopach::as.matrix(distancematrix(eset4plot,d=distRow)))
rowden <- as.dendrogram(hclust(dx,method=linkage))

op <- par(no.readonly=TRUE)

if (heat.kaplan %in% c('heat','both')) {
  if (missing(col)) {
browser()      
    mymax <- 2^max(abs(eset4plot))
    breaks <- log2(c(1,1.15,1.25,seq(1.5,4,length=20),5))
    if (mymax>4) breaks <- c(breaks,log2(mymax))
    breaks <- c(-breaks,breaks); breaks <- breaks[order(breaks)]; breaks <- unique(breaks)
    col <- greenred(length(breaks)-1)
    heatmap_plus(eset4plot,Rowv=rowden,Colv=as.dendrogram(hc1),clus=clus1,addvar=addvar,col=col,breaks=breaks,cluscol=cluscol,equalize=equalize)
  } else {
    heatmap_plus(eset4plot,Rowv=rowden,Colv=as.dendrogram(hc1),clus=clus1,addvar=addvar,col=col,cluscol=cluscol,equalize=equalize)
  }
  par(op)
}
  
if (heat.kaplan %in% c('kaplan','both')) {  
  if (!is.null(vars2test$survival)) {
    for (i in 1:nrow(vars2test$survival)) {
      tname <- vars2test$survival[i,'time']
      ename <- vars2test$survival[i,'event']
      s <- Surv(pData(x)[,tname], pData(x)[,ename]) ~ factor(clus1)
      coxph1 <- coxph(s)
      surv <- survfit(s)
      if (missing(survCol)) survCol <- 1:length(unique(clus1))
      plot(surv, col=survCol, ylab=paste('Survival (',ename,')',sep=''), xlab='Time', ...)
      pval <- 1-pchisq(2*diff(summary(coxph1)$loglik),df=1)
      hr <- round(exp(abs(coef(coxph1)))*sign(coef(coxph1)),3)
      hr <- ifelse(hr>10^3,Inf,hr); hr <- ifelse(hr< -10^3,-Inf,hr)
      text(.1*par('usr')[2],par('usr')[3]+.1*(par('usr')[4]-par('usr')[3]),paste('HR=',paste(hr,collapse=',')),pos=4)
      text(.1*par('usr')[2],par('usr')[3]+.05*(par('usr')[4]-par('usr')[3]),paste('P=',round(pval,4)),pos=4)
    }
  }
  par(op)
}

return(pvals)  
}
)

setMethod("heatmapPhenoTest",signature(x="ExpressionSet",signatures="missing"), function(x, signatures, vars2test, probes2genes=FALSE, filterVar, filteralpha=.05, distCol='pearson', nClust=2, distRow='cor', p.adjust.method='none', simulate.p.value=FALSE, B=10^5, linkage='average',equalize=FALSE,center=FALSE,col,survCol,heat.kaplan='both',...) {
  signatures <- featureNames(x)
  heatmapPhenoTest(x=x, signatures=signatures, vars2test=vars2test, probes2genes=probes2genes, filterVar=filterVar, filteralpha=filteralpha, distCol=distCol, nClust=nClust, distRow=distRow, p.adjust.method=p.adjust.method, simulate.p.value=simulate.p.value, B=B, linkage=linkage,equalize=equalize,center=center,col=col,survCol=survCol,toPDF=TRUE,heat.kaplan=heat.kaplan,...)
}
)

setMethod("heatmapPhenoTest",signature(x="ExpressionSet",signatures="list"), function(x, signatures, vars2test, probes2genes, filterVar, filteralpha=.05, distCol='manhattan', nClust=2, distRow='abscor', p.adjust.method='none', simulate.p.value=FALSE, B=10^5, linkage='complete',equalize=FALSE,center=FALSE,col,survCol,heat.kaplan='both',...) {
  if(!all(sapply(signatures,function(x) is.character(x)))) stop('All elements in argument signatures should be character vectors')
  for (i in 1:length(signatures)) {
    heatmapPhenoTest(x=x, signatures=signatures[[i]], vars2test=vars2test, probes2genes=probes2genes, filterVar=filterVar, distCol=distCol, nClust=nClust, distRow=distRow, p.adjust.method=p.adjust.method, simulate.p.value=simulate.p.value, B=B, linkage=linkage,equalize=equalize,center=center,col=col,survCol=survCol,heat.kaplan=heat.kaplan,...)
  }
}
)

setMethod("heatmapPhenoTest",signature(x="ExpressionSet",signatures="GeneSetCollection"), function(x, signatures, vars2test, probes2genes, filterVar, filteralpha=.05, distCol='manhattan', nClust=2, distRow='abscor', p.adjust.method='none', simulate.p.value=FALSE, B=10^5, linkage='complete',equalize=FALSE,center=FALSE,col,survCol,heat.kaplan='both',...) {
  signatures <- geneIds(signatures)
  heatmapPhenoTest(x=x, signatures=signatures, vars2test=vars2test, probes2genes=probes2genes, filterVar=filterVar, distCol=distCol, nClust=nClust, distRow=distRow, p.adjust.method=p.adjust.method, simulate.p.value=simulate.p.value, B=B, linkage=linkage,equalize=equalize,center=center,col=col,survCol=survCol,heat.kaplan=heat.kaplan,...)
}
)

setMethod("heatmapPhenoTest",signature(x="ExpressionSet",signatures="GeneSet"), function(x, signatures, vars2test, probes2genes, filterVar, filteralpha=.05, distCol='manhattan', nClust=2, distRow='abscor', p.adjust.method='none', simulate.p.value=FALSE, B=10^5, linkage='complete',equalize=FALSE,center=FALSE,col,survCol,heat.kaplan='both',...) {
  signatures <- geneIds(signatures)
  heatmapPhenoTest(x=x, signatures=signatures, vars2test=vars2test, probes2genes=probes2genes, filterVar=filterVar, distCol=distCol, nClust=nClust, distRow=distRow, p.adjust.method=p.adjust.method, simulate.p.value=simulate.p.value, B=B, linkage=linkage,equalize=equalize,center=center,col=col,survCol=survCol,heat.kaplan=heat.kaplan,...)
}
)
