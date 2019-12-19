getChrWithMedVar <- function(x) {
  sel <- table(x$chr) > quantile(table(x$chr),c(0.10)) #remove small chr
  sel <- names(sel[sel])
  x <- x[as.character(x$chr) %in% sel,]
  ans <- by(x,as.character(x$chr),function(x) sum(x$es^2,na.rm=T)/nrow(x))
  mynames <- names(ans)
  ans <- as.numeric(ans)
  names(ans) <- mynames
  ans <- names(sort(ans)[floor(length(ans)/2)])
  ans
}

getEsPositions <- function(epheno,phenoName,organism='human',logEs=T,center=FALSE) {
  #get gene pos from biomaRt
  stopifnot(class(epheno) == 'epheno')
  stopifnot(phenoName %in% phenoNames(epheno))
  if (annotation(epheno) == 'entrezid') {
    stopifnot(organism %in% c('human','mouse','drosophila'))
    if (organism=='human') mymart <- useMart("ensembl", "hsapiens_gene_ensembl")
    if (organism=='mouse') mymart <- useMart("ensembl", "mmusculus_gene_ensembl")
    if (organism=='dmelanogaster') mymart <- useMart("ensembl", "dmelanogaster_gene_ensembl")      
    homol <- getLDS(attributes = c("entrezgene"), mart = mymart, attributesL = c("start_position",'end_position','chromosome_name'), martL = mymart)
    homol <- homol[apply(homol,1,function(x) !any(is.na(x))),]
    homol[,1] <- as.character(homol[,1])
    gene.pos <- data.frame(ids=homol[,1],start=apply(homol[,2:3],1,min),end=apply(homol[,2:3],1,max),chr=homol[,4])
  } else {
    eval(parse(text=paste('library(',annotation(epheno),'.db)',sep='')))
    chr <- eval(parse(text=paste('unlist(lapply(AnnotationDbi::mget(featureNames(epheno),',annotation(epheno),'.db::',annotation(epheno),'CHR),function(x) x[1]))',sep='')))
    st <- eval(parse(text=paste('unlist(lapply(AnnotationDbi::mget(featureNames(epheno),',annotation(epheno),'CHRLOC),function(x) min(abs(x))))',sep='')))
    en <- eval(parse(text=paste('unlist(lapply(AnnotationDbi::mget(featureNames(epheno),',annotation(epheno),'CHRLOCEND),function(x) max(abs(x))))',sep='')))
    gene.pos <- data.frame(ids=featureNames(epheno),start=st,end=en,chr=chr)
  }
  #obtain es and merge with gene positions
  es <- as.numeric(getSummaryDif(epheno[,phenoName]))
  names(es) <- featureNames(epheno)
  mygene.pos <- gene.pos[match(names(es),gene.pos[,'ids']),]
  x <- data.frame(es,chr=mygene.pos$chr,pos=rowMeans(mygene.pos[,c('start','end')]),row.names=names(es))
  x <- x[apply(x,1,function(x) !any(is.na(x))),]
  #some preprocessing
  if (logEs) x$es <- log(abs(x$es))*sign(x$es)
  if (center) x$es <- x$es-mean(x$es)
  x$pos <- x$pos / 1e+6
  x <- x[order(x$chr,x$pos),]
  x$chr <- factor(as.character(x$chr))
  if (any(table(x$chr)<30)) {
    warning('chromosomes with less than 30 genes were removed from the analysis.\n')
    x <- x[x$chr %in% names(table(x$chr)[table(x$chr)>=30]),]
  }
  #
  return(x)
}

preProcess4Ace <- function(x,logEs=T,center=FALSE) {
  if (logEs) x$es <- log(abs(x$es))*sign(x$es)
  if (center) x$es <- x$es-mean(x$es)
  x$pos <- x$pos / 1e+6
  x <- x[order(x$chr,x$pos),]
  x$chr <- factor(as.character(x$chr))
  if (any(table(x$chr)<30)) {
    warning('chromosomes with less than 30 genes were removed from the analysis.\n')
    x <- x[x$chr %in% names(table(x$chr)[table(x$chr)>=30]),]
  }
  return(x)
}

#getPred <- function(es,pos='1',sp,k=round(length(es)/10)) {
getPred <- function(es,pos='1',sp,k=20) {
  if (missing(sp)) {
    mygam <- gam(es ~ s(pos,k=k,bs='cr'))
    sp <- mygam$sp
  } else {
    mygam <- gam(es ~ s(pos,sp=sp,k=k,bs='cr'))
  }
  m.pred <- predict(mygam,as.data.frame(pos))
  return(list(m.pred=m.pred,sp=sp))
}

getDistToClosest <- function(pos,numGenes) {
  ans <- sapply(1:length(pos),function(i) {
    tmp <- ifelse(i-(1:numGenes)<=0,NA,i-(1:numGenes))
    selprev <- ifelse(any(!is.na(tmp)),min(tmp,na.rm=T),0)
    selnext <- max(ifelse(i+(1:numGenes)>length(pos),0,i+(1:numGenes)))
    max(pos[i]-pos[selprev],pos[selnext]-pos[i])
  })
  ans
}

getIndivPvals <- function(obs.pred,sim.pred,p.adjust.method='BY',mc.cores=mc.cores,useAllPerm=F,d,B,useIntegrate=TRUE) {
  if (!useIntegrate) {
    pvals <- pnorm(-abs(obs.pred),mean(sim.pred),sd(sim.pred))*2
  } else {
    mymin <- min(c(obs.pred,sim.pred))*1.5
    mymax <- max(c(obs.pred,sim.pred))*1.5
    if (useAllPerm) {
      myDens <- density(sim.pred,from=mymin,to=mymax,n=2056)
      myFun <- approxfun(myDens$x,myDens$y)
      intAll <- integrate(myFun,mymin,mymax)$value
    }
    myDens <- density(sim.pred,from=mymin,to=mymax,n=2056*2)
    myFun <- approxfun(myDens$x,myDens$y)
    getPval <- function(x,sim.pred,useAllPerm,useIntegrate=TRUE) {
      if(useIntegrate) {
        if (!useAllPerm) {
          mymin <- min(c(sim.pred,x))
          mymax <- max(c(sim.pred,x))
          myDens <- density(sim.pred,from=mymin,to=mymax,n=2056)
          myFun <- approxfun(myDens$x,myDens$y)
          intAll <- integrate(myFun,mymin,mymax)$value
        }
        err <- try(
                   if (x>mean(sim.pred)) {
                     tmp <- integrate(myFun,x,mymax)$value
                   } else {
                     tmp <- integrate(myFun,mymin,x)$value
                   }
                   ,silent=TRUE)
        if (inherits(err, "try-error")) {
          ans <- NA
        } else {
          ans <- tmp / intAll
        }
      } else {
        ans <- pnorm(-abs(x),mean(sim.pred),sd(sim.pred))
      }
      ans
    }
    if (useAllPerm) {
      if (mc.cores==1) {
        pvals <- sapply(obs.pred,function(x) getPval(x,sim.pred,useAllPerm,useIntegrate)) * 2
      } else {
        if ('parallel' %in% loadedNamespaces()) {
          pvals <- unlist(parallel::mclapply(as.list(obs.pred),function(x) getPval(x,sim.pred,useAllPerm,useIntegrate),mc.preschedule=T,mc.cores=mc.cores)) * 2
        } else stop('parallel library has not been loaded!')
      }
    } else {
      mypos <- cbind(1:length(d),order(d))
      d <- d[mypos[,2]]
      sel <- lapply(as.list(1:length(d)),function(i) {
        shiftStart <- sum(((i-9):(i+9))<=0)
        shiftEnd <- sum(((i-9):(i+9))>length(d))
        which(order(abs(d[i]-d[((i-9):(i+9))+shiftStart-shiftEnd])) %in% 1:10) + i + shiftStart - shiftEnd - 10
      })
      mypos <- mypos[order(mypos[,2]),]
      sel <- lapply(sel,function(x) mypos[x,1])
      sel <- lapply(sel,function(sel) as.numeric(sapply((0:(B-1))*length(d),function(x) x + sel)))
      if (mc.cores==1) {      
        pvals <- sapply(1:length(obs.pred),function(i) getPval(obs.pred[i],sim.pred[sel[[i]]],useAllPerm,useIntegrate)) * 2
      } else {
      if ('parallel' %in% loadedNamespaces()) {
        pvals <- unlist(parallel::mclapply(as.list(1:length(obs.pred)),function(i) getPval(obs.pred[i],sim.pred[sel[[i]]],useAllPerm,useIntegrate),mc.preschedule=T,mc.cores=mc.cores)) * 2
      } else stop('parallel library has not been loaded!')
    }
    }
  }
  pvals <- ifelse(pvals>1,1,ifelse(pvals<0,0,pvals))
  pvals <- p.adjust(pvals,p.adjust.method)
  pvals
}

getRegions <- function(x,pvals,minGenes,pvalcutoff=0.01,obs.pred,exprScorecutoff=NA) {
  if (is.na(exprScorecutoff)) myrle <- rle(pvals<pvalcutoff) else myrle <- rle(pvals<pvalcutoff & abs(as.numeric(obs.pred))>exprScorecutoff)
  start <- c(names(pvals[1]),names(myrle$lengths[-length(myrle$lengths)]))
  end <- names(myrle$values) 
  tmp <- cbind(start,end)[myrle$values & myrle$lengths>minGenes,]
  ans <- matrix(x[as.character(tmp),'pos'],ncol=2)
  colnames(ans) <- c('start','end')
  ans
}

plotCopyNumber <- function(es,chr,pos,obs.pred,ssr,es.mean,genome,chrLengths) {
  xlim <- c(0,chrLengths[[chr]])
  chr <- gsub('CHR','',gsub('chr','',gsub('CHR0','',gsub('chr0','',chr))))
  lenNames <- gsub('CHR','',gsub('chr','',gsub('CHR0','',gsub('chr0','',names(chrLengths)))))
  plotCyto <- chr %in% lenNames & genome=='hg18' & chr %in% c(as.character(1:22),'X','Y','M')
  if (plotCyto) {
    def.par <- par(no.readonly = TRUE)    
    par(mar=c(1,3,1.5,1))
    layout(c(1,2),heights=c(2.5,0.5))
  }
  plot(pos,es,pch=20,main=paste("chr",chr),col=densCols(pos,es),xlab='Position in chromosome (Mb)',ylab='Expression Score',xlim=xlim/1e6)
  if (!missing(obs.pred)) lines(pos[order(pos)],obs.pred,col='darkred',lwd=2)
  abline(h=0,lty=2,col='grey')
  if (!missing(es.mean)) abline(h=es.mean,lty=3,col='red')
  abline(v=c(seq(0,xlim[2]/1e6,xlim[2]/1e6/10)),col='green',lty=3)
  if (!missing(ssr)) {
    if (nrow(ssr)>0) apply(ssr,1,function(y) rect(y[1],min(es)*2,y[2],max(es)*2,col=paste(rgb(165/255,42/255,42/255),'15',sep='')))
  }
  ## if (plotCyto) {
  ##   par(mar=c(4,3,1.5,1))
  ##   SNPchip::plotCytoband(chr,cex=0.6,build=genome,xlim=xlim,taper=T)
  ##   abline(v=c(seq(0,xlim[2],xlim[2]/10)),col='green',lty=1)
  ##   par(def.par)
  ## }
}

qcPlot <- function(ids,x,B,sim.pred,minGenes=10) {
  d <- unlist(lapply(ids,function(i) getDistToClosest(x[x$chr==i,'pos'],minGenes)))
  d <- rep(d,B)
  pdf('~/Desktop/qc.pdf')
  plot(d,sim.pred,col=densCols(d,sim.pred),pch=20,xlab='distance to the 10th neighbour',ylab='neighbourhood scores')
  mylm <- lm(sim.pred ~ d)
  summary(mylm)
  abline(mylm,col=2)
  boxplot(sim.pred ~ cut(d,50),xlab='distance to the 10th neighbour',ylab='neighbourhood scores' ,main='Hole genome. 50 boxplots')
  boxplot(sim.pred ~ cut(d,100),xlab='distance to the 10th neighbour',ylab='neighbourhood scores',main='Hole genome. 100 boxplots')
  boxplot(sim.pred ~ cut(d,200),xlab='distance to the 10th neighbour',ylab='neighbourhood scores',main='Hole genome. 200 boxplots')
  tmp <- cut(d,50)
  plot(density(sim.pred[tmp==levels(tmp)[1]]),xlim=c(-0.3,0.1),ylim=c(0,50))
  for (i in seq(2,50,5)) if (sum(tmp==levels(tmp)[i])>5) lines(density(sim.pred[tmp==levels(tmp)[i]]),xlim=c(-0.3,0.1),ylim=c(0,50),col=i)
  dev.off()
}

findCopyNumber <- function(x,minGenes=15,B=100,p.adjust.method='BH',pvalcutoff=0.05,exprScorecutoff=NA,mc.cores=1,useAllPerm=F,genome='hg19',chrLengths,sampleGenome=TRUE,useOneChr=FALSE,useIntegrate=TRUE,plot=TRUE,minGenesPerChr=100) {
  #Input has to be a data frame with gene ids as rownames and 3 columns: es (enrichment score), chr (chromosome) and pos (position in chromosome)
  stopifnot(identical(colnames(x),c('es','chr','pos')))
  #remove chr with less than 100 es
  if (any(table(x$chr)<100)) {
    sel.chr <- names(table(x$chr))[table(x$chr)>minGenesPerChr]
    x <- x[x$chr %in% sel.chr,]
    warning('Chromosomes with less than 100 scores have been removed from the analysis.')
  }
  #order
  x <- x[order(x$chr,x$pos),]
  #get smoothed observed values
  cat('get smoothed observed values...')  
  ids <- as.list(levels(factor(as.character(x$chr))))
  names(ids) <- unlist(ids)
  if (mc.cores==1) {
    tmp <- lapply(ids,function(i) getPred(es=x[x$chr==i,'es'],pos=x[x$chr==i,'pos']))
  } else {
    if ('parallel' %in% loadedNamespaces()) {
      tmp <- parallel::mclapply(ids,function(i) getPred(es=x[x$chr==i,'es'],pos=x[x$chr==i,'pos']),mc.preschedule=F,mc.cores=mc.cores)
    } else stop('parallel library has not been loaded!')      
  }
  obs.pred <- lapply(tmp,function(x) x[['m.pred']])
  mysp <- lapply(tmp,function(x) x[['sp']])
  #get smoothed simulated values
  cat('Ok\nget smoothed simulated values...')
  if (useOneChr) {
    selChr <- getChrWithMedVar(x)
    sim <- sapply(1:(B),function(i) sample(x[x$chr==selChr,'es'])) #sample over all genome
    if (mc.cores==1) {
      sim.pred <- unlist(lapply(as.list(1:(B)),function(j) getPred(es=sim[,j],pos=x[x$chr==selChr,'pos'],sp=mysp[[selChr]])$m.pred))
    } else {
      if ('parallel' %in% loadedNamespaces()) {
        sim.pred <- unlist(parallel::mclapply(as.list(1:(B)),function(j) getPred(es=sim[,j],pos=x[x$chr==selChr,'pos'],sp=mysp[[selChr]])$m.pred,mc.preschedule=T,mc.cores=mc.cores))
      } else stop('parallel library has not been loaded!')      
    }
  } else {
    if (sampleGenome) {
      sim <- sapply(1:(B),function(i) sample(x$es)) #sample over all genome
    } else {
      sim <- sapply(1:(B),function(i) unlist(sapply(as.character(unique(x$chr)),function(chr) sample(x[x$chr==chr,'es'])))) #sample each chromosome
    }
    if (mc.cores==1) {
      sim.pred <- unlist(lapply(as.list(1:(B)),function(j) lapply(ids,function(i) getPred(es=sim[x$chr==i,j],pos=x[x$chr==i,'pos'],sp=mysp[[i]])$m.pred)))
    } else {
      if ('parallel' %in% loadedNamespaces()) {
        sim.pred <- unlist(lapply(as.list(1:(B)),function(j) parallel::mclapply(ids,function(i) getPred(es=sim[x$chr==i,j],pos=x[x$chr==i,'pos'],sp=mysp[[i]])$m.pred,mc.preschedule=F,mc.cores=mc.cores)))
      } else stop('parallel library has not been loaded!')      
    }
  }
  #qc plot
#  qcPlot(ids,x,B,sim.pred,minGenes=10)
  #get pvals
  cat('Ok\nget pvals...')
  if (useOneChr) {
    pvals <- getIndivPvals(obs.pred=unlist(obs.pred),sim.pred=sim.pred,p.adjust.method=p.adjust.method,mc.cores=mc.cores,useAllPerm=T,d=d,B=B,useIntegrate=useIntegrate)
  } else {
    if (useAllPerm | !sampleGenome) {
      d <- 0
    } else {
      if (mc.cores==1) {
        d <- unlist(lapply(ids,function(i) getDistToClosest(x[x$chr==i,'pos'],minGenes)))
      } else {
        if ('parallel' %in% loadedNamespaces()) {
          d <- unlist(parallel::mclapply(ids,function(i) getDistToClosest(x[x$chr==i,'pos'],minGenes),mc.preschedule=F,mc.cores=mc.cores))
        } else stop('parallel library has not been loaded!')
      }
    }
    if (sampleGenome) {
      pvals <- getIndivPvals(obs.pred=unlist(obs.pred),sim.pred=sim.pred,p.adjust.method=p.adjust.method,mc.cores=mc.cores,useAllPerm=useAllPerm,d=d,B=B,useIntegrate=useIntegrate)
    } else {
      pvals <- unlist(sapply(as.character(unique(x$chr)),function(i) {
        sel <- x$chr==i
        ans <- getIndivPvals(obs.pred=obs.pred[[i]],sim.pred=sim.pred[sel],p.adjust.method='none',mc.cores=1,useAllPerm=T,d=d,B=B,useIntegrate=useIntegrate)
        ans
      }))
    }
  }
  names(pvals) <- rownames(x)
  #get enriched regions

##   i <- ids[[1]]
##   getRegions(x[x$chr==i,],pvals[x$chr==i],minGenes=minGenes,pvalcutoff=pvalcutoff,obs.pred=obs.pred[[i]],minScore=1)
  
  cat('Ok\nget enriched regions...')  
  ssr <- lapply(ids,function(i) getRegions(x[x$chr==i,],pvals[x$chr==i],minGenes=minGenes,pvalcutoff=pvalcutoff,obs.pred=obs.pred[[i]],exprScorecutoff=exprScorecutoff))
  #make plot
  if (plot) {
    cat('Ok\nmake plot...')  
    if (missing(chrLengths)) {
      chrLengths <- unlist(as.list(by(x,x$chr,function(x) max(x['pos'],na.rm=T)))) * 1e6
      chrLengths <- chrLengths[!is.na(chrLengths)]
    }
    is.num.ids <- !grepl('[a-zA-Z]',names(ids))
    ids.num <- names(ids)[is.num.ids][order(as.numeric(names(ids)[is.num.ids]))]
    ids.chr <- names(ids)[!is.num.ids][order(names(ids)[!is.num.ids])]
    ids <- ids[c(ids.num,ids.chr)]
    chrLengths <- chrLengths[c(ids.num,ids.chr)]
    lapply(ids,function(i) plotCopyNumber(es=x[x$chr==i,'es'],chr=i,pos=x[x$chr==i,'pos'],obs.pred=obs.pred[[i]],ssr=ssr[[i]],es.mean=mean(x$es),genome=genome,chrLengths=chrLengths))
  }
  #return regions
  cat('Ok\nreturn regions...')  
  ssr <- ssr[lapply(ssr,nrow)>0]
  if (length(ssr)>0) {
    ssr <- do.call(rbind,lapply(1:length(ssr),function(i) data.frame(chr=names(ssr)[i],regionInChr=1:nrow(ssr[[i]]),ssr[[i]])))
    rownames(ssr) <- paste('Region',1:nrow(ssr),sep='')
  }
  cat('Ok\n')
  return(ssr)
}

genesInArea <- function(x,regions) {
  if (length(regions)>0) {
    sel <- apply(regions,1,function(i) which(x$chr==as.character(i['chr']) & x$pos>=as.numeric(i['start']) & x$pos<=as.numeric(i['end'])))
    ans <- lapply(sel,function(i) data.frame(rownames(x[i,]),x[i,'es']))
    if (nrow(regions)==1) names(ans) <- rep('Region1',length(ans))
    ans <- lapply(1:length(ans),function(i) data.frame(rep(names(ans[i]),nrow(ans[[i]])),ans[[i]]))
    ans <- do.call(rbind,ans)
    colnames(ans) <- c('Region','gene','es')
  } else {
    ans <- NULL
  }
  ans
}
