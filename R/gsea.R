#####################
## 1 GET ES ES.SIM ##
#####################

preProcessX <- function(x,logScale=TRUE,absVals=FALSE,averageRepeats=FALSE) {
  if (absVals) x <- abs(x)
  if (averageRepeats) {
    x <- tapply(x,names(x),mean)
  }
  x <- x[complete.cases(x) & sign(x)!=0]
  x <- x[!is.na(as.character(names(x)))]
  x <- x[order(x,decreasing=TRUE)]  
  if (logScale) x <- ifelse(x==0,0,log(abs(x)) * sign(x)) #else x <- ifelse(x==0,0,x)
  return(x)
}

selSignatures <- function(x,signatures) {
  sel <- unlist(lapply(signatures,function(y) any(y %in% names(x))))
  if (any(sel)) {
    signatures <- signatures[sel]
  } else {
    stop('No genes of the provided signatures are in names(x).')
  }
  signNotInX <- any(unlist(lapply(signatures, function(y) any(!(y %in% names(x))))))
  if (signNotInX) {
    warning('There is at least one signature that contains elements that are not in x!')
  }
  return(signatures)
}

getSourceData <- function(x,s,B=1000,mc.cores=1,test='perm') {
  s <- names(x) %in% unique(as.character(s))
  enrichmentScore <- function(fchr,sign) {
    nfchr <- length(fchr)
    p <- numeric(nfchr)
    p[] <- -1 / (nfchr - length(sign))
    absfchr <- abs(fchr[sign])
    p[sign] <- absfchr / sum(absfchr)
    cumsum(p)
  }
  getEs <- function(x,s) {
#    es <- .Call('getEs',x,as.integer(which(s)),PACKAGE='phenoTest') #this line can be used instead of the next one to use C instead of R
    es <- enrichmentScore(x,as.integer(which(s))) 
    return(es)
  }
  getEsSim <- function(x,s,B,mc.cores) {
    myFun1 <- function(dummy) { return(which(sample(s,replace=FALSE))) } 
    myFun2 <- function(s.rand) {
#      ans <- .Call('getEs',x,as.integer(s.rand),PACKAGE='phenoTest') #this line can be used instead of the next one to use C instead of R
      ans <- enrichmentScore(x,as.integer(s.rand))
      ans <- ans[abs(ans)==max(abs(ans))]
      return(ans)
    }
    s.rand <- lapply(vector("list", B),myFun1)
    if (mc.cores==1) {
      es.sim <- unlist(lapply(s.rand,myFun2))      
    } else {
      if ('multicore' %in% loadedNamespaces()) {
        es.sim <- unlist(multicore::mclapply(s.rand,myFun2,mc.set.seed=FALSE,mc.preschedule=TRUE,mc.cores=mc.cores))
      } else stop('multicore library has not been loaded!')
    }
    return(es.sim)
  }
  es <- getEs(x,s)
  if (test=='perm') es.sim <- getEsSim(x,s,B=B,mc.cores=mc.cores) else es.sim <- NA
  ans <- list(es=es,es.sim=es.sim,signature=which(s))
  return(ans)
}



###############
## 2 SUMMARY ##
###############

getSummary <- function(es,es.sim,fchr,p.adjust.method='none',pval.comp.method='original',pval.smooth.tail=TRUE,signatures,test) {
  getNesSim <- function(es.sim) {
    nes.sim <- es.sim / ifelse(es.sim>0,mean(es.sim[es.sim>0]),-mean(es.sim[es.sim<0]))
    return(nes.sim)
  }
  nes.sim <- unlist(lapply(es.sim,getNesSim))
  #
  getEsScore <- function(es) {
    escore <- range(es,na.rm=TRUE)
    escore <- escore[abs(escore)==max(abs(escore))][1]
  }
  gseaSignificance <- function(es,es.sim,fchr.avg,pval.comp.method,pval.smooth.tail) {
    getNes <- function(escore,es.sim) {
      if (any(sign(es.sim)==sign(escore))) {
        ms <- abs(mean(es.sim[sign(es.sim)==sign(escore)],na.rm=TRUE))
        nes <- escore / ms
      } else {
        nes <- NA
      }
      return(nes)
    }
    getEsPval <- function(escore,es.sim,fchr.avg,pval.comp.method,pval.smooth.tail) {
      if (pval.comp.method=='original') {
        if (sum(sign(es.sim)==sign(escore),na.rm=TRUE)>1) {
          if (sign(fchr.avg)!=sign(escore)) {
            mypval <- 1
          } else {
            es.sim.sign <- abs(es.sim[sign(es.sim)==sign(escore)])
            if (pval.smooth.tail) {
              myDens <- density(es.sim.sign,from=0,to=max(abs(escore),es.sim.sign)*2,n=2056)
              myFun <- approxfun(myDens$x,myDens$y)
              mypval <- try((1 - (integrate(myFun,0,abs(escore))$value / integrate(myFun,0,max(abs(escore),es.sim.sign)*2)$value)),silent=TRUE)
              if (inherits(mypval, "try-error")) {
                mypval <- NA
              }
            } else {
              mypval <- sum(abs(es.sim.sign) > abs(escore)) / length(es.sim.sign)
            }
          }
        } else {
          mypval <- NA
        }
      } else if (pval.comp.method=='signed') {
          if (pval.smooth.tail) {
            myDens <- density(es,from=min(escore,es.sim)*2,to=max(escore,es.sim)*2,n=2056)
            myFun <- approxfun(myDens$x,myDens$y)
            if (sign(escore)==1) {
              tmp <- integrate(myFun,min(escore,es.sim)*2,escore)$value
              mypval <- (1 - (tmp / (tmp+integrate(myFun,escore,max(escore,es.sim)*2)$value)))
            } else {
              tmp <- integrate(myFun,escore,max(escore,es.sim)*2)$value
              mypval <- (1 - (tmp / (tmp+integrate(myFun,min(escore,es.sim)*2,escore)$value)))
            }
          } else {
            if (sign(escore)==1) mypval <- sum(es.sim > escore) / length(es.sim) else mypval <- sum(es.sim < escore) / length(es.sim)
          }
      } else {
        stop("pval.comp.method has to be 'original' or 'signed'")
      }
      mypval <- ifelse(mypval<0,0,ifelse(mypval>1,1,mypval))
      return(mypval)
    }
    getNesPval <- function(nes,fchr.avg,pval.comp.method,pval.smooth.tail) { # nes.sim has been computed in getNesSim function
      if (!is.na(nes)) {
        if (pval.comp.method=='original') {
          if (sum(sign(nes.sim)==sign(nes),na.rm=TRUE)>1) {
            if (sign(fchr.avg)!=sign(escore)) {
              mypval <- 1
            } else {
              nes.sim.sign <- abs(nes.sim[sign(nes.sim)==sign(nes)])
              if (pval.smooth.tail) {
                myDens <- density(nes.sim.sign,from=0,to=max(abs(nes),nes.sim.sign)*2,n=2^13)
                myFun <- approxfun(myDens$x,myDens$y)
                mypval <- (1 - (integrate(myFun,0,abs(nes))$value / integrate(myFun,0,max(abs(nes),nes.sim.sign)*2)$value))
              } else {
                mypval <- sum(abs(nes.sim.sign) > abs(nes)) / length(nes.sim.sign)
              }
            }
          } else {
            mypval <- NA
          }
        } else if (pval.comp.method=='signed') {
          if (pval.smooth.tail) {
            myDens <- density(nes.sim,from=min(nes,nes.sim)*2,to=max(nes,nes.sim)*2,n=2^13)
            myFun <- approxfun(myDens$x,myDens$y)
            if (sign(escore)==1) {
              tmp <- integrate(myFun,min(nes,nes.sim)*2,nes)$value
              mypval <- (1 - (tmp / (tmp+integrate(myFun,nes,max(nes,nes.sim)*2)$value)))
            } else {
              tmp <- integrate(myFun,nes,max(nes,nes.sim)*2)$value
              mypval <- (1 - (tmp / (tmp+integrate(myFun,min(nes,nes.sim)*2,nes)$value)))
            }
          } else {
            if (sign(escore)==1) mypval <- sum(nes.sim > nes) / length(nes.sim) else mypval <- sum(nes.sim < nes) / length(nes.sim)
          }
        } else {
          stop("pval.comp.method has to be 'original' or 'signed'")
        }
      } else {
        mypval <- NA
      }
      mypval <- ifelse(mypval<0,0,ifelse(mypval>1,1,mypval))
      return(mypval)
    }
    escore <- getEsScore(es)
    nes <- getNes(escore,es.sim)
    pval.nes <- getNesPval(nes,fchr.avg,pval.comp.method,pval.smooth.tail)
    pval.es <- getEsPval(escore,es.sim,fchr.avg,pval.comp.method,pval.smooth.tail)
    return(data.frame(es=escore,nes=nes,pval.es=pval.es,pval.nes=pval.nes))
  }
  getFdr <- function(nes,pval.nes) {
    fdr <- vector('numeric',length=length(nes)); names(fdr) <- names(nes)
    nes.numeric <- as.numeric(nes)
    nes.numeric <- nes.numeric[!is.na(nes.numeric)]
    myFun <- function(x,y) { x / (sum(abs(nes.numeric[sign(nes.numeric)==sign(y)])>=abs(y)) / length(nes.numeric[sign(nes.numeric)==sign(y)])) }
    fdr <- mapply(function(x,y) myFun(x,y),x=pval.nes,y=nes)
    fdr <- ifelse(fdr<0,0,ifelse(fdr>1,1,fdr))
    return(fdr)
  }
  #
  if (test!='perm') {
    if (test=='wilcox') {
#      pval <- mapply(function(x,y) wilcox.test(x[y],x[-y])$p.value,x=es,y=signatures)
#      pval <- mapply(function(x,y) wilcox.test(x[y],mu=0)$p.value,x=fchr,y=signatures)
      escore <- unlist(lapply(signatures,function(x) log2(mean(2^fchr[x]))))
      fchr.avg <- mean(fchr)
      pval <- unlist(lapply(signatures,function(x) max(c(wilcox.test(fchr[x],mu=0)$p.value,wilcox.test(fchr[x],mu=fchr.avg)$p.value))))
    } else if (test=='ttperm') {    
      myttperm <- function(x,y) {
        perm <- ttperm(matrix(x,nrow=1),factor(1:length(x) %in% y),B=1000)
        t.mine <- perm[[1]][,'statistic']
        t.perm <- unlist(lapply(perm[[2]], function(x) x[,'statistic']))
        ans <- sum((-abs(t.mine))<t.perm | abs(t.mine)>t.perm) / length(t.perm)
        return(ans)
      }
      escore <- unlist(lapply(es, getEsScore))
      pval <- mapply(function(x,y) myttperm(x,y),x=es,y=signatures)
    } else {
      stop('test has to be perm, wilcox or ttperm')
    }
    ans <- cbind(es=escore,pval)
    if (p.adjust.method!='none') ans[,'pval'] <- p.adjust(as.numeric(ans[,'pval']),p.adjust.method)
  } else {
    fchr.avg <- lapply(signatures,function(x) mean(fchr[x]))    
    ans <- mapply(function(x,y,z) gseaSignificance(x,y,z,pval.comp.method,pval.smooth.tail),x=es,y=es.sim,z=fchr.avg)
    if (class(ans)=='list') ans <- do.call(cbind,ans) #this line solves an error that appears on windows Windows Server 2003 R2 (32-bit) / x64 
    ans.tmp <- matrix(as.numeric(ans),ncol=ncol(ans)); colnames(ans.tmp) <- colnames(ans); rownames(ans.tmp) <- rownames(ans); ans <- ans.tmp
    fdr <- getFdr(ans['nes',],ans['pval.nes',])
    ans <- rbind(ans,fdr=fdr)
    if (p.adjust.method!='none') {
      ans['pval.es',] <- p.adjust(as.numeric(ans['pval.es',]),p.adjust.method)
      ans['pval.nes',] <- p.adjust(as.numeric(ans['pval.nes',]),p.adjust.method)
    }
    ans <- t(ans)
  }
  return(ans)
}



##################
## 3 MAKE PLOTS ##
##################

plotGSEA <- function(es.nes,fc.hr,s,mainTitle='',variable='',pvalfdr,p.adjust.method='',EsOrNes='ES',es.nes.ylim,test) {
  #functions
  plot1.perm <- function(x,myYlim,mainTitle,EsOrNes) {
    myDat <- approx(1:length(x),x,n=500)    
    plot(myDat$x,myDat$y,type='l',axes=FALSE,ylab=EsOrNes,xlab='',oma=c(0,0,0,0),main=mainTitle,ylim=myYlim,col='darkgreen',lwd=2)
    abline(h=0)
    axis(2)
  }
  plot1.wilcox <- function(fc.hr,s,mainTitle) {
    if (sum(s)>1) {
      xlim <- range(density(fc.hr[s])$x)[c(2,1)] #we invert the plot to show negatives first (like broad's gsea)
      ylim <- range(density(fc.hr[s])$y)
      plot(density(fc.hr[s]),xlim=xlim,ylim=ylim,col=3,main=mainTitle) 
    } else {
      plot(fc.hr[s],1,col=0)
    }
    abline(v=mean(fc.hr[s]),lty=2,col=3); abline(v=0,col=1,lty=2)
  }
  plot2 <- function(x,varLabel) {
    par(mar=c(1,4,1,1))
    x <- sort(x,decreasing=TRUE)
    myDat <- approx(1:length(x),x,n=300)
    plot(myDat$x,myDat$y,type='l',axes=FALSE,ylim=c(-1*max(abs(x)),max(abs(x))),ylab=paste(varLabel,'(log)'),xlab=''); abline(h=0); axis(2)
    if (any(myDat$y<0)) {
      abline(v=mean(myDat$x[myDat$y==min(myDat$y[myDat$y>0])],myDat$x[myDat$y==max(myDat$y[myDat$y<0])]),lty=2)
    }
    polygon(c(myDat$x,max(myDat$x),0),c(myDat$y,0,0), col=4, border = 4)
  }
  plot3 <- function(x) {
    par(mar=c(1,5.6,1,2.6))
    x.abs <- abs(x)
    x.abs.qua <- quantile(x.abs,(0:11)/11)
    x.abs.qua[1] <- 0
    x.qua <- c(-x.abs.qua[-1][(length(x.abs.qua)-1):1], x.abs.qua)
    x.idx <- ceiling((1:length(x))/50)
    x.small <- tapply(x,x.idx,mean)
    image(as.matrix(x.small),col=c(rgb(0,0.25,0.75,seq(0.1,1,length.out=11)[11:1]),rgb(1,0,0,seq(0.1,1,length.out=11))),breaks=x.qua,axes=FALSE)
  }
  plot4 <- function(s,subTitle) {
    par(mar=c(5,4,1,1))
    plot(NA,NA,xlim=c(1,length(s)),ylim=c(0,1),axes=FALSE,ylab='',xlab='Gene list rank',sub=subTitle,col='green',lwd=2)
    myTicks <- axTicks(1); myTicks[myTicks==0] <- 1; axis(1,lty=0,at=myTicks)
    if (sum(s)<50) myCols <- 1 else myCols <- densCols(which(s[1:length(s)]))
    abline(v=which(s[1:length(s)]),col=myCols)
  }
  #vars
  if (!missing(es.nes.ylim)) { myYlim <- es.nes.ylim } else { myYlim <- max(abs(range(es.nes))); myYlim <- c(-myYlim,myYlim) }
  if (EsOrNes=='ES') {
    mainTitle.es <- paste(EsOrNes,' plot / ',mainTitle,ifelse(is.na(pvalfdr),'(pval=NA)',ifelse(pvalfdr<0.001,'(pval<0.001)',paste('(pval=',round(pvalfdr,3),')',sep=''))),'*',sep=' ')
  } else {
    mainTitle.es <- paste(EsOrNes,' plot / ',mainTitle,ifelse(is.na(pvalfdr),'(fdr=NA)',ifelse(pvalfdr<0.001,'(fdr<0.001)',paste('(fdr=',round(pvalfdr,3),')',sep=''))),sep=' ')
  }
  if (EsOrNes=='ES') {
    subTitle <- ifelse(p.adjust.method=='none','(*) No pvalue adjustment method was used',paste('(*) pvalue adjustment method:',p.adjust.method))
  } else {
    subTitle <- ''
  }
  #plot
  def.par <- par(no.readonly = TRUE)
  par(mar=c(1,4,2,1))
  layout(c(1,2,3,4),height=c(4,2,1,2))
  if (test=='perm') {
    plot1.perm(es.nes,myYlim,mainTitle=mainTitle.es,EsOrNes=EsOrNes)
  } else {
    plot1.wilcox(fc.hr,s,mainTitle.es)
  }
  plot2(fc.hr,variable)
  plot3(fc.hr)
  plot4(s,subTitle)
  par(def.par)
}

plotGseaPreprocess <- function(x,y,z,variable='',es.ylim,nes.ylim,test,es.nes) {
  if (test=='wilcox') es.nes <- 'es'
  if (class(z)=='list') {
    gsl <- x[[1]]
    for (i in 1:length(z)) {
      es <- gsl[[i]][['es']]
      es.sim <- gsl[[i]][['es.sim']]
      fc.hr <- x[[2]]
      s <- names(fc.hr) %in% z[[i]]      
      myTitle <- paste(ifelse(missing(variable),'',paste('variable:',variable,' / ',sep='')),'signature:', names(z)[i],sep='')
      if (es.nes %in% c('es','both')) {
        #es plot
        if (test!='perm') pvalfdr <- y$summary[i,'pval'] else pvalfdr <- y$summary[i,][['pval.es']]
        plotGSEA(es.nes=es,fc.hr=fc.hr,s=s,mainTitle=myTitle,variable=ifelse(missing(variable),'',variable),pvalfdr=pvalfdr,p.adjust.method=y[[2]],EsOrNes='ES',es.nes.ylim=es.ylim,test=test)
      }
      if (es.nes %in% c('nes','both')) {
        #nes plot
        if (test=='perm') {
          es.range <- range(es)
          escore <- es.range[abs(es.range)==max(abs(es.range))]
          if (sum(sign(es.sim)==sign(escore))>2) {
            nes <- es / abs(mean(es.sim[sign(es.sim)==sign(escore)],na.rm=TRUE))
            plotGSEA(es.nes=nes,fc.hr=fc.hr,s=s,mainTitle=myTitle,variable=ifelse(missing(variable),'',variable),pvalfdr=y[[1]][i,][['fdr']],p.adjust.method=y[[2]],EsOrNes='NES',es.nes.ylim=nes.ylim,test=test)
          } else {
            warning('Too few simulations with the same sign as ES were obtained in some signature(s). NAs will be assigned to nes.')
          }
        }
      }
    }
  } else if (class(z)=='character') {
      if (es.nes %in% c('es','both')) {
        #es plot
        es <- x[[1]][[1]][['es']]
        es.sim <- x[[1]][[1]][['es.sim']]
        fc.hr <- x[['fc.hr']]
        s <- names(fc.hr) %in% z
        myTitle <- paste(ifelse(missing(variable),'',paste('variable:',variable,sep='')))
        if (test!='perm') pvalfdr <- summary(y)[i,'pval'] else pvalfdr <- y[[1]][,'pval.es'][[1]]
        plotGSEA(es.nes=es,fc.hr=fc.hr,s=s,mainTitle=myTitle,variable=ifelse(missing(variable),'',variable),pvalfdr=y[[1]][,'pval.es'][[1]],p.adjust.method=y[[2]],EsOrNes='ES',es.nes.ylim=es.ylim,test=test)
      }
      if (es.nes %in% c('nes','both')) {      
        #nes plot
        if (test=='perm') {
          es.range <- range(es)
          escore <- es.range[abs(es.range)==max(abs(es.range))]
          if (sum(sign(es.sim)==sign(escore))>2) {
            nes <- es / abs(mean(es.sim[sign(es.sim)==sign(escore)],na.rm=TRUE))
            plotGSEA(es.nes=nes,fc.hr=fc.hr,s=s,mainTitle=myTitle,variable=ifelse(missing(variable),'',variable),pvalfdr=y[[1]][,'fdr'][[1]],p.adjust.method=y[[2]],EsOrNes='NES',es.nes.ylim=nes.ylim,test=test)
          } else {
            warning('No simulations with the same sign as ES were obtained. NES could not be computed!')
          }
        }
      }
  } else {
    stop('signatures has to be of class list or character')
  }
}



###############
## 4 METHODS ##
###############

# 4.1 methods for gseaSignatures #
##################################

setGeneric("gseaSignatures",function (x,gsets,logScale=TRUE,absVals=FALSE,averageRepeats=FALSE,B=1000,mc.cores=1,test='perm') standardGeneric("gseaSignatures"))
setGeneric("getEs",function (x) standardGeneric("getEs"))
setGeneric("getEsSim",function (x) standardGeneric("getEsSim"))
setGeneric("getNes",function (x) standardGeneric("getNes"))
setGeneric("getFcHr",function (x) standardGeneric("getFcHr"))

setMethod("gseaSignatures",signature(x="numeric",gsets="list"),
 function(x,gsets,logScale=TRUE,absVals=FALSE,averageRepeats=FALSE,B=1000,mc.cores=1,test='perm') {
  x <- preProcessX(x,logScale=logScale,absVals=absVals,averageRepeats=averageRepeats)
  gsets <- selSignatures(x,gsets)
  if (test=='perm') {
    numPerm <- ceiling(B/length(gsets))
    if (numPerm<50) { #check that each gene set has minimum number of permutations
      warning('You are using less than 50 permutations per gene set. Results will be very inaccurate. You should increse the number of permutations!\n')
    } else {
      cat(paste(length(gsets),'gene sets were provided and',B,'permutations were assigned,\n therefore ',numPerm,'permutations will be computed on each gene set.\n'))
    }
    gsets.len <- unlist(lapply(gsets,length)); names(gsets.len) <- names(gsets)
    probNotCentered <- (1-pnorm(0,mean(x),sd(x) * sqrt((1/gsets.len) * (1-gsets.len/length(x)) * (length(x)/(length(x)-1)))))
    if (any(probNotCentered>0.99)) { #check if x is centered
      maxSigLen <- sum(1-pnorm(0,mean(x),sd(x) * sqrt((1/(1:length(x))) * (1-(1:length(x))/length(x)) * (length(x)/(length(x)-1)))) < 0.99)
      cat('Some gene sets are too big to compute permutation tests. Wilcoxon test will be used instead.\n')
      cat('Remove gene sets with more than',maxSigLen,'genes and run again to use permutation tests.\n')    
      test <- 'wilcox'
    }
  }
  ans <- lapply(gsets,function(y) getSourceData(x=x,s=y,B=numPerm,mc.cores=mc.cores,test=test))
  ans <- new("gseaSignaturesSign",list(es.esSim=new("gseaSignatures",ans),fc.hr=x,signatures=gsets,test=test))
  return(ans)
 }
)

setMethod("gseaSignatures",signature(x="epheno",gsets="list"),
 function(x,gsets,logScale=TRUE,absVals=FALSE,averageRepeats=FALSE,B=1000,mc.cores=1,test='perm') {
   phenoTest2list <- function(x) {
      ans <- getSummaryDif(x)
      ans.list <- vector('list',ncol(ans)); names(ans.list) <- colnames(ans)
      for (i in 1:length(ans.list)) ans.list[[i]] <- ans[,i]
      return(ans.list)
   }
   x <- phenoTest2list(x)
   ans <- lapply(x,function(x) gseaSignatures(x,gsets,logScale,absVals=absVals,averageRepeats=averageRepeats,B,mc.cores,test)) # call method(numeric,list)
   return(new("gseaSignaturesVar",ans))
 }
)

setMethod("gseaSignatures",signature(x="matrix",gsets="list"),
 function(x,gsets,logScale=TRUE,absVals=FALSE,averageRepeats=FALSE,B=1000,mc.cores=1,test) {
   matrix2list <- function(x) {
      ans <- x[,grep('\\.HR$|\\.fc$',colnames(x),perl=TRUE),drop=FALSE]
      ans.list <- vector('list',ncol(ans)); names(ans.list) <- colnames(ans)
      for (i in 1:length(ans.list)) ans.list[[i]] <- ans[,i]
      return(ans.list)
   }
   x <- matrix2list(x)
   B <- ceiling(B/length(x))
   ans <- lapply(x,function(x) gseaSignatures(x,gsets,logScale,absVals=absVals,averageRepeats=averageRepeats,B,mc.cores,test)) # call method(numeric,list)
   return(new("gseaSignaturesVar",ans))
 }
)

setMethod("gseaSignatures",signature(gsets="character"),
 function(x,gsets,logScale=TRUE,absVals=FALSE,averageRepeats=FALSE,B=1000,mc.cores=1,test) {
   gsets <- list(signature=gsets)
   ans <- gseaSignatures(x,gsets,absVals=absVals,averageRepeats=averageRepeats,logScale,B,mc.cores,test) #call method( ,list)
   return(ans)
 }
)

#method for signature=geneset has to be defined for every possible class of x (there seems to be an error with R identifying GeneSetCollection objects)
setMethod("gseaSignatures",signature(x="numeric",gsets="GeneSetCollection"),
 function(x,gsets,logScale=TRUE,absVals=FALSE,averageRepeats=FALSE,B=1000,mc.cores=1,test='perm') {
   gsets <- geneIds(gsets)
   ans <- gseaSignatures(x,gsets,absVals=absVals,averageRepeats=averageRepeats,logScale,B,mc.cores,test) #call method( ,list)
   return(ans)
 }
)
setMethod("gseaSignatures",signature(x="epheno",gsets="GeneSetCollection"),
 function(x,gsets,logScale=TRUE,absVals=FALSE,averageRepeats=FALSE,B=1000,mc.cores=1,test='perm') {
   tmpgsets <- geneIds(gsets); names(tmpgsets) <- names(gsets); gsets <- tmpgsets
   ans <- gseaSignatures(x,gsets,absVals=absVals,averageRepeats=averageRepeats,logScale,B,mc.cores,test) #call method( ,list)
   return(ans)
 }
)
setMethod("gseaSignatures",signature(x="matrix",gsets="GeneSetCollection"),
 function(x,gsets,logScale=TRUE,absVals=FALSE,averageRepeats=FALSE,B=1000,mc.cores=1,test='perm') {
   tmpgsets <- geneIds(gsets); names(tmpgsets) <- names(gsets); gsets <- tmpgsets   
   ans <- gseaSignatures(x,gsets,absVals=absVals,averageRepeats=averageRepeats,logScale,B,mc.cores,test) #call method( ,list)
   return(ans)
 }
)

setMethod("gseaSignatures",signature(gsets="GeneSet"),
 function(x,gsets,logScale=TRUE,absVals=FALSE,averageRepeats=FALSE,B=1000,mc.cores=1,test) {
   gsets <- geneIds(gsets)
   ans <- gseaSignatures(x,gsets,absVals=absVals,averageRepeats=averageRepeats,logScale,B,mc.cores,test) #call method( ,character)
   return(ans)
 }
)

setMethod("getEs",signature(x="gseaSignaturesSign"),
  function (x) {
    ans <- lapply(x[[1]],function(x) x[['es']])
    return(ans)
  }
)

setMethod("getEs",signature(x="gseaSignaturesVar"),
  function (x) {
    ans <- lapply(x,getEs)
    return(ans)
  }
)

setMethod("getEsSim",signature(x="gseaSignaturesSign"),
  function (x) {
    ans <- lapply(x[[1]],function(x) x[['es.sim']])
    return(ans)
  }
)

setMethod("getEsSim",signature(x="gseaSignaturesVar"),
  function (x) {
    ans <- lapply(x,getEsSim)
    return(ans)
  }
)

setMethod("getNes",signature(x="gseaSignaturesSign"),
  function (x) {
    myFun <- function(es,es.sim) {
      es.range <- range(es)
      escore <- es.range[abs(es.range)==max(abs(es.range))]
      ans <- es / abs(mean(es.sim[sign(es.sim)==sign(escore)],na.rm=TRUE))
      ans
    }
    es <- getEs(x)
    es.sim <- getEsSim(x)
    ans <- mapply(myFun, es=es,es.sim=es.sim)
    mynames <- colnames(ans)
    ans <- split(ans, 1:ncol(ans))
    names(ans) <- mynames
    return(ans)
  }
)

setMethod("getNes",signature(x="gseaSignaturesVar"),
  function (x) {
    ans <- lapply(x,getNes)
    return(ans)
  }
)
setMethod("getFcHr",signature(x="gseaSignaturesSign"),
  function (x) {
    return(x[['fc.hr']])
  }
)

setMethod("getFcHr",signature(x="gseaSignaturesVar"),
  function (x) {
    ans <- lapply(x,getFcHr)
    return(ans)
  }
)

setMethod("show",signature(object="gseaSignaturesSign"),
  function (object) {
    cat("Object of class 'gseaSignaturesSign'\n")
    cat("You can use the getEs, getNes, getEsSim and getFcHr methods to easily acces its data\n")
    cat("-The tested gene sets are:\n ",names(object[[1]]),'\n')
  }
)

setMethod("show",signature(object="gseaSignaturesVar"),
  function (object) {
    cat("Object of class 'gseaSignaturesVar'\n")
    cat("You can use the getEs, getNes, getEsSim and getFcHr methods to easily acces its data\n")
    cat("-The tested variables are:\n ",names(object),'\n')
    cat("-The tested gene sets (for each variable) are\n :",names(object[[1]][[1]]),'\n')
  }
)


# 4.2 methods for gseaSignificanceSign #
########################################

setGeneric("gseaSignificance",function (x,p.adjust.method='none',pval.comp.method='original',pval.smooth.tail=TRUE) standardGeneric("gseaSignificance"))

setMethod("gseaSignificance",signature(x="gseaSignaturesSign"),
 function(x,p.adjust.method='none',pval.comp.method='original',pval.smooth.tail=TRUE) {
  y <- x[[1]]
  es <- lapply(y,function(x) x[['es']])
  es.sim <- lapply(y,function(x) x[['es.sim']])
  signatures <- lapply(y,function(x) x[['signature']])
  test <- x[['test']]
  fchr <- x[['fc.hr']]
  ans <- cbind(n=unlist(lapply(x$s,length)),getSummary(es,es.sim,fchr,p.adjust.method=p.adjust.method,pval.comp.method,pval.smooth.tail,signatures,test))
  if (any(is.na(ans))) warning('Some NAs were produced due to not having enough simulations with the same sign as the observed ES.')
  return(new("gseaSignificanceSign",list(summary=as.matrix(ans),p.adjust.method=p.adjust.method)))
 }
)

setMethod("gseaSignificance",signature(x="gseaSignaturesVar"),
 function(x,p.adjust.method='none',pval.comp.method='original',pval.smooth.tail=TRUE) {
  ans <- lapply(x,function(x) gseaSignificance(x,p.adjust.method,pval.comp.method,pval.smooth.tail))
  return(new("gseaSignificanceVar",ans))
 }
)

summary.gseaSignificanceSign <- function (object,...) {
  return(object[[1]])
}

summary.gseaSignificanceVar <- function (object,...) {
  ans <- lapply(object,function(object) object[['summary']])
  idx <- rep(seq(1,length(ans),1),nrow(ans[[1]]))
  idx <- idx[order(idx)]
  ans <- data.frame(variable=names(ans)[idx],geneSet=rownames(ans[[1]]),do.call('rbind',ans),row.names=NULL)
  return(ans)
}

setMethod("show",signature(object="gseaSignificanceSign"),
  function (object) {
    cat("Object of class 'gseaSignificanceSign'\n")
    cat("You can use the summary method to easily acces its data\n")
    cat("-The tested gene sets are:\n",rownames(object[[1]]),'\n')
    cat("-P adjust method:",object[[2]],'\n')
  }
)

setMethod("show",signature(object="gseaSignificanceVar"),
  function (object) {
    cat("Object of class 'gseaSignificanceVar'\n")
    cat("You can use the summary method to easily acces its data\n")
    cat("-The tested variables are:\n",names(object),'\n')
    cat("-The tested gene sets (for each variable) are:\n",rownames(object[[1]][[1]]),'\n')
    cat("-P adjust method:",object[[1]][[2]],'\n')
  }
)


# 4.3 plot methods #
####################

plot.gseaSignaturesSign <- function(x,gseaSignificance,es.ylim,nes.ylim,es.nes='both',...) {
  if (class(gseaSignificance)!="gseaSignificanceSign") stop("if x is off class gseaSignaturesSign gseaSignificance has to be of class gseaSignificanceSign")
  if (missing(es.ylim)) {
    es.abs.max <- max(abs(summary(gseaSignificance)[,'es']),na.rm=TRUE)
    es.ylim <- c(-es.abs.max,es.abs.max)
  }
  test <- x['test']
  if (test!='perm') es.nes <- 'es'
  if (missing(nes.ylim) & missing(test)) {
    nes.abs.max <- max(abs(summary(gseaSignificance)[,'nes']),na.rm=TRUE)
    nes.ylim <- c(-nes.abs.max,nes.abs.max)
  }
  mysign <- x$signatures
  plotGseaPreprocess(x=x,y=gseaSignificance,z=mysign,es.ylim=es.ylim,nes.ylim=nes.ylim,test=test,es.nes=es.nes)
}

plot.gseaSignaturesVar <- function(x,gseaSignificance,es.ylim,nes.ylim,es.nes='both',...) {
  if (class(gseaSignificance)!="gseaSignificanceVar") stop("if x is off class gseaSignaturesVar gseaSignificance has to be of class gseaSignificanceVar")
  if (missing(es.ylim)) {
    es.abs.max <- max(abs(summary(gseaSignificance)[,'es']),na.rm=TRUE)
    es.ylim <- c(-es.abs.max,es.abs.max)
  }
  test <- unique(unlist(lapply(x,function(x) x['test'])))
  if (missing(nes.ylim) & missing(test)) {
    nes.abs.max <-max(abs(summary(gseaSignificance)[,'nes']),na.rm=TRUE)
    nes.ylim <- c(-nes.abs.max,nes.abs.max)
  }
  mysign <- x[[1]]$signatures
  for (i in 1:length(x)) plotGseaPreprocess(x=x[[i]],y=gseaSignificance[[i]],z=mysign,variable=names(x)[i],es.ylim=es.ylim,nes.ylim=nes.ylim,test=test,es.nes=es.nes)
}


# 4.4 Methods for gseaData #
############################

checkGsetLen <- function(gsets,minGenes,maxGenes) {
  gsets.len <- unlist(lapply(gsets,length))
  if (any(gsets.len<minGenes)) warning(paste(sum(gsets.len<minGenes),' gene sets were removed due to having less than minGenes (',minGenes,') genes.',sep=''))
  if (any(gsets.len>maxGenes)) warning(paste(sum(gsets.len>maxGenes),' gene sets were removed due to having more than maxGenes (',maxGenes,') genes.',sep=''))
  gsets <- gsets[gsets.len>=minGenes & gsets.len<=maxGenes] 
  gsets
}

gsea <- function(x,gsets,logScale=TRUE, absVals=FALSE, averageRepeats=FALSE, B=1000, mc.cores=1, test="perm",
                    p.adjust.method="none", pval.comp.method="original", pval.smooth.tail=TRUE, minGenes=10,maxGenes=500) {
  gsets <- checkGsetLen(gsets,minGenes,maxGenes)
  sim <- gseaSignatures(x=x, gsets=gsets, logScale=logScale, absVals=absVals, averageRepeats=averageRepeats, B=B, mc.cores=mc.cores, test=test)
  pva <- gseaSignificance(sim, p.adjust.method=p.adjust.method, pval.comp.method=pval.comp.method, pval.smooth.tail=pval.smooth.tail)
  ans <- list(simulations=sim, significance=pva,gsetOrigin='User')
  ans <- new('gseaData',ans)
  return(ans)
}

gsea.go <- function(x,species='Hs', ontologies='MF', logScale=TRUE, absVals=FALSE, averageRepeats=FALSE, B=1000, mc.cores=1, test="perm",
                    p.adjust.method="none", pval.comp.method="original", pval.smooth.tail=TRUE,minGenes=10,maxGenes=500) {
  if (class(x)=='numeric') {
    x.names <- names(x)
  } else if (class(x)=='matrix') {
    x.names <- rownames(x)
  } else if (class(x)=='epheno') {
    x.names <- featureNames(x)
  }
  geo <- getGo(species, ontologies)
  geo <- lapply(geo,function(x) x[x %in% x.names])
  geo <- geo[unlist(lapply(geo,length))>0]
  if (!any(x.names %in% unique(unlist(geo)))) stop('None of the identifiers in x exist in kegg:\nMake sure your identifiers are entrezids and that your data is at gene level!')
  ans <- gsea(x=x, gsets=geo,logScale=logScale, absVals=absVals, averageRepeats=averageRepeats, B=B, mc.cores=mc.cores, test=test, p.adjust.method=p.adjust.method, pval.comp.method=pval.comp.method, pval.smooth.tail=pval.smooth.tail, minGenes=minGenes, maxGenes=maxGenes)
  ans[['gsetOrigin']] <- 'GO'
  return(ans)
}

gsea.kegg <- function(x,species='Hs',logScale=TRUE, absVals=FALSE, averageRepeats=FALSE, B=1000, mc.cores=1, test="perm",
                    p.adjust.method="none", pval.comp.method="original", pval.smooth.tail=TRUE,minGenes=10,maxGenes=500) {
  stopifnot(class(x) %in% c('numeric','matrix','epheno'))
  if (class(x)=='numeric') {
    x.names <- names(x)
  } else if (class(x)=='matrix') {
    x.names <- rownames(x)
  } else if (class(x)=='epheno') {
    x.names <- featureNames(x)
  }
  kegg <- getKegg(species)
  kegg <- lapply(kegg,function(x) x[x %in% x.names])
  kegg <- kegg[unlist(lapply(kegg,length))>0]
  if (!any(x.names %in% unique(unlist(kegg)))) stop('None of the identifiers in x exist in kegg:\nMake sure your identifiers are entrezids and that your data is at gene level!')
  ans <- gsea(x=x, gsets=kegg,logScale=logScale, absVals=absVals, averageRepeats=averageRepeats, B=B, mc.cores=mc.cores, test=test, p.adjust.method=p.adjust.method, pval.comp.method=pval.comp.method, pval.smooth.tail=pval.smooth.tail, minGenes=minGenes, maxGenes=maxGenes)
  ans[['gsetOrigin']] <- 'KEGG'
  return(ans)
}

setMethod("show",signature(object="gseaData"),
  function (object) {
    cat("Object of class 'gseaData'\n")
    cat("You can use the summary method to produce result summaries\n")
    cat("You can use the getEs, getNes, getEsSim and getFcHr methods to easily acces its data\n")
    cat("-The tested variables are:\n ",paste(names(object[[1]]),collapse=', '),'\n')
    cat("-The tested gene sets (for each variable) are:\n",paste(names(object[[1]][[1]][[1]]),collapse=', '),'\n')
  }
)

setMethod("getEs",signature(x="gseaData"),
  function (x) {
    getEs(x[[1]])
  }
)

setMethod("getEsSim",signature(x="gseaData"),
  function (x) {
    getEsSim(x[[1]])
  }
)

setMethod("getNes",signature(x="gseaData"),
  function (x) {
    getNes(x[[1]])
  }
)

setMethod("getFcHr",signature(x="gseaData"),
  function (x) {
    getFcHr(x[[1]])
  }
)

summary.gseaData <- function (object,...) {
  summary(object[[2]])
}

gsea.selGsets <- function(x,selGsets) {
  sim <- x[['simulations']]
  pva <- x[['significance']]
  if (class(sim)=='gseaSignaturesVar') {
    sim <- lapply(sim,function(x) list(es.esSim=x$es.esSim[selGsets], fc.hr=x$fc.hr, signatures=x$signatures[selGsets], test=x$test))
    sim <- new("gseaSignaturesVar",sim)
  } else {
    sim <- list(es.esSim=sim$es.esSim[selGsets], fc.hr=sim$fc.hr, signatures=sim$signatures[selGsets], test=sim$test)
    sim <- new("gseaSignaturesSign",sim)
  }
  if (class(pva)=='gseaSignificanceVar') {
    pva <- lapply(pva,function(x) list(summary=x$summary[selGsets,,drop=F], p.adjust.method=x$p.adjust.method))
    pva <- new("gseaSignificanceVar",pva)
  } else {
    pva <- list(summary=pva$summary[selGsets,,drop=F], p.adjust.method=pva$p.adjust.method)
    pva <- new("gseaSignificanceSign",pva)
  }
  ans <- new('gseaData',list(simulations=sim, significance=pva))
  return(ans)
}

gsea.selVars <- function(x,selVars) {
  sim <- x[['simulations']]
  pva <- x[['significance']]
  sim <- sim[selVars]
  sim <- new("gseaSignaturesVar",sim)
  pva <- pva[selVars]
  pva <- new("gseaSignificanceVar",pva)
  ans <- new('gseaData',list(simulations=sim, significance=pva))
  return(ans)
}

plot.gseaData <- function(x,selGsets,selVars,...) {
  if (!missing(selGsets)) {
    if (class(x[[1]])=='gseaSignaturesVar') signatureNames <- names(x[[1]][[1]][['es.esSim']]) else signatureNames <- names(x[[1]][['es.esSim']])
    stopifnot(!any(!selGsets %in% signatureNames))
    x <- gsea.selGsets(x,selGsets)
  }
  if (!missing(selVars) & class(x[[1]])=='gseaSignaturesVar') {
    stopifnot(!any(!selVars %in% names(x[[1]])))
    x <- gsea.selVars(x,selVars)
  }
  plot(x[[1]],x[[2]],...)
}
