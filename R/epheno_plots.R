barplotCI <- function(x,groups,alpha=.05,exp=FALSE,base='nat',refgroup,plot.ci=TRUE,ci.l=NA,ci.u=NA,order,...) {
  if (!is.factor(groups)) {
    warning('groups have been converted to factor')
    groups <- factor(groups)
  }
  if (length(unique(groups))==1) {
    lm1 <- lm(x ~ 1)
    b <- c(coef(lm1),confint(lm1,level=1-alpha))
    names(b) <- as.character(groups)[1]
    if (exp & base=='nat') { b <- exp(b) } else if (exp & base!='nat') { b <- as.numeric(base)^b }
    barplot2(b[1],names=names(b)[1],plot.ci=TRUE,ci.l=b[2],ci.u=b[3],...)
  } else {
    lm1 <- lm(x ~ -1 + groups)
    b <- cbind(coef(lm1),confint(lm1,level=1-alpha))
    rownames(b) <- sub('groups','',rownames(b))
    if (exp & base=='nat') { b <- exp(b) } else if (exp & base!='nat') { b <- as.numeric(base)^b }
    if (!missing(refgroup)) { b <- b/b[refgroup,1] }
    if (!missing(order)) { b <- b[order,] }
    barplot2(b[,1],names=rownames(b),plot.ci=TRUE,ci.l=b[,2],ci.u=b[,3],...)
  }
}



setGeneric("barplotSignifSignatures",function (x, signatures, referenceSignature, testUpDown=FALSE, simulate.p.value = FALSE, B = 10^4, p.adjust.method='none', alpha=.05, ylab, ylim=ylim, cex.text=1, ...) standardGeneric("barplotSignifSignatures"))

setMethod("barplotSignifSignatures",signature(x="epheno",signatures="list"),
  function (x, signatures, referenceSignature, testUpDown=FALSE, simulate.p.value = FALSE, B = 10^4, p.adjust.method='none', alpha=.05, ylab, ylim=ylim, cex.text=1, ...) {
    if (!missing(referenceSignature)) if (!(referenceSignature %in% names(signatures))) stop('referenceSignature not in names(signatures)')
    lensig <- unlist(lapply(signatures,function(x) length(x)))
    if (all(lensig==0)) stop('All signatures are empty. Check the argument signatures.\n')
    if (any(lensig==0)) {
      warning('Removed some signatures with no genes\n')
      signatures <- signatures[lensig>0]
    }
    badsignature <- sapply(signatures,function(y) !(any(y %in% featureNames(x))))
    if (any(badsignature)) stop(paste('Identifiers in signatures',paste(names(signatures)[which(badsignature)],collapse=','),'not in featureNames'))
    if (!missing(referenceSignature)) { if (!(referenceSignature %in% names(signatures))) stop('referenceSignature has to be one of names(signatures) !') }

    #Compute percentages of DE genes
    myFun1 <- function(i) {
      signde <- sign(rowMeans(getSummaryDif(x[,i]))) * (getPvals(x[,i])<alpha)
      ans <- lapply(signatures, function(y) table(factor(signde[featureNames(x) %in% y],levels=-1:1)))
      ans <- do.call(rbind,ans)
      return(ans)
    }
    x2plot <- sapply(phenoNames(x), myFun1, simplify=FALSE)
    if (!missing(referenceSignature)) {
      myorder <- c(which(names(signatures) %in% referenceSignature),which(!(names(signatures) %in% referenceSignature)))
      x2plot <- lapply(x2plot,function(x) x[myorder,])
    }

    #Test for significance
    if (missing(referenceSignature)) {
      myFun2 <- function(x) ifelse(sum(x[c(1,3)])==0,1,binom.test(x=x[1],n=sum(x[c(1,3)]),alternative='two.sided')$p.value)
      pvals <- lapply(x2plot, function(x) apply(x,1,myFun2))
    } else {
      sel <- which(names(signatures) != referenceSignature)
      if (!testUpDown) {
        myFun3 <- function(x) {
          xtab <- sapply(sel,function(y) rbind(x[y,c(1,3)],x[referenceSignature,c(1,3)]),simplify=FALSE)
          names(xtab) <- names(signatures)[sel]
          ans <- unlist(lapply(xtab, function(y) ifelse((any(rowSums(y)==0) | any(colSums(y)==0)),1,chisq.test(y, simulate.p.value=simulate.p.value, B=B)$p.value)))
          return(ans)
        }
        pvals <- lapply(x2plot,myFun3)
      } else {
        myFun4 <- function(x) {
          xtab <- sapply(sel, function(y) matrix(c(x[y,1],sum(x[y,-1]),x[referenceSignature,1],sum(x[referenceSignature,-1])),nrow=2),simplify=FALSE)
          names(xtab) <- names(signatures)[sel]
          pvals.do <- unlist(lapply(xtab,function(x) ifelse((any(rowSums(x)==0) | any(colSums(x)==0)),1,chisq.test(x, simulate.p.value=simulate.p.value, B=B)$p.value)))
          xtab <- sapply(sel, function(y) matrix(c(x[y,3],sum(x[y,-3]),x[referenceSignature,3],sum(x[referenceSignature,-3])),nrow=2),simplify=FALSE)
          names(xtab) <- names(signatures)[sel]
          pvals.up <- unlist(lapply(xtab,function(x) ifelse((any(rowSums(x)==0) | any(colSums(x)==0)),1,chisq.test(x, simulate.p.value=simulate.p.value, B=B)$p.value)))
          ans <- cbind(Down=pvals.do,Up=pvals.up)
          return(ans)
        }
      pvals <- lapply(x2plot,myFun4)
      }
    }

    if (length(signatures)==1) {
      pvals <- as.list(p.adjust(unlist(pvals),method=p.adjust.method))
    } else {
      if (p.adjust.method!='none') pvals <- lapply(pvals,function(x) p.adjust(x,p.adjust.method))
    }
    
    #Plot
    z <- lapply(x2plot,function(x) as.matrix(100*t(x/rowSums(x))[-2,]))
    if (missing(ylim)) ylim <- c(0,1.2*max(sapply(z,'max')))
    if (missing(ylab)) ylab <- "% differentially expressed genes"
    cex <- ifelse(testUpDown,.8,1)
    if (length(signatures)==1) {
      pval2plot <- ifelse(unlist(pvals) < 1e-04, "P<0.0001", paste("P=", round(unlist(pvals), 4), sep = ""))
      myData <- matrix(unlist(z), ncol = length(z)); colnames(myData) <- names(z); rownames(myData) <- c("Down-regulated", "Up-regulated")
      barplot(myData,beside=TRUE,legend=TRUE,ylim=ylim,ylab=ylab,args.legend=list(cex=0.7,box.lty=0),...)
      text((1:length(pval2plot) * 3) - 1,apply(myData,2,max)+1,pval2plot,cex=cex)
    } else {
      for (i in 1:length(x2plot)) {
        ndec <- ifelse(testUpDown,3,4)
        pval2plot <- ifelse(pvals[[i]]<1/10^ndec,paste('P<',1/10^ndec,sep=''),paste('P=',round(pvals[[i]],ndec),sep=''))
        rownames(z[[i]]) <- c('Down-regulated','Up-regulated')
        barplot(z[[i]],beside=TRUE,legend=TRUE,ylim=ylim,ylab=ylab,args.legend=list(cex=0.7,box.lty=0),...)
        if (missing(referenceSignature)) {
          text((1:length(pval2plot) * 3) - 1,apply(z[[i]],2,max)+1,pval2plot,cex=cex)
        } else {
          if (!testUpDown) {
            text((1:length(pval2plot) * 3) + 2,apply(z[[i]],2,max)[-1]+1,pval2plot,cex=cex)
          } else {
            xpos <- (1:length(pval2plot) * 3)
            xpos <- rep(xpos,2); xpos <- xpos[order(xpos)]; xpos <- xpos - c(1.5,.5); xpos <- xpos[-1:-2]
            text(xpos,apply(z[[i]],2,max)[-1]+1,pval2plot,cex=cex)
          }
        }
      }
    }
  }
)

setMethod("barplotSignifSignatures",signature(x="epheno",signatures="character"),
  function (x, signatures, referenceSignature, testUpDown=FALSE, simulate.p.value = FALSE, B = 10^4, p.adjust.method='none', alpha=.05, ylab, ylim=ylim, cex.text=1, ...) {
    signatures <- list(mySign=signatures)
    barplotSignifSignatures(x, signatures, referenceSignature, testUpDown, simulate.p.value, B, p.adjust.method, alpha, ylab, ylim, cex.text,...)
  }
)

setMethod("barplotSignifSignatures",signature(x="epheno",signatures="GeneSetCollection"),
  function (x, signatures, referenceSignature, testUpDown=FALSE, simulate.p.value = FALSE, B = 10^4, p.adjust.method='none', alpha=.05, ylab, ylim=ylim, cex.text=1, ...) {
    signatures <- geneIds(signatures)
    barplotSignifSignatures(x, signatures, referenceSignature, testUpDown, simulate.p.value, B, p.adjust.method, alpha, ylab, ylim, cex.text,...)
  }
)

setMethod("barplotSignifSignatures",signature(x="epheno",signatures="GeneSet"),
  function (x, signatures, referenceSignature, testUpDown=FALSE, simulate.p.value = FALSE, B = 10^4, p.adjust.method='none', alpha=.05, ylab, ylim=ylim, cex.text=1, ...) {
    signatures <- geneIds(signatures)
    barplotSignifSignatures(x, signatures, referenceSignature, testUpDown, simulate.p.value, B, p.adjust.method, alpha, ylab, ylim, cex.text,...)
  }
)


setGeneric("barplotSignatures",function (x, signatures, referenceSignature, alpha=.05, p.adjust.method='none', ylab, cex.text=1, ...) standardGeneric("barplotSignatures"))

setMethod("barplotSignatures",signature(x="epheno",signatures="list"),
  function (x, signatures, referenceSignature, alpha=.05, p.adjust.method='none', ylab, cex.text=1, ...) {
    if (!missing(referenceSignature)) if (!(referenceSignature %in% names(signatures))) stop('referenceSignature not in names(signatures)')
    lensig <- unlist(lapply(signatures,function(x) length(x)))
    if (all(lensig==0)) stop('All signatures are empty. Check the argument signatures.\n')
    if (any(lensig==0)) {
      warning('Removed some signatures with no genes\n')
      signatures <- signatures[lensig>0]
    }
    badsignature <- sapply(signatures,function(y) !(any(y %in% featureNames(x))))
    if (any(badsignature)) stop(paste('Identifiers in signatures',paste(names(signatures)[which(badsignature)],collapse=','),'not in featureNames'))
    if (missing(ylab)) {
      ylab <- ifelse(length(signatures)==1,'log2 Fold Change / Hazard Ratio',c(rep('log2 Fold Change',ncol(getFc(x))),rep('log Hazard Ratio',ncol(getHr(x)))))
    }
    fc2plot <- lapply(phenoNames(x), function(y) lapply(signatures, function(z) logFcHr(x)[featureNames(x) %in% z,phenoNames(x) %in% y]))
    names(fc2plot) <- phenoNames(x)
    if (length(signatures)==1) {
      groups <- factor(rep(names(fc2plot),each=length(fc2plot[[1]][[1]])))
      o <- names(fc2plot)[unlist(lapply(fc2plot,function(x) length(x)))>0]
      barplotCI(unlist(fc2plot),groups=groups,order=o,alpha=alpha,ylab=ylab,...); abline(h=0,lty=2)
      abline(h=0,lty=2)
    } else {
      for (i in 1:ncol(x)) {
        groups <- factor(rep(names(fc2plot[[i]]),unlist(lapply(fc2plot[[i]],function(x) length(x)))))
        o <- names(fc2plot[[i]])[unlist(lapply(fc2plot[[i]],function(x) length(x)))>0]
        myxpos <- barplotCI(unlist(fc2plot[[i]]),groups=groups,order=o,alpha=alpha,ylab=ylab[i],main=names(fc2plot)[i], ...)
        abline(h=0,lty=2)
        if (!missing(referenceSignature)) {
          sigs2test <- which(names(signatures)!=referenceSignature)
          pvals <- sapply(sigs2test, function(y) wilcox.test(fc2plot[[i]][[y]],fc2plot[[i]][[referenceSignature]])$p.value)
          pvals <- p.adjust(pvals,method=p.adjust.method)
          pvals <- ifelse(pvals<.0001,'P<0.0001',paste('P=',round(pvals,4),sep=''))
          text(myxpos[sigs2test],rep(par("usr")[4],length(sigs2test)),pvals,pos=1,cex=cex.text)
        }
      }
    }
  }
)

setMethod("barplotSignatures",signature(x="epheno",signatures="character"),
  function (x, signatures, referenceSignature, alpha=.05, p.adjust.method='none', ylab, cex.text=1, ...) {
    signatures <- list(mySign=signatures)
    barplotSignatures(x, signatures, referenceSignature, alpha, p.adjust.method, ylab, cex.text, ...)
  }
)

setMethod("barplotSignatures",signature(x="epheno",signatures="GeneSetCollection"),
  function (x, signatures, referenceSignature, alpha=.05, p.adjust.method='none', ylab, cex.text=1, ...) {
    tmpsignatures <- geneIds(signatures); names(tmpsignatures) <- names(signatures); signatures <- tmpsignatures
    barplotSignatures(x, signatures, referenceSignature, alpha, p.adjust.method, ylab, cex.text, ...)
  }
)

setMethod("barplotSignatures",signature(x="epheno",signatures="GeneSet"),
  function (x, signatures, referenceSignature, alpha=.05, p.adjust.method='none', ylab, cex.text=1, ...) {
    signatures <- geneIds(signatures)
    barplotSignatures(x, signatures, referenceSignature, alpha, p.adjust.method, ylab, cex.text, ...)
  }
)


setGeneric("boxplotSignatures",function (x, signatures, referenceSignature, outline=FALSE, cex.text=1, ...) standardGeneric("boxplotSignatures"))

setMethod("boxplotSignatures",signature(x="epheno",signatures="list"),
  function (x, signatures, referenceSignature, outline=FALSE, cex.text=1, ...) {
    if (!missing(referenceSignature)) if (!(referenceSignature %in% names(signatures))) stop('referenceSignature not in names(signatures)')
    lensig <- unlist(lapply(signatures,function(x) length(x)))
    if (all(lensig==0)) stop('All signatures are empty. Check the argument signatures.\n')
    if (any(lensig==0)) {
      warning('Removed some signatures with no genes\n')
      signatures <- signatures[lensig>0]
    }
    badsignature <- sapply(signatures,function(y) !(any(y %in% featureNames(x))))
    if (any(badsignature)) stop(paste('Identifiers in signatures',paste(names(signatures)[which(badsignature)],collapse=','),'not in featureNames'))
    for (i in 1:ncol(x)) {
      logfc <- logFcHr(x)
      fc2plot <- vector("list",length(signatures))
      fc2plot <- lapply(signatures,function(y) logfc[featureNames(x) %in% y])
      boxplot(fc2plot, outline=outline, ...); abline(h=0,lty=2)
      if (!missing(referenceSignature)) {
        sigs2test <- which(names(signatures)!=referenceSignature)
        pvals <- double(length(sigs2test))
        pvals <- unlist(lapply(sigs2test, function(x) wilcox.test(fc2plot[[x]],fc2plot[[referenceSignature]])$p.value))
        pvals <- ifelse(pvals<.0001,'P<0.0001',paste('P=',round(pvals,4),sep=''))
        text(sigs2test,rep(par("usr")[4],length(sigs2test)),pvals,pos=1,cex=cex.text)
      }
    }
  }
)
