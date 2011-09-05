rowVar <- function(x,...) { return((rowMeans(x^2,...)-rowMeans(x,...)^2)*ncol(x)/(ncol(x)-1)) }
rowSD <- function(x,...) { return(sqrt(rowVar(x,...))) }

#
# FUNCTION THAT CREATES THE EPHENO OBJECT
#
ExpressionPhenoTest <- function (x,vars2test,adjustVars,p.adjust.method='BH',continuousCategories=3,mc.cores=1) {

mycoxph <- function(x,eset,coxSurvEvent,coxSurvTime,adjustVarsTxt) {
  val <- try({
    if (adjustVarsTxt=='') myFormula <- "Surv(coxSurvTime, coxSurvEvent) ~ I(as.numeric(x))" else myFormula <- paste("Surv(coxSurvTime, coxSurvEvent) ~ I(as.numeric(x))",adjustVarsTxt,sep=' + ')
    coxph1 <- coxph(eval(parse(text=myFormula)))
    eCoef <- exp((abs(coef(coxph1)[1])))*sign(coef(coxph1)[1])
    summaryDif <- eCoef
    pval <- summary(coxph1)$logtest["pvalue"]
  },silent=TRUE)
  if (inherits(val, "try-error")) {
    summaryDif <- NA; pval <- NA
#    warning('Some errors happened analyzing survival vars. NAs were produced.')
  }
  return(list(summaryDif,pval))
}

if (class(x) != 'ExpressionSet') stop('x must be of class ExpressionSet')
if (!is.numeric(exprs(x))) stop('exprs(x) must be numeric. Maybe it is stored as a data.frame?')
if (!is.numeric(mc.cores)) stop('mc.cores must be numeric')

if (!is.null(vars2test$continuous)) {
  if (sum(vars2test$continuous %in% names(pData(x))) != length(vars2test$continuous)) stop("Invalid variable names in vars2test$continuous")
  if (!any(!is.na(pData(x)[,vars2test$continuous]))) stop('Some variable in vars2test$continuous contains only missing values.')
}
if (!is.null(vars2test$ordinal)) {
  if (sum(vars2test$ordinal %in% names(pData(x))) != length(vars2test$ordinal)) stop("Invalid variable names in vars2test$ordinal")
  if (!any(!is.na(pData(x)[,vars2test$ordinal]))) stop('Some variable in vars2test$ordinal contains only missing values.')
}
if (!is.null(vars2test$categorical)) {
  if (sum(vars2test$categorical %in% names(pData(x))) != length(vars2test$categorical)) stop("Invalid variable names in vars2test$categorical")
  if (!any(!is.na(pData(x)[,vars2test$categorical]))) stop('Some variable in vars2test$categorical contains only missing values.')
}
if (!is.null(vars2test$survival)) {
  if (any(!(as.character(vars2test$survival) %in% names(pData(x))))) stop('Invalid variable names in vars2test$survival')  
  if (ncol(vars2test$survival)!=2) stop("vars2test$survival must have 2 columns")
  if (sum(colnames(vars2test$survival) %in% c('event','time'))<2) stop("colnames(vars2test$survival) must be 'event' and 'time' (in any order)")
}

if (continuousCategories<2) stop('continuousCategories can not be less than 2')
if (continuousCategories>4) stop('continuousCategories can not be bigger than 4')
if (!missing(adjustVars)) {
  if (class(adjustVars)!='character') stop('adjustVars has to be of class character')
  if (!all(adjustVars %in% colnames(pData(x)))) stop('adjustVars is not in colnames(pData(x))')
}

phenoName <- NULL; phenoClass <- NULL; phenoClass.pval <- NULL; phenoType <- NULL; meanLabel <- NULL; survTime <- NULL

summaryDif <- pval <- vector("list", 4)
names(summaryDif) <- names(pval) <- c("cont", "categ", "surv", "ordi")
if (!is.null(vars2test$continuous)) {
  numContVars <- (continuousCategories*2)-1
  summaryDif$cont <- matrix(NA,nrow=nrow(x),ncol=(length(vars2test$continuous)*numContVars))
  pval$cont <- matrix(NA,nrow=nrow(x),ncol=length(vars2test$continuous))

  #epheno names
  phenoName <- c(phenoName,as.character(sapply(vars2test$continuous, function(x) rep(x,numContVars))))
  phenoClass <- c(phenoClass,rep('continuous',numContVars*length(vars2test$continuous)))
  phenoType <- c(phenoType,as.character(unlist(sapply(vars2test$continuous, function(x) c(rep('mean',continuousCategories),rep('summaryDif',continuousCategories-1))))))
  survTime <- c(survTime,rep(NA,numContVars*length(vars2test$continuous)))

  colnames(pval$cont) <- vars2test$continuous
  mynames <- character(0)
  for (i in 1:length(vars2test$continuous)) {
    cat(paste('\nPerforming analysis for continuous variable ',vars2test$continuous[i],' ',sep=''))
    colsel <- !is.na(pData(x)[,vars2test$continuous[i]])
    tmpVar <- pData(x)[colsel,vars2test$continuous[i]]
    
    if (missing(adjustVars)) {
      myFormula <- "~ tmpVar"
    } else {
      adjustVarsTxt <- paste(paste("pData(x)[colsel,'",adjustVars[!(adjustVars %in% vars2test$continuous[i])],"']",sep='',collapse=' + '))
      myFormula <- paste("~ tmpVar",adjustVarsTxt,sep=' + ')
    }
    design <- model.matrix(eval(parse(text=myFormula)))
    colnames(design)[1:2] <- c('Intercept',vars2test$continuous[i])
    lm1 <- lmFit(x[,colsel],design=design) #select arrays in the design matrix (e.g. NAs excluded)
    eb <- eBayes(contrasts.fit(lm1,contrasts=c(0,1,rep(0,ncol(design)-2))))
    ncateg <- length(unique(pData(x)[,vars2test$continuous[i]]))
    if (any(is.na(pData(x)[,vars2test$continuous[i]]))) ncateg <- ncateg-1
    if (ncateg>continuousCategories) {
      categories <- cut2(pData(x)[,vars2test$continuous[i]],g=continuousCategories)
    } else {
      categories <- factor(pData(x)[,vars2test$continuous[i]])
    }

    #epheno names
    meanLabel <- c(meanLabel,levels(categories),rep(NA,length(levels(categories))-1))
    phenoClass.pval <- c(phenoClass.pval,'continuous')
    
    mynames <- c(mynames,paste(vars2test$continuous[i],levels(categories),sep='.'))
    mynames <- c(mynames,paste(vars2test$continuous[i],'fc',levels(categories)[2:length(levels(categories))],sep='.'))
    sel.mean <- ((i-1)*numContVars) + 1:continuousCategories
    summaryDif$cont[,sel.mean] <- t(apply(exprs(x),1,function(x) sapply(levels(categories),function(y) mean(x[categories %in% y]))))
    numitems <- (2*continuousCategories)-1
    sel.fc <- ((i-1)*numContVars+continuousCategories) + 1:(continuousCategories-1)
    summaryDif$cont[,sel.fc] <- t(apply(summaryDif$cont[,sel.mean],1,function(x) 2^(abs(x[-1]-x[1]))*ifelse(x[-1]-x[1]<0,-1,1)))
    if (length(featureNames(x))==1) { pval$cont[,i] <- summary(eb)$coefficients[2,'Pr(>|t|)'] } #2 if there is an intercept variable
    else { pval$cont[,i] <- eb$p.value[,1] }
  }
  colnames(summaryDif$cont) <- mynames
}

if (!is.null(vars2test$ordinal)) {
  nc <- unlist(apply(data.frame(pData(x)[,vars2test$ordinal]),2,function(x) length(unique(x[!is.na(x)]))))
  names(nc) <- vars2test$ordinal
  fc <- nc-1
  ncCoef <- c(1,1+cumsum(nc))
  fcCoef <- c(0,cumsum(fc))

  #epheno names
  phenoName <- c(phenoName,as.character(unlist(sapply(vars2test$ordinal,function(y) rep(y,2* length(unique(pData(x)[,y][!is.na(pData(x)[,y])])) -1)))))
  phenoType <- c(phenoType,as.character(unlist(sapply(nc, function(x) c(rep('mean',x),rep('summaryDif',x-1))))))
  phenoClass <- c(phenoClass,rep('ordinal',sum(2*nc+-1)))
  survTime <- c(survTime,rep(NA,sum(2*nc+-1)))
  
  summaryDif$ordi <- matrix(NA,nrow=nrow(x),ncol=(sum(nc)+sum(fc)))
  pval$ordi <- matrix(NA,nrow=nrow(x),ncol=length(vars2test$ordinal))
  colnames(pval$ordi) <- vars2test$ordinal
  mynames <- character(0)
  for (i in 1:length(vars2test$ordinal)) {
    cat(paste('\nPerforming analysis for ordinal variable ',vars2test$ordinal[i],' ',sep=''))

    #epheno names
    myLevels <- levels(as.factor(as.character(pData(x)[,vars2test$ordinal[i]])))
    meanLabel <- c(meanLabel,myLevels,rep(NA,length(myLevels)-1))
    phenoClass.pval <- c(phenoClass.pval,'ordinal')

    colsel <- !is.na(pData(x)[,vars2test$ordinal[i]])
    design <- model.matrix(~ -1 + as.factor(pData(x)[colsel,vars2test$ordinal[i]]))
    colnames(design) <- levels(factor(pData(x)[,vars2test$ordinal[i]]))
    lm1 <- lmFit(x[,colsel],design=design) #to get means
    mycontrasts <- rbind(rep(-1,nlevels(tmpVar)-1),diag(nlevels(tmpVar)-1),matrix(0,nrow=ncol(design)-nlevels(tmpVar),ncol=nlevels(tmpVar)-1))
    lm2 <- contrasts.fit(lm1,contrasts=mycontrasts) #to get fc
    design4p <- model.matrix(~ as.numeric(as.factor(pData(x)[colsel,vars2test$ordinal[i]])))
    lm4p <- lmFit(x[,colsel],design=design4p) 
    mycontrasts <- c(0,1)
    eb <- eBayes(contrasts.fit(lm4p,contrasts=mycontrasts)) #to get pvals
    categories <- levels(as.factor(pData(x)[,vars2test$ordinal[i]][colsel]))
    mynames <- c(mynames,paste(vars2test$ordinal[i],categories,sep='.'))
    mynames <- c(mynames,paste(vars2test$ordinal[i],categories[2:length(categories)],'fc',sep='.'))
    summaryDif$ordi[,(ncCoef[i]+fcCoef[i]):(ncCoef[i]+fcCoef[i]+fc[i])] <- coef(lm1) #means
    summaryDif$ordi[,(ncCoef[i]+fcCoef[i]+nc[i]):(ncCoef[i]+fcCoef[i]+nc[i]+(fc[i]-1))] <- (2^abs(coef(lm2)))*ifelse(coef(lm2)<0,-1,1) #fc
    if (length(featureNames(x))==1) {
      pval$ordi[,i] <- summary(eb)$coefficients[2,'Pr(>|t|)']  #2 if there is an intercept variable
    } else {
      pval$ordi[,i] <- as.numeric(eb$p.value)
    }
  }
  colnames(summaryDif$ordi) <- mynames
}

if (!is.null(vars2test$categorical)) {
  nc <- unlist(apply(data.frame(pData(x)[,vars2test$categorical]),2,function(x) return(length(unique(x[!is.na(x) & x!=""])))))
  names(nc) <- vars2test$categorical
  fc <- nc-1
  ncCoef <- c(1,1+cumsum(nc))
  fcCoef <- c(0,cumsum(fc))

  #epheno names
  phenoName <- c(phenoName,as.character(unlist(sapply(vars2test$categorical,function(y) rep(y,2* length(unique(pData(x)[,y][!is.na(pData(x)[,y]) & pData(x)[,y]!=""])) -1)))))
  phenoType <- c(phenoType,as.character(unlist(sapply(nc, function(x) c(rep('mean',x),rep('summaryDif',x-1))))))
  phenoClass <- c(phenoClass,rep('categorical',sum(2*nc-1)))
  survTime <- c(survTime,rep(NA,sum(2*nc+-1)))
    
  summaryDif$categ <- matrix(NA,nrow=nrow(x),ncol=sum(nc)+sum(fc))
  pval$categ <- matrix(NA,nrow=nrow(x),ncol=length(vars2test$categorical)); colnames(pval$categ) <- vars2test$categorical
  mynames <- character(0)
  for (i in 1:length(vars2test$categorical)) {
    cat(paste('\nPerforming analysis for categorical variable ',vars2test$categorical[i],' ',sep=''))

    #epheno names
    myLevels <- levels(as.factor(as.character(pData(x)[,vars2test$categorical[i]])))
    myLevels <- myLevels[myLevels!='']
    meanLabel <- c(meanLabel,myLevels,rep(NA,length(myLevels)-1))
    phenoClass.pval <- c(phenoClass.pval,'categorical')

    colsel <- !is.na(pData(x)[,vars2test$categorical[i]]) & pData(x)[,vars2test$categorical[i]]!=''
    tmpVar <- as.factor(as.character(pData(x)[colsel,vars2test$categorical[i]]))
    if (missing(adjustVars)) {
      myFormula <- "~ -1 + tmpVar"
    } else {
      adjustVarsTxt <- paste(paste("pData(x)[colsel,'",adjustVars[!(adjustVars %in% vars2test$categorical[i])],"']",sep='',collapse=' + '))
      myFormula <- paste("~ -1 + tmpVar",adjustVarsTxt,sep=' + ')
    }
    design <- model.matrix(eval(parse(text=myFormula)))
    lm1 <- lmFit(x[,colsel],design=design) #select arrays in the design matrix (e.g. NAs excluded)
    mycontrasts <- rbind(rep(-1,nlevels(tmpVar)-1),diag(nlevels(tmpVar)-1),matrix(0,nrow=ncol(design)-nlevels(tmpVar),ncol=nlevels(tmpVar)-1))
    eb <- eBayes(contrasts.fit(lm1,contrasts=mycontrasts))
    tmpCoef <- coef(lm1)[,1:(length(levels(tmpVar)))]
    colnames(tmpCoef) <- paste(vars2test$categorical[i],levels(tmpVar),sep='.')
    summaryDif$categ[,(ncCoef[i]+fcCoef[i]):(ncCoef[i]+fcCoef[i]+fc[i])] <- tmpCoef
    pval$categ[,i] <- eb$F.p.value
    mynames <- c(mynames,colnames(tmpCoef))
    mynames <- c(mynames,paste(colnames(tmpCoef)[-1],'fc',sep='.')) 
    refCol <- summaryDif$categ[,ncCoef[i]+fcCoef[i]]
    selCol <- summaryDif$categ[,(ncCoef[i]+fcCoef[i]+1):(ncCoef[i]+fcCoef[i]+ncol(tmpCoef)-1)]
    summaryDif$categ[,(ncCoef[i]+fcCoef[i]+nc[i]):(ncCoef[i]+fcCoef[i]+nc[i]+(fc[i]-1))] <- 2^(abs(selCol-refCol))*ifelse(selCol-refCol<0,-1,1)
  }
  colnames(summaryDif$categ) <- mynames
}

if (!is.null(vars2test$survival)) {
  summaryDif$surv <- pval$surv <- matrix(NA,nrow=nrow(x),ncol=nrow(vars2test$survival))
  colnames(summaryDif$surv) <- paste(as.character(vars2test$survival[,'event']),'HR',sep='.')
  colnames(pval$surv) <- as.character(vars2test$survival[,'event'])
  for (i in 1:nrow(vars2test$survival)) {
    phenoName <- c(phenoName,vars2test$survival[i,"event"])
    phenoClass <- c(phenoClass,"survival")
    phenoClass.pval <- c(phenoClass.pval,'survival')
    phenoType <- c(phenoType,"summaryDif")
    meanLabel <- c(meanLabel,NA)
    survTime <-  c(survTime,vars2test$survival[i,"time"])

    cat(paste('\nPerforming analysis for survival variable ',vars2test$survival[i],' ',sep=''))
    exprsX <- exprs(x)/rowSD(exprs(x),na.rm=TRUE)
    coxSurvEvent <- pData(x)[,vars2test$survival[i,"event"]]
    coxSurvTime <- pData(x)[,vars2test$survival[i,"time"]]
    if (mc.cores==1) {
      exprs.list <- list(exprsX)
    } else {
      exprs.list <- by(exprsX,as.numeric(cut(1:nrow(x),b=mc.cores)),function(x) x,simplify=FALSE)
    }
    if (missing(adjustVars)) adjustVarsTxt <- '' else adjustVarsTxt <- paste(paste("pData(eset)[,'",adjustVars[!(adjustVars %in% vars2test$survival[i])],"']",sep='',collapse=' + '))
    if (mc.cores==1) {
      result.list <- lapply(exprs.list,function(y) do.call('rbind',apply(y,1,function(z) mycoxph(x=z,eset=x,coxSurvEvent,coxSurvTime,adjustVarsTxt))))
    } else {
      if ('multicore' %in% loadedNamespaces()) {
        result.list <- multicore::mclapply(exprs.list,function(y) do.call('rbind',apply(y,1,function(z) mycoxph(x=z,eset=x,coxSurvEvent,coxSurvTime,adjustVarsTxt))),mc.cores=mc.cores)
      } else stop('multicore library has not been loaded!')
    }
    result.matrix <- do.call('rbind',result.list)
    summaryDif$surv[,i] <- as.numeric(result.matrix[,1])
    pval$surv[,i] <- as.numeric(result.matrix[,2])
  }
}

allsummaryDif <- cbind(summaryDif$cont,summaryDif$ordi,summaryDif$categ,summaryDif$surv)
allpvals <- cbind(pval$cont, pval$ordi, pval$categ, pval$surv)
rownames(allsummaryDif) <- rownames(allpvals) <- featureNames(x)
if (p.adjust.method!='none') allpvals <- apply(allpvals,2,function(x) p.adjust(x,p.adjust.method))
cat('\n')

phenoName <- c(phenoName,colnames(allpvals))
phenoClass <- c(phenoClass,phenoClass.pval)
if (class(allpvals)=='numeric') allpvals <- t(as.matrix(allpvals))
phenoType <- c(phenoType,rep('pval',ncol(allpvals)))
meanLabel <- c(meanLabel,rep(NA,ncol(allpvals)))
survTime <-  c(survTime,rep(NA,ncol(allpvals)))

pdata <- data.frame(phenoName=as.factor(phenoName),phenoClass=as.factor(phenoClass),phenoType=as.factor(phenoType),meanLabel=as.factor(meanLabel),survTime=as.character(survTime))

pdata$survTime <- as.character(pdata$survTime)
myName <- pdata[pdata$phenoClass=='survival' & is.na(pdata$survTime),'phenoName']
pdata[pdata$phenoName %in% myName & is.na(pdata$survTime),'survTime'] <- pdata[pdata$phenoName %in% myName & !is.na(pdata$survTime),'survTime']
pdata$survTime <- as.factor(pdata$survTime)

rownames(pdata) <- c(colnames(allsummaryDif),colnames(allpvals))
pdata <- new('AnnotatedDataFrame',pdata)
ans <- new('epheno',exprs=as.matrix(cbind(allsummaryDif,allpvals)),phenoData=pdata,annotation=annotation(x))
ans@p.adjust.method <- p.adjust.method

return(ans)
}


#
# METHODS 4 THE EPHENO OBJECT (THERE ARE MORE OUT OF THIS FILE)
#

setMethod("dim",signature(x="epheno"),
  function (x) {
    ans <- c(nrow(exprs(x)),length(unique(pData(x)$phenoName)))
    names(ans) <- c('Features','Phenotypes')
    return(ans)
  }
)

setGeneric("phenoNames",function (x) standardGeneric("phenoNames"))
setMethod("phenoNames",signature(x="epheno"),
  function (x) {
    return(as.character(unique(pData(x)$phenoName)))
  }
)

setGeneric("p.adjust.method",function (x) standardGeneric("p.adjust.method"))
setMethod("p.adjust.method",signature(x="epheno"),
  function (x) {
    return(x@p.adjust.method)
  }
)

setGeneric("phenoClass",function (x) standardGeneric("phenoClass"))
setMethod("phenoClass",signature(x="epheno"),
  function (x) {
    tmp <- unique(pData(x)[,c('phenoName','phenoClass')])
    ans <- as.character(tmp[,2])
    names(ans) <- tmp[,1]
    return(ans)
  }
)

setMethod("show",signature(object="epheno"),
  function (object) {
    cat("Object of class 'epheno'\n")
    if (nrow(object)>3) {
      cat(paste("featureNames: ",paste(featureNames(object)[1:3],collapse=", ")," ... (",nrow(object),") feature(s)\n",sep=""))
    } else {
      cat(paste("featureNames: ",paste(head(featureNames(object)),collapse=", "),".  (",nrow(object),") feature(s)\n",sep=""))
    }
    if (ncol(object)>3) {
      cat(paste("phenoNames: ",paste(phenoNames(object)[1:3],collapse=', ')," ... (",ncol(object),") phenotype(s)\n",sep=''))
    } else {
      cat(paste("phenoNames: ",paste(head(phenoNames(object)),collapse=', '),". (",ncol(object),") phenotype(s)\n",sep=''))
    }
    cat("P-value adjustment method: ",p.adjust.method(object),'\n')
    cat("Annotation:",annotation(object),'\n')
#    cat("phenotypes variable class(es):\n")
#    show(phenoClass(object))
    cat("\nType \"showMethods(classes='epheno')\" for a list of ALL methods\n")
  }
)

setMethod("[",signature(x="epheno"),
  function (x, i, j) {
    if (!missing(i)) {
      if (class(i)=="character") {
        sel <- featureNames(x) %in% i
        exprs(x) <- exprs(x)[sel,,drop=FALSE]
      } else if (class(i)=="numeric" | class(i)=="integer") {
        if (any(c(i<1, i>nrow(x)))) stop('Error: subscript out of bounds')
        exprs(x) <- exprs(x)[i,,drop=FALSE]
      } else if (class(i)=="logical") {
        exprs(x) <- x[which(i),,drop=FALSE]
      }
    }
    if (!missing(j)) {
      if (class(j)=="character") {
        sel <- pData(x)$phenoName %in% j
      } else if (class(j)=="numeric" | class(j)=="integer") {
        if (any(c(j<1, j>length(phenoNames(x))))) stop('Error: subscript out of bounds')
        sel <- pData(x)$phenoName %in% phenoNames(x)[j]
      } else if (class(j)=="logical") {
        if (length(j)!=length(phenoNames(x))) stop('Error: subscript out of bounds')
        sel <- pData(x)$phenoName %in% phenoNames(x)[j]
      }
      pData(x) <- pData(x)[sel,,drop=FALSE]
      exprs(x) <- exprs(x)[,sel,drop=FALSE]
    }
    return(x)
  }
)

setGeneric("getMeans",function (x) standardGeneric("getMeans"))
setMethod("getMeans",signature(x="epheno"),
  function (x) {
    sel <- pData(x)$phenoType=="mean"
    return(exprs(x)[,sel,drop=FALSE])
  }
)

setGeneric("getPvals",function (x) standardGeneric("getPvals"))
setMethod("getPvals",signature(x="epheno"),
  function (x) {
    sel <- pData(x)$phenoType=="pval"
    return(exprs(x)[,sel,drop=FALSE])
  }
)

setGeneric("getSummaryDif",function (x) standardGeneric("getSummaryDif"))
setMethod("getSummaryDif",signature(x="epheno"),
  function (x) {
    sel <- pData(x)$phenoType=="summaryDif"
    return(exprs(x)[,sel,drop=FALSE])
  }
)

setGeneric("getFc",function (x) standardGeneric("getFc"))
setMethod("getFc",signature(x="epheno"),
  function (x) {
    sel <- pData(x)$phenoType=="summaryDif" & pData(x)$phenoClass!="survival"
    return(exprs(x)[,sel,drop=FALSE])
  }
)

setGeneric("getHr",function (x) standardGeneric("getHr"))
setMethod("getHr",signature(x="epheno"),
  function (x) {
    sel <- pData(x)$phenoType=="summaryDif" & pData(x)$phenoClass=="survival"
    return(exprs(x)[,sel,drop=FALSE])
  }
)

setGeneric("logFcHr",function (x) standardGeneric("logFcHr"))
setMethod("logFcHr",signature(x="epheno"),
  function (x) {
    fc <- getFc(x); fc <- log(abs(fc),2)*sign(fc)
    hr <- getHr(x); hr <- log(abs(hr),exp(1))*sign(hr) 
    ans <- cbind(fc,hr)
    return(ans)
  }
)

setGeneric("export2CSV",function(x, file, row.names=FALSE, ...) standardGeneric("export2CSV"))
setMethod("export2CSV",signature(x="epheno"), function(x, file, row.names=FALSE, ...) {
pvals <- getPvals(x)
colnames(pvals) <- paste(colnames(pvals),'.Pval',sep='')
if (is.null(fData(x))) {
  ans <- data.frame(featureNames=featureNames(x), getMeans(x), getSummaryDif(x), pvals)
} else {
  ans <- data.frame(fData(x), getMeans(x), getSummaryDif(x), pvals)
}
write.csv(ans, file=file, row.names=row.names, ...)
}
)

setGeneric("pAdjust",function(x, method='BH') standardGeneric("pAdjust"))
setMethod("pAdjust",signature(x="epheno"), function(x, method='BH') {
sel <- which(pData(x)$phenoType=="pval")
exprs(x)[,sel] <- p.adjust(exprs(x)[,sel],method=method)
x@p.adjust.method <- method
return(x)
}
)

setGeneric("getVars2test",function(x) standardGeneric("getVars2test"))
setMethod("getVars2test",signature(x="epheno"), function(x) {
ans <- list()
if (any(pData(x)$phenoClass=='continuous')) ans$continuous <- unique(as.character(pData(x)[pData(x)$phenoClass=='continuous','phenoName']))
if (any(pData(x)$phenoClass=='ordinal')) ans$ordinal <- unique(as.character(pData(x)[pData(x)$phenoClass=='ordinal','phenoName']))
if (any(pData(x)$phenoClass=='categorical')) ans$categorical <- unique(as.character(pData(x)[pData(x)$phenoClass=='categorical','phenoName']))
if (any(pData(x)$phenoClass=='survival')) {
  survival <- as.matrix(unique(pData(x)[pData(x)$phenoClass=='survival',c('phenoName','survTime')]))
  colnames(survival) <- colnames(survival) <- c('event','time')
  rownames(survival) <- NULL
  ans$survival <- survival
}
return(ans)
}
)
