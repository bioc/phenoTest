sortDragHtmlTable <- function(filename) {
#Adds javascript to html to enable sorting.
# INPUT
# - filename: html file to modify.

#get outputdir from filename
lastSlashPos <- gregexpr(.Platform$file.sep,filename)[[1]][length(gregexpr(.Platform$file.sep,filename)[[1]])]
outputdir <- ifelse(lastSlashPos== -1, getwd(), substr(filename,0,lastSlashPos))
  
#import javascript files
st <- system.file("javascript", "sorttable.js", package = "phenoTest")
dt <- system.file("javascript", "dragtable.js", package = "phenoTest")
file.copy(st,outputdir)
file.copy(dt,outputdir)

#rewrite html file
tmpTxt <- readLines(filename)
tmpTxt[1] <- paste('<script src="sorttable.js"></script>\n',tmpTxt[1])
tmpTxt[1] <- paste('<script src="dragtable.js"></script>\n',tmpTxt[1])
tmpTxt <- sub('TABLE','TABLE class="draggable sortable"',tmpTxt)
writeLines(tmpTxt,con=filename)
}

epheno2html <- function(x,epheno,outputdir,prefix='',genelimit=50,categories=3,withPlots=TRUE,mc.cores=1) {
  require(paste(annotation(x),'.db',sep=''),character.only = TRUE)
  stopifnot(annotation(x)==annotation(epheno))
  stopifnot(identical(featureNames(x),featureNames(epheno))  )
  if (!(categories %in% c(2,3,4))) stop('categories has to be 2, 3 or 4!')
  if (genelimit>nrow(x)) genelimit <- nrow(x)
  stopifnot(file.exists(outputdir))
  id.entrezid <- 'org.' %in% substr(annotation(x),0,4)
    
  export2html <- function(varName,x,epheno,varType,categories,outputdir) {
    #keep only columns that contain our variable
    epheno <-  epheno[,varName]

    #make xls file
    if (id.entrezid) {
      probeid <- rep(NA,nrow(x))
      entrezid <- featureNames(x)
    } else {
      probeid <- featureNames(x)
      envir <- eval(parse(text=paste(annotation(x),'ENTREZID',sep='')))
      entrezid <- as.character(AnnotationDbi::mget(featureNames(x),envir,ifnotfound=NA))
    }
    envir <- eval(parse(text=paste(annotation(x),'SYMBOL',sep='')))
    symbol <- as.character(AnnotationDbi::mget(featureNames(x),envir,ifnotfound=NA))
    envir <- eval(parse(text=paste(annotation(x),'GENENAME',sep='')))
    genename <- as.character(AnnotationDbi::mget(featureNames(x),envir,ifnotfound=NA))
    if (varType=='survival') {
      xout <- data.frame(probeid,entrezid,symbol,genename,getSummaryDif(epheno)[featureNames(x),,drop=FALSE],p.value=as.numeric(getSignif(epheno)[featureNames(x),,drop=FALSE]))
    } else {
      xout <- data.frame(probeid,entrezid,symbol,genename,getMeans(epheno)[featureNames(x),,drop=FALSE],getSummaryDif(epheno)[featureNames(x),,drop=FALSE],p.value=as.numeric(getSignif(epheno)[featureNames(x),,drop=FALSE]))
    }
    if (id.entrezid) xout <- xout[,-1]
    if (approach(epheno)=='bayesian') colnames(xout)[grepl('p.value',colnames(xout))] <- 'postProb'
    write.csv(xout,paste(outputdir,'/',ifelse(prefix=='','',paste(prefix,'_',sep='')),varType,'_',varName,'.csv',sep=''),row.names=FALSE)
    if (approach(epheno)=='frequentist') xout <- xout[order(xout$p.value),] else xout <- xout[order(xout$postProb,decreasing=TRUE),]
    
    #keep amount of genes equal to genelimit
    epheno <- epheno[order(getSignif(epheno),decreasing=(approach(epheno)=='bayesian')),]
    epheno <- epheno[1:genelimit,]
    x <- x[featureNames(epheno),]
    if (id.entrezid) sel <- match(featureNames(epheno),entrezid) else sel <- match(featureNames(epheno),probeid)
    probeid <- probeid[sel]
    entrezid <- entrezid[sel]
    symbol <- symbol[sel]
    genename <- genename[sel]
    if (id.entrezid) xout <- xout[xout$entrezid %in% entrezid,] else xout <- xout[xout$probeid %in% probeid,]
    
    #remove samples with nulls on pData
    x <- x[,!is.na(pData(x)[,varName])]
    
    #make plots
    if (withPlots) {
      if (!file.exists(paste(outputdir,'phenoPlots',sep='/'))) dir.create(paste(outputdir,'phenoPlots',sep='/')) 
      geneCategories <- t(apply(exprs(x), 1, function(x) as.numeric(cut2(x,g=categories)))); colnames(geneCategories) <- colnames(exprs(x))
      #make plots
      dirOut <- paste(outputdir,'/phenoPlots/',varName,sep=''); if (!file.exists(dirOut)) dir.create(dirOut)
      if (varType=='survival') {
        varName.time <- as.character(unique(pData(epheno)[varName,'survTime']))
        varName.event <- varName
        s <- Surv(as.numeric(pData(x)[,varName.time]),pData(x)[,varName.event])
        fit <- apply(geneCategories, 1, function(x) survfit(s~x))
        if (categories==2) uniqueCat <- c('Low','High') else if (categories==3) uniqueCat <- c('Low','Medium','High') else if (categories==4) uniqueCat <- c('Low','Medium low','Medium high','High')
        myFun <- function(i) {
          png(paste(dirOut,'/survival_',gsub('/','_',featureNames(x)[i]),'.png',sep=''))
          if (categories==2) {
            plot(fit[[i]],ylab='Survival',xlab=paste('Time to ',varName,sep=''), col=c('green','red'))
            legend("bottomleft",uniqueCat,col=c('green','red'),lty=1)
          } else if (categories==3) {
            plot(fit[[i]],ylab='Survival',xlab=paste('Time to ',varName,sep=''), col=c('green','black','red'))
            legend("bottomleft",uniqueCat,col=c('green','black','red'),lty=1)
          } else if (categories==4) {
            plot(fit[[i]],ylab='Survival',xlab=paste('Time to ',varName,sep=''), col=c('green','blue','yellow','red'))
            legend("bottomleft",uniqueCat,col=c('green','blue','yellow','red'),lty=1)
          }
          title(paste(ifelse(approach(epheno)=='frequentist','p-value','postProb'),ifelse(getSignif(epheno)[i,]<.0001,'<0.0001',paste('=',round(getSignif(epheno)[i,],4),sep='')),sep=''))
          dev.off()
          png(paste(dirOut,'/smooth_',gsub('/','_',featureNames(x)[i]),'.png',sep=''))
          dummy <- try(smoothCoxph(pData(x)[,varName.time],pData(x)[,varName.event],as.numeric(exprs(x[i,]))),silent=T)
          dev.off()
        }
        lapply(as.list(1:nrow(x)),myFun)
      }
      myFun <- function(i) {
        png(paste(dirOut,'/plot_',gsub('/','_',featureNames(x)[i]),'.png',sep=''))
        if (varType=='survival') {
          if (ncol(x)>20) {
            boxplot(as.numeric(exprs(x[i,])) ~ as.factor(pData(x)[,varName.event]),ylab='Expression level',xlab=varName.time)
          } else {
            dotchart(as.numeric(exprs(x[i,])),as.factor(pData(x)[,varName.event]),xlab='Expression level',ylab=varName.time,color=as.numeric(as.factor(pData(x)[,varName.event]))+1,groups=as.factor(pData(x)[,varName.event]))
          }
        } else if (varType=='continuous') {
          mylevels <- levels(pData(epheno)[pData(epheno)$phenoName==varName,'meanLabel'])
          mylevels <- gsub('\\[','',gsub('\\]','',gsub('\\)','',mylevels)))
          mylevels <- strsplit(mylevels,',')
          mylevels <- lapply(mylevels,as.numeric)
          myContPheno <- pData(x)[,varName]
          myCateg <- vector('numeric',length=length(myContPheno))
          for (j in 1:length(mylevels)) {
            sel <- myContPheno >= mylevels[[j]][1] & myContPheno <= mylevels[[j]][2]
            myCateg[sel] <- j
          }
          myCateg <- factor(myCateg,labels=levels(pData(epheno)[pData(epheno)$phenoName==varName,'meanLabel']))
          if (ncol(x)>20) {
            boxplot(as.numeric(exprs(x[i,])) ~ myCateg,ylab='Expression level',xlab=varName)
          } else {
            dotchart(as.numeric(exprs(x[i,])),as.factor(myCateg),xlab='Expression level',ylab=varName,color=as.numeric(as.factor(myCateg))+1,groups=as.factor(myCateg))
          }
        } else {
          if (ncol(x)>20) {
            boxplot(as.numeric(exprs(x[i,])) ~ as.factor(pData(x)[,varName]),ylab='Expression level',xlab=varName)
          } else {
            dotchart(as.numeric(exprs(x[i,])),as.factor(pData(x)[,varName]),xlab='Expression level',ylab=varName,color=as.numeric(as.factor(pData(x)[,varName]))+1,groups=as.factor(pData(x)[,varName]))
          }
        }
        title(paste(ifelse(approach(epheno)=='frequentist','p-value','postProb'),ifelse(getSignif(epheno)[i,]<.0001,'<0.0001',paste('=',round(getSignif(epheno)[i,],4),sep='')),sep=''))
        dev.off()
      }
      lapply(as.list(1:nrow(xout)),myFun)
    }

    #make html file
    xout[is.na(xout)] <- '---'
    if (withPlots) {
      if (varType=='survival') {
        txt <- rep('Open',nrow(xout))
        xout <- data.frame(xout[,c(1:(ifelse(id.entrezid,3,4)))],Plot=txt,Kaplan=txt,SmoothedHR=txt,xout[,c((ifelse(id.entrezid,4,5)):ncol(xout))])
      } else {
        xout <- data.frame(xout[,c(1:(ifelse(id.entrezid,3,4)))],plot=rep('Open',nrow(xout)),xout[,c((ifelse(id.entrezid,4,5)):ncol(xout))])
      }
    }
    links <- tiny.pic <- vector('list',length=ncol(xout))
    if (!id.entrezid) links[[1]] <- paste('https://www.affymetrix.com/LinkServlet?probeset=',xout$probeid,sep='')
    links[[ifelse(id.entrezid,1,2)]] <- paste('http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids=',xout$entrezid,sep='')
    links[[ifelse(id.entrezid,2,3)]] <- paste('http://www.genecards.org/index.php?path=/Search/keyword/',xout$symbol,sep='')
    if (withPlots) {
      tiny.pic[[ifelse(id.entrezid,4,5)]] <- links[[ifelse(id.entrezid,4,5)]] <- paste('phenoPlots/',varName,'/plot_',gsub('/','_',featureNames(x)),'.png',sep='')
      if (varType=='survival') {
        tiny.pic[[ifelse(id.entrezid,5,6)]] <- links[[ifelse(id.entrezid,5,6)]] <- paste('phenoPlots/',varName,'/survival_',gsub('/','_',featureNames(x)),'.png',sep='')
        tiny.pic[[ifelse(id.entrezid,6,7)]] <- links[[ifelse(id.entrezid,6,7)]] <- paste('phenoPlots/',varName,'/smooth_',gsub('/','_',featureNames(x)),'.png',sep='')
      }
    }
    file <- paste(outputdir,'/',prefix,'_',varType,'_',varName,'.html',sep='')
    title <- paste(varName,' (',varType,' variable)',sep='')
    write.html(xout,file=file,links=links,tiny.pic=tiny.pic,title=title)
  }
   
  myFun <- function(varName) {
    varType <- as.character(phenoClass(epheno[,varName]))
    export2html(varName=varName,x=x,epheno=epheno,varType=varType,categories=categories,outputdir=outputdir)
  }
  varName <- as.list(phenoNames(epheno))
  if (mc.cores==1) {
    dummy <- lapply(varName, myFun)
  } else {
    if ('parallel' %in% loadedNamespaces()) {
      dummy <- parallel::mclapply(varName, myFun, mc.cores=mc.cores)
    } else stop('parallel library has not been loaded!')
  }

  ##TO DEBUG USE for INSTEAD OF mclapply
  #for (i in 1:ncol(epheno)) {
  #  varName <- phenoNames(epheno)
  #  varType <- as.character(phenoClass(epheno)[phenoClass(epheno)==varName,'phenoClass'])
  #  export2html(varName=varName,x=x,epheno=epheno,varType=varType,categories=categories,outputdir=outputdir)
  #}
}
