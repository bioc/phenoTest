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

getQuery4LOCALFILE <- function (ids) {
#Get directory for LOCAL FILE
    blanks <- ids == "&nbsp;"
    out <- paste("LOCALFILE_",ids,".html",sep = "")
    out <- sub(":","",out) #remove strange characters from file name
    out[blanks] <- "&nbsp;"
    return(out)
}
getQuery4miRBase <- function (ids) {
#Get directory for miRBase
    blanks <- ids == "&nbsp;"
    out <- paste("http://microrna.sanger.ac.uk/cgi-bin/sequences/query.pl?terms=",ids,sep='')
    out[blanks] <- "&nbsp"
    return(out)
}
getQuery4BROAD <- function (ids) {
#Get URL for BROAD gene set
    blanks <- ids == "&nbsp;"
    out <- paste("http://www.broad.mit.edu/gsea/msigdb/geneset_page.jsp?geneSetName=", ids, sep = "")
    out[blanks] <- "&nbsp;"
    return(out)
}
getQuery4MIRANDA <- function (ids) {
#Get directory for miRBase
    blanks <- ids == "&nbsp;"
    out <- paste("http://microrna.sanger.ac.uk/cgi-bin/targets/v5/hit_list.pl?genome_id=native&mirna_id=",ids,sep='')
    out[blanks] <- "&nbsp"
    return(out)
}
getQuery4KEGG <- function (ids) {
#Get URL for KEGG gene set
    blanks <- ids == "&nbsp;"
    ids <- gsub(" ","0",format(ids,width=5,justify='right'))  #zero-fill
    out <- paste("http://www.genome.jp/dbget-bin/www_bget?pathway+ko",ids,sep="")
    out[blanks] <- "&nbsp;"
    return(out)
}

epheno2html <- function(x,epheno,outputdir,prefix='',genelimit=50,categories=3,withPlots=TRUE,id.entrezid=FALSE,organism='human',homol.symbol,homol.genename,mc.cores=1) {
  if (!id.entrezid) require(paste(annotation(x),'.db',sep=''),character.only = TRUE)
  if (!(categories %in% c(2,3,4))) stop('categories has to be 2, 3 or 4!')
  if (genelimit>nrow(x)) genelimit <- nrow(x)
   
  getHomolSymbol <- function(organism='human',entrezid) {
    if (organism=='human') mart <- useMart("ensembl","hsapiens_gene_ensembl") else if (organism=='mouse') mart <- useMart("ensembl","mmusculus_gene_ensembl") else stop('Only human and mouse organisms are implemented!')
    homol.symbol <- getLDS(attributes = c("entrezgene"),mart=mart,attributesL=c("external_gene_id"),martL=mart,filters="entrezgene",values=entrezid)
    return(homol.symbol)
  }
   
  getHomolGenename <- function(organism='human',entrezid) {
    if (organism=='human') mart <- useMart("ensembl","hsapiens_gene_ensembl") else if (organism=='mouse') mart <- useMart("ensembl","mmusculus_gene_ensembl") else stop('Only human and mouse organisms are implemented!')
    homol.genename <- getLDS(attributes = c("entrezgene"),mart=mart,attributesL=c("description"),martL=mart,filters="entrezgene",values=entrezid)
    return(homol.genename)
  }
   
  entrezid.anotation <- function(entrezid,homol.symbol,homol.genename) {
    symbol <- homol.symbol[match(entrezid,homol.symbol[,1]),2]
    genename <- homol.genename[match(entrezid,homol.genename[,1]),2]
    genename <- substr(genename,0,unlist(gregexpr('\\[',genename))-1)
    ans <- data.frame(symbol,genename); rownames(ans) <- entrezid
    return(ans)
  }
   
  export2html <- function(varName,x,epheno,varType,categories,outputdir,homol.symbol,homol.genename) {
    #keep only columns that contain our variable
    epheno <-  epheno[,varName]
    
    #make xls file
    if (id.entrezid) {
      tmp <- entrezid.anotation(featureNames(x),homol.symbol=homol.symbol,homol.genename=homol.genename)
      if (varType=='survival') {
        xout <- data.frame(entrezid=featureNames(x),tmp,getSummaryDif(epheno),p.value=as.numeric(getPvals(epheno)))
      } else {
        xout <- data.frame(entrezid=featureNames(x),tmp,getMeans(epheno),getSummaryDif(epheno),p.value=as.numeric(getPvals(epheno)))
      }
    } else {
      entrezid <- unlist(AnnotationDbi::mget(featureNames(x),AnnotationDbi::get(paste(annotation(x),"ENTREZID", sep = "")),ifnotfound=NA))
      symbol <- unlist(AnnotationDbi::mget(featureNames(x),AnnotationDbi::get(paste(annotation(x),"SYMBOL", sep = "")),ifnotfound=NA))
      genename <- unlist(AnnotationDbi::mget(featureNames(x),AnnotationDbi::get(paste(annotation(x),"GENENAME", sep = "")),ifnotfound=NA))
      if (varType=='survival') {
        xout <- data.frame(probeid=featureNames(x),entrezid,symbol,genename,getSummaryDif(epheno),p.value=as.numeric(getPvals(epheno)))
      } else {
        xout <- data.frame(probeid=featureNames(x),entrezid,symbol,genename,getMeans(epheno),getSummaryDif(epheno),p.value=as.numeric(getPvals(epheno)))
      }
    }
    xout <- xout[order(xout$p.value),]
    write.csv(xout,paste(outputdir,'/',prefix,'_',varType,'_',varName,'.csv',sep=''),row.names=FALSE)
    
    #keep amount of genes equal to genelimit
    epheno <- epheno[order(getPvals(epheno)),];  epheno <- epheno[1:genelimit,]
    x <- x[featureNames(epheno),]
    if (id.entrezid) tmp <- tmp[featureNames(epheno),]
    
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
          pdf(paste(dirOut,'/survival_',gsub('/','_',featureNames(x)[i]),'.pdf',sep=''))
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
          title(paste('p-value',ifelse(getPvals(epheno)[i,]<.0001,'<0.0001',paste('=',round(getPvals(epheno)[i,],4),sep='')),sep=''))
          dev.off()
          pdf(paste(dirOut,'/smooth_',gsub('/','_',featureNames(x)[i]),'.pdf',sep=''))
          smoothCoxph(pData(x)[,varName.time],pData(x)[,varName.event],as.numeric(exprs(x[i,])))
          dev.off()
        }
        lapply(as.list(1:nrow(x)),myFun)
      }
      myFun <- function(i) {
        pdf(paste(dirOut,'/boxplot_',gsub('/','_',featureNames(x)[i]),'.pdf',sep=''))
        if (varType=='survival') {
          boxplot(as.numeric(exprs(x[i,])) ~ as.factor(pData(x)[,varName.event]),ylab='Expression level',xlab=varName.time)
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
          boxplot(as.numeric(exprs(x[i,])) ~ myCateg,ylab='Expression level',xlab=varName)
        } else {
          boxplot(as.numeric(exprs(x[i,])) ~ as.factor(pData(x)[,varName]),ylab='Expression level',xlab=varName)
        }
        title(paste('p-value',ifelse(getPvals(epheno)[i,]<.0001,'<0.0001',paste('=',round(getPvals(epheno)[i,],4),sep='')),sep=''))
        dev.off()
      }
      lapply(as.list(1:nrow(x)),myFun)
    }
   
    #make html file
    if (id.entrezid) {
      symbol <- tmp[,'symbol']; genename <- tmp[,'genename']
      #columns with links
      genelist <- list(entrezid=featureNames(x)); head <- c('Entrez ID'); repository <- list("en")
    } else {
      entrezid <- unlist(AnnotationDbi::mget(featureNames(x),AnnotationDbi::get(paste(annotation(x),"ENTREZID", sep = "")),ifnotfound=NA))
      symbol <- unlist(AnnotationDbi::mget(featureNames(x),AnnotationDbi::get(paste(annotation(x),"SYMBOL", sep = "")),ifnotfound=NA))
      genename <- unlist(AnnotationDbi::mget(featureNames(x),AnnotationDbi::get(paste(annotation(x),"GENENAME", sep = "")),ifnotfound=NA))
      #columns with links
      genelist <- list(probeids=featureNames(x)); head <- c('Probe ID'); repository <- list("affy")
      genelist$entrez <- entrezid; head <- c(head,'Entrez ID'); repository <- c(repository,"en")
    }
    if (withPlots) {
      genelist$linkText  <- paste(varName,'/boxplot_',gsub('/','_',featureNames(x)),sep=''); head <- c(head,'Box plot'); repository <- c(repository,"phenoplot")
      if (varType=='survival') {
        genelist$linkText2  <- paste(varName,'/survival_',gsub('/','_',featureNames(x)),sep=''); head <- c(head,'Survival plot'); repository <- c(repository,"phenoplot")
        genelist$linkText3  <- paste(varName,'/smooth_',gsub('/','_',featureNames(x)),sep=''); head <- c(head,'Smoothed HR'); repository <- c(repository,"phenoplot")
      }
    }
    #columns without live links
    if (varType=='survival') {
      othernames <- data.frame(symbol,genename,getSummaryDif(epheno),p.value=getPvals(epheno))
      head <- c(head,'Symbol','Gene Name',colnames(getSummaryDif(epheno)),'p.value')
    } else {
      othernames <- data.frame(symbol,genename,getMeans(epheno),getSummaryDif(epheno),p.value=getPvals(epheno))
      head <- c(head,'Symbol','Gene Name',colnames(getMeans(epheno)),colnames(getSummaryDif(epheno)),'p.value')
    }
    #create html file
    phenoTest::htmlpage(genelist=genelist,filename=paste(outputdir,'/',prefix,varType,'_',varName,'.html',sep=''),title=paste(varName,' (',varType,' variable)',sep=''),othernames=othernames,table.head=head,repository=repository)
    sortDragHtmlTable(paste(outputdir,'/',prefix,varType,'_',varName,'.html',sep=''))
  }
   
  #get homol.symbol and homol.genename if necessary
  if (id.entrezid & (missing(homol.symbol) | missing(homol.genename))) {
    homol.symbol <- getHomolSymbol(organism=organism,entrezid=featureNames(x))
    homol.genename <- getHomolGenename(organism=organism,entrezid=featureNames(x))
  }
   
  myFun <- function(varName) {
    varType <- as.character(phenoClass(epheno[,varName]))
    export2html(varName=varName,x=x,epheno=epheno,varType=varType,categories=categories,outputdir=outputdir,homol.symbol=homol.symbol,homol.genename=homol.genename)
  }
  varName <- as.list(phenoNames(epheno))
  if (mc.cores==1) {
    dummy <- lapply(varName, myFun)
  } else {
    if ('multicore' %in% loadedNamespaces()) {
      dummy <- multicore::mclapply(varName, myFun, mc.cores=mc.cores)
    } else stop('multicore library has not been loaded!')
  }

  ##TO DEBUG USE for INSTEAD OF mclapply
  #for (i in 1:ncol(epheno)) {
  #  varName <- phenoNames(epheno)
  #  varType <- as.character(phenoClass(epheno)[phenoClass(epheno)==varName,'phenoClass'])
  #  export2html(varName=varName,x=x,epheno=epheno,varType=varType,categories=categories,outputdir=outputdir,homol.symbol=homol.symbol,homol.genename=homol.genename)
  #}
}
