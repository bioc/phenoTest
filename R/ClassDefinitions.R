setClass('epheno',contains='ExpressionSet',representation=representation(p.adjust.method="character"))
setValidity("epheno", function(object){
  msg <- NULL
  pd <- pData(object)
  if (ncol(pData(object))==0) {
    msg <- 'You can not create an epheno object without pData'
  } else {
    if (sum(c('phenoName','phenoClass','phenoType','meanLabel','survTime') %in% colnames(pd))<5) {
      msg <- 'pData colnames have to be phenoName,phenoClass,phenoType,meanLabel and survTime'
    } else {
      if (any(pd[,'phenoClass']=='survival' & is.na(pd[,'survTime']))) {
        msg <- 'Survival variables must have time variables'
      } else {
        if (any(pd[,'phenoType']=='mean' & is.na(pd[,'meanLabel']))) msg <- 'Survival variables must have time variables'
      }
    }
  }
  if (is.null(msg)) TRUE else msg
})

setClass("gseaSignatures",contains="list",representation(es="numeric",es.sim="numeric",signature="numeric"))
setClass("gseaSignaturesSign",contains="list",representation(gseaSignatures="gseaSignatures",fc.hr="character",s="logical",test='character'))
setClass("gseaSignaturesVar",contains="list",representation(gseaSignatures="gseaSignaturesSign"))
setClass("gseaSignificanceSign",contains="list",representation(gseaSignificance="matrix",p.adjust.method="character"))
setClass("gseaSignificanceVar",contains="list",representation(gseaSignificance="gseaSignificanceSign"))
