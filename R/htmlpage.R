getQuery4PHENOPLOT <- function (ids) {
#Get directory for LOCAL FILE
    blanks <- ids == "&nbsp;"
    path <- paste("phenoPlots",ids,sep=.Platform$file.sep)
    out <- paste(path,".pdf",sep = "")
    out <- sub(":","",out) #remove strange characters from file name
    out[blanks] <- "&nbsp;"
    return(out)
}

getQueryLink <- function (ids, repository = "ug") {
#Redefines getQueryLink function from 'annotate' to add localfile repository
    switch(tolower(repository), ug = return(annotate:::getQuery4UG(ids)), 
        ll = return(annotate:::getQuery4LL(ids)), affy = return(annotate:::getQuery4Affy(ids)), 
        gb = return(annotate:::getQuery4GB(ids)), sp = return(annotate:::getQuery4SP(ids)), 
        omim = return(annotate:::getQuery4OMIM(ids)), fb = return(annotate:::getQuery4FB(ids)), 
        en = return(annotate:::getQuery4EN(ids)), tr = return(annotate:::getQuery4TR(ids)), 
        go = return(annotate:::getQuery4GO(ids)),
        phenoplot = return(getQuery4PHENOPLOT(ids)), #line added           
        localfile = return(getQuery4LOCALFILE(ids)), #line added
        mirbase = return(getQuery4miRBase(ids)), #line added
        broad = return(getQuery4BROAD(ids)), #line added
        kegg = return(getQuery4KEGG(ids)), #line added
        miranda = return(getQuery4MIRANDA(ids)), #line added                      
        stop("Unknown repository name"))
}

getCells <- function (ids, repository = "ug") {
#introduce exception 4 localfile repository  
    if (is.list(ids)) {
        out <- vector()
        temp <- lapply(ids, getQueryLink, repository = repository)
        for (i in seq(along = ids)) {
            if (temp[i] != "&nbsp;") 
                out[i] <- paste("<P><A HREF=\"", temp[[i]], "\">", 
                  ids[[i]], "</A></P>", sep = "", collapse = "")
            else out[i] <- temp[i]
        }
    }
    else {
        temp <- getQueryLink(ids, repository)
        blanks <- temp == "&nbsp;"
        #epl: introduced an exception 4 localfile repository
        if (repository=='localfile') {out <- paste(" <A HREF=\"", temp, "\">", 'View', "</A>", sep = "")}
          else {out <- paste(" <A HREF=\"", temp, "\">", ids, "</A>", sep = "")}
        out[blanks] <- "&nbsp;"
    }
    return(out)
}

htmlpage <- function (genelist, filename, title, othernames, table.head,
    table.center = TRUE, repository = list("en"), ...) {
    if (!is.list(repository)) 
        stop("The repository argument must be a list!", call. = FALSE)
    chklen <- function(x) {
        if (is.data.frame(x) || is.matrix(x)) 
            dim(x)[1]
        else length(x)
    }
    getRows <- function(x) {
        paste("<P>", x, "</P>", collapse = "", sep = "")
    }
    if (is.data.frame(genelist)) 
        len.vec <- chklen(genelist)
    else if (is.list(genelist)) 
        len.vec <- sapply(genelist, chklen)
    else stop("The 'genelist' should be either a data.frame or a list", 
        call. = FALSE)
    if (!missing(othernames)) {
        if (is.data.frame(othernames)) 
            len.vec <- c(len.vec, chklen(othernames))
        else if (is.list(othernames)) 
            len.vec <- c(len.vec, sapply(othernames, chklen))
        else stop("The 'othernames' should be either a data.frame or a list", 
            call. = FALSE)
    }
    if (any(len.vec != len.vec[1])) 
        stop(paste("Some items in either", genelist, "or", othernames, 
            "have mis-matched lengths.\nPlease check this", "discrepancy and re-run.\n"), 
            .call = FALSE)
    if (is.list(repository)) {
        out <- NULL
        for (i in seq(along = repository)) {
            out <- cbind(out, getCells(genelist[[i]], repository[[i]]))
        }
    }
    else out <- getCells(genelist, repository)
    if (!missing(othernames)) {
        if (is.data.frame(othernames)) 
            out <- data.frame(out, othernames)
        else if (is.list(othernames)) {
            others <- vector("list", length(othernames))
            for (i in seq(along = othernames)) {
                if (is.data.frame(othernames[[i]])) 
                  others[[i]] <- othernames[[i]]
                else if (is.list(othernames[[i]])) {
                  others[[i]] <- sapply(othernames[[i]], getRows)
                }
                else {
                  others[[i]] <- othernames[[i]]
                }
            }
            out <- data.frame(out, as.data.frame(others))
        }
    }
    colnames(out) <- table.head
    out <- xtable(out, caption = if (!missing(title)) 
        title, ...)
    print(out, type = "html", file = filename, caption.placement = "top", 
        include.rownames = FALSE, sanitize.text.function = function(x) x, 
        ...)
}
