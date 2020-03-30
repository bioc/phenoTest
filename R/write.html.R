write.html <- function(x, links, tiny.pic, tiny.pic.size=100, title='', file, digits=3) {
  stopifnot(is(x, 'data.frame'))
  if (missing(links)) links <- vector('list',ncol(x))
  if (missing(tiny.pic)) tiny.pic <- vector('list',ncol(x))  
  stopifnot(is(links, 'list'))
  stopifnot(is(tiny.pic, 'list'))
  stopifnot(length(links)==ncol(x))
  stopifnot(length(tiny.pic)==ncol(x))  
  stopifnot(!missing(file))

  #set factors as characters
  column.class <- unlist(lapply(x,class))
  for (j in 1:ncol(x)) {
    if (column.class[j]=='factor') x[,j] <- as.character(x[,j])
    if (column.class[j]=='numeric') x[,j] <- round(x[,j],digits=digits)
  }

  #header
  cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">\n',sep='',file=file,append=F)
  cat('<html>\n',file=file,append=T)
  cat('<body>\n',file=file,append=T)  
  cat(paste('<CAPTION ALIGN="top"><center><B>',title,'</B></center></CAPTION><BR>\n'),sep='',file=file,append=T)
  cat('<TABLE border=1>\n',file=file,append=T)

  #first row with column names
  cat('<TR>\n',file=file,append=T)
  for (j in 1:ncol(x)) {
    cat('<TH>',file=file,append=T)
    cat(colnames(x)[j],file=file,append=T)
    cat('</TH>\n',file=file,append=T)
  }
  cat('</TR>\n',file=file,append=T)     

  #rest of the data
  for (i in 1:nrow(x)) {
    cat('<TR>\n',file=file,append=T)
    for (j in 1:ncol(x)) {
      cat('<TD>',file=file,append=T)
      if (is.null(links[[j]]) & is.null(tiny.pic[[j]])) {
        cat(x[i,j],file=file,append=T)
      } else if (is.null(links[[j]]) & !is.null(tiny.pic[[j]])) { #
        cat(paste('<A HREF="',links[[j]][[i]],'"><img src="',tiny.pic[[j]][[i]],'" height="',tiny.pic.size,'" width="',tiny.pic.size,'" /></A>',sep=''),file=file,append=T)
      } else if (!is.null(links[[j]]) & is.null(tiny.pic[[j]])) {
        cat(paste('<A HREF="',links[[j]][[i]],'">',x[i,j],'</A>',sep=''),file=file,append=T)
      } else if (!is.null(links[[j]]) & !is.null(tiny.pic[[j]])) { #
        cat(paste('<A HREF="',links[[j]][[i]],'"><img src="',tiny.pic[[j]][[i]],'" height="',tiny.pic.size,'" width="',tiny.pic.size,'" /></A>',sep=''),file=file,append=T)
      }
      cat('</TD>\n',file=file,append=T)
    }
    cat('</TR>\n',file=file,append=T)
  }

  #close table
  cat('</TABLE>\n',file=file,append=T)
  cat('</body>\n',file=file,append=T)
  cat('</html>\n',file=file,append=T)  

  #add javascript to order & drag columns
  phenoTest:::sortDragHtmlTable(filename=file)
}
