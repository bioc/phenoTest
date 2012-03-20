\name{write.html}
\alias{write.html}
\title{
Write a data.frame to an html file.
}
\description{
Creates an html file with links and plots from a table.
}
\usage{
write.html(x, links, tiny.pic, tiny.pic.size = 100, title = "", file, digits = 3)
}
\arguments{
  \item{x}{
    Object of class \code{data.frame}.
  }
  \item{links}{
    Object of class \code{list} with one item per column in x.
    If we want the ith column of x to have links to a site or local file
    we will have to write those links into the ith element of links.
  }
  \item{tiny.pic}{
    Object of class \code{list} with one item per column in x.
    If we want the ith column of x to show plots instead of text we will
    have to write the path to the plots into the ith element of links.
  }
  \item{tiny.pic.size}{
    size of the pictures if any.
  }
  \item{title}{
    Title that will be shown on top of the html file.
  }
  \item{file}{
    path and name of the file that will be created.
  }
  \item{digits}{
    number of digits that will be shown in numeric columns of x.
  }
}
\author{
Evarist Planet
}
\seealso{
  write.csv, write.table, htmlpage
}
\examples{
##
##Code has been commented to avoid the creation of files
##
#(x <- data.frame(gene.symbol=c('AARS','ABCF1','ABLIM1'),value=c(2.054,30.024,5.0221),plot=rep('Open',3)))
#tiny.pic <- links <- vector('list',length=ncol(x))
#links[[1]] <- paste('http://www.genecards.org/index.php?path=/Search/keyword/',x[,1])
#for (i in 1:nrow(x)) {
#  png(paste('~/Desktop/',x[i,1],'.png',sep=''))
#  plot(1:3,log(1:3))
#  dev.off()
#}
#tiny.pic[[3]] <- links[[3]] <- paste(x[,1],'.png',sep='')
#write.html(x,links=links,tiny.pic=tiny.pic,file='~/Desktop/x.html',title='My html test')
}
\keyword{ ~kwd1 }
