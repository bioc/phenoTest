pca.2d <- function(x, pc, pc.lab, group, group2, pair, names, ellipse=FALSE, main='', components= c(1, 2), legend=TRUE) {
  #control errors
  stopifnot(is(x, 'ExpressionSet'))
  if (!missing(group)) {
    stopifnot(group %in% colnames(pData(x)))
  }
  if (!missing(group2)) {
    stopifnot(group2 %in% colnames(pData(x)))
  }
  if (!missing(pair)) {
    stopifnot(pair %in% colnames(pData(x)))
  }
  if (!missing(names)) {
    stopifnot(names %in% colnames(pData(x)))
  }
  stopifnot(is(components, 'numeric') | is(components, 'integer'))
  stopifnot(all(!is.na(components)))
#
  colnames(pc) <- c('pc1', 'pc2')
  #get plot limits
  lim <- range(pc)
  d <- dist(lim) *.05
  lim <- c(lim[1]-d, lim[2]+d)  
#
  #parameters for plot
  add <- c()
  if (!missing(group)) {
    colour <- x[[group]]
    if (legend) legend.position <- 'top' else legend.position <- 'none'
  } else {
    colour <- 'black'
    legend.position <- 'none'
  }
  if (!missing(group2)) {
    shape <- x[[group2]]
  } else {
    shape <- ''
    add <- c(add, 'guides(shape=FALSE)')
  }
  if (!missing(pair)) {
    line <- x[[pair]]
    add <- c(add, "geom_line(aes(group=line), colour='grey')")
  } else {
    line <- ''
  }
  if (!missing(names)) {
    label <- x[[names]]
    add <- c(add, 'geom_text(aes(label=label), size=3)')
  } else {
    label <- ''
  }
  dat <- data.frame(pc, colour, shape, line, label)
#
  #ellipse with 95% CI of the mean of each group
  if (ellipse) {
    df.ell <- data.frame()
    for (i in 1:length(levels(colour))) {
      sel <- colour %in% levels(colour)[i]
      df.ell <- rbind(df.ell, cbind(as.data.frame(ellipse(x=cov(data.frame(pc[sel, ])/sqrt(nrow(pc))),
                                                             centre=colMeans(pc[sel, ]))),
                      colour=levels(colour)[i], shape=''))
    }
    add <- c(add, 'geom_path(data=df.ell, aes(x=pc1, y=pc2, group=colour))')
  }
#
  #plot
  tmp <- qplot(x=pc1, y=pc2, data=dat, xlab=pc.lab[1], ylab=pc.lab[2], colour=colour, shape=shape, main=main) + 
    geom_point() + coord_cartesian(xlim=lim, ylim=lim) + theme(legend.position=legend.position) +
      ggtitle(main)
  if (length(add)>0) for (i in 1:length(add)) tmp <- tmp + eval(parse(text=add[i]))
  tmp
}

## pca.multipleDim <- function(x, pc, pc.lab, group, group2, pair, names, ellipse, main, components) {
##   idx <- combn(1:length(components), 2)
##   grid.newpage()
##   pushViewport(viewport(layout = grid.layout(length(components), length(components)-1,
##                           heights = unit(c(0.5, rep(10/(length(components)-1), length(components)-1)), "null"))))
##   for (i in 1:ncol(idx)) {
##     pc.sel <- pc[, idx[, i]]
##     tmp <- pca.2d(x, pc.sel, pc.lab[idx[, i]], group, group2, pair, names, ellipse, main='', components[idx[, i]],
##                   legend=F)
##     print(tmp, vp=viewport(layout.pos.row = idx[1, i]+1, layout.pos.col=idx[2, i]-1))
##   }
##   #add title (we plot a blank plot)
##   vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
##   dummy <- ggplot(data.frame(main=main), aes(x=1,y=1,label=main)) + geom_text(size=5) +
##     theme(panel.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
##           axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
##           axis.title.x=element_blank(), axis.title.y=element_blank())
##   print(dummy,vp=vplayout(1,1:(length(components)-1)))
##   #add legend
##   add <- c()
##   if (!missing(group2)) {
##     shape <- x[[group2]]
##   } else {
##     shape <- ''
##     add <- c(add, 'guides(shape=FALSE)')
##   }
##   dat <- data.frame(colour=x[[group]], shape=shape, dummy1=rep(.5, ncol(x)), dummy2=rep(.5, ncol(x)))
##   tmp <- ggplot(aes(x=dummy1, y=dummy2, colour=colour, shape=shape), data=dat) + geom_point()
##   tmp <- tmp + theme(panel.background=element_blank(), panel.grid.minor=element_blank(),
##                      panel.grid.major=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
##                      axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
##                      legend.position=c(.5, .5))
##   if (length(add)>0) for (i in 1:length(add)) tmp <- tmp + eval(parse(text=add[i]))
##   print(tmp, vp=viewport(layout.pos.row = length(components), layout.pos.col=1))
## }

pca <- function(x, group, group2, pair, names, ellipse=FALSE, main='', components= c(1, 2)) {
  #libs
  require(ggplot2)
#  require(gridExtra)
  #error control
  stopifnot(length(components)<=4)
  #calculate principal components
  pcdat <- prcomp(t(exprs(x)))
  pc <- data.frame(pcdat$x[, components])
  #pc.lab <- round(pcdat$sdev/sum(pcdat$sdev)*100,1)[components]
  pc.lab <- round(pcdat$sdev^2/sum(pcdat$sdev^2)*100,1)[components]
  pc.lab <- paste('PC', components, ' (', pc.lab, '%)', sep='')
  #choose 2d or 3d
  if (length(components)==2) {
    pca.2d(x, pc, pc.lab, group, group2, pair, names, ellipse, main, components)
  } else {
#    pca.multipleDim(x, pc, pc.lab, group, group2, pair, names, ellipse, main, components)
    stop('multiple dimension pca will not work until problems in gridExtra package get solved')
  }
}
