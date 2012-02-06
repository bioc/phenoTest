smoothCoxph <- function(time, event, x, xlim, ylim, xlab, ylab, logrisk=TRUE, ...) {
  z <- data.frame(time = time, event = event, x = x)
  z <- z[!is.na(z$event) & !is.na(z$time) & !is.na(z$x),]
  z <- z[order(z$x), ]
  coxph1 <- coxph(Surv(time, event) ~ pspline(x, df = 4), data = z)
  if (logrisk) {
    pred <- predict(coxph1, type = "lp", se.fit = TRUE)
    y <- cbind(pred$fit,pred$fit-1.96*pred$se.fit,pred$fit+1.96*pred$se.fit)
  } else {
    pred <- predict(coxph1, type='risk', se.fit=TRUE)
    y <- cbind(log(pred$fit),log(pred$fit-1.96*pred$se.fit),log(pred$fit+1.96*pred$se.fit))
  }
  if (missing(xlim)) xlim <- range(x)
  if (missing(ylim)) ylim <- c(min(y[,1],median(y[,2]),na.rm=T),max(y[,1],median(y[,3]),na.rm=T)) * 1.5
  if (missing(xlab)) xlab <- 'Gene expression'
  if (missing(ylab)) ylab <- 'log Hazard Ratio'  
  if (any(is.infinite(ylim)))
    plot(z$x[1:length(pred$fit)], y[, 1], type = "l", xlim = xlim,xlab=xlab,ylab=ylab)
  else
  plot(z$x[1:length(pred$fit)], y[, 1], type = "l", ylim = ylim, xlim = xlim,xlab=xlab,ylab=ylab,...)
  lines(z$x[1:length(pred$fit)], y[, 2], lty = 2,...)
  lines(z$x[1:length(pred$fit)], y[, 3], lty = 2,...)
  abline(h=0,lty=2,col='grey')
}
