smoothCoxph <- function(time, event, x, xlim, ylim, ...) {
  z <- data.frame(time = time, event = event, x = x)
  z <- z[order(z$x), ]
  coxph1 <- coxph(Surv(time, event) ~ pspline(x, df = 4), data = z)
  pred <- predict(coxph1, type = "risk", se.fit = TRUE)
  tmp <- pred$fit - 1.96 * pred$se.fit
  sel <- tmp>0; tmp[sel] <- log(tmp[sel]); tmp[!sel] <- NA
  y <- cbind(log(pred$fit), tmp, log(pred$fit + 1.96 * pred$se.fit))
  if (missing(xlim)) xlim <- range(x)
  if (missing(ylim)) ylim <- 2*log(range(pred$fit))
  if (any(is.infinite(ylim)))
    plot(z$x[1:length(pred$fit)], y[, 1], type = "l", xlim = xlim,xlab='Gene expression',ylab='Hazard')
  else
    plot(z$x[1:length(pred$fit)], y[, 1], type = "l", ylim = 2*ylim, xlim = xlim,xlab='Gene expression',ylab='Hazard')
  lines(z$x[1:length(pred$fit)], y[, 2], lty = 2)
  lines(z$x[1:length(pred$fit)], y[, 3], lty = 2)
  abline(h=0,lty=2,col='grey')
}
