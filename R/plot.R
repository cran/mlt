
plot.ctm <- function(x, newdata, type = c(
    "distribution", "logdistribution",
    "survivor", "logsurvivor", 
    "density", "logdensity", 
    "hazard", "loghazard", 
    "cumhazard", "logcumhazard", 
    "odds", "logodds", 
    "quantile", "trafo"),
    q = NULL, prob = 1:(K - 1) / K, K = 50, col = rgb(.1, .1, .1, .1), lty = 1, add = FALSE, ...) {

    args <- list(...)
    y <- variable.names(x, "response")

    if (is.null(q))
        q <- mkgrid(x, n = K)[[y]]
    type <- match.arg(type)

    if (type == "quantile") {
        pr <- predict(x, newdata = newdata, type = "distribution", q = q)
        pmin <- max(apply(pr, 2, function(p) min(p, na.rm = TRUE)))
        pmax <- min(apply(pr, 2, function(p) max(p, na.rm = TRUE)))
        prob <- prob[prob > pmin & prob < pmax]
    }
    pr <- predict(x, newdata = newdata, type = type, q = q, prob = prob)
    pr[!is.finite(pr)] <- NA
    rpr <- range(pr, na.rm = TRUE)
    if (is.null(dim(pr))) pr <- matrix(pr, ncol = 1)
    ylim <- switch(type, "distribution" = c(0, 1),
                         "survivor" = c(0, 1),
                         "density" = c(0, rpr[2]),
                         "hazard" = c(0, rpr[2]),
                         "cumhazard" = c(0, rpr[2]),
                         rpr)
    if (!is.null(args$ylim)) ylim <- args$ylim
    if (type == "quantile")  {
        q <- prob
        if (is.null(args$xlab)) args$xlab <- type
        if (is.null(args$ylab)) args$ylab <- y
    }
    if (length(col) == 1) col <- rep_len(col, ncol(pr))
    if (length(lty) == 1) lty <- rep_len(lty, ncol(pr))
    lwd <- 1
    if (!is.null(args$lwd)) lwd <- args$lwd
    if (length(lwd) == 1) lwd <- rep_len(lwd, ncol(pr))
    
    if (!add) {
        args$x <- unclass(q)
        args$y <- rep_len(rpr[1], length(q))
        args$ylim <- ylim
        args$xlab <- ifelse(is.null(args$xlab), y, args$xlab)
        args$ylab <- ifelse(is.null(args$ylab), type, args$ylab)
        args$type <- "n"
        args$axes <- FALSE
        do.call("plot", args)
        if (is.factor(q)) {
            axis(1, at = unclass(q), labels = levels(q))
        } else {
            axis(1)
        }
        axis(2)
        box()
    }
    y <- as.vars(x)[[y]]
    if (inherits(y, "continuous_var")) {
        for (i in 1:ncol(pr)) 
            lines(q, pr[,i], col = col[i], lty = lty[i], lwd = lwd[i])
    } else {
        for (i in 1:ncol(pr)) 
            lines(stepfun(q, c(ylim[1], pr[,i])), 
                  col = col[i], lty = lty[i], lwd = lwd[i])
    }
    invisible(pr)
}

plot.mlt <- function(x, ...) {
    ctmobj <- x$model
    coef(ctmobj) <- coef(x)
    plot(ctmobj, ...)
}
