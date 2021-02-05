    
### simulate from model object with data newdata
simulate.ctm <- function(object, nsim = 1, seed = NULL, 
                         newdata, K = 50, q = NULL, interpolate = FALSE,
                         bysim = TRUE, ...) {

    if (interpolate)
        warning("Argument interpolate ignored in mlt >= 1.2-1")

    ### from stats:::simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }

    if (!is.data.frame(newdata))
        stop("not yet implemented")


    y <- variable.names(object, "response")
    if (is.null(q))
        q <- mkgrid(object, n = K)[[y]]
    newdata <- newdata[colnames(newdata) %in% variable.names(object)]
    newdata[y] <- NULL
    if (NCOL(newdata) == 0L) newdata <- data.frame(1)
    U <- matrix(runif(nsim * NROW(newdata)), ncol = nsim)
    pr <- predict(object, newdata = newdata, q = q, type = "distribution")
    if (!is.matrix(pr)) pr <- matrix(pr, ncol = 1)

    ret <- .invf(object, f = t(pr), q = q, z = U)

    if (nsim > 1) {
        if (bysim) {
            tmp <- vector(mode = "list", length = nsim)
            for (i in 1:nsim) {
                idx <- 1:NROW(newdata) + (i - 1) * NROW(newdata)
                if (is.data.frame(ret)) 
                    tmp[[i]] <- ret[idx,]
                else 
                    tmp[[i]] <- ret[idx]
            }
        } else {
            tmp <- vector(mode = "list", length = NROW(newdata))
            for (i in 1:NROW(newdata)) {
                idx <- NROW(newdata) * 0:(nsim - 1) + i
                if (is.data.frame(ret)) 
                    tmp[[i]] <- ret[idx,]
                else 
                    tmp[[i]] <- ret[idx]
            }
        }
        ret <- tmp
    }
    return(ret)
}

simulate.mlt <- function(object, nsim = 1, seed = NULL, newdata = object$data, bysim = TRUE, ...) {
    ctmobj <- object$model
    coef(ctmobj) <- coef(object)
    simulate(ctmobj, nsim = nsim, seed = seed, newdata = newdata, bysim = bysim, ...)
}

paraboot <- function(object, ...)
    UseMethod("paraboot")

paraboot.mlt <- function(object, parm = coef(object), B = 100, ...) {

    coefs <- matrix(0, nrow = B, ncol = length(parm))
    colnames(coefs) <- names(parm)
    logLR <- numeric(B)

    ndf <- object$data
    y <- variable.names(object, "response")

    for (b in 1:B) {
        ndf[[y]] <- simulate(object, ...)
        m <- mlt(model = object$model, data = ndf, theta = parm, weights = object$weights, 
                 scale = object$scale, offset = object$offset, check = FALSE, checkGrad = FALSE)
        coefs[b,] <- coef(m)
        logLR[b] <- -2 * (logLik(m, parm = parm) - logLik(m))
    }

    object$coefs <- coefs
    object$logLR <- logLR
    class(object) <- c("mlt_paraboot", class(object))
    object
}
