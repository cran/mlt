
.frmt <- function(q) {
    if (is.factor(q)) 
        return(as.character(q))
    return(formatC(q, digits = 3, width = 5))
}

### transformation function
tmlt <- function(object, newdata = NULL, q = NULL, ...) {

    vn <- unlist(variable.names(object))
    y <- variable.names(object, "response")
    vnx <- vn[!(vn %in% y)]
    model <- object$model

    stopifnot(!is.null(newdata[[y]]) || !is.null(q))

    if (is.data.frame(newdata)) {

        ### unconditional
        if (length(vnx) == 0 & !is.null(q)) {
            newdata <- data.frame(q)
            names(newdata) <- y
            q <- NULL
        }

        ### in sample predictions
        ### this will _not_ work for censored responses
        if (!is.null(newdata[[y]]) & is.null(q)) {
            stopifnot(is.atomic(newdata[[y]]))
            ret <- c(predict(model, newdata = newdata, 
                             coef = coef(object), ...))
            names(ret) <- rownames(newdata)
            ### P(Y \le y_K) = 1 but trafo can be < Inf
            ### depending on parameterisation
            if (is.factor(f <- newdata[[y]])) {
                i <- f == levels(f)[nlevels(f)]
                if (any(i)) 
                    ret[i] <- Inf
            } else {
                ### Y in (b[1], b[2]) => P(Y \le b[1]) = 0, P(Y \le b[2]) = 1
                ### <FIXME> what happens with deriv in ...? </FIXME>
                b <- bounds(as.vars(object)[[y]])[[1]]
                ret[f < (b[1] + .Machine$double.eps)] <- -Inf
                ret[f > (b[2] - .Machine$double.eps)] <- Inf
            }
            return(ret)
        }

        ### extra quantiles, compute transformation
        ### for each q and each row of newdata
        stopifnot(is.atomic(q))
#        stopifnot(length(unique(q)) == length(q))
        dim <- c(length(q), nrow(newdata))

        ### <FIXME> this triggers a trick in 
        ### basefun:::predict.basis; better checks needed </FIXME>
        names(dim) <- c(y, vnx[1])
        newdata <- as.list(newdata)
        newdata[[y]] <- q
        ret <- predict(object$model, newdata = newdata, 
                       coef = coef(object), dim = dim, ...)
        dn <- vector(mode = "list", length = 2)
        names(dn) <- c(y, "newdata") ### deparse(substitute(newdata))) ?
        dn[[y]] <- .frmt(q)
        dn[["newdata"]] <- rownames(newdata)
        dimnames(ret) <- dn

        ### trafo of last level is always Inf, see above
        if (is.factor(f <- newdata[[y]])) {
            i <- f == levels(f)[nlevels(f)]
            if (any(i))
                ret[i,] <- Inf
        } else {
            ### Y in (b[1], b[2]) => P(Y \le b[1]) = 0, P(Y \le b[2]) = 1
                ### <FIXME> what happens with deriv in ...? </FIXME>
            b <- bounds(as.vars(object)[[y]])[[1]]
            ret[f < (b[1] + .Machine$double.eps)] <- -Inf
            ret[f > (b[2] - .Machine$double.eps)] <- Inf
        }
        return(ret)
    }

    ### need to generate newdata outside tmlt such that
    ### the rows of expand.grid(newdata) match the elements of
    ### the return value
    stopifnot(is.atomic(newdata[[y]]))
    stopifnot(is.null(q))
    stopifnot(y %in% names(newdata))
    ret <- predict(object$model, newdata = newdata, 
                   coef = coef(object), dim = TRUE, ...)
    dn <- lapply(newdata, .frmt)
    dimnames(ret) <- dn

    ### trafo of last level is always Inf, see above
    if (is.factor(f <- newdata[[y]])) {
        i <- f == levels(f)[nlevels(f)]
        if (any(i)) {
            args <- lapply(names(dn), function(d) {
                if (d == y)
                    return(i)
                return(1:(dim(ret)[which(names(dn) == d)]))
            })
            ret <- do.call("[<-", c(list(i = ret), args, 
                                    list(value = Inf)))
        }
    } else {
        ### Y in (b[1], b[2]) => P(Y \le b[1]) = 0, P(Y \le b[2]) = 1
                ### <FIXME> what happens with deriv in ...? </FIXME>
        b <- bounds(as.vars(object)[[y]])[[1]]
        i <- (f < (b[1] + .Machine$double.eps))
        if (any(i)) {
            args <- lapply(names(dn), function(d) {
                if (d == y)
                    return(i)
                return(1:(dim(ret)[which(names(dn) == d)]))
            })
            ret <- do.call("[<-", c(list(i = ret), args, 
                                    list(value = -Inf)))
        }
        i <- (f > (b[2] - .Machine$double.eps))
        if (any(i)) {
            args <- lapply(names(dn), function(d) {
                if (d == y)
                    return(i)
                return(1:(dim(ret)[which(names(dn) == d)]))
            })
            ret <- do.call("[<-", c(list(i = ret), args, 
                                    list(value = Inf)))
        }
    }
    return(ret)
}

### distribution function
pmlt <- function(object, newdata = NULL, q = NULL, lower.tail = TRUE, log = FALSE, ...)
    object$todistr$p(tmlt(object = object, newdata = newdata,
                          q = q, ...), lower.tail = lower.tail, log.p = log)

### survivor function
smlt <- function(object, newdata = NULL, q = NULL, log = FALSE, ...)
    pmlt(object = object, newdata = newdata, q = q, lower.tail = FALSE, 
         log = log, ...)

### cumulative hazard function
Hmlt <- function(object, newdata = NULL, q = NULL, log = FALSE, ...) {
    if (object$todistr$name == "minimum extreme value") {
        ### Cox model
        logH <- tmlt(object = object, newdata = newdata, q = q, ...)
        if (log) return(logH)
        return(exp(logH))
    }
    if (object$todistr$name == "exponential") {
        ### Aalen model
        H <- tmlt(object = object, newdata = newdata, q = q, ...)
    } else {
        ### generic
        H <- -smlt(object = object, newdata = newdata, q = q, log = TRUE, 
                   ...)
    }
    if (log) return(log(H))
    return(H)
}

### odds function
Omlt <- function(object, newdata = NULL, q = NULL, log = FALSE, ...) {
    if (object$todistr$name == "logistic") {
        logO <- tmlt(object = object, newdata = newdata, q = q, ...)
        if (log) return(logO)
        return(exp(logO))
    }
    F <- pmlt(object = object, newdata = newdata, q = q, log = log, ...)
    S <- smlt(object = object, newdata = newdata, q = q, log = log, ...)
    if (log) return(F - S)
    return(F / S)
}

### numerically invert z = f(q)
.invf <- function(object, f, q, z) {

    if (!is.matrix(f)) f <- matrix(f, nrow = 1)
    N <- nrow(f)
    K <- length(q)
    stopifnot(ncol(f) == K)

    if (!is.matrix(z)) z <- matrix(z, nrow = N)
    stopifnot(nrow(z) == N)
    nsim <- ncol(z)

    y <- variable.names(object, "response")
    bounds <- bounds(object)[[y]]

    discrete <- !inherits(as.vars(object)[[y]],
                          "continuous_var")
    if (discrete) {
        ### use "old" code
        f <- f[rep(1:N, nsim), , drop = FALSE]
        f <- cbind(-Inf, f, Inf)
        i <- rowSums(f < as.vector(z))
        return(q[i])
    }

    ### use spline/approx to evaluate quantile function
    nleft <- nexact <- nright <- 
        matrix(NA, nrow = N, ncol = nsim)

    for (i in 1:N) {
             
        pr <- f[i,]
        pr0 <- which(pr < min(pr) + sqrt(.Machine$double.eps))
        pr1 <- which(pr > max(pr) - sqrt(.Machine$double.eps))
        rmi <- which(!is.finite(pr))
        if (length(pr0) > 0L)
            rmi <- c(rmi, pr0[-length(pr0)])
        if (length(pr1) > 0L)
            rmi <- c(rmi, pr1[-1L])
        if (length(rmi) > 0L) {
            qq <- q[-rmi]
            pr <- pr[-rmi]
        } else {
            qq <- q
        }
        s <- spline(x = qq, y = pr, method = "hyman")
        ynew <- approx(x = s$y, y = s$x, xout = z[i, ], 
                       yleft = -Inf, yright = Inf)$y
        if (bounds[1L] == min(q))
            ynew[ynew == -Inf] <- min(q)
        if (bounds[2L] == max(q))
            ynew[ynew == Inf] <- max(q)
        nleft[i, ynew == -Inf] <- bounds[1L]
        nright[i, ynew == -Inf] <- min(q)
        nleft[i, ynew == Inf] <- max(q)
        nright[i, ynew == Inf] <- bounds[2L]
        ynew[!is.finite(ynew)] <- NA
        nexact[i, ] <- ynew
    }

    if (all(is.na(nleft)) && all(is.na(nright)))
        return(as.vector(nexact))
    return(R(as.vector(nexact), 
             cleft = as.vector(nleft), 
             cright = as.vector(nright)))
}


### quantile function
qmlt <- function(object, newdata = NULL, q = NULL, prob = .5, n = 50, 
                 interpolate = FALSE, ...) {

   if (interpolate)
        warning("Argument interpolate ignored in mlt >= 1.2-1")

    stopifnot(all(prob > sqrt(.Machine$double.eps)) && 
              all(prob < 1 - sqrt(.Machine$double.eps)))

    y <- variable.names(object, "response")
    if (is.null(q))
        q <- mkgrid(object, n = n)[[y]]
    if (!is.null(newdata) & !is.data.frame(newdata)) {
        newdata[[y]] <- NULL
        nm <- names(newdata)
        newdata[[y]] <- q
        newdata <- newdata[c(y, nm)]
        tr <- pmlt(object, newdata, ...)
    } else {
        tr <- pmlt(object, newdata = newdata, q = q, ...)
    } 

    ### convert potential array-valued distribution function
    ### to matrix where rows correspond to observations newdata 
    ### and columns to quantiles q
    trm <- matrix(tr, ncol = length(q), byrow = TRUE)
    ### Note: we invert the cdf direct; 
    ### inverting transformation functions using approx didn't work well
    qu <- matrix(prob, 
                 nrow = nrow(trm), ncol = length(prob), byrow = TRUE)

    ret <- .invf(object, f = trm, q = q, z = qu)
    i <- as.vector(matrix(1:length(ret), ncol = nrow(trm), byrow = TRUE))
    ret <- if (is.matrix(ret)) ret[i,] else ret[i]

    ### arrays of factors are not allowed
    if (is.factor(q)) return(ret)

    ### return "response" object
    if (inherits(ret, "response")) 
        return(ret)

    dim <- dim(tr)
    dim[1] <- length(prob)
    dn <- c(list(prob = .frmt(prob)), dimnames(tr)[-1L])
    return(array(ret, dim = dim, dimnames = dn))
}

### density
dmlt <- function(object, newdata = NULL, q = NULL, log = FALSE, ...) {

    y <- variable.names(object, "response")
    response <- mkgrid(object, n = 10)[[y]]

    ### Lebesgue density only for double
    if (.type_of_response(response) %in% c("double", "survival")) {
        trafo <- tmlt(object, newdata = newdata, q = q, ...)
        deriv <- 1
        names(deriv) <- y
        trafoprime <- tmlt(object, newdata = newdata, q = q, 
                           deriv = deriv, ...)
        ### <FIXME> trafoprime is +/-Inf at boundaries, so use 0 density
        trafoprime[!is.finite(trafoprime)] <- .Machine$double.eps
        trafoprime <- pmax(.Machine$double.eps, trafoprime)
        if (log)
            return(object$todistr$d(trafo, log = TRUE) + log(trafoprime))
        return(object$todistr$d(trafo) * trafoprime)
    }

    stopifnot(!is.null(newdata[[y]]) || !is.null(q))

    ### for factors and integers compute density as F(y) - F(y - 1)
    lev <- levels(response)
    if (is.data.frame(newdata)) {

        ### in sample density
        if (!is.null(newdata[[y]]) & is.null(q)) {
            stopifnot(is.atomic(newdata[[y]]))
            q <- newdata[[y]]
            first <- q == lev[1]
            qwoK <- factor(lev[pmax(unclass(q) - 1, 1)], 
                           levels = lev, labels = lev, ordered = is.ordered(q))
            p <- pmlt(object, newdata = newdata, ...)
            newdata[[y]] <- qwoK
            pwoK <- pmlt(object, newdata = newdata, ...)
            pwoK[first] <- 0
            ret <- p - pwoK
        } else {
            ### extra quantiles, compute density
            ### for each q and each row of newdata 
            stopifnot(is.atomic(q))
            first <- q == lev[1]
            qfirst <- q[first]
            qwoK <- q[q != lev[length(lev)]]
            qwo1 <- q[q != lev[1]]

            pfirst <- pmlt(object, newdata = newdata, q = qfirst, ...)
            pwo1 <- pmlt(object, newdata = newdata, q = qwo1, ...)
            pwoK <- pmlt(object, newdata = newdata, q = qwoK, ...)
            ret <- matrix(0, nrow = length(first), ncol = NCOL(pfirst))
            ret[!first,] <- pwo1 - pwoK
            ret[first,] <- pfirst
            rownames(ret) <- as.character(q)
       }
    } else {

        ### need to generate newdata outside tmlt such that
        ### the rows of expand.grid(newdata) match the elements of
        ### the return value
        stopifnot(is.atomic(newdata[[y]]))
        stopifnot(is.null(q))
        stopifnot(y %in% names(newdata))
        dim <- sapply(newdata, NROW)
        q <- newdata[[y]]

        first <- q == lev[1]
        qfirst <- q[first]
        qwoK <- q[q != lev[length(lev)]]
        qwo1 <- q[q != lev[1]]

        newdata[[y]] <- qfirst
        pfirst <- pmlt(object, newdata = newdata, ...)
        newdata[[y]] <- qwo1
        pwo1 <- pmlt(object, newdata = newdata, ...)
        newdata[[y]] <- qwoK
        pwoK <- pmlt(object, newdata = newdata, ...)

        dn <- dim(pfirst)
        names(dn) <- names(dimnames(pfirst))

        frst <- lapply(names(dn), function(d) {
            if (d != y) return(1:dn[d])
            return(first)
        })
        ntfrst <- lapply(names(dn), function(d) {
            if (d != y) return(1:dn[d])
            return(!first)
        })

        dn <- lapply(names(dn), function(d) {
            if (d != y) return(dimnames(pfirst)[[d]])
            return(as.character(q))
        })

        ret <- array(0, dim = dim, dimnames = dn)
        ret <- do.call("[<-", 
            c(list(i = ret), frst, list(value = pfirst)))
        ret <- do.call("[<-", 
            c(list(i = ret), ntfrst, list(value = pwo1 - pwoK)))
    }

    if (log) return(log(ret))
    return(ret)
}

### hazard function
hmlt <- function(object, newdata = object$data, q = NULL, log = FALSE, ...) {

    y <- variable.names(object, "response")
    response <- mkgrid(object, n = 10)[[y]]

    if (.type_of_response(response) %in% c("double", "survival")) {
        d <- dmlt(object, newdata = newdata, q = q, log = log, ...)
        s <- smlt(object, newdata = newdata, q = q, log = log, ...)
        if (log) 
            return(d - s)
        return(d / s)
    }
    ### discrete hazard: Prob(Y = y | Y >= y)
    d <- dmlt(object, newdata = newdata, q = q, log = log, ...)
    p <- pmlt(object, newdata = newdata, q = q, ...)
    if (log)
        return(d - log1p(-(p - exp(d))))
    return(d / (1 - (p - d)))
}
