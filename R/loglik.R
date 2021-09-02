
.drop <- function(x) {
    if (inherits(x, "Matrix")) x <- x@x
    drop(x)
}

.pmax <- function(x, y) {
    y[y < x] <- x
    return(y)
}

.xmb <- function(X, beta) {
    if (is.matrix(beta)) {
        stopifnot(NROW(beta) == NROW(X))
        return(rowSums(X * beta))
    }
    return(X %*% beta)
}

### if (finite) fun(X %*% beta + offset) * Xmult else value
.dealinf <- function(X, beta, offset, fun, value, Xmult = FALSE) {
    if (is.null(X)) return(value)
    lpr <- .drop(offset + .xmb(X, beta))
    OK <- is.finite(lpr)
    if (all(OK)) {
        ret <- fun(lpr)
    } else {
        ret <- rep(value, length(lpr))
        if (all(!OK)) return(ret)
        ret[OK] <- fun(lpr[OK])
    }
    if (Xmult) {
        ### this is costly, can we avoid this?
        if (!all(OK))
            X[!OK,] <- 0
        ret <- ret * X
    }
    return(ret)
}

.log <- function(x) {
    return(log(.pmax(.Machine$double.eps, x)))
    pos <- (x > .Machine$double.eps)
    if (all(pos)) return(log(x))
#    mx <- min(x)
#    if (mx < -.Machine$double.eps)
#        warning(paste("negative contribution to likelihood (", mx, "); 
#                      constraints violated!"))
    ret <- rep(log(.Machine$double.eps), length(x))
    ret[pos] <- log(x[pos])
    return(ret)
}

..mlt_loglik_interval <- function(d, mml, mmr, offset = 0, beta) {
    RC <- !is.finite(mmr[,1])
    LC <- !is.finite(mml[,1])
    ret <- numeric(nrow(mmr))
    ### right-censored: ret = log(1 - F(...))
    if (any(irc <- RC & !LC)) {
        if (is.matrix(beta)) 
            ret[irc] <- .dealinf(mml[irc,,drop = FALSE], beta[irc,,drop = FALSE], 
                offset[irc], function(x) d$p(x, lower.tail = FALSE, log.p = TRUE), 0)
        else
            ret[irc] <- .dealinf(mml[irc,,drop = FALSE], beta, offset[irc], 
                function(x) d$p(x, lower.tail = FALSE, log.p = TRUE), 0)
    }
    ### left-censored:  ret = log(F(...))
    if (any(ilc <- !RC & LC)) {
        if (is.matrix(beta))
            ret[ilc] <- .dealinf(mmr[ilc,,drop = FALSE], beta[ilc,,drop = FALSE],
                offset[ilc], function(x) d$p(x, lower.tail = TRUE, log.p = TRUE), 1)
        else
            ret[ilc] <- .dealinf(mmr[ilc,,drop = FALSE], beta, offset[ilc], 
                function(x) d$p(x, lower.tail = TRUE, log.p = TRUE), 1)
    }
    ### interval-censored: log(F(...) - F(...))
    iic <- !RC & !LC
    if (is.matrix(beta)) 
        ret[iic] <- .log(.dealinf(mmr[iic,,drop = FALSE], beta[iic,,drop = FALSE], 
                                  offset[iic], d$p, 1) - 
                         .dealinf(mml[iic,,drop = FALSE], beta[iic,,drop = FALSE], 
                                  offset[iic], d$p, 0))
    else
        ret[iic] <- .log(.dealinf(mmr[iic,,drop = FALSE], beta, offset[iic], d$p, 1) - 
                         .dealinf(mml[iic,,drop = FALSE], beta, offset[iic], d$p, 0))
    ret
}

.mlt_loglik_interval <- function(d, mml, mmr, offset = 0, mmtrunc = NULL) {
    if (is.null(mmtrunc))
        return(function(beta) 
                   ..mlt_loglik_interval(d, mml, mmr, offset, beta))
    return(function(beta)
               ..mlt_loglik_interval(d, mml, mmr, offset, beta) - 
               ..mlt_loglik_interval(d, mmtrunc$left, 
                                     mmtrunc$right, offset, beta))
}

..mlt_score_interval <- function(d, mml, mmr, offset = 0, beta, Xmult = TRUE)
    (.dealinf(mmr, beta, offset, d$d, 0, Xmult = Xmult) -
     .dealinf(mml, beta, offset, d$d, 0, Xmult = Xmult)) / 
    .pmax(.Machine$double.eps^(1/3), (.dealinf(mmr, beta, offset, d$p, 1) - 
     .dealinf(mml, beta, offset, d$p, 0)))

.mlt_score_interval <- function(d, mml, mmr, offset = 0, mmtrunc = NULL) {
    if (is.null(mmtrunc))
        return(function(beta, Xmult = TRUE)
                   ..mlt_score_interval(d, mml, mmr, offset, beta, Xmult))
    return(function(beta, Xmult = TRUE)
               ..mlt_score_interval(d, mml, mmr, offset, beta, Xmult) - 
               ..mlt_score_interval(d, mmtrunc$left, 
                                    mmtrunc$right, offset, beta, Xmult))
}

..mlt_hessian_interval <- function(d, mml, mmr, offset = 0, w = 1, beta) {
    Fr <- .dealinf(mmr, beta, offset, d$p, 1)
    Fl <- .dealinf(mml, beta, offset, d$p, 0)
    fr <- .dealinf(mmr, beta, offset, d$d, 0)
    fl <- .dealinf(mml, beta, offset, d$d, 0)
    dfr <- .dealinf(mmr, beta, offset, d$dd, 0)
    dfl <- .dealinf(mml, beta, offset, d$dd, 0)
    if (length(w) != length(Fr)) w <- rep(w, length(Fr))
    Frl <- Fr - Fl
    w1 <- dfr / Frl * w
    w2 <- dfl / Frl * w
    w3 <- fr / Frl * sqrt(w)
    w4 <- fl / Frl * sqrt(w)
    if (is.null(mmr)) 
        mmr <- matrix(0, nrow = nrow(mml), ncol = ncol(mml))
    if (is.null(mml)) 
        mml <- matrix(0, nrow = nrow(mmr), ncol = ncol(mmr))
    mmr[!is.finite(mmr)] <- 0
    mml[!is.finite(mml)] <- 0
    W3 <- mmr * w3
    W4 <- mml * w4
    return(-(crossprod(mmr * w1, mmr) - crossprod(mml * w2, mml) - 
            (crossprod(W3) - crossprod(W3, W4) - 
                crossprod(W4, W3) + crossprod(W4))))
}

.mlt_hessian_interval <- function(d, mml, mmr, offset = 0, mmtrunc = NULL, 
                                  w = 1) {
    if (is.null(mmtrunc))
        return(function(beta)
                   ..mlt_hessian_interval(d, mml, mmr, offset, w, beta))
    return(function(beta) 
               ..mlt_hessian_interval(d, mml, mmr, offset, w, beta) - 
               ..mlt_hessian_interval(d, mmtrunc$left, 
                                      mmtrunc$right, offset, w, beta))
}

.mlt_loglik_exact <- function(d, mm, mmprime, offset = 0, mmtrunc = NULL) {
    if (is.null(mmtrunc)) 
        return(function(beta)
            d$d(offset + .xmb(mm, beta), log = TRUE) + .log(.xmb(mmprime, beta)))
    return(function(beta)
               d$d(offset + .xmb(mm, beta), log = TRUE) + .log(.xmb(mmprime, beta)) - 
               ..mlt_loglik_interval(d, mmtrunc$left, 
                                     mmtrunc$right, offset, beta))
}

.mlt_score_exact <- function(d, mm, mmprime, offset = 0, mmtrunc = NULL) {                          
    function(beta, Xmult = TRUE) {
        mmb <- drop(.xmb(mm, beta)) + offset
        ret <- d$dd2d(mmb) 
        if (Xmult)
            ret <- ret * mm + mmprime / drop(.xmb(mmprime, beta))
        if (!is.null(mmtrunc))
            ret <- ret - ..mlt_score_interval(d, mmtrunc$left, 
                                              mmtrunc$right, offset, beta, Xmult) 
        return(ret)
    }
}

.mlt_hessian_exact <- function(d, mm, mmprime, offset = 0, mmtrunc = NULL, 
                               w = 1) {
    function(beta) {
        mmb <- drop(.xmb(mm, beta)) + offset
        if (length(w) != length(mmb)) w <- rep(w, length(mmb))
        w1 <- -(d$ddd(mmb) / d$d(mmb) - d$dd2d(mmb)^2) * w
        w2 <- w / (drop(.xmb(mmprime, beta))^2)
        ret <- crossprod(mm * w1, mm) + crossprod(mmprime * w2, mmprime) 
        if (!is.null(mmtrunc))
            ret <- ret - ..mlt_hessian_interval(d, mmtrunc$left, 
                                                mmtrunc$right, offset, w, beta)
        return(ret)
    }
}
