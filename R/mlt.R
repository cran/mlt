
### <FIXME> rename fixed to coef and allow for specification of coefs, 
###         ie fitted models? </FIXME>
.mlt_setup <- function(model, data, y, offset = NULL, 
                       fixed = coef(model)) {

    response <- variable.names(model, "response")
    stopifnot(length(response) == 1)
    todistr <- model$todistr

    eY <- .mm_exact(model, data = data, response = response, object = y)
    iY <- .mm_interval(model, data = data, response = response, object = y)

    if (is.null(eY)) {
        Y <- iY$Yleft
    } else {
        Y <- eY$Y
    }

    ui <- attr(Y, "constraint")$ui
    ci <- attr(Y, "constraint")$ci

    if (!is.null(fixed)) {
        fixed <- fixed[!is.na(fixed)]
        if (length(fixed) == 0) fixed <- NULL
    }

    if (!is.null(fixed)) {
        stopifnot(all(names(fixed) %in% colnames(Y)))
        fix <- colnames(Y) %in% names(fixed)
        fixed <- fixed[colnames(Y)[fix]]

        ### adjust contrasts a fixed parameter contributes to
        ci <- ci - c(as(ui[, fix, drop = FALSE], "matrix") %*% fixed)

        ### remove columns corresponding to fixed parameters
        ui <- ui[,!fix,drop = FALSE]
        .parm <- function(beta) {
            nm <- names(beta)
            if (is.matrix(beta)) nm <- colnames(beta)
            ### check if beta already contains fix 
            ### parameters
            if (all(names(fixed) %in% nm))
                return(beta)
            ret <- numeric(ncol(Y))
            ret[fix] <- fixed
            if (is.matrix(beta)) {
                ret <- matrix(ret, nrow = 1)
                ret <- ret[rep(1L, NROW(beta)),,drop = FALSE]
                ret[,!fix] <- beta
            } else {
                ret[!fix] <- beta
            }
            ret
        }
    } else {
        .parm <- function(beta) beta
        fix <- rep(FALSE, ncol(Y))
    } 

    .ofuns <- function(weights, subset = NULL, offset = NULL, 
                       perm = NULL, permutation = NULL, 
                       distr = todistr) ### <- change todistr via update(, distr =)
    {
        if (is.null(offset)) offset <- rep(0, nrow(data))
        es <- .exact_subset(.exact(y), subset)
        exY <- NULL
        iYleft <- NULL

        if (!is.null(perm)) {
            if (length(unique(weights)) > 1 || !is.null(subset))
                stop("permutations not implemented for weights or subsets")
        }

        nm <- colnames(Y)
        if (!is.null(perm)) {
            stopifnot(all(perm %in% nm))
            if (is.null(permutation)) 
                permutation <- sample(1:NROW(Y))
            X <- matrix(0, nrow = NROW(y), ncol = ncol(Y))
            if (!is.null(eY)) {
                X[eY$which,] <- eY$Y
                colnames(X) <- colnames(eY$Y)
            }
            if (!is.null(iY)) {
                tmpl <- iY$Yleft
                cc <- complete.cases(tmpl)
                tmpl[!cc,] <- iY$Yright[!cc,]
                X[iY$which,] <- tmpl
                colnames(X) <- colnames(iY$Yleft)
            }
            Xperm <- X[permutation, perm, drop = FALSE]
        }

        if (!is.null(es$full_ex)) {
            exY <- eY$Y
            exYprime <- eY$Yprime
            exoffset <- offset[.exact(y)]
            exweights <- weights[.exact(y)]
            extrunc <- eY$trunc
            if (!is.null(perm)) {
                stopifnot(is.null(extrunc))
                exY[, perm] <- Xperm[eY$which,]
            }
            if (!is.null(es$redu_ex)) {
                exY <- exY[es$redu_ex,,drop = FALSE]
                exYprime <- exYprime[es$redu_ex,,drop = FALSE]
                exoffset <- exoffset[es$redu_ex]
                exweights <- exweights[es$redu_ex]
                if (!is.null(extrunc)) {
                    extrunc$left <- extrunc$left[es$redu_ex,,drop = FALSE]
                    extrunc$right <- extrunc$right[es$redu_ex,,drop = FALSE]
                }
            }
        }
        if (!is.null(es$full_nex)) {
            iYleft <- iY$Yleft
            iYright <- iY$Yright
            ioffset <- offset[!.exact(y)]
            iweights <- weights[!.exact(y)]
            itrunc <- iY$trunc
            if (!is.null(perm)) {
                stopifnot(is.null(itrunc))
                iYleft[, perm] <- Xperm[iY$which,]
                iYright[, perm] <- Xperm[iY$which,]
            }
            if (!is.null(es$redu_nex)) {
                iYleft <- iYleft[es$redu_nex,,drop = FALSE]
                iYright <- iYright[es$redu_nex,,drop = FALSE]
                ioffset <- ioffset[es$redu_nex]
                iweights <- iweights[es$redu_nex]
                if (!is.null(itrunc)) {
                    itrunc$left <- itrunc$left[es$redu_nex,,drop = FALSE]
                    itrunc$right <- itrunc$right[es$redu_nex,,drop = FALSE]
                }
            }
        }
        ret_ll <- numeric(nrow(data))
        .matrix <- matrix
        if (inherits(exY, "Matrix") || inherits(iYleft, "Matrix")) 
            .matrix <- function(...) Matrix(..., sparse = TRUE)
        ret_scM <- .matrix(0, nrow = nrow(data), ncol = length(fix))
        ret_sc <- matrix(0, nrow = nrow(data), ncol = 1L)
        rownames(ret_scM) <- rownames(ret_sc) <- rownames(data)
        EX_ONLY <- isTRUE(all.equal(es$full_ex, 1:nrow(data)))
        IN_ONLY <- isTRUE(all.equal(es$full_nex, 1:nrow(data)))
        return(list(
            ll = function(beta) {
                ret <- ret_ll 
                if (is.matrix(beta)) {
                    beta_ex <- beta[es$full_ex,,drop = FALSE]
                    beta_nex <- beta[es$full_nex,,drop = FALSE]
                } else {
                    beta_ex <- beta_nex <- beta
                }
                if (!is.null(es$full_ex))
                    ret[es$full_ex] <- .mlt_loglik_exact(distr, 
                        exY, exYprime, exoffset, extrunc)(.parm(beta_ex))
                if (!is.null(es$full_nex))
                    ret[es$full_nex] <- .mlt_loglik_interval(distr, 
                        iYleft, iYright, ioffset, itrunc)(.parm(beta_nex))
                return(ret)
            },
            sc = function(beta, Xmult = TRUE) {
                ### Xmult = FALSE means
                ### score wrt to an intercept term. This avoids
                ### multiplication with the whole design matrix.
                ### Don't use on your own. 
                if (Xmult) {
                    ret <- ret_scM 
                } else {
                    ret <- ret_sc
                }
                if (is.matrix(beta)) {
                    beta_ex <- beta[es$full_ex,,drop = FALSE]
                    beta_nex <- beta[es$full_nex,,drop = FALSE]
                    nm <- colnames(beta)
                } else {
                    beta_ex <- beta_nex <- beta
                    nm <- names(beta)
                }
                if (!is.null(es$full_ex)) {
                    scr <- .mlt_score_exact(distr, 
                        exY, exYprime, exoffset, extrunc)(.parm(beta_ex), Xmult)
                    if (EX_ONLY) {
                        ret <- scr
                    } else {
                        ret[es$full_ex,] <- scr
                    }
                }
                if (!is.null(es$full_nex)) {
                    scr <- .mlt_score_interval(distr, 
                        iYleft, iYright, ioffset, itrunc)(.parm(beta_nex), Xmult)
                    if (IN_ONLY) {
                        ret <- scr
                    } else {
                        ret[es$full_nex,] <- scr
                    }
                }
                if (!Xmult) return(ret)
                colnames(ret) <- colnames(Y)
                ### in case beta contains fix parameters,
                ### return all scores
                if (!is.null(fixed)) {
                    if (all(names(fixed) %in% nm))
                        return(ret)
                }
                return(ret[, !fix, drop = FALSE])
            },
            he = function(beta) {
                ret <- 0
                if (is.matrix(beta)) {
                    beta_ex <- beta[es$full_ex,,drop = FALSE]
                    beta_nex <- beta[es$full_nex,,drop = FALSE]
                    nm <- colnames(beta)
                } else {
                    beta_ex <- beta_nex <- beta
                    nm <- names(beta)
                }
                if (!is.null(es$full_ex))
                    ret <- ret + .mlt_hessian_exact(distr, 
                        exY, exYprime, exoffset, extrunc, 
                        exweights)(.parm(beta_ex))
                if (!is.null(es$full_nex))
                    ret <- ret + .mlt_hessian_interval(distr, 
                        iYleft, iYright, ioffset, itrunc, 
                        iweights)(.parm(beta_nex))
                colnames(ret) <- rownames(ret) <- colnames(Y)
                ### in case beta contains fix parameters,
                ### return all scores
                if (!is.null(fixed)) {
                    if (all(names(fixed) %in% nm))
                        return(ret)
                }

                return(ret[!fix, !fix, drop = FALSE])
            })
        )
    }

    if (all(!is.finite(ci))) {
        ui <- ci <- NULL
    } else {
        ui <- as(ui[is.finite(ci),,drop = FALSE], "matrix")
        ci <- ci[is.finite(ci)]
        r0 <- rowSums(abs(ui)) == 0
        ui <- ui[!r0,,drop = FALSE]
        ci <- ci[!r0]
        if (nrow(ui) == 0) ui <- ci <- NULL
#        ci <- ci + sqrt(.Machine$double.eps) ### we need ui %*% theta > ci, not >= ci
    }

    optimfct <- function(theta, weights, subset = NULL, offset = NULL, 
                         scale = FALSE, optim, ...) {
        of <- .ofuns(weights = weights, subset = subset, 
                     offset = offset, ...)
        loglikfct <- function(beta, weights)  
            -sum(weights * of$ll(beta))
        score <- function(beta, weights, Xmult = TRUE) 
            weights * of$sc(beta, Xmult = Xmult)
        hessian <- function(beta, weights) 
            of$he(beta)
        scorefct <- function(beta, weights) 
            -colSums(score(beta, weights), na.rm = TRUE)
        logliki <- function(beta, weights)
            of$ll(beta)

        if (scale) {
            Ytmp <- Y
            Ytmp[!is.finite(Ytmp)] <- NA
            sc <- apply(abs(Ytmp[, !fix, drop = FALSE]), 2, max, na.rm = TRUE)
            lt1 <- sc < 1.1
            gt1 <- sc >= 1.1
            sc[gt1] <- 1 / sc[gt1]
            sc[lt1] <- 1
            f <- function(gamma) loglikfct(sc * gamma, weights)
            g <- function(gamma) scorefct(sc * gamma, weights) * sc
            theta <- theta / sc
            if (!is.null(ui))
                ui <- t(t(ui) * sc)
        } else {
            f <- function(gamma) loglikfct(gamma, weights)
            g <- function(gamma) scorefct(gamma, weights)
        }
        for (i in 1:length(optim)) {
            ret <- optim[[i]](theta, f, g, ui, ci)
            if (ret$convergence == 0) break()
        }
        if (ret$convergence != 0)
            warning("Optimisation did not converge")
        ### degrees of freedom: number of free parameters, ie #parm NOT meeting the constraints
        ret$df <- length(ret$par)
        ### <FIXME> check on alternative degrees of freedom
#        if (!is.null(ui)) 
#            ret$df <- ret$df - sum(ui %*% ret$par - ci < .Machine$double.eps)
        ### </FIXME>
        if (scale) ret$par <- ret$par * sc

        ### NOTE: this overwrites the functions taking the possibly updated
        ### offset into account
        ret$score <- score
        ret$hessian <- hessian
        ret$loglik <- loglikfct
        ret$logliki <- logliki

        return(ret)
    }

    coef <- rep(NA, length(fix))
    coef[fix] <- fixed
    names(coef) <- colnames(Y)

    ret <- list()
    ret$parm <- .parm
    ret$coef <- coef
    ret$fixed <- fixed
    ret$model <- model
    ret$data <- data
    ret$offset <- offset
    ret$todistr <- todistr
    ret$optimfct <- optimfct

    ret$loglik <- function(beta, weights)  
        -sum(weights * .ofuns(weights = weights, offset = offset)$ll(beta))
    ret$logliki <- function(beta, weights)
        .ofuns(weights = weights, offset = offset)$ll(beta)
    ret$score <- function(beta, weights, Xmult = TRUE) 
        weights * .ofuns(weights = weights, offset = offset)$sc(beta, Xmult = Xmult)
    ret$hessian <- function(beta, weights) 
        .ofuns(weights = weights, offset = offset)$he(beta)

    class(ret) <- c("mlt_setup", "mlt")
    return(ret)
}

.mlt_start <- function(model, data, y, pstart, offset = NULL, fixed = NULL, weights = 1) {

    stopifnot(length(pstart) == nrow(data))

    response <- variable.names(model, "response")
    stopifnot(length(response) == 1)

    eY <- .mm_exact(model, data = data, response = response, object = y)
    iY <- .mm_interval(model, data = data, response = response, object = y)

    if (is.null(eY)) {
        Y <- iY$Yleft
    } else {
        Y <- eY$Y
    }

    ui <- as.matrix(attr(Y, "constraint")$ui)
    ci <- attr(Y, "constraint")$ci

    if (!is.null(fixed)) {
        stopifnot(all(names(fixed) %in% colnames(Y)))
        fix <- colnames(Y) %in% names(fixed)
        fixed <- fixed[colnames(Y)[fix]]
        ui <- ui[,!fix,drop = FALSE]
    } else {
        fix <- rep(FALSE, ncol(Y))
    } 

    X <- matrix(0, nrow = NROW(y), ncol = ncol(Y))
    if (!is.null(eY))
        X[eY$which,] <- as(eY$Y, "matrix")
    if (!is.null(iY))
        X[iY$which,] <- as(iY$Yright, "matrix")
    X[!is.finite(X[,1]),] <- 0

    if (any(fix)) {
        offset <- X[, fix, drop = FALSE] %*% fixed
        X <- X[, !fix, drop = FALSE]
    }

    todistr <- model$todistr
    Z <- todistr$q(pmax(.01, pmin(pstart, .99))) - offset

    X <- X * sqrt(weights)
    Z <- Z * sqrt(weights)

    dvec <- crossprod(X, Z)
    Dmat <- crossprod(X)
    diag(Dmat) <- diag(Dmat) + 1e-08

    if (all(!is.finite(ci))) {
        ui <- ci <- NULL
    } else {
        ui <- as(ui[is.finite(ci),,drop = FALSE], "matrix")
        ci <- ci[is.finite(ci)]
        r0 <- rowSums(abs(ui)) == 0
        ui <- ui[!r0,,drop = FALSE]
        ci <- ci[!r0]
        if (nrow(ui) == 0) ui <- ci <- NULL
        ci <- ci + sqrt(.Machine$double.eps) ### we need ui %*% theta > ci, not >= ci
    }

    if (!is.null(ui)) {    
        ret <- suppressWarnings(try(c(coneproj::qprog(Dmat, dvec, ui, ci, msg = FALSE)$thetahat)))
        if (inherits(ret, "try-error")) {
            diag(Dmat) <- diag(Dmat) + 1e-3
            ret <- c(coneproj::qprog(Dmat, dvec, ui, ci, msg = FALSE)$thetahat)
        }
        ### was: solve.QP(Dmat, dvec, t(ui), ci, meq = 0)$solution
    } else {
        ret <- lm.fit(x = X, y = Z)$coef
    }
    names(ret) <- names(coef(model))[!fix]
    ret
}

.mlt_fit <- function(object, weights, subset = NULL, offset = NULL, 
                     theta = NULL, scale = FALSE, optim, ...) {

    if (is.null(theta))
        stop(sQuote("mlt"), "needs suitable starting values")

    ### BBoptim issues a warning in case of unsuccessful convergence
    ret <- try(object$optimfct(theta, weights = weights, 
        subset = subset, offset = offset, scale = scale, 
        optim = optim, ...))    

    cls <- class(object)
    object[names(ret)] <- NULL
    object <- c(object, ret)
    object$coef[] <- object$parm(ret$par) ### [] preserves names
    object$theta <- theta ### starting value
    object$subset <- subset
    object$scale <- scale ### scaling yes/no
    object$weights <- weights
    object$offset <- offset
    object$optim <- optim
    class(object) <- c("mlt_fit", cls)
    
    return(object)
}

mlt <- function(model, data, weights = NULL, offset = NULL, fixed = NULL,
                theta = NULL, pstart = NULL, scale = FALSE,
                dofit = TRUE, optim = mltoptim(), ...) {

    vars <- as.vars(model)
    response <- variable.names(model, "response")
    responsevar <- vars[[response]]
    ### <FIXME>: how can we check Surv objects?
    if (!inherits(data[[response]], "Surv"))
        stopifnot(check(responsevar, data))
    ### </FIXME>
    bounds <- bounds(responsevar)
    stopifnot(length(response) == 1)
    ### <FIXME> add bounds?
    y <- R(object = data[[response]])
    ### </FIXME>

    if (!.checkR(y, weights)) dofit <- FALSE

    if (is.null(weights)) weights <- rep(1, nrow(data))
    if (is.null(offset)) offset <- rep(0, nrow(data))
    stopifnot(nrow(data) == length(weights))
    stopifnot(nrow(data) == length(offset))

    s <- .mlt_setup(model = model, data = data, y = y, 
                    offset = offset, fixed = fixed) 
    s$convergence <- 1
    if (!dofit) return(s)

    if (is.null(theta)) {
        ### unconditional ECDF, essentially
###        if (is.null(pstart)) pstart <- y$rank / max(y$rank)
        if (is.null(pstart)) pstart <- attr(y, "prob")(weights)(y$approxy) ### y$rank / max(y$rank)
        theta <- .mlt_start(model = model, data = data, y = y, 
                            pstart = pstart, offset = offset, fixed = fixed, weights = weights)
    }

    args <- list()
    args$object <- s
    args$weights <- weights
    args$offset <- offset
    args$theta <- theta
    args$subset <- NULL ### only available in update()
    args$scale <- scale
    args$optim <- optim
    ret <- do.call(".mlt_fit", args)
    ret$call <- match.call()
    ret$bounds <- bounds
    ret$response <- y
    ret
}

update.mlt_fit <- function(object, weights = stats::weights(object), 
                           subset = NULL, offset = object$offset,
                           theta = coef(object, fixed = FALSE), 
                           ...) {

    stopifnot(length(weights) == NROW(object$data))
    if (!is.null(subset))
        stopifnot(is.integer(subset) && 
                  min(subset) >= 1L &&
                  max(subset) <= NROW(object$data))
    args <- list(...)
    if (inherits(object, "mlt_fit")) 
        class(object) <- class(object)[-1L]
    args$object <- object
    if (missing(weights)) {
        args$weights <- stats::weights(object)
    } else {
        args$weights <- weights
    }
    args$subset <- subset
    args$offset <- offset
    args$theta <- theta
    args$scale <- object$scale
    args$optim <- object$optim
    ret <- do.call(".mlt_fit", args)
    ret$call <- match.call()
    ret
}

### frailty models; ie F_Z with additional scale parameter
fmlt <- function(object, frailty = c("Gamma", "InvGauss", "PositiveStable"), 
                 interval = fr$support, ...) {
    frailtyfun <- paste0(".", frailty <- match.arg(frailty), "Frailty")
    fr <- do.call(frailtyfun, list())
    object <- as.mlt(object)
    model <- object$model
    ll <- function(fparm)
        -logLik(update(object, distr = do.call(frailtyfun, list(fparm))))
    om <- optimise(ll, interval = interval, maximum = FALSE)
    model$todistr <- do.call(frailtyfun, list(om$minimum))
    ret <- mlt(model = model, data = object$data, 
               weights = weights(object),
               subset = object$subset, offset = object$offset, 
               theta = coef(object), fixed = object$fixed, 
               scale = object$scale, optim = object$optim)
    class(ret) <- c("fmlt", class(ret))
    ret
}

### cure mixture models (for a constant cure probability)
### e.g. 10.1002/sim.687
cmlt <- function(object, interval = fr$support, ...) {
    object <- as.mlt(object)
    model <- object$model
    ### take F_Z from fitted model
    distr <- get(object$model$todistr$call)
    ### no unknown parameters in frailties as of now
    stopifnot(length(distr$parm()) == 0L)
    fr <- .CureRate(distr = distr)
    ### <FIXME> allow for additional parameters in frailty
    ll <- function(logitrho)
        -logLik(update(object, distr = .CureRate(logitrho, distr = distr)))
    om <- optimise(ll, interval = interval, maximum = FALSE)
    model$todistr <- .CureRate(om$minimum, distr = distr)
    ### </FIXME>
    ret <- mlt(model = model, data = object$data, weights = weights(object),
               subset = object$subset, offset = object$offset, 
               theta = coef(object), fixed = object$fixed, 
               scale = object$scale, optim = object$optim)
    class(ret) <- c("fmlt", class(ret))
    ret
}
