
.rbind <- function(x) {
    if (!is.list(x)) return(matrix(x, nrow = 1))
    return(do.call("rbind", x))
}

.models <- function(..., strict = TRUE) {

    ### ... may contain junk when arguments were not fully matched
    ismlt <- sapply(list(...), function(x) inherits(x, "mlt"))
    m <- lapply(list(...)[ismlt], function(x) as.mlt(x))
    # nm <- abbreviate(sapply(m, function(x) x$model$response), 4)
    nm <- sapply(m, function(x) x$model$response)
    J <- length(m)
    Jp <- J * (J - 1) / 2
    normal <- sapply(m, function(x) x$todistr$name == "normal")
  
    w <- lapply(m, weights)
    wlength <- sapply(w, function(x) isTRUE(all.equal(x, w[[1]])))
    if (!all(wlength))
        stop("Number of (non-missing) observations differs between models supplied")
    w <- w[[1L]]
  
    ### determine if response is numeric and was measured exactly
    ### (censoring and discreteness are treated the same here)
    mm <- lapply(m, function(mod) {
      eY <- get("eY", environment(mod$parm))
      iY <- get("iY", environment(mod$parm))
      list(eY = eY, iY = iY)
    })

    cmod <- sapply(mm, function(x) !is.null(x$eY))  
    dmod <- sapply(mm, function(x) !is.null(x$iY))  
    if (strict)
        stopifnot(all(xor(cmod, dmod)))

    ### determine if response is conceptually numeric
    cresp <- sapply(m, function(x) 
        inherits(attr(x$model$bases$response, "variables"), 
                 "continuous_var"))

    ### nobs() gives length(weights[weights > 0])
    nobs <- unique(sapply(m, function(x) NROW(x$data)))
    stopifnot(length(nobs) == 1L)
    nobs <- nobs[[1L]]

    mcf <- lapply(m, function(x) coef(x, fixed = TRUE))
    P <- sapply(mcf, length)
    fpar <- factor(rep(1:J, times = P))

    ### par always includes marginally fixed parameters
    parm <- function(par) {
        mpar <- par[1:sum(P)]
        split(mpar, fpar)
    }

    constr <- lapply(mm, function(m) {
        if (is.null(m$eY)) return(attr(m$iY$Yleft, "constraint"))
        return(attr(m$eY$Y, "constraint"))
    })

    ui <- do.call("bdiag", lapply(constr, function(x) x$ui))
    ci <- do.call("c", lapply(constr, function(x) x$ci))
    ui <- as(ui[is.finite(ci),,drop = FALSE], "matrix")
    ci <- ci[is.finite(ci)]

    mf <- lapply(1:J, function(j) {
        mf <- m[[j]]$data ###model.frame(m[[j]])
        if (cmod[j]) return(mf)
        yl <- m[[j]]$response$cleft
        yr <- m[[j]]$response$cright
        rp <- m[[j]]$model$response
        ml <- mr <- mf
        ml[[rp]] <- yl
        mr[[rp]] <- yr
        return(list(left = ml, right = mr))
    })

    fixed <- vector(mode = "list", length = J)
    for (j in 1:J) {
        if (!is.null(m[[j]]$fixed)) {
            fj <- m[[j]]$fixed
            names(fj) <- paste(nm[j], names(fj), sep = ".")
            fixed[[j]] <- fj
        }
    }

    return(list(models = m, mf = mf, cont = cmod, cresp = cresp, 
                normal = normal, nobs = nobs, weights = w, 
                nparm = P, mcoef = mcf, parm = parm, ui = ui, ci = ci, mm = mm, 
                names = nm, fixed = fixed))
}

.mget <- function(models, j = 1, parm, newdata = NULL,
                  what = c("trafo", "dtrafo", "z", "zleft", 
                           "dzleft", "zright", "dzright", "zprime", 
                           "trafoprime", "estfun", "scaleparm"), ...) {

    what <- match.arg(what)

    if (length(j) > 1) {
        ret <- lapply(j, .mget, models = models, parm = parm, 
                      newdata = newdata,
                      what = what, ...)
        return(ret)
    }

    if (what == "scaleparm") {
        cf <- coef(models$models[[j]], fixed = TRUE)
        parsc <- models$models[[j]]$parsc
        cf[] <- 1
        cf[names(parsc)] <- parsc
        return(cf)
    }

    prm <- models$parm(parm)[[j]]
    ### remove marginally fix parameters
    if (!is.null(models$fixed[[j]]))
        prm <- prm[!names(prm) %in% names(models$fixed[[j]])]
    tmp <- as.mlt(models$models[[j]])
    if (!is.null(newdata)) {
        tmp <- mlt(tmp$model, data = newdata,
                   fixed = tmp$fixed, theta = prm,
                   scaleparm = tmp$scaleparm, dofit = FALSE)
    }


    if (what == "trafoprime") {
        ret <- models$models[[j]]$trafoprime(prm)
        if (models$cont[j])
            return(ret[c("exY", "exYprime")])
        ret$iYleft[!is.finite(ret$iYleft[,1])] <- 0
        ret$iYright[!is.finite(ret$iYright[,1])] <- 0
        return(ret[c("iYleft", "iYright")])
    }

    ret <- tmp$trafo(prm)
    ### extract both exact and interval (former might be needed for
    ### predictions)
    tr <- ret$trex
    trp <- ret$trexprime
    if (!models$normal[j])
        trd <- tmp$todistr$d(tr) * trp
    trl <- ret$trleft
    trr <- ret$trright

    if (what == "trafo") {
        return(tr)
    }
    if (what == "dtrafo") {
        return(tmp$todistr$d(tr))
    }
    if (what == "z") {
        if (models$normal[j]) 
            return(tr)
        return(qnorm(tmp$todistr$p(tr, log.p = TRUE), log.p = TRUE))
    }
    if (what == "zleft") {
        if (models$normal[[j]])
            return(trl)
        return(qnorm(tmp$todistr$p(trl, log.p = TRUE), log.p = TRUE))
    }
    if (what == "dzleft") {
        if (models$normal[[j]])
            return(rep_len(1, length(trl)))
        qn <- qnorm(tmp$todistr$p(trl, log.p = TRUE), log.p = TRUE)
        dn <- dnorm(qn)
        dn[!is.finite(dn)] <- 1
        return(tmp$todistr$d(trl) / dn)
    }
   if (what == "zright") {
        if (models$normal[[j]])
            return(trr)
        return(qnorm(tmp$todistr$p(trr, log.p = TRUE), log.p = TRUE))
    }
    if (what == "dzright") {
        if (models$normal[[j]])
            return(rep_len(1, length(trr)))
        qn <- qnorm(tmp$todistr$p(trr, log.p = TRUE), log.p = TRUE)
        dn <- dnorm(qn)
        dn[!is.finite(dn)] <- 1
        return(tmp$todistr$d(trr) / dn)
    }
    if (what == "zprime") {
        if (models$normal[[j]])
            return(trp)
        qn <- qnorm(tmp$todistr$p(tr, log.p = TRUE), log.p = TRUE)
        return(trd / dnorm(qn))
    }
    if (what == "estfun") {
        return(tmp$scorei(prm, Xmult = TRUE))
    }
}

.MCw <- function(J, M, seed, type = c("MC", "ghalton")) {

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

    type <- match.arg(type)
    if (type == "ghalton") {
        if (requireNamespace("qrng", quietly = TRUE))
            return(t(qrng::ghalton(d = J - 1, n = M)))
        warning("package qrng is not available for ghalton(); using type = MC")
        type <- "MC"
    }
    return(matrix(runif((J - 1) * M), ncol = M))
}

.start <- function(m, xnames) {

    J <- length(m$models)
    Jp <- J * (J - 1) / 2
    Jnames <- m$names
    ### fixed = FALSE?
    margin_par <- do.call("c", lapply(m$models, 
                                      function(mod) coef(as.mlt(mod))))
    names(margin_par) <- paste(rep(Jnames, times = m$nparm), 
                               names(margin_par), sep = ".")

    rn <- rownames(unclass(ltMatrices(1:Jp, names = Jnames, byrow = TRUE)))
    lnames <- paste(rep(rn, each = length(xnames)),
                    rep(xnames, times = length(rn)), sep = ".")
    lambda_par <- rep_len(0, length(lnames))
    names(lambda_par) <- lnames

    start <- c(margin_par, lambda_par)
    return(start)
}


.ll <- function(dim, standardize = TRUE, args = list()) {

    if (length(dim) == 1L)
        dim <- c(dim, 0L)

    cJ <- dim[1L]
    dJ <- dim[2L]

    if (!dJ) {

        cll <- function(obs, Lambda) {

            if (dim(Lambda)[2L] > 1)
                stopifnot(!attr(Lambda, "diag"))

            return(logLik(mvnorm(invchol = Lambda), obs = obs, 
                          standardize = standardize, logLik = FALSE))
        }

        csc <- function(obs, Lambda) {

            if (dim(Lambda)[2L] > 1)
                stopifnot(!attr(Lambda, "diag"))

            ret <- lLgrad(mvnorm(invchol = Lambda), obs = obs, standardize = standardize)
            return(list(Lambda = ret$scale, obs = ret$obs))
       }

       return(list(logLik = cll, score = csc))
    }

    ll <- function(obs = NULL, lower, upper, Lambda) {

        if (dim(Lambda)[2L] > 1)
            stopifnot(!attr(Lambda, "diag"))

        a <- args
        a$object <- mvnorm(invchol = Lambda)
        a$obs <- obs
        a$lower <- lower
        a$upper <- upper
        a$standardize <- standardize
        a$logLik <- FALSE
        return(do.call("logLik", a))
    }

    sc <- function(obs = NULL, lower, upper, Lambda) {

        a <- args
        a$object <- mvnorm(invchol = Lambda)
        a$obs <- obs
        a$lower <- lower
        a$upper <- upper
        a$standardize <- standardize
        ret <- do.call("lLgrad", a)
        ret <- list(Lambda = ret$scale,
                    obs = ret$obs,
                    mean = ret$mean, 
                    lower = ret$lower, 
                    upper = ret$upper)
        return(ret)
    }

    return(list(logLik = ll, score = sc))
}

.mmlt_setup <- function(models, formula = ~ 1, data, conditional = FALSE, dofit = TRUE, 
                        args = list(seed = 1, M = 1000))
{
  
    if (isTRUE(all.equal(formula, ~ 1))) {
        lX <- matrix(1)
        colnames(lX) <- "(Intercept)"
        bx <- NULL
    } else {
        bx <- formula
        if (inherits(formula, "formula")) {
            bx <- as.basis(formula, data, drop.unused.levels = FALSE)
            ### this function is called on small subsets of data
            ### and formula ~ x might not contain all levels
        } 
        lX <- model.matrix(bx, data = data)
        if (conditional)
            warning("Conditional models with covariate-dependent", " ",
                    "correlations are order-dependent")
    }

    parnames <- names(.start(models, colnames(lX)))

    if (conditional && !all(models$normal))
        stop("Conditional models only available", 
             "for marginal probit-type models.")

    cJ <- sum(models$cont)
    dJ <- sum(!models$cont)
    J <- cJ + dJ
    Jp <- J * (J - 1) / 2
    llsc <- .ll(c(cJ, dJ), standardize = !conditional, args)

    if (dJ > 1L) {
        if (is.null(args$w)) {
            args$w <- .MCw(J = dJ, M = args$M, seed = args$seed, args$type)
            args$type <- NULL
        } else {
            if (!is.matrix(args$w)) args$w <- matrix(args$w, nrow = 1)
            if (nrow(args$w) < dJ - 1) stop("incorrect dimension of w")
            if (nrow(args$w) > dJ - 1)
                ### make sure only dJ - 1 columns are present
                args$w <- args$w[-(dJ:nrow(args$w)),,drop = FALSE]
        }
    } else {
        args$type <- NULL
    }

    .Xparm <- function(parm) {
        parm <- parm[-(1:sum(models$nparm))]
        return(matrix(parm, nrow = ncol(lX)))
    }

    LAMBDA <- ltMatrices(matrix(0, nrow = Jp, ncol = nrow(lX)),
                         byrow = TRUE, diag = FALSE, names = models$names) #names(models$models))

    ll <- function(parm, newdata = NULL) {

        if (!is.null(newdata) && !isTRUE(all.equal(formula, ~ 1))) 
            lX <- model.matrix(bx, data = newdata)

        # Lambda <- ltMatrices(t(lX %*% .Xparm(parm)), byrow = TRUE, 
        #                      diag = FALSE, names = models$names) ##names(models$models))
        # saves time in ltMatrices
        Lambda <- LAMBDA
        Lambda[] <- t(lX %*% .Xparm(parm))
        ret <- 0
        if (cJ) {
            z <- .rbind(.mget(models, j = which(models$cont), parm = parm, what = "z", 
                              newdata = newdata))
            rownames(z) <- models$names[which(models$cont)]
            zp <- .rbind(.mget(models, j = which(models$cont), parm = parm, 
                               what = "zprime", newdata = newdata))
            ret <- colSums(.log(zp))
            if (!dJ) return(ret + llsc$logLik(obs = z, Lambda = Lambda))
        }
        if (dJ) {
            lower <- .rbind(.mget(models, j = which(!models$cont), parm = parm, 
                                  what = "zleft", newdata = newdata))
            upper <- .rbind(.mget(models, j = which(!models$cont), parm = parm, 
                                  what = "zright", newdata = newdata))
            rownames(lower) <- rownames(upper) <- models$names[which(!models$cont)]
            if (!cJ)
                return(llsc$logLik(lower = lower, upper = upper, 
                                   Lambda = Lambda))
        }
        return(ret + llsc$logLik(obs = z, lower = lower, upper = upper, 
                                 Lambda = Lambda))
    }

    sc <- function(parm, newdata = NULL) {

        if (!is.null(newdata) && !isTRUE(all.equal(formula, ~ 1))) 
            lX <- model.matrix(bx, data = newdata)

        # Lambda <- ltMatrices(t(lX %*% .Xparm(parm)), byrow = TRUE, 
        #                      diag = FALSE, names = models$names) # names(models$models))
        # saves time in ltMatrices
        Lambda <- LAMBDA
        Lambda[] <- t(lX %*% .Xparm(parm))

        if (cJ) {
            z <- .rbind(.mget(models, j = which(models$cont), parm = parm, what = "z", 
                              newdata = newdata))
            rownames(z) <- models$names[which(models$cont)]
            if (!dJ)
                sc <- llsc$score(obs = z, Lambda = Lambda)
        }
        if (dJ) {
            lower <- .rbind(.mget(models, j = which(!models$cont), parm = parm, 
                                  what = "zleft", newdata = newdata))
            upper <- .rbind(.mget(models, j = which(!models$cont), parm = parm, 
                                  what = "zright", newdata = newdata))
            rownames(lower) <- rownames(upper) <- models$names[which(!models$cont)]
            if (!cJ)
                sc <- llsc$score(lower = lower, upper = upper, 
                                 Lambda = Lambda)
        }
        if (cJ && dJ)
            sc <- llsc$score(obs = z, lower = lower, upper = upper, 
                             Lambda = Lambda)

        ### <FIXME> explain subset </FIXME>
        Lmat <- Lower_tri(sc$Lambda)[rep(1:Jp, each = ncol(lX)), , drop = FALSE]
        if (identical(c(lX), 1)) {
            scL <- Lmat ### NaN might appear in scores
        } else {
            scL <- Lmat * t(lX[,rep(1:ncol(lX), times = Jp), drop = FALSE])
        }
      
        scp <- vector(mode = "list", length = cJ + dJ)

        if (cJ) {
            mm <- lapply(which(models$cont), 
                function(j) .mget(models, j = j, parm = parm, what = "trafoprime", 
                                  newdata = newdata))
            if (all(models$normal)) {
                zp <- .rbind(.mget(models, j = which(models$cont), parm = parm, 
                                   what = "zprime", newdata = newdata))
                scp[which(models$cont)] <- lapply(1:cJ, function(j) {
                    mm[[j]]$exY * c(sc$obs[j,]) + 
                        mm[[j]]$exYprime / c(zp[j,])
                })
            } else {
                dz <- .rbind(.mget(models, j = which(models$cont), parm = parm, 
                                   what = "dtrafo", newdata = newdata))
                ### these are the unweighted score contributions
                ef <- lapply(which(models$cont), 
                             function(j) 
                                 .mget(models, j = j, parm = parm, what = "estfun", 
                                       newdata = newdata))
                ### note: estfun() gives negative weighted gradient of loglik
                scp[which(models$cont)] <- lapply(1:cJ, function(j) {
                    (mm[[j]]$exY * c(sc$obs[j,] + z[j,]) / 
                        c(dnorm(z[j,])) * c(dz[j,])) + ef[[j]]
                })
            }
        }

        if (dJ) {
            mm <- lapply(which(!models$cont), 
                function(j) .mget(models, j = j, parm = parm, what = "trafoprime", 
                                  newdata = newdata))
            if (all(models$normal)) {
                scp[which(!models$cont)] <- lapply(1:dJ, function(j) {
                    mm[[j]]$iYleft * c(sc$lower[j,]) +
                    mm[[j]]$iYright * c(sc$upper[j,])
                })
            } else {
                dzl <- .rbind(.mget(models, j = which(!models$cont), parm = parm, 
                                    what = "dzleft", newdata = newdata))
                dzl[!is.finite(dzl)] <- 0
                dzr <- .rbind(.mget(models, j = which(!models$cont), parm = parm, 
                                    what = "dzright", newdata = newdata))
                dzr[!is.finite(dzr)] <- 0
                scp[which(!models$cont)] <- lapply(1:dJ, function(j) {
                    return((mm[[j]]$iYleft * c(dzl[j,]) * c(sc$lower[j,])) +
                           (mm[[j]]$iYright * c(dzr[j,]) * c(sc$upper[j,])))
                })
            }
        }
        
        ret <- cbind(do.call("cbind", scp), t(scL))
        colnames(ret) <- parnames
        return(ret)
    }


    ret <- list()

    ### N contributions to the log-likelihood, UNWEIGHTED
    ret$logliki <- function(parm, newdata = NULL) ll(parm, newdata = newdata)

    ### sum of log-likelihood contributions, WEIGHTED
    ret$loglik <- function(parm, weights = NULL, newdata = NULL) { 
        if (is.null(weights)) return(sum(ll(parm, newdata = newdata)))
        sum(weights * ll(parm, newdata = newdata))
    }

    ### N contributions to the score function, UNWEIGHTED
    ret$scorei <- function(parm, newdata = NULL, Xmult = TRUE) {
        if (!Xmult) stop("Xmult not implemented")
        sc(parm, newdata = newdata)
    }

    ### N contributions to score function, WEIGHTED
    ret$score <- function(parm, weights = NULL, newdata = NULL, Xmult = TRUE) {
        if (!Xmult) stop("Xmult not implemented")
        if (is.null(weights)) return(sc(parm, newdata = newdata))
        weights * sc(parm, newdata = newdata)
    }

    scl <- rep(apply(abs(lX), 2, max, na.rm = TRUE), times = Jp)
    lt1 <- scl < 1.1
    gt1 <- scl >= 1.1
    scl[gt1] <- 1 / scl[gt1]
    scl[lt1] <- 1
    scl <- c(do.call("c", .mget(models, j = 1:J, parm = NULL, what = "scaleparm")), 
             scl)
    names(scl) <- parnames
    ret$scl <- scl

    ret$dofit <- dofit
    ret$ui <- models$ui
    ret$ci <- models$ci

    ### FIXME: old interface, deprecate
    ret$ll <- ret$loglik
    ret$sc <- function(...) colSums(ret$score(...))
   
    ret$args <- args
    ret$models <- models
    ret$formula <- formula
    ret$bx <- bx
    ret$parnames <- parnames
    ret$parm <- function(par, flat = FALSE) {
        par <- par[parnames]
        if (flat) return(par)
        return(c(models$parm(par), list(.Xparm(par))))
    }
    if (!missing(data))
        ret$data <- data
    ret$names <- models$names
    class(ret) <- c("mmlt_setup", "mmlt")
    ret
}

.mmlt_optimfct <- function(loglik, score, scl, ui, ci, parnames, dofit) {

    optimfct <- function(theta, weights, subset, scaleparm = FALSE, optim, fixed = NULL, ...) {

        eparnames <- parnames
        if (!is.null(fixed)) eparnames <- parnames[!(parnames %in% names(fixed))]
    
        ### negative log-likelihood
        f <- function(par, scl, ...) {
            if (!is.null(fixed)) {
                p <- par
                names(p) <- eparnames
                p <- c(p, fixed)
                par <- p[parnames]
            }
            return(-loglik(par * scl, weights = weights, ...))
        }

        ### gradient of negative log-likelihood
        g <- function(par, scl, ...) {
            if (!is.null(fixed)) {
                p <- par
                names(p) <- eparnames
                p <- c(p, fixed)
                par <- p[parnames]
            }
            ret <- - colSums(score(par * scl, weights = weights, ...) * scl)
            if (is.null(fixed)) return(ret)
            if (is.matrix(ret))
                return(ret[, !parnames %in% names(fixed)])
            return(ret[!parnames %in% names(fixed)])
        }

        ui <- cbind(ui, matrix(0, nrow = nrow(ui), 
                               ncol = length(parnames) - ncol(ui)))
        if (!is.null(fixed)) {
            d <- ui[, parnames %in% names(fixed), drop = FALSE] %*% fixed
            ui <- ui[, !parnames %in% names(fixed), drop = FALSE]
            ci <- ci - d
        }
        if (!scaleparm) scl[] <- 1
        start <- theta / scl[eparnames]
        ### check if any non-fixed parameters come with constraints
        if (any(abs(ui) > 0) && any(ci > -Inf)) {
            ui <- t(t(ui) * scl[eparnames])
        } else {
            ui <- ci <- NULL
        }

        if (dofit) {
            for (i in 1:length(optim)) {
                ret <- optim[[i]](theta = start, 
                                  f = function(par) f(par, scl = scl), 
                                  g = function(par) g(par, scl = scl), 
                                  ui = ui, ci = ci)
                if (ret$convergence == 0) break()
            }
            if (ret$convergence != 0)
                warning("Optimisation did not converge")
        } else {
            ret <- list(par = theta, value = NA, ### f(theta, scl = 1), (problem with dofit = FALSE)
                        convergence = NA, optim_hessian = NA)
        }
        names(ret$par) <- eparnames
        ret$par[eparnames] <- ret$par[eparnames] * scl[eparnames]
        return(ret)
    }
    return(optimfct)
}


.mmlt_fit <- function(object, weights, subset = NULL, theta = NULL, 
                      scaleparm = FALSE, optim, fixed = NULL, ...)
{

    if (!is.null(fixed)) 
        stopifnot(all(names(fixed) %in% object$parnames))

    if (is.null(theta))
        stop(sQuote("mlt"), "needs suitable starting values")

    optimfct <- .mmlt_optimfct(loglik = object$loglik, 
                               score = object$score, scl = object$scl, 
                               ui = object$ui, ci = object$ci, parnames = object$parnames, 
                               dofit = object$dofit)

    ### BBoptim issues a warning in case of unsuccessful convergence
    ret <- try(optimfct(theta, weights = weights, fixed = fixed,
        subset = subset, scaleparm = scaleparm, optim = optim, ...))

    cls <- class(object)
    object[names(ret)] <- NULL
    object <- c(object, ret)
    object$coef <- ret$par
    object$fixed <- fixed
    object$theta <- theta ### starting value
    object$subset <- subset
    object$scaleparm <- scaleparm ### scaling yes/no
    object$weights <- weights
    object$optim <- optim
    ll <- object$loglik
    object$loglik <- function(parm, ...) {
        if (!is.null(fixed)) {
            eparnames <- object$parnames[!object$parnames %in% names(fixed)]
            stopifnot(length(parm) == length(eparnames))
            p <- parm
            names(p) <- eparnames
            p <- c(p, fixed)
            parm <- p[object$parnames]
        }
        ll(parm = parm, ...)
    }
    sc <- object$score
    object$score <- function(parm, ...) {
        if (!is.null(fixed)) {
            eparnames <- object$parnames[!object$parnames %in% names(fixed)]
            stopifnot(length(parm) == length(eparnames))
            p <- parm
            names(p) <- eparnames
            p <- c(p, fixed)
            parm <- p[object$parnames]
        }
        sc(parm = parm, ...)
    }
    class(object) <- c("mmlt_fit", cls)

    return(object)
}


mmlt <- function(..., formula = ~ 1, data, conditional = FALSE, 
                 theta = NULL, fixed = NULL, scaleparm = FALSE,
                 optim = mltoptim(hessian = TRUE),  ### provides hessian
                 args = list(seed = 1, M = 1000), 
                 dofit = TRUE, domargins = TRUE)
{
  
    call <- match.call()

    if (conditional && !domargins)
        stop("Conditional models must fit marginal and joint parameters.")

    models <- .models(..., strict = FALSE)
    x <- numeric(models$nobs)
    xb <- lapply(models$mm, function(m) {
        if (is.null(m$iY)) return(x)
        x[m$iY$which] <- 1
        return(x)
    })
    cdpat <- do.call("interaction", xb)[,drop = TRUE]

    if (nlevels(cdpat) > 1L) {
        mm <- vector(mode = "list", length = nlevels(cdpat))
        names(mm) <- levels(cdpat)
        for (j in names(mm)) {
            idx <- which(j == cdpat)
            tmp <- data[idx,,drop = FALSE]
            nm <- lapply(models$models, function(mod) {
                ret <- mlt(mod$model, data = tmp, theta = coef(as.mlt(mod), fixed = FALSE), 
                           fixed = mod$fixed, 
                           scaleparm = FALSE, ### because dofit = FALSE
                           weights = mod$weights[idx],
                           offset = mod$offset[idx], dofit = FALSE)
                ret
            })
            sargs <- list(models = do.call(".models", nm))
            sargs$formula <- formula
            sargs$data <- tmp
            sargs$conditional <- conditional
            sargs$dofit <- dofit
            sargs$args <- args
            mm[[j]] <- do.call(".mmlt_setup", sargs)
        }
        idx <- match(1:length(cdpat), do.call("c", split(1:length(cdpat), cdpat)))
        ret <- mm[[1L]]
        ret$data <- data
        ret$logliki <- function(parm, newdata = NULL) {
            if (!is.null(newdata))
                stop("newdata not implemented")
            ll <- do.call("c", lapply(mm, function(m) m$logliki(parm)))
            ll <- ll[idx]
            return(ll)
        }
        ret$scorei <- function(parm, newdata = NULL) {
            if (!is.null(newdata))
                stop("newdata not implemented")
            sc <- do.call("rbind", lapply(mm, function(m) m$scorei(parm)))
            sc <- sc[idx,,drop = FALSE]
            return(sc)
        }
        ret$loglik <- function(parm, weights, ...)
            sum(weights * ret$logliki(parm, ...))
        ret$score <- function(parm, weights, ...)
            weights * ret$scorei(parm, ...)
    } else {
        ret <- .mmlt_setup(models = models, formula = formula, dofit = dofit,
                           data = data, conditional = conditional, 
                           args = args)
    }

    mfixed <- NULL
    if (!is.null(models$fixed)) {
        mfixed <- do.call("c", models$fixed)
    }

    if (!is.null(fixed)) {
        stopifnot(all(!names(fixed) %in% names(mfixed)))
    }
    fixed <- c(mfixed, fixed)

    ### compute starting values for lambda
    if (is.null(theta) && conditional) {
        cl <- match.call()
        cl$conditional <- FALSE
        cl$domargins <- FALSE
        sm <- eval(cl, parent.frame())
        theta <- coef(sm, fixed = TRUE)
        if (conditional) {
            ### theta are conditional parameters, scale with sigma
            class(sm)[1] <- "cmmlt" ### do NOT standardize Lambda
            d <- rowMeans(mvtnorm::diagonals(coef(sm, newdata = data, 
                                             type = "Sigma")))
            theta[1:sum(models$nparm)] <- 
                theta[1:sum(models$nparm)] * rep(sqrt(d), times = models$nparm)
        }
    } 

    if (is.null(theta) || !domargins) {
        ### starting values for lambda
        if (is.null(theta))
            theta <- numeric(length(ret$parnames))
        names(theta) <- ret$parnames
        mpar <- do.call("c", models$mcoef)
        theta[ret$parnames[1:length(mpar)]] <- mpar
        ### if data is continuous for all models and
        ### lambda is constant, compute starting values analytically without
        ### paying attention to fixed parameters
        av <- 1
        if (inherits(formula, "formula"))
            av <- length(all.vars(formula))
        if (nlevels(cdpat) == 1 && all(models$cont) 
            && av == 0L) {
            J <- length(models$models)
            Z <- do.call("cbind", .mget(models, j = 1:J, parm = mpar, what = "trafo"))
            N <- NROW(Z)
            S <- var(Z) * (N - 1) / N
            L <- try(solve(t(chol(S))))
            if (!inherits(L, "try-error")) {
                L <- L %*% diag((1 / diag(L)))
                Lambda <- ltMatrices(L[lower.tri(L, diag = FALSE)], byrow = FALSE, diag = FALSE)
                Lambda <- ltMatrices(Lambda, byrow = TRUE)
                ### starting values for the starting values computation ;-)
                theta[-(1:length(mpar))] <- c(Lower_tri(Lambda, diag = FALSE))
            }
        } 
        mtheta <- theta
        mfixed <- fixed
        if (!is.null(fixed))
            theta[names(fixed)] <- fixed
        nfixed <- unique(c(ret$parnames[1:length(mpar)], names(fixed)))
        sfixed <- theta[nfixed]
        stheta <- theta[!(names(theta) %in% nfixed)]
        if (length(theta)) {
            mret <- .mmlt_fit(ret, weights = models$weights, 
                             subset = subset, fixed = sfixed,  optim = optim, theta = stheta)
            class(mret) <- c(ifelse(conditional, "cmmlt", "mmmlt"), "mmlt")
            mret$mmlt <- "Multivariate Conditional Transformation Model"
            mret$call <- call
            if (!domargins) return(mret)
            theta <- coef(mret, fixed = TRUE)
        } else {
            ### all parameters were fixed (= all lambda parameters)
            theta <- mtheta
            fixed <- mfixed
        }       
    }

    ret <- .mmlt_fit(ret, weights = models$weights, subset = subset, optim = optim,
                     theta = theta[!names(theta) %in% names(fixed)], fixed = fixed)
    class(ret) <- c(ifelse(conditional, "cmmlt", "mmmlt"), "mmlt")
    ret$mmlt <- "Multivariate Conditional Transformation Model"
    ret$call <- call
    return(ret)
}


.coef.mmlt <- function(object, newdata,
                       type = c("all", "Lambdapar", "Lambda", "Lambdainv", 
                                "Precision", "PartialCorr", "Sigma", "Corr", 
                                "Spearman", "Kendall"), 
                       fixed = TRUE, ...)
{
  
    type <- match.arg(type)
    cf <- c(object$par, object$fixed)[object$parnames]
    if (type == "all") {
        if (!fixed) return(object$par)
        return(cf)
    }

    if (type == "Spearman")
        return(6 * asin(coef(object, newdata = newdata, type = "Cor") / 2) / pi)
  
    if (type == "Kendall")
        return(2 * asin(coef(object, newdata = newdata, type = "Cor")) / pi)

    prm <- object$parm(cf)
    prm <- prm[[length(prm)]]

    if (missing(newdata) || is.null(object$bx)) {
        if (NROW(prm) > 1L && type != "Lambda")
            stop("newdata not specified")
        ret <- ltMatrices(t(prm), byrow = TRUE, diag = FALSE, 
                          names = object$names)
    } else {
        X <- model.matrix(object$bx, data = newdata)
        ret <- ltMatrices(t(X %*% prm), byrow = TRUE, diag = FALSE, 
                          names = object$names)
    }

    if (inherits(object, "mmmlt")) {
        ret0 <- ret
        ret <- mvtnorm::standardize(invchol = ret)
    }

    ret <- switch(type, "Lambdapar" = ret0, 
                        "Lambda" = ret,
                        "Lambdainv" = solve(ret),
                        "Precision" = invchol2pre(ret),
                        "PartialCorr" = invchol2pc(ret),
                        "Sigma" = invchol2cov(ret),
                        "Corr" = invchol2cor(ret))
    return(ret)
}

coef.cmmlt <- function(object, newdata,
                       type = c("all", "conditional", "Lambdapar", "Lambda", 
                                "Lambdainv", "Precision", "PartialCorr", 
                                "Sigma", "Corr", "Spearman", "Kendall"), 
                       fixed = TRUE,
                       ...)
{

    type <- match.arg(type)
    if (type == "conditional") {
        cf <- c(object$par, object$fixed)[object$parnames]
        if (!fixed) return(object$par)
        prm <- object$parm(cf)
        ret <- prm[-length(prm)]
        names(ret) <- object$names
        return(ret)
    }
    return(.coef.mmlt(object = object, newdata = newdata, type = type, 
                      fixed = fixed, ...))
}

coef.mmmlt <- function(object, newdata,
                       type = c("all", "marginal", "Lambdapar", "Lambda", 
                                "Lambdainv", "Precision", "PartialCorr", 
                                "Sigma", "Corr", "Spearman", "Kendall"), 
                       fixed = TRUE,
                       ...)
{

    type <- match.arg(type)
    if (type == "marginal") {
        cf <- c(object$par, object$fixed)[object$parnames]
        if (!fixed) return(object$par)
        prm <- object$parm(cf)
        ret <- prm[-length(prm)]
        names(ret) <- object$names
        return(ret)
    }
    return(.coef.mmlt(object = object, newdata = newdata, type = type, 
                      fixed = fixed, ...))
}

"coef<-.mmlt" <- function(object, value) {
    cf <- coef(object, fixed = TRUE)
    stopifnot(length(cf) == length(value))
    if (!is.null(names(value)))
        stopifnot(all.equal(names(cf), names(value)))
    object$par <- object$parm(value, flat = TRUE)
    object
}

Hessian.mmlt <- function(object, parm = coef(object, fixed = FALSE), ...) {
    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", names(args), " are ignored")
    H <- object$optim_hessian
    if (is.null(H) || !isTRUE(all.equal(parm, coef(object, fixed = FALSE)))) {
        if (requireNamespace("numDeriv")) {
            H <- numDeriv::hessian(function(par) -logLik(object, parm = par), parm)
        } else {
            stop("Hessian not available")
        }
    }
    return(H)
}

Gradient.mmlt <- Gradient.mlt

vcov.mmlt <- function(object, parm = coef(object, fixed = FALSE), 
                      complete = FALSE, ...) {
    step <- 0
    lam <- 1e-6

    if (complete) stop("argument complete not implemented")

    H <- Hessian(object, parm = parm, ...)

    ### <NOTE> add an option to compute vcov for selected 
    ### parameters (eg marginal effects) only and use Schur
    while((step <- step + 1) <= 3) {
          ret <- try(solve(H + (step - 1) * lam * diag(nrow(H))))
          if (!inherits(ret, "try-error")) break
    }
    if (inherits(ret, "try-error"))
        stop("Hessian is not invertible")
    if (step > 1)
        warning("Hessian is not invertible, an approximation is used")
    rownames(ret) <- colnames(ret) <- names(parm)
    ret
}

logLik.mmlt <- function (object, parm = coef(object, fixed = FALSE), w = NULL, newdata = NULL, ...) 
{
    args <- list(...)
    if (length(args) > 0) 
        warning("Arguments ", names(args), " are ignored")
    if (is.null(w) && is.null(newdata))
        w <- weights(object)
    ret <- object$loglik(parm, newdata = newdata, weights = w)
    attr(ret, "df") <- length(object$par)
    class(ret) <- "logLik"
    ret
}

estfun.mmlt <- function(x, parm = coef(x, fixed = FALSE), 
                        w = NULL, newdata = NULL, ...) {
    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", names(args), " are ignored")
    if (is.null(w) && is.null(newdata))
        w <- weights(x)
    return(-x$score(parm, newdata = newdata, weights = w))
}

weights.mmlt <- weights.mlt

summary.mmlt <- function(object, ...) {

    ret <- list(call = object$call,
                convergence = object$convergence,
                type = "Multivariate transformation model",
                logLik = logLik(object),
#                AIC = AIC(object),
                coef = coef(object))
    class(ret) <- "summary.mlt"
    ret
}

print.mmlt <- function(x, ...) {
    cat("\n", "Multivariate conditional transformation model", "\n")
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(coef(x))
    invisible(x)
}

predict.mmlt <- function (object, newdata, margins = 1:J, 
    type = c("trafo", "distribution", "survivor", "density", "hazard"), 
    log = FALSE, args = object$args, ...) 
{
    J <- length(object$models$models)
    margins <- sort(margins)
    stopifnot(all(margins %in% 1:J))

    if (length(margins) == 1L) {
        ### ... may carry q = something
        tmp <- object$models$models[[margins]]
        cf <- coef(tmp, fixed = TRUE)
        ncf <- names(cf)
        names(cf) <- paste(variable.names(tmp)[1L], names(cf), sep = ".")
        cfm <- object$models$parm(coef(object, fixed = TRUE))[[margins]]
        cf[names(cfm)] <- cfm
        names(cf) <- ncf
        coef(tmp) <- cf
        ### marginal models
        if (!inherits(object, "cmmlt")) {
            ret <- predict(tmp, newdata = newdata, type = type, log = log, ...)
            return(ret)
        }
        ### conditional models
        mcov <- coef(object, newdata = newdata, type = "Sigma")
        msd <- sqrt(mvtnorm::diagonals(mcov)[margins,])
        if (length(unique(msd)) == 1L && 
            !"bscaling" %in% names(tmp$model$model)) { ### no stram model
            cf <- cf / unique(msd)
            coef(tmp) <- cf
            ret <- predict(tmp, newdata = newdata, type = type, log = log, ...)
            return(ret)
        }
        type <- match.arg(type)
        tr <- predict(tmp, newdata = newdata, type = "trafo", ...) 
        msd <- matrix(msd, nrow = NROW(tr), ncol = NCOL(tr), byrow = TRUE)
        tr <- tr / msd
        switch(type, "trafo" = return(tr),
                     "distribution" = return(pnorm(tr, log.p = log)),
                     "survivor" = return(pnorm(tr, log.p = log, 
                                               lower.tail = FALSE)),
                     "density" = {
                         dx <- 1
                         names(dx) <- variable.names(tmp)[1L]
                         dtr <- predict(tmp, newdata = newdata, 
                                        type = "trafo", deriv = dx, ...)
                         ret <- dnorm(tr, log = TRUE) - .log(msd) + .log(dtr)
                         if (log) return(ret)
                         return(exp(ret))
                     },
                     stop("not yet implemented"))
    }

    type <- match.arg(type)
    ### don't feed ...
    z <- .mget(object$models, margins, parm = coef(object, fixed = TRUE),
               newdata = newdata, what = "z")
    z <- .rbind(z)

    if (type == "trafo") {
        stopifnot(!log)
        L <- coef(object, newdata = newdata, type = "Lambda")
        if (length(margins) != J) 
            L <- marg_mvnorm(invchol = L, which = margins)$invchol
        return(Mult(L, z))
    }
    if (type == "distribution") {
        lower <- matrix(-Inf, ncol = ncol(z), nrow = nrow(z))
        upper <- z
        Linv <- coef(object, newdata = newdata, type = "Lambdainv")
        if (length(margins) != J) 
            Linv <- marg_mvnorm(chol = Linv, which = margins)$chol
        a <- args
        a$lower <- lower
        a$upper <- upper
        a$logLik <- FALSE
        a$chol <- Linv
        ret <- do.call("lpmvnorm", a)
        if (log) return(ret)
        return(exp(ret))
    }
    if (type == "survivor") {
        lower <- z 
        upper <- matrix(Inf, ncol = ncol(z), nrow = nrow(z))
        Linv <- coef(object, newdata = newdata, type = "Lambdainv")
        if (length(margins) != J) 
            Linv <- marg_mvnorm(chol = Linv, which = margins)$chol
        a <- args
        a$lower <- lower
        a$upper <- upper
        a$logLik <- FALSE
        a$chol <- Linv
        ret <- do.call("lpmvnorm", a)
        if (log) return(ret)
        return(exp(ret))
    }
    stopifnot(type == "density")
    stopifnot(all(object$models$cresp))
    zprime <- .mget(object$models, margins, parm = coef(object, fixed = TRUE),
                    newdata = newdata, what = "zprime")
    if (length(margins) > 1L) {
        zprime <- .rbind(zprime)
    } else {
        zprime <- matrix(zprime, nrow = 1)
    }
    L <- coef(object, newdata = newdata, type = "Lambda")
        if (length(margins) != J) 
            L <- marg_mvnorm(invchol = L, which = margins)$invchol
    ret <- ldmvnorm(obs = z, invchol = L, logLik = FALSE)
    ret <- ret + colSums(.log(zprime))
    if (log) return(ret)
    return(exp(ret))
}

simulate.mmlt <- function(object, nsim = 1L, seed = NULL, newdata, K = 50, 
                          ...) {

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

    args <- list(...)
    if (length(args) > 0L)
        stop("argument(s)", paste(names(args), collapse = ", "), "ignored")

    if (nsim > 1L) 
        return(replicate(nsim, simulate(object, newdata = newdata, K = K, ...), 
                         simplify = FALSE))

    J <- length(object$models$models)
    L <- coef(object, newdata = newdata, type = "Lambda")
    N <- nrow(newdata)

    Z <- matrix(rnorm(J * N), ncol = N)
    Ztilde <- solve(L, Z)

    ret <- matrix(0.0, nrow = N, ncol = J)

    if (inherits(object, "cmmlt")) {
        for (j in 1:J) {
            tmp <- object$models$models[[j]]
            q <- mkgrid(tmp, n = K)[[1L]]
            cf <- coef(tmp, fixed = TRUE)
            ncf <- names(cf)
            names(cf) <- paste(variable.names(tmp)[1L], names(cf), sep = ".")
            cfm <- object$models$parm(coef(object, fixed = TRUE))[[j]]
            cf[names(cfm)] <- cfm
            names(cf) <- ncf
            coef(tmp) <- cf
            pr <- predict(tmp, newdata = newdata, type = "trafo", q = q)
            if (!is.matrix(pr)) 
                pr <- matrix(pr, nrow = length(pr), ncol = NROW(newdata))
            ret[,j] <- as.double(.invf(tmp, f = t(pr), q = q, 
                                       z = t(Ztilde[j,,drop = FALSE])))
        }
    } else {
        Ztilde <- pnorm(Ztilde, log.p = TRUE)
        for (j in 1:J) {
            tmp <- object$models$models[[j]]
            q <- mkgrid(tmp, n = K)[[1L]]
            cf <- coef(tmp, fixed = TRUE)
            ncf <- names(cf)
            names(cf) <- paste(variable.names(tmp)[1L], names(cf), sep = ".")
            cfm <- object$models$parm(coef(object, fixed = TRUE))[[j]]
            cf[names(cfm)] <- cfm
            names(cf) <- ncf
            coef(tmp) <- cf
            pr <- predict(tmp, newdata = newdata, type = "logdistribution", q = q)
            if (!is.matrix(pr)) 
                pr <- matrix(pr, nrow = length(pr), ncol = NROW(newdata))
            ret[,j] <- as.double(.invf(tmp, f = t(pr), q = q, 
                                       z = t(Ztilde[j,,drop = FALSE])))
        }
    }
    colnames(ret) <- variable.names(object, response_only = TRUE)
    return(ret)
}

variable.names.mmlt <- function(object, response_only = FALSE, ...) {

    if (response_only)
        return(sapply(object$models$models, function(x) variable.names(x)[1L]))
    vn <- unique(c(sapply(object$models$models, function(x) variable.names(x)), 
                 all.vars(object$formula)))
    return(vn)
}
    
mkgrid.mmlt <- function(object, ...) {

    lx <- mkgrid(as.basis(object$formula, data = object$data), ...)
    grd <- do.call("c", lapply(object$models$models, mkgrid, ...))
    grd <- c(grd, lx)
    do.call("expand.grid", grd[unique(names(grd))])
}
