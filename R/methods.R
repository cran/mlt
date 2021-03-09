
as.mlt <- function(object)
    UseMethod("as.mlt")

as.mlt.mlt <- function(object)
    object

weights.mlt <- function(object, ...) {
    if (!is.null(object$weights))
        return(object$weights)
    rep(1, NROW(object$data))
}

coef.mlt <- function(object, fixed = TRUE, ...) {
    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", names(args), " are ignored")
    if (fixed) 
        return(object$coef)
    return(object$par)
}

"coef<-" <- function(object, value)
    UseMethod("coef<-")

"coef<-.mlt" <- function(object, value) {
    cf <- coef(object, fixed = TRUE)
    stopifnot(length(cf) == length(value))
    if (!is.null(names(value)))
        stopifnot(all.equal(names(cf), names(value)))
    object$par <- object$parm(value)
    object$coef[] <- value ### preserve names
    object
}

### use maxLik::hessian?
Hessian <- function(object, ...)
    UseMethod("Hessian")

Hessian.mlt <- function(object, parm = coef(object, fixed = FALSE), ...) {
    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", names(args), " are ignored")
    w <- weights(object)
    if (!is.null(object$subset)) 
        w[-object$subset] <- 0
    object$hessian(parm, weights = w)
}
    
Gradient <- function(object, ...)
    UseMethod("Gradient")

Gradient.mlt <- function(object, parm = coef(object, fixed = FALSE), ...) {
    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", names(args), " are ignored")
    as.vector(colSums(estfun(object, parm = parm)))
}

vcov.mlt <- function(object, parm = coef(object, fixed = FALSE), complete = FALSE, ...) {
    ### <FIXME> implement complete argument </FIXME>
    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", names(args), " are ignored")
    step <- 0
    lam <- 1e-6
    H <- Hessian(object, parm = parm)
    while((step <- step + 1) <= 3) {
        ret <- try(solve(H + (step - 1) * lam * diag(nrow(H))))
        if (!inherits(ret, "try-error")) break
    }
    if (inherits(ret, "try-error"))
        stop("Hessian is not invertible")
    if (step > 1)
        warning("Hessian is not invertible, an approximation is used")
    ret
}

logLik.mlt <- function(object, parm = coef(object, fixed = FALSE), 
                       w = NULL, newdata, ...) {
    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", names(args), " are ignored")
    if (!missing(newdata)) {
        tmpmod <- mlt(object$model, data = newdata, dofit = FALSE)
        coef(tmpmod) <- coef(object)
        return(logLik(tmpmod, parm = parm, w = w))
    }
    if (is.null(w))
        w <- weights(object)
    if (!is.null(object$subset)) {
        ret <- sum(object$logliki(parm, weights = w)[object$subset] * 
                   w[object$subset])
    } else {
        ret <- -object$loglik(parm, weights = w)
    }
    ###    attr(ret, "df") <- length(coef(object, fixed = FALSE))
    attr(ret, "df") <- object$df
    class(ret) <- "logLik"
    ret
}

estfun.mlt <- function(object, parm = coef(object, fixed = FALSE), 
                       w = NULL, newdata, ...) {
    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", names(args), " are ignored")
    if (!missing(newdata)) {
        tmpmod <- mlt(object$model, data = newdata, dofit = FALSE)
        coef(tmpmod) <- coef(object)
        return(estfun(tmpmod, parm = parm, weights = w))
    }
    if (is.null(w))
        w <- weights(object)
    sc <- -object$score(parm, weights = w)
    if (!is.null(object$subset))
        sc <- sc[object$subset,,drop = FALSE]
    return(sc)
}

residuals.mlt <- function(object, parm = coef(object, fixed = FALSE), 
                          w = NULL, newdata, ...) {
    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", names(args), " are ignored")
    if (!missing(newdata)) {
        tmpmod <- mlt(object$model, data = newdata, dofit = FALSE)
        coef(tmpmod) <- coef(object)
        return(resid(tmpmod, parm = parm, w = w))
    }
        if (is.null(w))
        w <- weights(object)
    sc <- -object$score(parm, weights = w, Xmult = FALSE)
    if (!is.null(object$subset))
        sc <- sc[object$subset,,drop = TRUE]
    if (inherits(sc, "matrix"))
        return(sc[, 1L, drop = TRUE])
    return(sc)
}

mkgrid.mlt <- function(object, n, ...)
    mkgrid(object$model, n = n, ...)

mkgrid.ctm <- function(object, n = n, ...)
    mkgrid(object$model, n = n, ...)

variable.names.mlt <- function(object, ...)
    variable.names(object$model, ...)

model <- function(object)
    UseMethod("model")

.one_factor_only <- function(object) {
    f <- inherits(object, "formula_basis")
    v <- as.vars(object)
    f && (length(v) == 1 && inherits(v[[1L]], "factor_var"))
}
    
model.ctm <- function(object) {
    x <- object$bases
    ret <- list(response_trafo = c("continuous", "discrete")[.one_factor_only(x$response) + 1L],
                response_type = sapply(x$response, class),
                response_var = as.vars(x$response),
                interaction_trafo = !is.null(x$interacting),
                shift_trafo = !is.null(x$shifting))
    if (ret$interaction_trafo) {
        ret$interaction_vars = as.vars(x$interacting)
        ret$interaction_type = c("continuous", "discrete")[.one_factor_only(x$interacting) + 1L]
    }
    if (ret$shift_trafo) 
        ret$shift_vars = as.vars(x$shifting) 
    return(ret)
}

model.mlt <- function(object) {
    c(model(object$model), 
      list(todistr = object$todistr$name,
           fixed = object$fixed))
}

### take estimated "logrho" parameter into account in coef and vcov 
coef.fmlt <- function(object, addparm = FALSE, ...) {
    if (addparm)
        return(object$model$todistr$parm())
    class(object) <- class(object)[-1L]
    coef(object, ...)
}

### hessian, partially analytically / numerically
Hessian.fmlt <- function(object, parm = coef(object, fixed = FALSE), ...) {
    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", names(args), " are ignored")
    w <- weights(object)
    if (!is.null(object$subset)) 
        w[-object$subset] <- 0
    addparm <- coef(object, addparm = TRUE)
    parm <- c(addparm, parm)
    model <- object$model
    ll <- function(addparm, parm) {
        model$todistr <- do.call(model$todistr$call, as.list(addparm))
        m <- mlt(model = model, data = object$data, weights = w,
                 subset = object$subset, offset = object$offset, dofit = FALSE,
                 theta = parm,
                 fixed = object$fixed, scale = object$scale, optim = object$optim)
        -logLik(m, parm = parm, w = w)
    }
    sc <- function(addparm, parm, which = 1L) {
        model$todistr <- do.call(model$todistr$call, as.list(addparm))
        m <- mlt(model = model, data = object$data, weights = w,
                 subset = object$subset, offset = object$offset, dofit = FALSE,
                 theta = parm,
                 fixed = object$fixed, scale = object$scale, optim = object$optim)
        ### note: weights(m) are used by estfun
        Gradient(m, parm = parm)[which]
    }

    ### compute the analytical hessian for coef(object) for fixed optimal
    ### logrho analytically
    mltobj <- object
    class(mltobj) <- class(mltobj)[-1L]
    Htheta <- Hessian(mltobj, parm = parm[-1L])
    ret <- matrix(NA, nrow = length(parm), ncol = length(parm))
    rownames(ret) <- colnames(ret) <- names(parm)
    idx <- 1:length(addparm)
    ret[-idx,-idx] <- Htheta

    ### compute the hessian for logrho for fixed coef(object) 
    ### numerically (we don't know the form of the likelihood
    ### here and this is just fast and accurate enough)
    ret[idx, idx] <- numDeriv::hessian(ll, addparm, parm = parm[-idx])

    ### compute the hessian for logrho and coef(object):
    ### first compute the analytical gradient wrt coef(object) and then the
    ### numerical gradient wrt logrho for each component
    ### same reason as above
    Hlogrho_theta <- sapply(1:(length(parm) - length(idx)), 
        function(i) numDeriv::grad(sc, addparm, parm = parm[-idx], which = i))
    ret[idx, -idx] <- ret[-idx, idx] <- Hlogrho_theta
    ret
}

vcov.fmlt <- function(object, parm = coef(object, fixed = FALSE), 
                          complete = FALSE, addparm = FALSE, ...) {
    ### <FIXME> implement complete argument </FIXME>
    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", names(args), " are ignored")
    step <- 0
    lam <- 1e-6
    H <- Hessian(object, parm = parm)
    while((step <- step + 1) <= 3) {
        ret <- try(solve(H + (step - 1) * lam * diag(nrow(H))))
        if (!inherits(ret, "try-error")) break
    }
    if (inherits(ret, "try-error"))
        stop("Hessian is not invertible")
    if (step > 1)
        warning("Hessian is not invertible, an approximation is used")
    if (addparm)
        return(ret)
    idx <- 1:length(coef(object, addparm = TRUE))
    return(ret[-idx, -idx])
}

description <- function(object) {
    stopifnot(inherits(object, "mlt"))
    m <- model(object)
    cond <- m$interaction_trafo || m$shift_trafo
    strat <- m$interaction_trafo && m$interaction_type == "discrete"
    lin <- cond && (!m$interaction_trafo || strat)
    if (lin) pm <- switch(m$todistr, "logistic" = "odds",
                                     "minimum extreme value" = "hazards", "")

    ret <- paste(m$response_trafo, 
       if (!cond) "unconditional",
       if (strat) "stratified",
       if (lin) "linear", 
       if (!lin && cond) "conditional",
       "transformation model", 
       "(transformed", m$todistr, "distribution)")
    if (cond && lin && pm != "")
        ret <- c(ret, paste(m$response_trafo, if (strat) "stratified", "proportional", pm, "model"))
    ret <- gsub("\\s\\s*", " ", ret)
    return(ret)
}

summary.mlt <- function(object, ...) {

    ret <- list(call = object$call,
                convergence = object$convergence,
                type = paste(description(object), collapse = "\n\t"),
                logLik = logLik(object),
#                AIC = AIC(object),
                coef = coef(object))
    class(ret) <- "summary.mlt"
    ret
}

print.summary.mlt <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

    cat("\nCall:\n")
    print(x$call)
    if (x$convergence != 0L)
    cat("\nCould not estimate parameters; optimisation did not converge!\n")
    cat("\nType: ", x$type)
#    cat("\nAIC: ", x$AIC)
    cat("\nLog-Likelihood: ", x$logLik, " (df = ", attr(x$logLik, "df"), ")", sep = "")
    cat("\n")
    cat("\nCoefficients:", x$coef)
    cat("\n\n")
    invisible(x)
}

print.mlt <- function(x, ...)
    print(summary(x, ...))

as.vars.ctm <- function(object)
    as.vars(object$model)

as.vars.mlt <- function(object)
    as.vars(object$model)

bounds.ctm <- function(object)
    bounds(as.vars(object))

bounds.mlt <- function(object)
    bounds(as.vars(object))

print.response <- function(x, digits = getOption("digits"), ...) {

    ac <- function(x) {
        if (inherits(x, "factor")) return(levels(x)[x])
        format(c(x), digits = digits)
    }
    obs <- paste(ifelse(!is.na(x$exact), ac(x$exact), 
                 paste("(", ac(x$cleft), ", ", ac(x$cright), "]", sep = "")))

    if (all(is.na(x$tleft) & is.na(x$tright))) {
        print(obs, quote = FALSE, ...)
        return(invisible(obs))
    }

    trc <- character(length(obs))
    i <- (!is.na(x$tleft) & is.na(x$tright))
    if (sum(i) > 0)
        trc[i] <- paste("| >", ac(x$tleft[i]))
    i <- (is.na(x$tleft) & !is.na(x$tright))
    if (sum(i) > 0)
        trc[i] <- paste("| <", ac(x$tright[i]))
    i <- (!is.na(x$tleft) & !is.na(x$tright))
    if (sum(i) > 0)
        trc[i] <- paste("| (", ac(x$tleft[i]), ", ", ac(x$tright[i]), "]", sep = "")
    ret <- paste("{", obs, trc, "}", sep = "")
    print(ret, quote = FALSE, ...)
}

bread.mlt_fit <- function(x) vcov(x) * nrow(x$data) # sandwich estimator
bread.mlt <- function(x) vcov(x) * nrow(x$data) # sandwich estimator
