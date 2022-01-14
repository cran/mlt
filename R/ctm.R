
ctm <- function(response, interacting = NULL, shifting = NULL, 
                scaling = NULL, scale_shift = FALSE, data = NULL,
                todistr = c("Normal", "Logistic", "MinExtrVal", 
                            "MaxExtrVal", "Exponential"),
                sumconstr = inherits(interacting, c("formula", "formula_basis")), 
                ...) {

    ### mkgrid() will not work if data is missing
    if (.is.formula(response)) 
        response <- as.basis(response, data = data)
    if (.is.formula(interacting)) 
        interacting <- as.basis(interacting, data = data)
    if (.is.formula(shifting)) 
        shifting <- as.basis(shifting, data = data, remove_intercept = TRUE,
                             prefix = ifelse(is.null(scaling), "", "sft_"), ...)
    if (.is.formula(scaling)) 
        scaling <- as.basis(scaling, data = data, remove_intercept = TRUE, 
                            prefix = "scl_", ...)

    if (is.character(todistr))
        todistr <- .distr(todistr)

    bases <- list(response = response, interacting = interacting, 
                  shifting = shifting, scaling = scaling)

    if (!is.null(interacting))
        interacting <- b(iresponse = response, iinteracting = interacting, 
                         sumconstr = sumconstr)

    args <- bases
    names(args) <- paste0("b", names(args))
    args$binteracting <- interacting
    if (!is.null(interacting))
        args$bresponse <- NULL
    args <- args[!sapply(args, is.null)]
    mod <- do.call("c", args)

    ret <- list(model = mod, response = variable.names(response), 
                todistr = todistr, bases = bases, scale_shift = scale_shift)
    class(ret) <- "ctm"
    nd <- lapply(mkgrid(ret, n = 1), function(x) x[1]) ### integer may have more values
    nd <- do.call("expand.grid", nd)
    X <- model.matrix(ret, data = nd)
    cn <- colnames(X)
    if (any(duplicated(cn)))
        warning("Non-unique coefficient names; this might cause problems.")
    cf <- numeric(NCOL(X))
    names(cf) <- colnames(X)
    cf[] <- NA
    ret$coef <- cf
    return(ret)
}

model.matrix.ctm <- function(object, data, ...)
    return(model.matrix(object$model, data = data, ...))

variable.names.ctm <- function(object, 
    which = c("all", "response", "interacting", "shifting", "scaling"), 
    ...) {

    which <- match.arg(which)
    m <- object$bases
    if (which == "all")
        return(variable.names(object$model))
    if (which == "response")
        return(variable.names(m$response))
    if (which == "interacting") {
        if (!is.null(m$interacting))
            return(variable.names(m$interacting))
        return(NULL)
    }
    if (which == "shifting") {
        if (!is.null(m$shifting))
            return(variable.names(m$shifting))
        return(NULL)
    }
    if (which == "scaling") {
        if (!is.null(m$scaling))
            return(variable.names(m$scaling))
        return(NULL)
    }
}

coef.ctm <- function(object, ...)
    object$coef

"coef<-.ctm" <- function(object, value) {
    cf <- coef(object)
    stopifnot(length(cf) == length(value))
    if (!is.null(names(value))) {
        stopifnot(all(names(value) %in% names(cf)))
    } else {
        stopifnot(length(value) == length(cf))
        names(value) <- names(cf)
    }
    object$coef[names(value)] <- value
    object
}
