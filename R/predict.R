
predict.ctm <- function(object, newdata, type = c("trafo", 
                                                  "distribution", "logdistribution", 
                                                  "survivor", "logsurvivor",
                                                  "density", "logdensity", 
                                                  "hazard", "loghazard", 
                                                  "cumhazard", "logcumhazard", 
                                                  "odds", "logodds", 
                                                  "quantile"), 
    terms = c("bresponse", "binteracting", "bshifting"), q = NULL, prob = NULL, K = 50,
    interpolate = FALSE, ...) {

    type <- match.arg(type)
    terms <- match.arg(terms, several.ok = TRUE)

    if (type == "quantile")
        stopifnot(!is.null(prob) && (min(prob) > 0 & max(prob) < 1))

    y <- variable.names(object, "response")
    if (!missing(newdata)) {
        if (!is.data.frame(newdata)) { ### newdata is list with _all_ variables
            stopifnot(is.null(q))
            if (type != "quantile") ### hm, why???
                stopifnot(variable.names(object, "response") %in% names(newdata))
        }
        if (!(y %in% names(newdata)) && is.null(q))
            q <- mkgrid(object, n = K)[[y]]
    } else {
        newdata <- NULL
        if (is.null(q))
            q <- mkgrid(object, n = K)[[y]]
    }

    ret <- switch(type, 
        "trafo" = tmlt(object = object, newdata = newdata, q = q, terms = terms, ...),
        "distribution" = pmlt(object = object, newdata = newdata, q = q, ...),
        "logdistribution" = pmlt(object = object, newdata = newdata, q = q, log = TRUE, ...),
        "survivor" = smlt(object = object, newdata = newdata, q = q, ...),
        "logsurvivor" = smlt(object = object, newdata = newdata, q = q, log = TRUE, ...),
        "density" = dmlt(object = object, newdata = newdata, q = q, ...),
        "logdensity" = dmlt(object = object, newdata = newdata, q = q, log = TRUE, ...),
        "hazard" = hmlt(object = object, newdata = newdata, q = q, ...),
        "loghazard" = hmlt(object = object, newdata = newdata, q = q, log = TRUE, ...),
        "cumhazard" = Hmlt(object = object, newdata = newdata, q = q, ...),
        "logcumhazard" = Hmlt(object = object, newdata = newdata, q = q, log = TRUE, ...),
        "odds" = Omlt(object = object, newdata = newdata, q = q, ...),
        "logodds" = Omlt(object = object, newdata = newdata, q = q, log = TRUE, ...),
        "quantile" = qmlt(object = object, newdata = newdata, q = q, n = K,
                          prob = prob, interpolate = interpolate, ...))

    return(ret)
}

predict.mlt <- function(object, newdata = object$data, ...) {
    ctmobj <- object$model
    coef(ctmobj) <- coef(object)
    predict(ctmobj, newdata = newdata, ...)
}
