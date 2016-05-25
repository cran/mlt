
predict.ctm <- function(object, newdata, type = c("trafo", "distribution", "survivor", 
    "density", "logdensity", "hazard", "loghazard", "cumhazard", "quantile"), 
    terms = c("bresponse", "binteracting", "bshifting"), q = NULL, p = NULL, K = 50,
    interpolate = TRUE, ...) {

    type <- match.arg(type)
    terms <- match.arg(terms, several.ok = TRUE)

    if (type == "quantile")
        stopifnot(!is.null(p) && (min(p) > 0 & max(p) < 1))

    if (!is.null(newdata)) {
        if (!is.data.frame(newdata)) {
            if (type != "quantile")
                stopifnot(variable.names(object, "response") %in% names(newdata))
        }
    }

    ret <- switch(type, 
        "trafo" = tmlt(object = object, newdata = newdata, q = q, terms = terms, ...),
        "distribution" = pmlt(object = object, newdata = newdata, q = q),
        "survivor" = smlt(object = object, newdata = newdata, q = q),
        "density" = dmlt(object = object, newdata = newdata, q = q, log = FALSE),
        "logdensity" = dmlt(object = object, newdata = newdata, q = q, log = TRUE),
        "hazard" = hmlt(object = object, newdata = newdata, q = q, log = FALSE),
        "loghazard" = hmlt(object = object, newdata = newdata, q = q, log = TRUE),
        "cumhazard" = Hmlt(object = object, newdata = newdata, q = q),
        "quantile" = qmlt(object = object, newdata = newdata, n = K,
                          p = p, interpolate = interpolate))

    return(ret)
}

predict.mlt <- function(object, newdata = object$data, ...) {
    ctmobj <- object$model
    coef(ctmobj) <- coef(object)
    predict(ctmobj, newdata = newdata, ...)
}
