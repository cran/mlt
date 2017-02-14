
.is.formula <- function(x) {
    if (is.null(x)) return(FALSE)
    inherits(x, "formula")
}

.is.Surv <- function(x)
    inherits(x, "Surv")

.is.R <- function(x)
    inherits(x, "response")

.type_of_response <- function(y) {
    if (.is.Surv(y)) return("survival")
    if (.is.R(y)) {
        if (any(.exact(y))) {
            y <- y$exact[.exact(y)]
        } else {
            y <- y$cleft[.cleft(y)]
        }
    }
    if (storage.mode(y) == "double") return("double")
    if (is.integer(y)) return("integer")
    if (is.ordered(y)) return("ordered")
    if (is.factor(y)) return("unordered")
    return(NA)
}

.checkR <- function(x, weights = NULL) {
    if (!.is.R(x)) return(FALSE)
    if (is.null(weights)) weights <- 1
    if (all(.cleft(x) & weights > 0)) {
        warning("response contains left-censored observations only")
        return(FALSE)
    }
    if (all(.cright(x) & weights > 0)) {
        warning("response contains right-censored observations only")
        return(FALSE)
    }
    return(TRUE)
}
