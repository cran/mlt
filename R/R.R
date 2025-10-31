
R <- function(object, ...)
    UseMethod("R")

R.Surv <- function(object, as.R.ordered = FALSE, as.R.interval = FALSE, ...) {

    type <- attr(object, "type")
    stopifnot(type %in% c("left", "right", "interval", 
                          "interval2", "counting"))

    ### missing values in Surv() objects
    ### recode such that the log-likelihood contribution is zero
    if (any(ona <- is.na(object[, "status"]))) {
        ### attr(Surv(..., type = "interval2"), "type") == "interval"
        stopifnot(type %in% c("interval", "right", "left"))
        if (type == "interval") {
            ### recode as -Inf, Inf
            object[ona, "status"] <- 3 ### "interval"
            object[ona, "time1"] <- -Inf
            object[ona, "time2"] <- Inf
        }
        if (type == "right") {
            ### recode as -Inf, Inf
            object[ona, "status"] <- 0 ### "censoring"
            object[ona, "time"] <- -Inf
        }
        if (type == "left") {
            ### recode as -Inf, Inf
            object[ona, "status"] <- 0 ### "censoring"
            object[ona, "time"] <- Inf
        }
    }        
    status <- object[, "status"]

    if (as.R.interval) {
        if (as.R.ordered) warning("argument as.R.ordered is ignored")
        ret <- R(object, as.R.ordered = TRUE)
        utm <- attr(ret, "unique_obs")
        ### if utm contains 0, stay with 0 (users can chance data if they
        ### want to). Otherwise, use sqrt(.Machine$double.eps) as smallest 
        ### lower bound to facilitate log_first = TRUE in models
        ### if min(utm) < 0 (not really a survival time), use -Inf
        ### as in R(<numeric>, as.R.interval = TRUE)
        if (min(utm) > 0) {
            utm <- c(sqrt(.Machine$double.eps), utm, Inf)
        } else {
            if (min(utm) < 0) {
                utm <- c(-Inf, utm, Inf)
            } else {
                utm <- c(0, utm, Inf)
            }
        }
        left <- unclass(ret$cleft) + 1L
        left[is.na(left)] <- 1L
        right <- unclass(ret$cright) + 1L
        right[is.na(right)] <- nlevels(ret$cright) + 2L
        obj <- numeric(length(left))
        obj[] <- NA
        if (!is.ordered(ret$tleft)) {
            ret <- R(object = obj, cleft = utm[left], 
                     cright = utm[right])
        } else {
            tleft <- unclass(ret$tleft) + 1L
            tleft[is.na(tleft)] <- 1L
            ret <- R(object = obj, cleft = utm[left], 
                     cright = utm[right], tleft = utm[tleft])
        }
        return(ret)
    }

    if (as.R.ordered && !type %in% c("right", "counting"))
        stop("as.R.ordered only implemented for right-censored observations")
    if (as.R.ordered && type %in% c("right", "counting")) {
         ### code response as ordered factor with right-censoring
         ### this defines the nonparametric likelihood
         ### for right-censored data in terms of the observed event times
         tm <- if(type == "right") object[,"time"] else object[, "stop"]
         ### observed event times
         utm <- sort(unique(tm[status == 1]))
         if (all(tm[status == 0] < utm[length(utm)]))
             utm <- utm[-length(utm)]
         ### convert to ordered factor
         ct <- cut(tm, breaks = c(-Inf, utm, Inf), ordered_result = TRUE)
         ### events in category k contribute
         ### Prob(k) - Prob(k - 1)
         lf <- rg <- ct
         lf[status == 1] <- rg[status == 1] <- NA
         ### censored obs in category k contribute
         ### 1 - Prob(k - 1)
         rg[status == 0] <- levels(ct)[nlevels(ct)]
         ### Note: Censoring before first event contributes
         ### 1 - 0 = 1 (coded as interval with cleft = NA, cright = NA)
         ### handle censoring times tied to event times separately
         idx <- which(status != 1)
         idx <- idx[!(tm[idx] %in% utm)]	### censoring tied to event
         lf[idx] <- c(NA, levels(ct))[lf[idx]]
         ### left truncation
         tl <- NA
         if (type == "counting")
             tl <- cut(object[, "start"], breaks = c(-Inf, utm, Inf), 
                       ordered_result = TRUE)
         ### is this a "response" representation of an ordered factor now
         ret <- R(object = ct, cleft = lf, cright = rg, tleft = tl)
         attr(ret, "unique_obs") <- utm
         return(ret)
    }

    ret <- switch(type,
        "right" = R(object = ifelse(status == 1, object[, "time"], NA),
                    cleft = ifelse(status != 1, object[, "time"], NA),
                    cright = ifelse(status != 1, Inf, NA)),
        "left" =  R(object = ifelse(status == 1, object[, "time"], NA),
                    cleft = ifelse(status != 1, -Inf, NA),
                    cright = ifelse(status != 1, object[, "time"], NA)),

        "interval2" = {
            ret <- cbind(left = object[, "time1"], 
                         right = object[, "time2"])
            ret$left[is.na(ret$left)] <- -Inf
            ret$right[is.na(ret$right)] <- Inf
            R(cleft = ret$left, cright = ret$right)
        },
        "interval" = {
            status <- factor(status, levels = 0:3, 
                             labels = c("right", "exact", "left", "interval"))
            tmp <- matrix(NA, nrow = nrow(object), ncol = 2)
            colnames(tmp) <- c("left", "right")
            for (s in levels(status)) {
                idx <- which(status == s)
                tmp[idx, ] <- switch(s, 
                    "right" = cbind(object[idx, "time1"], Inf),
                    "exact" = cbind(object[idx, "time1"], NA),
                    "left" = cbind(-Inf, object[idx, "time1"]),
                    "interval" = object[idx, c("time1", "time2")])
            }
            R(object = ifelse(is.na(tmp[, "right"]), tmp[, "left"], NA), 
              cleft = ifelse(is.na(tmp[, "right"]), NA, tmp[, "left"]), 
              cright = tmp[, "right"])
        },
        ### left truncation, right censoring
        "counting" = R(object = ifelse(status == 1, object[, "stop"], NA),
                       cleft = ifelse(status != 1, object[, "stop"], NA),
                       cright = ifelse(status != 1, Inf, NA),
                       tleft = object[, "start"])
    )
    attr(ret, "prob") <- function(weights) {
        if (attr(object, "type") == "interval") {
            ### Turnbull estimation takes tooo looong
            ### maybe fall back to Kaplan-Meier?
            if (FALSE) {
            time1 <- object[,"time1"]
            time2 <- object[,"time2"]
            status <- object[,"status"]
            ### status = 3 means interval
            time1[status == 3L] <- time1[status == 3L] + 
                (time2[status == 3] - time1[status == 3]) / 2
            status[status == 3] <- 1L
            ### status = 2 is left censoring: FIXME
            status[status == 2] <- 0L
            object <- Surv(time = time1, event = status)
            }
        }
        ### FIXME: subsetting doesn't change this fct
        if (length(weights) == NROW(object)) {
            sf <- survival::survfit(object ~ 1, subset = weights > 0, weights = weights)
        } else {
            sf <- survival::survfit(object ~ 1)
        }
        function(y) {
            uy <- sort(unique(y))
            s <- summary(sf, times = uy)$surv
            if (length(s) < length(uy))
                s <- c(s, rep_len(0, length(uy) - length(s)))
            s[is.na(s)] <- 0
            p <- 1 - s
            p[match(y, uy)]
        }
    }
    ret
}

R.factor <- function(object, ...) {

    if (nlevels(object) > 2)
        warning("response is unordered factor;
                 results may depend on order of levels")
    return(R(ordered(object, levels = levels(object), 
                     labels = levels(object)), ...))
}

R.ordered <- function(object, cleft = NA, cright = NA, ...) {

    lev <- levels(object)
    ret <- .mkR(exact = object, cleft = cleft, cright = cright, ...)
    ret[is.na(ret$cright), "cright"] <- ret$exact[is.na(ret$cright)]
    ret[is.na(ret$cleft), "cleft"] <- factor(unclass(object)[is.na(ret$cleft)] - 1,
        levels = 1:length(lev), labels = lev, exclude = 0, ordered = TRUE)
    ret$exact <- NA
    ret[which(ret$cright == lev[nlevels(object)]), "cright"] <- NA
    attr(ret, "prob") <- function(weights) {
        ### FIXME: subsetting doesn't change this fct
        if (length(weights) != length(object))
            weights <- rep_len(1, length(object))
        prt <- cumsum(prop.table(xtabs(weights ~ object)))
        function(y) prt[y]
    }
    ret
}

### <FIXME> is bounds = c(min(object), Inf) the right thing, always?
### </FIXME>
R.integer <- function(object, cleft = NA, cright = NA, bounds = c(min(object), Inf), ...) {

    ret <- .mkR(exact = object, cleft = cleft, cright = cright, ...)
    ret$cright[is.na(ret$cright)] <- ret$exact[is.na(ret$cright)]
    ret$cright[which(ret$cright == bounds[2])] <- NA
    ret$cleft[is.na(ret$cleft)] <- ret$exact[is.na(ret$cleft)] - 1
    ret$cleft[which(ret$cleft < bounds[1])] <- NA
    ret$exact <- NA
    ret$approxy <- object
    attr(ret, "prob") <- function(weights)
        .wecdf(object, weights)
    ret
}

R.interval <- function(object, ...) {
    breaks <- attr(object, "levels")
    cleft <- breaks[-length(breaks)]
    cright <- breaks[-1L]
    R(cleft = cleft[object], cright = cright[object], ...)
}

### handle exact integer / factor as interval censored
R.numeric <- function(object = NA, cleft = NA, cright = NA, 
                      tleft = NA, tright = NA, tol = sqrt(.Machine$double.eps), 
                      as.R.ordered = FALSE, as.R.interval = FALSE, ...) {

    ### treat extremely small intervals as `exact' observations
    d <- cright - cleft
    if (length(object) == 1L && all(is.na(object)))
        object <- rep_len(NA, length(d))
    if (any(!is.na(d) | is.finite(d))) {
        if (any(d < 0, na.rm = TRUE)) stop("cleft > cright")
        i <- (d < tol)
        if (any(i, na.rm = TRUE)) {
            i <- which(i)
            object[i] <- cleft[i]
            cleft[i] <- cright[i] <- NA
        }
    }

    if (as.R.ordered && any(!is.na(cleft) || !is.na(cright)))
        warning("as.R.ordered only implemented for exact observations")

    if (as.R.ordered) {
      ### code response as ordered factor
      ### this defines the nonparametric likelihood
      ### in terms of the observed event times
      utm <- sort(unique(object))
      utm <- utm[-length(utm)]
      ### convert to ordered factor
      ct <- cut(object, breaks = c(-Inf, utm, Inf), ordered_result = TRUE)
      tl <- tr <- NA
      if (!all(is.na(tleft)))
         tl <- cut(tleft, breaks = c(-Inf, utm, Inf), ordered_result = TRUE)
      if (!all(is.na(tright)))
         tr <- cut(tright, breaks = c(-Inf, utm, Inf), ordered_result = TRUE)
      ret <- R(object = ct, tleft = tl, tright = tr)
      attr(ret, "unique_obs") <- utm
      return(ret)
    }

    if (as.R.interval) {
      ### code response as interval-censored defining the nonparametric
      ### likelihood BUT keep data as numeric
      utm <- sort(unique(object))
      utm <- utm[-length(utm)]     
      ct <- cut(object, breaks = c(-Inf, utm, Inf))
      return(R(Surv(time = c(-Inf, utm)[ct], time2 = c(utm, Inf)[ct], type = "interval2")))
    }

    ret <- .mkR(exact = object, cleft = cleft, cright = cright,
                tleft = tleft, tright = tright)
    ### if exact, left and right were NA, recode as -Inf, Inf
    ### such that the log-likelihood contribution is zero
    if (any(rna <- rowSums(is.na(ret[, c("exact", "cleft", "cright")])) == 3L)) {
        ret[rna, "cleft"] <- -Inf
        ret[rna, "cright"] <- Inf
    }
    ### <FIXME>
    ### this fails if is.na(object) and only cleft/cright are given
    # attr(ret, "prob") <- function(weights)
    #     .wecdf(object, weights)
    attr(ret, "prob") <- function(weights)
        .wecdf(ret$approxy, weights)
    ### we want something like survfit(Surv(... type = "interval")
    ### with adjustment to min(obs) = 0
    ### </FIXME>
    ret
}

### for object = NA, ie censored observations only
R.logical <- R.numeric

R.response <- function(object, ...)
    return(object)

R.default <- function(object, ...)
    stop("cannot deal with response class", class(object))

.mkR <- function(...) {

    args <- list(...)
    cls <- unique(sapply(args, function(a) class(a)[1]))
    cls <- cls[cls != "logical"]
    stopifnot(length(cls) <= 1)

    n <- unique(sapply(args, length))
    stopifnot(length(n) <= 2)
    if (length(n) == 2) stopifnot(min(n) == 1)

    if (all(sapply(args, function(x) all(is.na(x))))) {
        args$approxy <- NA
        ret <- do.call("as.data.frame", list(x = args))
        class(ret) <- c("response", class(ret))
        return(ret[, c("tleft", "cleft", "exact", "cright", "tright", "approxy"), drop = FALSE])
    }

    ret <- do.call("as.data.frame", list(x = args))
    if (is.null(ret$exact)) ret$exact <- NA
    if (is.null(ret$cleft) || all(is.na(ret$cleft))) {
        ret$cleft <- NA
        if (is.ordered(ret$exact)) 
            ret$cleft <- factor(ret$cleft, levels = 1:nlevels(ret$exact), 
                                labels = levels(ret$exact), ordered = TRUE)
    }
    if (is.null(ret$cright) || all(is.na(ret$cright))) {
        ret$cright <- NA
        if (is.ordered(ret$exact)) 
            ret$cright <- factor(ret$cright, levels = 1:nlevels(ret$exact), 
                                 labels = levels(ret$exact), ordered = TRUE)
    }
    if (is.null(ret$tleft)) ret$tleft <- NA
    if (is.null(ret$tright)) ret$tright <- NA

    if (all(is.finite(ret$exact))) {
#        ret$rank <- rank(ret$exact, ties.method = "max")
        ret$approxy <- ret$exact
    } else {
        ### some meaningful ordering of observations
        tmpexact <- as.numeric(ret$exact)
        tmpleft <- as.numeric(ret$cleft)
        tmpright <- as.numeric(ret$cright)
        tmpler <- c(tmpleft, tmpexact, tmpright)
        tmpleft[!is.finite(tmpleft)] <- min(tmpler[is.finite(tmpler)])
        tmpright[!is.finite(tmpright)] <- max(tmpler[is.finite(tmpler)])
        tmpexact[is.na(tmpexact)] <-
            (tmpleft + ((tmpright - tmpleft) / 2))[is.na(tmpexact)]
#        ret$rank <- rank(tmpexact, ties.method = "max")   
        ret$approxy <- tmpexact
    }
    class(ret) <- c("response", class(ret))
    ret[, c("tleft", "cleft", "exact", "cright", "tright", "approxy"), drop = FALSE]
}

.exact <- function(object)
    !is.na(object$exact)

.cleft <- function(object)
    is.finite(object$cleft) 

.cright <- function(object)
    is.finite(object$cright)

.cinterval <- function(object)
    !.exact(object)
#    .cleft(object) | .cright(object)

.tleft <- function(object)
    is.finite(object$tleft) 

.tright <- function(object)
   is.finite(object$tright)

.tinterval <- function(object)
    .tleft(object) | .tright(object)

.mm_exact <- function(model, data, response, object) {

    e <- .exact(object)
    if (!any(e)) return(NULL)
    tmp <- data[e,,drop = FALSE]
    tmp[[response]] <- object$exact[e]
    Y <- model.matrix(model, data = tmp)
    deriv <- 1
    names(deriv) <- response
    Yprime <- model.matrix(model, data = tmp, deriv = deriv)

    .matrix <- matrix
    if (inherits(Y, "Matrix")) 
        .matrix <- function(...) Matrix(..., sparse = TRUE)

    trunc <- NULL
    if (any(.tinterval(object) & e)) {
        Ytleft <- .matrix(-Inf, nrow = nrow(Y), ncol = ncol(Y))
        Ytright <- .matrix(Inf, nrow = nrow(Y), ncol = ncol(Y))
        if (any(il <- (.tleft(object) & e))) {
            tmp <- data[il,]
            tmp[[response]] <- object$tleft[il]
            Ytleft[.tleft(object)[e],] <- model.matrix(model, data = tmp)
        }
        if (any(ir <- (.tright(object) & e))) {
            tmp <- data[ir,,drop = FALSE]
            tmp[[response]] <- object$tright[ir]
            Ytright[.tright(object)[e],] <- model.matrix(model, data = tmp)
        }
        trunc <- list(left = Ytleft, right = Ytright)
    }

    list(Y = Y, Yprime = Yprime, trunc = trunc, which = which(e))
}

.mm_interval <- function(model, data, response, object) {

    i <- .cinterval(object)
    if (!any(i)) return(NULL)
    tmpdata <- data[i,,drop = FALSE]
    object <- object[i,, drop = FALSE]

    Yleft <- NULL
    if (any(il <- .cleft(object))) {
        tmp <- tmpdata[il,,drop = FALSE]
        tmp[[response]] <- object$cleft[il]
        Ytmp <- model.matrix(model, data = tmp)
        .matrix <- matrix
        if (inherits(Ytmp, "Matrix")) 
            .matrix <- function(...) Matrix(..., sparse = TRUE)
        Yleft <- .matrix(-Inf, nrow = length(il), ncol = ncol(Ytmp))
        colnames(Yleft) <- colnames(Ytmp)
        rownames(Yleft) <- rownames(tmpdata)
        Yleft[il,] <- Ytmp
        attr(Yleft, "constraint") <- attr(Ytmp, "constraint")
        attr(Yleft, "Assign") <- attr(Ytmp, "Assign")
    }

    Yright <- NULL
    if (any(ir <- .cright(object))) {
        tmp <- tmpdata[ir,, drop = FALSE]
        tmp[[response]] <- object$cright[ir]
        Ytmp <- model.matrix(model, data = tmp)
        .matrix <- matrix
        if (inherits(Ytmp, "Matrix")) 
            .matrix <- function(...) Matrix(..., sparse = TRUE)
        Yright <- .matrix(Inf, nrow = length(ir), ncol = ncol(Ytmp))
        colnames(Yright) <- colnames(Ytmp)
        rownames(Yright) <- rownames(tmpdata)
        Yright[ir,] <- Ytmp
        attr(Yright, "constraint") <- attr(Ytmp, "constraint")
        attr(Yright, "Assign") <- attr(Ytmp, "Assign")
    }

    ### only -Inf, Inf intervals present
    if (is.null(Yleft) & is.null(Yright)) {
        ### evaluate model.matrix and obtain dimensions and
        ### constraints; no data from Ytmp is used!"
        tmp <- tmpdata
        yfake <- mkgrid(model, n = 3)[[response]][2]
        tmp[[response]] <- yfake
        Ytmp <- model.matrix(model, data = tmp)
        .matrix <- matrix
        if (inherits(Ytmp, "Matrix")) 
            .matrix <- function(...) Matrix(..., sparse = TRUE)
        Yleft <- .matrix(-Inf, nrow = nrow(Ytmp), ncol = ncol(Ytmp))
        colnames(Yleft) <- colnames(Ytmp)
        rownames(Yleft) <- rownames(tmpdata)
        attr(Yleft, "constraint") <- attr(Ytmp, "constraint")
        attr(Yleft, "Assign") <- attr(Ytmp, "Assign")
        Yright <- .matrix(Inf, nrow = nrow(Ytmp), ncol = ncol(Ytmp))
        colnames(Yright) <- colnames(Ytmp)
        rownames(Yright) <- rownames(tmpdata)
        attr(Yright, "constraint") <- attr(Ytmp, "constraint")
        attr(Yright, "Assign") <- attr(Ytmp, "Assign")
    }

    if (is.null(Yright)) { 
        Yright <- matrix(Inf, nrow = nrow(Yleft), ncol = ncol(Yleft))
        colnames(Yright) <- colnames(Yleft)
        attr(Yright, "constraint") <- attr(Yleft, "constraint")
        attr(Yright, "Assign") <- attr(Yleft, "Assign")
    }
    if (is.null(Yleft)) {
        Yleft <- matrix(-Inf, nrow = nrow(Yright), ncol = ncol(Yright))
        colnames(Yleft) <- colnames(Yright)
        attr(Yleft, "constraint") <- attr(Yright, "constraint")
        attr(Yleft, "Assign") <- attr(Yright, "Assign")
    }

    trunc <- NULL
    if (any(.tinterval(object))) {
        .matrix <- matrix
        if (inherits(Yleft, "Matrix")) 
            .matrix <- function(...) Matrix(..., sparse = TRUE)
        Ytleft <- .matrix(-Inf, nrow = nrow(Yleft), ncol = ncol(Yleft))
        Ytright <- .matrix(Inf, nrow = nrow(Yleft), ncol = ncol(Yleft))
        colnames(Ytleft) <- colnames(Ytright) <- colnames(Yleft)
        if (any(il <- (.tleft(object)))) {
            tmp <- tmpdata[il,,drop = FALSE]
            tmp[[response]] <- object$tleft[il]
            Ytleft[il,] <- model.matrix(model, data = tmp)
        }
        if (any(ir <- (.tright(object)))) {
            tmp <- tmpdata[ir,,drop = FALSE]
            tmp[[response]] <- object$tright[ir]
            Ytright[ir,] <- model.matrix(model, data = tmp)
        }
        trunc <- list(left = Ytleft, right = Ytright)
    }

    list(Yleft = Yleft, Yright = Yright, trunc = trunc, which = which(i))
}

.wecdf <- function(x, weights) {
    ### FIXME: subsetting doesn't change this fct
    if (length(weights) != length(x)) weights <- rep_len(1, length(x))
    ### from: spatstat::ewcdf
    ox <- order(x)
    x <- x[ox]
    w <- weights[ox]
    vals <- sort(unique(x))
    xmatch <- factor(match(x, vals), levels = seq_along(vals))
    wmatch <- tapply(w, xmatch, sum)
    wmatch[is.na(wmatch)] <- 0
    cumwt <- cumsum(wmatch) / sum(wmatch)
    approxfun(vals, cumwt, method = "constant", yleft = 0, 
              yright = sum(wmatch), f = 0, ties = "ordered")
}

as.Surv <- function(object)
    UseMethod("as.Surv")

as.Surv.response <- function(object) {

    stopifnot(all(!.tright(object)))

    exact <- .exact(object)
    cleft <- .cleft(object)
    cright <- .cright(object)
    tleft <- .tleft(object)
    stopifnot(all(!.tright(object)))

    if (any(tleft)) {
        stopifnot(all(!cright))
        tm <- ifelse(exact, object$exact, object$cleft)
        return(Surv(time = object$tleft, time2 = tm, event = exact, 
                    type = "counting"))
    }

    if (any(cleft & cright) || (any(cleft) && any(cright))) {
        stopifnot(all(!tleft))
        return(Surv(time = ifelse(exact, object$exact, object$cleft),
                    time2 = ifelse(exact, object$exact, object$cright),
                    type = "interval2"))
    }

    if (any(cleft) & all(!cright))
        return(Surv(time = ifelse(exact, object$exact, object$cleft),
                    event = exact, type = "right"))

    return(Surv(time = ifelse(exact, object$exact, object$cright),
                event = exact, type = "left"))
}

as.Surv.numeric <- function(object) Surv(object, rep_len(TRUE, length(object)))

R.list <- function(object, ...) {

    ret <- lapply(object, R)
    do.call(".mkR", do.call("rbind", ret))
}

"[.response" <- function(x, i, j, ..., drop = FALSE) {
    cls <- class(x) 
    class(x) <- "data.frame"
    ret <- x[i,j, drop = FALSE]
    class(ret) <- cls
    ret
}

"[<-.response" <- function(x, i, j, value) {
    cls <- class(x)
    class(x) <- "data.frame"
    x[i,j] <- value
    class(x) <- cls
    x
}

### coerse to double vector (get rid of censoring)
as.double.response <- function(x, ...) {
    ex <- x$exact
    le <- x$cleft
    ri <- x$cright
    rex <- ex
    rle <- le
    rri <- ri
    rex[is.na(ex)] <- 0
    rle[is.na(le) | !is.finite(le)] <- 0
    rri[is.na(ri) | !is.finite(ri)] <- 0
    ### (-Inf, x] -> x and (x, Inf) -> x
    rex + (rle + ifelse(is.finite(ri) & is.finite(le), (rri - rle)/2, rri))
}

pstart <- function(x, weights = rep.int(1, NROW(x)))
    UseMethod("pstart")

pstart.numeric <- function(x, weights = rep.int(1, NROW(x))) {
    if (is.null(weights)) weights <- rep.int(1, NROW(x))
    .wecdf(x, weights = weights)(x)
}

pstart.factor <- function(x, weights = rep.int(1, NROW(x))) {
    if (is.null(weights)) weights <- rep.int(1, NROW(x))
    .wecdf(x, weights = weights)(x)
}

pstart.Surv <- function(x, weights = rep.int(1, NROW(x)))
    ### always
    pstart(R(x), weights = weights)

pstart.response <- function(x, weights = rep.int(1, NROW(x))) {
    if (is.null(weights)) weights <- rep.int(1, NROW(x))
    if (is.ordered(x$approxy) || is.integer(x$approxy))
        return(pstart(x$approxy, weights = weights))
    lwr <- x$cleft
    upr <- x$cright
    ex <- x$exact
    if (all(is.na(lwr)) && all(is.na(upr)))
        return(pstart(ex, weights = weights))
    if (any(nex <- !is.na(ex)))
        lwr[nex] <- upr[nex] <- ex[nex]
    lwr <- unclass(lwr)
    upr <- unclass(upr)
    lwr[is.na(lwr)] <- -Inf
    upr[is.na(upr)] <- Inf
    d <- data.frame(lwr = lwr, upr = upr, grp = gl(1, 1))
    w0 <- weights > 0
    TB <- icenReg::ic_np(cbind(lwr, upr) ~ grp, data = d[w0,,drop = FALSE], 
                         weights = weights[w0])
    sc <- TB$scurves
    xint <- sc[[1]]$Tbull_ints
    Prb <- 1 - sc[[1]]$S_curves$baseline
    Plwr <- stepfun(xint[-1,"lower"], Prb)
    Pupr <- stepfun(xint[-1,"upper"], Prb)
    ret <- rowMeans(cbind(Plwr(lwr), Plwr(ex), Pupr(ex), Pupr(upr)), 
                    na.rm = TRUE)
    return(ret)
}

findsupport <- function(x, weights = rep.int(1, NROW(x)), probs = c(.1, .9))
    UseMethod("findsupport")

findsupport.numeric <- function(x, weights = rep.int(1, NROW(x)), probs = c(.1, .9)) {
   if (max(abs(probs - c(0, 1))) < .Machine$double.eps) {
        return(range(x, na.rm = TRUE))
    }
    if (is.null(weights)) weights <- rep.int(1, NROW(x))    
    x <- rep(x, times = weights)
    quantile(x[is.finite(x)], probs = probs, na.rm = TRUE)
}

findsupport.Surv <- function(x, weights = rep.int(1, NROW(x)), probs = c(.1, .9)) {
    if (is.null(weights)) weights <- rep.int(1, NROW(x))
    if (attr(x, "type") %in% c("left", "right")) {
        if (max(abs(probs - c(0, 1))) < .Machine$double.eps) {
            return(range(x[,"time"], na.rm = TRUE))
        }
        sf <- survfit(x ~ 1, data = data.frame(x = x),
                      subset = weights > 0, weights = weights)
        probs <- max(1 - sf$surv) * probs
        return(quantile(sf, prob = probs)$quantile)
    }
    return(findsupport(R(x), weights = weights, probs = probs))
}

findsupport.response <- function(x, weights = rep.int(1, NROW(x)), probs = c(.1, .9)) {

    if (is.null(weights)) weights <- rep.int(1, NROW(x))
    lwr <- x$cleft
    upr <- x$cright
    ex <- x$exact
    if (max(abs(probs - c(0, 1))) < .Machine$double.eps) {
        x <- c(lwr, ex, ex)
        return(range(x[is.finite(x)], na.rm = TRUE))
    }
    if (all(is.na(lwr)) && all(is.na(upr)))
        return(findsupport(ex, weights = weights, probs = probs))
    if (any(nex <- !is.na(ex)))
        lwr[nex] <- upr[nex] <- ex[nex]
    lwr <- unclass(lwr)
    upr <- unclass(upr)
    lwr[is.na(lwr)] <- -Inf
    upr[is.na(upr)] <- Inf
    d <- data.frame(lwr = lwr, upr = upr, grp = gl(1, 1))
    w0 <- weights > 0
    TB <- icenReg::ic_np(cbind(lwr, upr) ~ grp, data = d[w0,,drop = FALSE], 
                         weights = weights[w0])
    sc <- TB$scurves
    xint <- sc[[1]]$Tbull_ints
    Prb <- 1 - sc[[1]]$S_curves$baseline
    isf <- rowSums(is.finite(xint)) == 2
    xint <- xint[isf,,drop = FALSE]
    Prb <- Prb[isf]
    xint <- xint[!duplicated(Prb),,drop = FALSE]
    Prb <- Prb[!duplicated(Prb)]
    ### this always gives a valid support
    probs <- max(Prb) * probs
    Plwr <- splinefun(Prb, xint[,"lower"])
    Pupr <- splinefun(Prb, xint[,"upper"])
    retlwr <- Plwr(probs) 
    retupr <- Pupr(probs)
    rowMeans(cbind(retlwr, retupr))
}

pretrafo <- function(y, weights)
    UseMethod("pretrafo")

.trafo <- function(support, order, coeff) {

    yv <- numeric_var("y", support = support)
    by <- Bernstein_basis(yv, order = order, ui = "increasing")

    ret <- list(
        trafo = function(y) {
            Ym <- model.matrix(by, data = data.frame(y = y))
            ### return h
            return(ifelse(is.finite(y), c(Ym %*% coeff), y))

            ### return plogis(h) in (0, 1)
            ret <- c(plogis(Ym %*% coeff))
            yinf <- !is.finite(y)
            if (any(yinf))
                ret[yinf] <- ifelse(y[yinf] < 0, 0, 1)
            return(ret)
        },
        dtrafo = function(y) {
            Ym1 <- model.matrix(by, data = data.frame(y = y), deriv = c("y" = 1))
            return(ifelse(is.finite(y), c(Ym1 %*% coeff), y))

            Ym <- model.matrix(by, data = data.frame(y = y))
            ret <- c(dlogis(Ym %*% coeff)) * c(Ym1 %*% coeff)
            yinf <- !is.finite(y)
            if (any(yinf))
                ret[yinf] <- 0
            return(ret)
        },
        ddtrafo = function(y) {
            Ym2 <- model.matrix(by, data = data.frame(y = y), deriv = c("y" = 2))
            return(ifelse(is.finite(y), c(Ym2 %*% coeff), y))

            Ym <- model.matrix(by, data = data.frame(y = y))
            Ym1 <- model.matrix(by, data = data.frame(y = y), deriv = c("y" = 1))

            ret <- c(.Logistic()$dd(Ym %*% coeff)) * c(Ym1 %*% coeff)^2
            ret <- ret + c(dlogis(Ym %*% coeff)) * c(Ym2 %*% coeff)
            yinf <- !is.finite(y)
            if (any(yinf))
                ret[yinf] <- 0
            return(ret)
        }
    )
    attr(ret, "support") <- support
    attr(ret, "add") <- c(0, 0) ### c(-1, 1) * diff(support) / 2
    return(ret)
}

pretrafo.numeric <- function(y, weights = NULL) {

    sup <- range(y, na.rm = TRUE)
    N <- max(c(length(y[!is.na(y)]), ifelse(!is.null(weights), 
                                            sum(weights[!is.na(y)]), 0)))
    order <- 6 + ceiling(sqrt(N))
    grd <- seq(from = sup[1L], to = sup[2L],
               length.out = order + 1)
    Fgrd <- .wecdf(y, weights)(grd)
    hgrd <- pmin(qlogis((N - 1) / N), qlogis(Fgrd))
    ### derivative must be > 0 everywhere
    hgrd <- hgrd + cumsum(c(0, (diff(hgrd) < sqrt(.Machine$double.eps)) * .1))
    return(.trafo(support = sup, order = order, coeff = hgrd))
}

pretrafo.ordered <- function(y, weights = 1)
    pretrafo(unclass(y), weights = weights)

pretrafo.integer <- function(y, weights = 1)
    pretrafo(as.numeric(y), weights = weights)

pretrafo.Surv <- function(y, weights = 1) {

    if (length(weights) == NROW(y)) {
        sf <- survival::survfit(y ~ 1, subset = weights > 0, weights = weights)
        N <- max(NROW(y[!is.na(y[,1]),]), sum(weights[!is.na(y[,1])])) ### case weights?
    } else {
        sf <- survival::survfit(y ~ 1)
        N <- NROW(y[!is.na(y[,1]),])
    }
    tm <- sf$time
    sv <- sf$surv
    tm <- tm[diff(c(1, sv)) < 0]
    sv <- sv[diff(c(1, sv)) < 0]

    sup <- range(tm)
    order <- 6 + ceiling(sqrt(N))
    grd <- seq(from = sup[1L], to = sup[2L],
               length.out = order + 1)
    Fgrd <- stepfun(tm[-1L], 1 - sv)(grd)
    hgrd <- pmax(qlogis(1 / N), pmin(qlogis((N - 1) / N), qlogis(Fgrd)))
    ### derivative must be > 0 everywhere
    hgrd <- hgrd + cumsum(c(0, (diff(hgrd) < sqrt(.Machine$double.eps)) * .1))

    return(.trafo(support = sup, order = order, coeff = hgrd))
}

pretrafo.response <- function(y, weights = 1) {

    if (is.ordered(y$approxy) || is.integer(y$approxy))
        return(pretrafo(y$approxy, weights = weights))

    lwr <- y$cleft
    upr <- y$cright
    ex <- y$exact
    if (all(is.na(lwr)) && all(is.na(upr)))
        return(pretrafo(ex, weights = weights))

    if (any(nex <- !is.na(ex)))
        lwr[nex] <- upr[nex] <- ex[nex]
    lwr <- unclass(lwr)
    upr <- unclass(upr)
    lwr[is.na(lwr)] <- -Inf
    upr[is.na(upr)] <- Inf
    d <- data.frame(lwr = lwr, upr = upr, grp = gl(1, 1))
    d <- d[complete.cases(d),,drop = FALSE]
    if (length(weights) == NROW(y)) {
        w0 <- weights > 0
        TB <- icenReg::ic_np(cbind(lwr, upr) ~ grp, data = d[w0,,drop = FALSE], 
                             weights = weights[w0])
        N <- max(nrow(d), sum(weights))
    } else {
        TB <- icenReg::ic_np(cbind(lwr, upr) ~ grp, data = d)
        N <- nrow(d)
    }

    sc <- TB$scurves
    tm <- sc[[1]]$Tbull_ints[, "upper"]
    sv <- sc[[1]]$S_curves$baseline
    sv <- sv[is.finite(tm)]
    tm <- tm[is.finite(tm)]

    sup <- range(tm)
    order <- 6 + ceiling(sqrt(N))
    grd <- seq(from = sup[1L], to = sup[2L],
               length.out = order + 1)
    Fgrd <- stepfun(tm[-1L], 1 - sv)(grd)
    hgrd <- pmax(qlogis(1 / N), pmin(qlogis((N - 1) / N), qlogis(Fgrd)))
    ### derivative must be > 0 everywhere
    hgrd <- hgrd + cumsum(c(0, (diff(hgrd) < sqrt(.Machine$double.eps)) * .1))

    return(.trafo(support = sup, order = order, coeff = hgrd))
}
