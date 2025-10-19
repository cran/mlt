
### <FIXME> rename fixed to coef and allow for specification of coefs, 
###         ie fitted models? </FIXME>
.mlt_setup <- function(model, data, y, offset = NULL, 
                       fixed = coef(model), dofit = TRUE) {

    response <- variable.names(model, "response")
    stopifnot(length(response) == 1)
    todistr <- model$todistr

    eY <- .mm_exact(model, data = data, response = response, object = y)
    iY <- .mm_interval(model, data = data, response = response, object = y)
    N <- length(eY$which) * length(iY$which)

    if (is.null(eY)) {
        Y <- Yt <- iY$Yleft
    } else {
        Y <- Yt <- eY$Y
        if (!is.null(iY)) Yt <- rbind(Yt, iY$Yleft)
    }
    Assign <- attr(Y, "Assign")
    SCALE <- "bscaling" %in% Assign[2,]
    if (SCALE) Z <- model.matrix(model$model$bscaling, data = data)
    ### <FIXME> fixed current does not apply to scale terms in Z </FIXME>

    ui <- attr(Y, "constraint")$ui
    ci <- attr(Y, "constraint")$ci

    if (!is.null(fixed)) {
        fixed <- fixed[!is.na(fixed)]
        if (length(fixed) == 0) fixed <- NULL
    }

    ### check for constant columns and fix corresponding parameters
    const <- NULL
    cYt <- which(colSums(abs(Yt)) < .Machine$double.eps)
    if (length(cYt)) {
        nf <- rep.int(0, length(cYt))
        names(nf) <- colnames(Yt)[cYt]
        fixed <- c(fixed, nf)
    }

    if (!is.null(fixed)) {
        if (!all(names(fixed) %in% colnames(Yt)))
            stop("Fixing parameters failed, maybe they apply to scale terms?")
        fix <- colnames(Yt) %in% names(fixed)
        fixed <- fixed[colnames(Yt)[fix]]

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
                ret <- matrix(ret, nrow = nrow(beta), 
                              ncol = length(ret), byrow = TRUE)
                ret[,!fix] <- beta
            } else {
                ret[!fix] <- beta
            }
            ret
        }
    } else {
        .parm <- function(beta) beta
        fix <- rep_len(FALSE, ncol(Y))
    } 

    .sparm <- function(beta, soffset = 0) {
        if (SCALE) {
            parm <- .parm(beta)
            if (is.matrix(parm)) {
#               slp <- base::rowSums(soffset + 
#                                    Z * parm[,Assign[2,] == "bscaling"])
                slp <- .Call("R_offrowSums", soffset, Z, 
                             parm[,Assign[2,] == "bscaling"])
            } else {
                slp <- c(soffset + Z %*% parm[Assign[2,] == "bscaling"])
            }
            sterm <- exp(.5 * slp)
            if (is.matrix(parm)) {
                parm[,Assign[2,] == "bscaling"] <- 0L
                Parm <- parm
            } else {
                parm[Assign[2,] == "bscaling"] <- 0L
                Parm <- matrix(parm, nrow = nrow(Z), ncol = length(parm), 
                               byrow = TRUE)
            }
            if (model$scale_shift) {
                idx <- !Assign[2,] %in% "bscaling"
            } else {
                idx <- !Assign[2,] %in% c("bshifting", "bscaling")
            }
            Parm[, idx] <- Parm[, idx] * sterm
            Parm
        } else {
            .parm(beta)
        }
    }

    .ofuns <- function(weights = NULL, subset = NULL, offset = NULL, 
                       perm = NULL, permutation = NULL, 
                       distr = todistr) ### <- 
                       ### change todistr via update(, distr =)
    {
        if (is.null(offset)) {
            offset <- rep_len(0, nrow(data))
            soffset <- rep_len(0, nrow(data))
        }
        ### weight are only necessary for computing the hessian
        ### as we only compute the weighted sum, not the 
        ### individual (to be weighted) contributions to the Fisher info
        if (is.null(weights))
            weights <- rep_len(1, nrow(data))

        if (SCALE) {
            if (is.matrix(offset)) {
                stopifnot(ncol(offset) == 2 && nrow(offset) == nrow(data))
                soffset <- offset[,2]
                offset <- offset[1]
            } else {
                soffset <- rep_len(0, nrow(data))
            }
        } else {
            soffset <- rep_len(0, nrow(data))
        }
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
                permutation <- sample(1:N)
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

        ### evaluate the transformation function 
        ### and it's derivative wrt its first argument (maybe externally
        ### tram::mmlt)
        trafo <- function(beta) {
            trex <- trexprime <- trleft <- trright <- ret_ll
            beta <- .sparm(beta, soffset)
            if (is.matrix(beta)) {
                beta_ex <- beta[es$full_ex,,drop = FALSE]
                beta_nex <- beta[es$full_nex,,drop = FALSE]
            } else {
                beta_ex <- beta_nex <- beta
            }
            id <- function(x) x
            if (!is.null(es$full_ex)) {
                trex[es$full_ex] <- .dealinf(exY, beta_ex, exoffset, id, 0)
                trexprime[es$full_ex] <- .dealinf(exYprime, beta_ex, 
                                                  exoffset, id, 0)
            }
            if (!is.null(es$full_nex)) {
                trleft[es$full_nex] <- .dealinf(iYleft, beta_nex, ioffset, 
                                                id, -Inf)
                trright[es$full_nex] <- .dealinf(iYright, beta_nex, ioffset, 
                                                 id, Inf)
            }
            ret <- list(trex = trex, trexprime = trexprime,
                        trleft = trleft, trright = trright)
            return(ret)
        }

        ### evaluate the derivative of the transformation function
        ### wrt to its parameters (this is just the design matrix
        ### for non-shift-scale models but more complex in this latter
        ### class)
        trafoprime <- function(beta) {
            ret <- vector(mode = "list", length = 4L)
            names(ret) <- c("exY", "exYprime", "iYleft", "iYright")
            idx <- !Assign[2,] %in% "bscaling"
            if (!is.null(es$full_ex)) {
                ret$exY <- exY[,idx,drop = FALSE]
                ret$exYprime <- exYprime[,idx,drop = FALSE]
            }
            if (!is.null(es$full_nex)) {
                ret$iYleft <- iYleft[,idx,drop = FALSE]
                ret$iYright <- iYright[,idx,drop = FALSE]
            }
            if (!SCALE) return(ret)

            sparm <- .parm(beta)
            slp <- c(soffset + Z %*% sparm[Assign[2,] == "bscaling"])
            sterm <- exp(.5 * slp)
            tr <- trafo(beta)
            names(tr) <- c("exY", "exYprime", "iYleft", "iYright")

            if (model$scale_shift) {
                if (!is.null(es$full_ex)) {
                    ret[c("exY", "exYprime")] <- lapply(c("exY", "exYprime"), 
                        function(j) 
                            cbind(sterm[es$full_ex] * 
                                    ret[[j]][es$full_ex,,drop = FALSE], 
                                  .5 * tr[[j]][es$full_ex] * 
                                    Z[es$full_ex,,drop = FALSE]))
                }
                if (!is.null(es$full_nex)) {
                    ret[c("iYleft", "iYright")] <- lapply(c("iYleft", "iYright"), 
                        function(j) 
                            cbind(sterm[es$full_nex] * 
                                      ret[[j]][es$full_nex,,drop = FALSE], 
                                  .5 * tr[[j]][es$full_nex] * 
                                      Z[es$full_nex,,drop = FALSE]))
                }
            } else {
                bbeta <- beta
                if (is.matrix(bbeta)) {
                    bbeta[,Assign[2,] %in% c("bshifting", "bscaling")] <- 0
                } else {
                    bbeta[Assign[2,] %in% c("bshifting", "bscaling")] <- 0
                }
                tr <- trafo(bbeta)
                names(tr) <- c("exY", "exYprime", "iYleft", "iYright")

                idx <- !Assign[2,idx] %in% "bshifting"
                
                if (!is.null(es$full_ex)) {
                    ret[c("exY", "exYprime")] <- lapply(c("exY", "exYprime"), 
                       function(j) 
                            cbind(sterm[es$full_ex] * 
                                      ret[[j]][es$full_ex,idx,drop = FALSE], 
                                  ret[[j]][es$full_ex,!idx,drop = FALSE], 
                                  sterm[es$full_ex] * .5 * 
                                      tr[[j]][es$full_ex] * 
                                       Z[es$full_ex, , drop = FALSE]))
                }
                if (!is.null(es$full_nex)) {
                    ret[c("iYleft", "iYright")]  <- lapply(c("iYleft", "iYright"), 
                        function(j) 
                            cbind(sterm[es$full_nex] * 
                                      ret[[j]][es$full_nex,idx,drop = FALSE], 
                                  ret[[j]][es$full_nex,!idx,drop = FALSE], 
                                  sterm[es$full_nex] * .5 * 
                                     tr[[j]][es$full_nex] * 
                                     Z[es$full_nex, , drop = FALSE]))
                }
            }
            return(ret)
        }

        ll <- function(beta) {
            ret <- ret_ll 
            beta <- .sparm(beta, soffset)
            if (is.matrix(beta)) {
                beta_ex <- beta[es$full_ex,,drop = FALSE]
                beta_nex <- beta[es$full_nex,,drop = FALSE]
            } else {
                beta_ex <- beta_nex <- beta
            }
            if (!is.null(es$full_ex))
                ret[es$full_ex] <- .mlt_loglik_exact(distr, 
                    exY, exYprime, exoffset, extrunc)(beta_ex)
            if (!is.null(es$full_nex))
                ret[es$full_nex] <- .mlt_loglik_interval(distr, 
                    iYleft, iYright, ioffset, itrunc)(beta_nex)
            return(ret)
        }

        ### score without shift-scale
        sc0 <- function(beta, Xmult = TRUE, ret_all = FALSE) {
            ### Xmult = FALSE means
            ### score wrt to an intercept term. This avoids
            ### multiplication with the whole design matrix.
            ### Don't use on your own. 
            if (Xmult) {
                ret <- ret_scM 
            } else {
                ret <- ret_sc
            }
            beta <- .sparm(beta, soffset)
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
                    exY, exYprime, exoffset, extrunc)(beta_ex, Xmult)
                if (EX_ONLY) {
                    ret <- matrix(scr, nrow = NROW(scr), 
                                  dimnames = dimnames(scr))
                } else {
                    ret[es$full_ex,] <- scr
                }
            }
            if (!is.null(es$full_nex)) {
                scr <- .mlt_score_interval(distr, 
                    iYleft, iYright, ioffset, itrunc)(beta_nex, Xmult)
                if (IN_ONLY) {
                    ret <- matrix(scr, nrow = NROW(scr), 
                                  dimnames = dimnames(scr))
                } else {
                    ret[es$full_nex,] <- scr
                 }
            }
            if (!Xmult) return(ret)
            colnames(ret) <- colnames(Y)
            ### return all scores, no questions asked
            if (ret_all) return(ret)
            ### in case beta contains fix parameters,
            ### return all scores
            if (!is.null(fixed)) {
                if (all(names(fixed) %in% nm))
                    return(ret)
            }
            return(ret[, !fix, drop = FALSE])
        }

        ### with potential shift-scale
        sc <- function(beta, Xmult = TRUE) {
            if (!"bscaling" %in% Assign[2,])
                return(sc0(beta, Xmult = Xmult))
            sc <- sc0(beta, Xmult = TRUE, ret_all = TRUE)
            sparm <- .parm(beta)
            slp <- c(soffset + Z %*% sparm[Assign[2,] == "bscaling"])
            sterm <- exp(.5 * slp)
            if (model$scale_shift) {
                idx <- (!Assign[2,] %in% "bscaling")
            } else {
                idx <- (!Assign[2,] %in% c("bshifting", "bscaling"))
            }
            if (!Xmult) {
                fct <- c(sc0(beta, Xmult = FALSE))
                return(cbind(shifting = if (model$scale_shift) sterm * fct 
                                        else fct,
                             scaling = sterm * c(sc[, idx, drop = FALSE] %*% 
                                                 .parm(beta)[idx]) * .5))
            }
            ret <- cbind(sterm * sc[, idx, drop = FALSE],
                  if (!model$scale_shift) sc[, Assign[2, ] == "bshifting", 
                                             drop = FALSE],
                  sterm * c(sc[, idx, drop = FALSE] %*% 
                      .parm(beta)[idx]) * .5 * Z)
            colnames(ret) <- colnames(sc)
            if (!is.null(fixed)) {
                if (all(names(fixed) %in% names(beta)))
                    return(ret)
            }
            return(ret[, !fix, drop = FALSE])
        }

        he = function(beta) {
            if ("bscaling" %in% Assign[2,])
                stop("Analytical Hessian not available for location-scale models")
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
                    exweights)(.sparm(beta_ex, soffset))
            if (!is.null(es$full_nex))
                ret <- ret + .mlt_hessian_interval(distr, 
                    iYleft, iYright, ioffset, itrunc, 
                    iweights)(.sparm(beta_nex, soffset))
            colnames(ret) <- rownames(ret) <- colnames(Y)
            ### in case beta contains fix parameters,
            ### return all scores
            if (!is.null(fixed)) {
                if (all(names(fixed) %in% nm))
                    return(ret)
            }
            return(ret[!fix, !fix, drop = FALSE])
        }

        return(list(
            offset = offset,
            soffset = soffset,
            ll = ll,
            sc = sc,
            he = he,
            trafo = trafo,
            trafoprime = trafoprime
        ))
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
#        ci <- ci + sqrt(.Machine$double.eps) ### 
         ### we need ui %*% theta > ci, not >= ci
    }

    optimfct <- function(theta, weights, subset = NULL, offset = NULL, 
                         scaleparm = TRUE, optim, ...) {
        of <- .ofuns(weights = weights, subset = subset, 
                     offset = offset, ...)
        inweights <- weights

        ### N contributions to the log-likelihood, UNWEIGHTED
        logliki <- function(beta, weights = NULL)
            of$ll(beta)
        ### NEGATIVE sum of log-likelihood contributions, WEIGHTED
        loglikfct <- function(beta, weights)  
            -sum(weights * of$ll(beta))
        ### N contributions to the score function, UNWEIGHTED
        scorei <- function(beta, weights = NULL, Xmult = TRUE)
            of$sc(beta, Xmult = Xmult)
        ### gradient of negative log-likelihood, WEIGHTED
        scorefct <- function(beta, weights)
            -colSums(weights * scorei(beta), na.rm = TRUE)
        ### hessian of negative log-lik, ALWAYS WEIGHTED
        hessian <- function(beta, weights) 
            .ofuns(weights = weights, offset = offset)$he(beta)

        if (scaleparm) {
            Ytmp <- Y
            Ytmp[!is.finite(Ytmp)] <- NA
            sc <- apply(abs(Ytmp[, !fix, drop = FALSE]), 2, max, na.rm = TRUE)
            lt1 <- sc < 1.1
            gt1 <- sc >= 1.1
            sc[gt1] <- 1 / sc[gt1]
            sc[lt1] <- 1
            f <- function(gamma) {
                ## nloptr sometimes forgets about names(gamma)
                ## but this name matching shouldn't be necessary anyhow
                ## gamma[names(sc)] <- gamma[names(sc)] * sc
                loglikfct(gamma * sc, weights)
            }
            g <- function(gamma) {
                ## gamma[names(sc)] <- gamma[names(sc)] * sc
                ret <- scorefct(gamma * sc, weights) * sc
                ## ret[names(sc)] <- ret[names(sc)] * sc
                ret
            }
            h <- function(gamma) {
                hessian(gamma * sc, weights) * sc^2
            }
            theta <- theta / sc
            if (!is.null(ui))
                ui <- t(t(ui) * sc)
        } else {
            f <- function(gamma) loglikfct(gamma, weights)
            g <- function(gamma) scorefct(gamma, weights)
            h <- function(gamma) hessian(gamma, weights)
        }

        if (dofit) {
            for (i in 1:length(optim)) {
                if (i > 1) {
                    msg <- paste(names(optim)[i - 1], "did not converge, trying", names(optim)[i], sep = " ")
                    warning(msg)
                }	
                ret <- optim[[i]](theta = theta, f = f, g = g, ui = ui, ci = ci, 
                                  h = h)
                if (ret$convergence == 0) break()
            }
        } else {
            ret <- list(par = theta, value = f(theta), convergence = 0)
        }
        if (ret$convergence != 0)
            warning("Optimisation did not converge")
        ### degrees of freedom: number of free parameters, ie
        ###  #parm NOT meeting the constraints
        ret$df <- length(ret$par)
        ### <FIXME> check on alternative degrees of freedom
#        if (!is.null(ui)) 
#            ret$df <- ret$df - sum(ui %*% ret$par - ci < .Machine$double.eps)
        ### </FIXME>
        if (scaleparm) ret$par <- ret$par * sc

        if (SCALE) {
            ret$hessian <- function(beta, weights) {
                H <- ret$optim_hessian
                if (!scaleparm && !is.null(H) && 
                    max(c(abs(beta - ret$par), 
                          abs(weights - inweights))) < .Machine$double.eps) {
                    ret <- H
                } else {
                    # warning("Analytical Hessian not available, using numerical approximation")
                    ret <- try(stats::optimHess(ret$par, fn = loglikfct, gr = scorefct, 
                                                weights = weights))
                    ### this may fail because numDeriv can't deal with -Inf
                    ### values of the target function
                    if (inherits(ret, "try-error"))
                        ret <- try(numDeriv::hessian(loglikfct, beta, weights = weights))
                    if (inherits(ret, "try-error"))
                        stop("Approximation of Hessian failed")
                }
                rownames(ret) <- colnames(ret) <- names(beta)
                return(ret)
            }
        } else {
            ret$hessian <- hessian
        }

        ret$logliki <- logliki

        ### sum of log-likelihood contributions, WEIGHTED
        ret$loglik <- function(beta, weights)  
            sum(weights * ret$logliki(beta))

        ret$scorei <- scorei

        ### N contributions to score function, WEIGHTED
        ret$score <- function(beta, weights, Xmult = TRUE)
            weights * ret$scorei(beta, Xmult = Xmult)

        ### FIXME: remove weight arguments (needed in tram::mmlt)
        ret$trafo <- function(beta, weights) of$trafo(beta)
        ret$trafoprime <- function(beta, weights) of$trafoprime(beta)

        if (scaleparm) ret$parsc <- sc

        return(ret)
    }

    coef <- rep_len(NA, length(fix))
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

    ret$logliki <- function(beta, weights = NULL)
        .ofuns(weights = weights, offset = offset)$ll(beta)

    ### sum of log-likelihood contributions, WEIGHTED
    ret$loglik <- function(beta, weights)  
        sum(weights * .ofuns(weights = weights, offset = offset)$ll(beta))

    ret$scorei <- function(beta, weights = NULL, Xmult = TRUE)
        .ofuns(weights = weights, offset = offset)$sc(beta, Xmult = Xmult)

    ### N contributions to score function, WEIGHTED
    ret$score <- function(beta, weights = NULL, Xmult = TRUE)
        weights * ret$scorei(beta, Xmult = Xmult)

    hessian <- function(beta, weights) 
        .ofuns(weights = weights, offset = offset)$he(beta)

    class(ret) <- c("mlt_setup", "mlt")
    return(ret)
}

.mlt_start <- function(model, data, y, pstart, offset = NULL, fixed = NULL, 
                       weights = 1) {

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
        fix <- rep_len(FALSE, ncol(Y))
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
    ### <TH> this makes sure we obtain the ML solution
    ###      under a completely random outcome, ie the correct
    ###      solution for an unconditional model
    Z <- todistr$q(pmin(0, log(pstart)), log.p = TRUE) - offset

    X <- X * sqrt(weights)
    Z <- Z * sqrt(weights)
    X <- X[is.finite(Z),,drop = FALSE]
    Z <- Z[is.finite(Z)]
    ### </TH>

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
        cipls <- ci + sqrt(.Machine$double.eps) 
            ### we need ui %*% theta > ci, not >= ci
    }

    if (!is.null(ui)) {    
        iter <- 0
        while(iter < 5) {
            iter <- iter + 1
            ### always produces solutions meeting the constraints;
            ### problems with non-pd Dmat
            ret <- try(c(quadprog::solve.QP(Dmat, dvec, 
                                          t(ui), cipls)$solution),  
                       silent = TRUE)
            ### takes non-pd Dmat but does not always meet constraints
            if (inherits(ret, "try-error"))
                ret <- try(c(coneproj::qprog(Dmat, dvec, 
                                             ui, cipls, msg = FALSE)$thetahat), 
                           silent = TRUE)
            ### make sure constraints apply before moving on
            if (!inherits(ret, "try-error")) {
                if (all(ui %*% ret > ci)) break
            }
            diag(Dmat) <- diag(Dmat) + 1e-3
        }
        if (!all(ui %*% ret > ci))
            warning("Starting values violate contraints")
    } else {
        ret <- lm.fit(x = X, y = Z)$coefficients
    }
    names(ret) <- names(coef(model))[!fix]
    ret[!is.finite(ret)] <- 0
    ret
}

.mlt_fit <- function(object, weights, subset = NULL, offset = NULL, 
                     theta = NULL, scaleparm = TRUE, optim, fixed = NULL, ...) {

    if (is.null(theta))
        stop(sQuote("mlt"), "needs suitable starting values")

    ### BBoptim issues a warning in case of unsuccessful convergence
    ret <- try(object$optimfct(theta, weights = weights, 
        subset = subset, offset = offset, scaleparm = scaleparm, 
        optim = optim, ...))    

    cls <- class(object)
    object[names(ret)] <- NULL
    object <- c(object, ret)
    object$coef[] <- object$parm(ret$par) ### [] preserves names
    object$theta <- theta ### starting value
    object$subset <- subset
    object$scaleparm <- scaleparm ### scaling yes/no
    object$weights <- weights
    object$offset <- offset
    object$optim <- optim
    class(object) <- c("mlt_fit", cls)
    
    return(object)
}

mlt <- function(model, data, weights = NULL, offset = NULL, fixed = NULL,
                theta = NULL, pstart = NULL, scaleparm = TRUE,
                dofit = TRUE, optim = mltoptim(hessian = has_scale(model))) {

    vars <- as.vars(model)
    response <- variable.names(model, "response")
    responsevar <- vars[[response]]
    ### <FIXME>: how can we check Surv or R objects (which may contain
    ### values outside bounds)?
    if (!(inherits(data[[response]], "Surv") ||
          inherits(data[[response]], "response")))
        stopifnot(check(responsevar, data))
    ### </FIXME>
    bounds <- bounds(responsevar)
    stopifnot(length(response) == 1)
    ### <FIXME> add bounds?
    y <- R(object = data[[response]])
    ### </FIXME>

    if (!.checkR(y, weights)) dofit <- FALSE

    if (is.null(weights)) weights <- rep_len(1, nrow(data))
    if (is.null(offset)) offset <- rep_len(0, nrow(data))
    stopifnot(nrow(data) == length(weights))
    stopifnot(nrow(data) == length(offset))

    s <- .mlt_setup(model = model, data = data, y = y, 
                    offset = offset, fixed = fixed, dofit = dofit) 
    s$convergence <- 1
    if (!dofit && is.null(theta)) return(s)

    if (is.null(theta)) {
        ### unconditional ECDF, essentially
        ### this doesn't really work for censored data, any alternative?
        if (is.null(pstart)) 
            pstart <- pstart(data[[response]], weights = weights)
        theta <- .mlt_start(model = model, data = data, y = y, 
                            pstart = pstart, offset = offset, 
                            fixed = s$fixed, weights = weights)
    }

    args <- list()
    args$object <- s
    args$weights <- weights
    args$offset <- offset
    args$theta <- theta
    args$subset <- NULL ### only available in update()
    args$scaleparm <- scaleparm
    args$optim <- optim
    args$fixed <- fixed
    ret <- do.call(".mlt_fit", args)
    ret$call <- match.call()
    ret$bounds <- bounds
    ret$response <- y
    ret
}

update.mlt_fit <- function(object, weights = stats::weights(object), 
                           subset = NULL, offset = object$offset,
                           theta = coef(object, fixed = FALSE), 
                           fixed = NULL,
                           ...) {

    stopifnot(is.null(subset) || is.null(fixed))
    if (!is.null(fixed)) {
        ### refit completely, because fixed needs to be handled
        ### by .mlt_setup
        strt <- theta
        strt <- strt[!names(strt) %in% names(fixed)]
        return(mlt(object$model, data = object$data, weights = weights,
            offset = offset, theta = strt, fixed = fixed,
            scaleparm = object$scaleparm, optim = object$optim))
    }

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
    args$scaleparm <- object$scaleparm
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
               offset = object$offset, 
               theta = coef(object), fixed = object$fixed, 
               scaleparm = object$scaleparm, optim = object$optim)
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
    stopifnot(length(distr()$parm()) == 0L)
    fr <- .CureRate(distr = distr)
    ### <FIXME> allow for additional parameters in frailty
    ll <- function(logitrho)
        -logLik(update(object, distr = .CureRate(logitrho, distr = distr)))
    om <- optimise(ll, interval = interval, maximum = FALSE)
    model$todistr <- .CureRate(om$minimum, distr = distr)
    ### </FIXME>
    ret <- mlt(model = model, data = object$data, weights = weights(object),
               offset = object$offset, 
               theta = coef(object), fixed = object$fixed, 
               scaleparm = object$scaleparm, optim = object$optim)
    class(ret) <- c("fmlt", class(ret))
    ret
}
