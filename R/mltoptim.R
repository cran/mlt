
mltoptim <- function(auglag = list(maxtry = 5, kkt2.check = FALSE), 
                     spg = list(maxit = 10000, quiet = TRUE, checkGrad = FALSE),
                     nloptr = NULL, ### list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1.0e-8),
                     trace = FALSE) 
{
    ret <- list()
    if (!is.null(auglag))
        ret$auglag <- function(theta, f, g, ui = NULL, ci = NULL) {
            control <- auglag
            maxtry <- control$maxtry
            control$maxtry <- NULL
            control$trace <- trace
            atheta <- theta
            if (!is.null(ui)) {
                for (tr in 1:maxtry) {
                    ret <- try(alabama::auglag(par = atheta, fn = f, gr = g, hin = function(par) ui %*% par - ci, hin.jac = function(par) ui,
                                               control.outer = control)[c("par", "convergence", "value")])
                    atheta <- runif(length(atheta))
                    names(atheta) <- names(theta)
                    if (inherits(ret, "try-error")) next()
                    if (ret$convergence == 0) break()
                }
            } else { 
                control <- spg
                control$trace <- trace
                quiet <- control$quiet
                control$quiet <- NULL
                ret <- try(BBoptim(par = theta, fn = f, gr = g, control = control, quiet = quiet))
            }
            if (inherits(ret, "try-error"))
                ret <- list(par = theta, convergence = 1)
            return(ret)
        }
    if (!is.null(spg)) 
        ret$spg <- function(theta, f, g, ui = NULL, ci = NULL) {
            control <- spg
            control$trace <- trace
            quiet <- control$quiet
            control$quiet <- NULL
            if (!is.null(ui)) {
                ret <- try(BB::BBoptim(par = theta, fn = f, gr = g, project = "projectLinear",
                                       projectArgs = list(A = ui, b = ci, meq = 0), control = control,
                                       quiet = quiet))
            } else { 
                ret <- try(BBoptim(par = theta, fn = f, gr = g, control = control, quiet = quiet))
            }
            if (inherits(ret, "try-error"))
                ret <- list(par = theta, convergence = 1)
            return(ret)
        }
    if (!is.null(nloptr))
        ### Note: This is still experimental (and switched off by default)
        ret$nloptr <- function(theta, f, g, ui = NULL, ci = NULL) {
            control <- nloptr
            if (trace) control$print_level <- 2
            atheta <- theta
            if (!is.null(ui)) {
                mui <- -ui
                ret <- try(nloptr::nloptr(
                    x0 = atheta, 
                    eval_f = f, 
                    eval_grad = g, opts = control,
                    eval_g_ineq = function(par) mui %*% par + ci,
                    eval_jac_g_ineq = function(par) mui))
                ret$convergence <- 0L # (0:1)[ret$status %in% c(1L, 4L) + 0L]
            } else { 
                control$algorithm <- "NLOPT_LD_LBFGS"
                ret <- try(nloptr::nloptr(
                    x0 = atheta, 
                    eval_f = f, 
                    eval_grad = g, opts = control))
                ret$convergence <- 0L # (0:1)[ret$status %in% c(1L, 4L) + 0L]
            }
            if (inherits(ret, "try-error")) {
                ret <- list(par = theta, convergence = 1)
            } else {
                ret$par <- ret$solution
            }
            return(ret)
        }

    return(ret)
}
