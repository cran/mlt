
mltoptim <- function(
    auglag = list(
        kkt2.check = hessian,         ### turn off/on numerical hessian
               eps = abstol,          ### absolute tolerance for _parameter_ updates 
             itmax = 1000L,           ### max number of outer iterations
            method = "BFGS",          ### inner algorithm
             maxit = 500L             ### max number of inner (BFGS) iterations    
        ), 
    spg = list(
              ftol = abstol,          ### absolute tolerance for _neg. logLik_
             quiet = TRUE,            ### don't talk
         checkGrad = FALSE            ### don't check analytical gradient
    ),
    nloptr = list(
         algorithm = "NLOPT_LD_MMA",  ### inner algorithm 
          ftol_rel = reltol,          ### relative change for _neg. logLik_
          ftol_abs = abstol,          ### absolute tolerance for _neg. logLik_
           maxeval = 1000L            ### max number of evaluations
    ),
    constrOptim = list(
            method = "BFGS",          ### inner algorithm 	
             maxit = 1000L,           ### max number of inner (BFGS) iterations 
  outer.iterations = 500L,            ### max number of outer iterations
         outer.eps = reltol          ### relative change for _neg. logLik_          
    ),        
    optim = list(
  checkconstraints = TRUE,            ### return -Inf if violated
            method = "BFGS",          ### inner algorithm 	
             maxit = 1000L,           ### max number of inner (BFGS) iterations 
            reltol = reltol          ### relative change for _neg. logLik_          
    ),        
    nlminb = list(
  checkconstraints = TRUE,            ### return -Inf if violated
          iter.max = 1000L,           ### max number of iterations 
          eval.max = 1500L,           ### max number of function evaluations
           rel.tol = reltol,          ### relative change for _neg. logLik_
           abs.tol = 0.0,             ### absolute tolerance (nll is not >= 0) 
            xf.tol = 1e-10
    ),  
  abstol = 1e-07, 
  reltol = 1e-6, 
   trace = FALSE, 
 hessian = FALSE) 
{
    ret <- list()

    if (!is.null(auglag))
        ret$auglag <- function(theta, f, g, ui = NULL, ci = NULL, ...) {
            control <- auglag
            maxtry <- control$maxtry
            control$maxtry <- NULL
            control$trace <- trace
            atheta <- theta
            if (!is.null(ui)) {
                ret <- try(alabama::auglag(par = atheta, fn = f, gr = g, 
                    hin = function(par) ui %*% par - ci, hin.jac = function(par) ui,
                    control.outer = control))
                if (inherits(ret, "try-error"))
                    return(list(par = theta, convergence = 1))
                rtn <- c("par", "convergence", "value")
                if ("hessian" %in% names(ret)) {
                    ret$optim_hessian <- ret$hessian
                    rtn <- c(rtn, "optim_hessian")
                }
                ret <- ret[rtn]
            } else { 
                control <- spg
                control$trace <- trace
                quiet <- control$quiet
                control$quiet <- NULL
                ret <- try(BB::BBoptim(par = theta, fn = f, gr = g, control = control, quiet = quiet))
                if (inherits(ret, "try-error"))
                    return(list(par = theta, convergence = 1))
                ### we often only use this part for generating starting
                ### values, so refrain from issuing a warning
                # if (hessian)
                #    warning("Cannot compute Hessian using BB::BBoptim")
            }
            return(ret)
        }
    if (!is.null(spg))
        ret$spg <- function(theta, f, g, ui = NULL, ci = NULL, ...) {
            control <- spg
            control$trace <- trace
            quiet <- control$quiet
            control$quiet <- NULL
            if (!is.null(ui)) {
                ret <- try(BB::BBoptim(par = theta, fn = f, gr = g, project = "projectLinear",
                                       projectArgs = list(A = ui, b = ci, meq = 0), control = control,
                                       quiet = quiet))
            } else { 
                ret <- try(BB::BBoptim(par = theta, fn = f, gr = g, control = control, quiet = quiet))
            }
            if (inherits(ret, "try-error"))
                ret <- list(par = theta, convergence = 1)
            if (hessian)
                warning("Cannot compute Hessian using BB::BBoptim")
            return(ret)
        }
    if (!is.null(nloptr))
        ret$nloptr <- function(theta, f, g, ui = NULL, ci = NULL, ...) {
            control <- nloptr
            if (trace) control$print_level <- 2
            atheta <- theta
            if (!is.null(ui)) {
                mui <- -ui
                ### <FIXME> replace with nloptr::auglag </FIXME>
                ret <- try(nloptr::nloptr(
                    x0 = atheta, 
                    eval_f = f, 
                    eval_grad_f = g, opts = control,
                    eval_g_ineq = function(par) mui %*% par + ci,
                    eval_jac_g_ineq = function(par) mui))
                if (ret$status < 0) {
                    ### some form of failure
                    ret$convergence <- 1
                    warning(ret$message)
                } else {
                    ret$convergence <- 0L
                    ### some tolerance criterion reached
                    # if (ret$status > 1)
                    #   warning(ret$message)
                }
            } else { 
                control$algorithm <- "NLOPT_LD_LBFGS"
                ret <- try(nloptr::nloptr(
                    x0 = atheta, 
                    eval_f = f, 
                    eval_grad_f = g, opts = control))
                if (ret$status < 0) {
                    ### some form of failure
                    ret$convergence <- 1
                    warning(ret$message)
                } else {
                    ret$convergence <- 0L
                    ### some tolerance criterion reached
                    # if (ret$status > 1)
                    #    warning(ret$message)
                }
            }
            if (inherits(ret, "try-error")) {
                ret <- list(par = theta, convergence = 1)
            } else {
                ret$par <- ret$solution
                ret$value <- ret$objective
            }
            if (hessian)
                warning("Cannot compute Hessian using nloptr::nloptr")
            return(ret)
        }
    if (!is.null(constrOptim))
        ret$constrOptim <- function(theta, f, g, ui = NULL, ci = NULL, ...) {
            control <- constrOptim$control
            control$trace <- trace
            control[c("mu", "method", "outer.iterations", "outer.eps")] <- NULL
            if (!is.null(ui)) {
                ret <- try(stats::constrOptim(theta = theta, f = f, grad = g,
                                              ui = ui, ci = ci, mu = constrOptim$mu, 
                                              control = control, method = constrOptim$method,
                                              outer.iterations = constrOptim$outer.iterations, 
                                              outer.eps = constrOptim$outer.eps, 
                                              hessian = hessian))
                rtn <- c("par", "convergence", "value")
                if ("hessian" %in% names(ret)) {
                    ret$optim_hessian <- ret$hessian
                    rtn <- c(rtn, "optim_hessian")
                }
            } else { 
                ret <- try(stats::optim(par = theta, fn = f, gr = g,
                                        control = control, 
                                        method = constrOptim$method,
                                        hessian = hessian))
                rtn <- c("par", "convergence", "value")
                if ("hessian" %in% names(ret)) {
                    ret$optim_hessian <- ret$hessian
                    rtn <- c(rtn, "optim_hessian")
                }
            }
            if (inherits(ret, "try-error"))
                ret <- list(par = theta, convergence = 1)
            return(ret)

        }
    if (!is.null(optim)) {
        ret$optim <- function(theta, f, g, h, ui = NULL, ci = NULL) {
            control <- optim
            control$trace <- trace
            objfun <- f
            if (!is.null(ui) && control$checkconstraints) {
                objfun <- function(theta) {
                    if (any(ui %*% theta - ci < -.Machine$double.eps^(1/3)))
                        return(-Inf)
                    return(f(theta))
                }
            }
            control$checkconstraints <- NULL
            method <- control$method
            control$method <- NULL
            ### _unconstraint optimisation_
            ret <- try(stats::optim(par = theta, fn = objfun, gr = g, method = method,
                                    control = control, hessian = hessian))
            if (inherits(ret, "try-error")) {
                return(list(par = theta, convergence = 1))
            } else {
                ### check constraints
                if (!is.null(ui)) {
                    if (any(ui %*% ret$par - ci < -.Machine$double.eps^(1/3)))
                        warning("optim result violates inequality constraints")
                }
            }
            rtn <- c("par", "convergence", "value")
            if ("hessian" %in% names(ret)) {
                ret$optim_hessian <- ret$hessian
                rtn <- c(rtn, "optim_hessian")
            }
            return(ret[rtn])
        }
    }
    if (!is.null(nlminb)) {
        ret$nlminb <- function(theta, f, g, h, ui = NULL, ci = NULL) {
            control <- nlminb
            control$trace <- trace
            objfun <- f
            if (!is.null(ui) && control$checkconstraints) {
                objfun <- function(theta) {
                    if (any(ui %*% theta - ci < -.Machine$double.eps^(1/3)))
                        return(Inf) ### we minimise
                    return(f(theta))
                }
            }
            control$checkconstraints <- NULL
            ### _unconstraint optimisation_
            ret <- try(stats::nlminb(start = theta, objective = objfun, gradient = g, hessian = h,
                                     control = control))
            if (inherits(ret, "try-error")) {
                return(list(par = theta, convergence = 1))
            } else {
                ### check constraints
                if (!is.null(ui)) {
                    if (any(ui %*% ret$par - ci < -.Machine$double.eps^(1/3)))
                        warning("nlminb result violates inequality constraints")
                }
            }
            ret$value <- ret$objective
            rtn <- c("par", "convergence", "value", "optim_hessian")
            return(ret[rtn])
        }
    }
    if (hessian)
        return(ret[c("auglag", "constrOptim", "optim")])
    return(ret[c("auglag", "spg", "nloptr", "constrOptim", "nlminb", "optim")])
}
