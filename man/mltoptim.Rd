\name{mltoptim}
\alias{mltoptim}
\title{
  Control Optimisation
}
\description{
  Define optimisers and their control parameters
}
\usage{
mltoptim(
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
}
\arguments{
  \item{auglag}{
  A list with control parameters for the \code{\link[alabama]{auglag}} optimiser.
  \code{maxtry} is the number of times the algorithm is re-started 
  in case it failed.
}
  \item{spg}{
  A list with control parameters for the \code{\link[BB]{BBoptim}} optimiser (calling
   \code{\link[BB]{spg}} internally).
}
  \item{nloptr}{
  A list with control parameters for the \code{\link[nloptr]{nloptr}}
  family of optimisers. 
}
  \item{constrOptim}{
  A list with control parameters for the \code{\link[stats]{constrOptim}}
  optimiser. 
}
  \item{optim}{
  A list with control parameters for the \code{\link[stats]{optim}}
  optimiser producing an unconstrained fit.
}
  \item{nlminb}{
  A list with control parameters for the \code{\link[stats]{nlminb}}
  optimiser producing an unconstrained fit.
}
  \item{abstol,reltol}{
  Absolute and relative tolerances used as stopping criterion by various
  algorithms.
}
  

  \item{trace}{
  A logical switching trace reports by the optimisers off.
}
  \item{hessian}{A logical indicating if a numerically differentiated Hessian matrix
          be returned.
}
}
\details{
  This function sets-up functions to be called in \code{\link{mlt}} internally.
}
\value{
  A list of functions with arguments \code{theta} (starting values), \code{f} (log-likelihood),
  \code{g} (scores), \code{h} (Hessian), \code{ui} and \code{ci} (linear inequality constraints).
  Adding further such functions is a way to add more optimisers to \code{\link{mlt}}.
  The first one in this list converging  defines the resulting model.

  All procedures except \code{optim} and \code{nlminb} perform constained optimisation. The 
  model is only defined for parameters meeting the constraints. For
  parameter configurations not meeting the constraints, the resulting
  log-likelihood is \code{-Inf}. This, however, does not mean that
  unconstraint optimisation will always produce parameter estimates 
  which lead to valid models, so one should always check the unconstrained
  (but probably faster obtained) result against the (slower) constrained
  solution.

}
\examples{

  ### set-up linear transformation model for conditional
  ### distribution of dist given speed
  dist <- numeric_var("dist", support = c(2.0, 100), bounds = c(0, Inf))
  ctmm <- ctm(response = Bernstein_basis(dist, order = 4, ui = "increasing"),
              shifting = ~ speed, data = cars)

  ### the numerically determined
  ### hessian is returned as "optim_hessian" slot
  op <- mltoptim(hessian = TRUE)
  mltm <- mlt(ctmm, data = cars, scale = FALSE, optim = op)

  ### compare analytical and numerical hessian
  all.equal(c(Hessian(mltm)), c(mltm$optim_hessian), tol = 1e-4)

}
\keyword{list}
