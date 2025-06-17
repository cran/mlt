\name{mltoptim}
\alias{mltoptim}
\title{
  Control Optimisation
}
\description{
  Define optimisers and their control parameters
}
\usage{
mltoptim(auglag = list(maxtry = 5, kkt2.check = hessian), 
         spg = list(maxit = 10000, quiet = TRUE, checkGrad = FALSE), 
         nloptr = list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1.0e-8, maxeval = 1000L), 
         constrOptim = list(method = "BFGS", control = list(), mu = 1e-04, 
                            outer.iterations = 100, outer.eps = 1e-05, hessian = hessian),
         trace = FALSE, hessian = FALSE)
}
\arguments{
  \item{auglag}{
  A list with control parameters for the \code{\link[alabama]{auglag}} optimiser.
  \code{maxtry} is the number of times the algorithm is started on random starting
  values in case it failed with the precomputed ones.
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
  \code{g} (scores), \code{ui} and \code{ci} (linear inequality constraints).
  Adding further such functions is a way to add more optimisers to \code{\link{mlt}}.
  The first one in this list converging  defines the resulting model.
}
\examples{

  ### set-up linear transformation model for conditional
  ### distribution of dist given speed
  dist <- numeric_var("dist", support = c(2.0, 100), bounds = c(0, Inf))
  ctmm <- ctm(response = Bernstein_basis(dist, order = 4, ui = "increasing"),
              shifting = ~ speed, data = cars)

  ### use auglag with kkt2.check = TRUE => the numerically determined
  ### hessian is returned as "optim_hessian" slot
  op <- mltoptim(auglag = list(maxtry = 5, kkt2.check = TRUE))[1]
  mltm <- mlt(ctmm, data = cars, scale = FALSE, optim = op)

  ### compare analytical and numerical hessian
  all.equal(c(Hessian(mltm)), c(mltm$optim_hessian), tol = 1e-4)

}
\keyword{list}
