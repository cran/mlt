\name{mlt}
\alias{mlt}
\title{
  Most Likely Transformations
}
\description{
  Likelihood-based model estimation in conditional transformation models
}
\usage{
mlt(model, data, weights = NULL, offset = NULL, fixed = NULL, theta = NULL, 
    pstart = NULL, scale = FALSE, checkGrad = FALSE, 
    trace = FALSE, quiet = TRUE, dofit = TRUE, ...)
}
\arguments{
  \item{model}{a conditional transformation model as specified by \code{\link{ctm}}}
  \item{data}{a \code{data.frame} containing all variables specified in \code{model}}
  \item{weights}{an optional vector of weights}
  \item{offset}{an optional vector of offset values}
  \item{fixed}{a named vector of fixed regression coefficients; the names
               need to correspond to column names of the design matrix}
  \item{theta}{optional starting values for the model parameters}
  \item{pstart}{optional starting values for the distribution function
                evaluated at the data}
  \item{scale}{a logical indicating if (internal) scaling shall be applied to the
               model coefficients}
  \item{checkGrad}{a logical triggering checks of the analytic gradient
                   during optimisation}
  \item{trace}{a logical triggering status messages of the optimiser}
  \item{quiet}{a logical turning off any messages from the optimiser}
  \item{dofit}{a logical indicating if the model shall be fitted to the
               data (\code{TRUE}) or not}
  \item{\dots}{additional arguments to \code{\link[BB]{spg}}}
}
\details{
  This function fits a conditional transformation model by searching for 
  the most likely transformation as described in Hothorn et al. (2016).
}
\value{
  An object of class \code{mlt} with corresponding methods.
}
\references{

  Torsten Hothorn, Lisa Moest, Peter Buehlmann (2016), Most Likely 
  Transformations, http://arxiv.org/abs/1508.06749

}
\examples{
 
  ### set-up conditional transformation model for conditional
  ### distribution of dist given speed
  dist <- numeric_var("dist", support = c(2.0, 100), bounds = c(0, Inf))
  speed <- numeric_var("speed", support = c(5.0, 23), bounds = c(0, Inf)) 
  ctmm <- ctm(response = Bernstein_basis(dist, order = 4, ui = "increasing"),
              interacting = Bernstein_basis(speed, order = 3))

  ### fit model
  (mltm <- mlt(ctmm, data = cars))

  ### plot data
  plot(cars)
  ### predict quantiles and overlay data with model via a "quantile sheet"
  q <- predict(mltm, newdata = data.frame(speed = 0:24), type = "quantile", 
               p = 2:8 / 10, K = 500)
  tmp <- apply(q, 1, function(x) lines(0:24, x, type = "l"))

}