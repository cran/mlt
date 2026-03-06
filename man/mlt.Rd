\name{mlt}
\alias{mlt}
\title{
  Most Likely Transformations
}
\description{
  Likelihood-based model estimation in conditional transformation models
}
\usage{
mlt(model, data, weights = NULL, offset = NULL, fixed = NULL,
    theta = NULL, pstart = NULL, scaleparm = TRUE,
    dofit = TRUE, optim = mltoptim(hessian = has_scale(model)))
}
\arguments{
  \item{model}{a conditional transformation model as specified by \code{\link{ctm}}}
  \item{data}{a \code{data.frame} containing all variables specified in \code{model}}
  \item{weights}{an optional vector of case weights}
  \item{offset}{an optional vector of offset values; offsets are not added
                to an optional \code{scaling} term (see \code{link{ctm}})}
  \item{fixed}{a named vector of fixed regression coefficients; the names
               need to correspond to column names of the design matrix}
  \item{theta}{optional starting values for the model parameters}
  \item{pstart}{optional starting values for the distribution function
                evaluated at the data}
  \item{scaleparm}{a logical indicating if (internal) scaling shall be applied to the
               model parameters; \code{TRUE} unless a location-scale model
               is fitted where the numerically approximated Hessian is only
               available on the original parameter scale}
  \item{dofit}{a logical indicating if the model shall be fitted to the
               data (\code{TRUE}) or not. If \code{theta} is given,
               a model of class \code{mlt} (a full "fitted" model)
               featuring these parameters is returned. Otherwise,
               an unfitted model of class \code{ctm} is returned}
  \item{optim}{a list of functions implementing suitable optimisers,
               requires numerical Hessian for location-scale models}
}
\details{
  This function fits a conditional transformation model by searching for 
  the most likely transformation as described in
  \bibcitet{mlt::Hothorn:Moest:Buehlmann:2017}, including location-scale models
  \bibcitep{Siegfried_Kook_Hothorn_2023}. Implementation details are given
  in \bibcitep{mlt::Hothorn:2018}.
}
\value{
  An object of class \code{mlt} with corresponding methods.
}
\references{\bibshow{*}}
\examples{
 
  ### set-up conditional transformation model for conditional
  ### distribution of dist given speed
  dist <- numeric_var("dist", support = c(2.0, 100), bounds = c(0, Inf))
  speed <- numeric_var("speed", support = c(5.0, 23), bounds = c(0, Inf)) 
  ctmm <- ctm(response = Bernstein_basis(dist, order = 4, ui = "increasing"),
              interacting = Bernstein_basis(speed, order = 3))

  ### fit model
  mltm <- mlt(ctmm, data = cars)

  ### plot data
  plot(cars)
  ### predict quantiles and overlay data with model via a "quantile sheet"
  q <- predict(mltm, newdata = data.frame(speed = 0:24), type = "quantile", 
               p = 2:8 / 10, K = 500)
  tmp <- apply(q, 1, function(x) lines(0:24, x, type = "l"))

}
