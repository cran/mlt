\name{mmlt-methods}
\alias{coef<-.mmlt}
\alias{weights.mmlt}
\alias{logLik.mmlt}
\alias{vcov.mmlt}
\alias{Hessian.mmlt}
\alias{Gradient.mmlt}
\alias{estfun.mmlt}
\alias{mkgrid.mmlt}
\alias{variable.names.mmlt}
\title{
  Methods for mmlt Objects
}
\description{
  Methods for objects of class mmlt
}
\usage{
\method{weights}{mmlt}(object, ...)
\method{logLik}{mmlt}(object, parm = coef(object, fixed = FALSE), w = NULL, newdata = NULL, ...)
\method{vcov}{mmlt}(object, parm = coef(object, fixed = FALSE), complete = FALSE, ...)
\method{Hessian}{mmlt}(object, parm = coef(object, fixed = FALSE), ...)
\method{Gradient}{mmlt}(object, parm = coef(object, fixed = FALSE), ...)
\method{estfun}{mmlt}(x, parm = coef(x, fixed = FALSE),
       w = NULL, newdata = NULL, ...)
\method{mkgrid}{mmlt}(object, ...)
\method{variable.names}{mmlt}(object, response_only = FALSE, ...)
}
\arguments{
  \item{object, x}{a fitted multivariate transformation model as returned by \code{\link{mmlt}}}
  \item{fixed}{a logical indicating if only estimated coefficients (\code{fixed = FALSE}) 
               should be returned OR (for \code{update})
               a named vector of fixed regression coefficients; the names
               need to correspond to column names of the design matrix}
  \item{parm}{model parameters}
  \item{w}{model weights}
  \item{weights}{model weights}
  \item{newdata}{an optional data frame of new observations. Allows
                 evaluation of the log-likelihood for a given
                 model \code{object} on these new observations. The
                 parameters \code{parm} and \code{w} are ignored in this situation.}
  \item{response_only}{only return the names of the response variables}
  \item{complete}{currently ignored}
  \item{\dots}{additional arguments}
}
\details{
  \code{coef} can be used to get and set model parameters, \code{weights} and 
  \code{logLik} extract weights and evaluate the log-likelihood (also for
  parameters other than the maximum likelihood estimate). \code{Hessian}
  returns the Hessian (of the \emph{negative} log-likelihood) and \code{vcov} the inverse thereof. \code{Gradient}
  gives the negative gradient (minus sum of the score contributions) 
  and \code{estfun} the \emph{negative} score contribution by each observation. \code{mkgrid}
  generates a grid of all variables (as returned by \code{variable.names}) in the model.
}
