\name{mmlt}
\alias{mmlt}
\alias{coef.cmmlt}
\alias{coef.mmmlt}
\alias{predict.mmlt}
\alias{simulate.mmlt}
\title{
  Multivariate Conditional Transformation Models
}
\description{
  Conditional transformation models for multivariate continuous, discrete,
  or a mix of continuous and discrete outcomes
}
\usage{
mmlt(..., formula = ~ 1, data, conditional = FALSE, 
     theta = NULL, fixed = NULL, scaleparm = FALSE,
     optim = mltoptim(hessian = TRUE),  ### provides hessian
     args = list(seed = 1, M = 1000), 
     dofit = TRUE, domargins = TRUE)
\method{coef}{cmmlt}(object, newdata, 
     type = c("all", "conditional", "Lambdapar", 
              "Lambda", "Omega", "Lambdainv", "Omegainv",
              "Precision", "PartialCorr", "Sigma", "Corr", 
              "Spearman", "Kendall"), fixed = TRUE, 
     ...)
\method{coef}{mmmlt}(object, newdata, 
     type = c("all", "marginal", "Lambdapar", 
              "Lambda", "Omega", "Lambdainv", "Omegainv",
              "Precision", "PartialCorr", "Sigma", "Corr", 
              "Spearman", "Kendall"), fixed = TRUE,
     ...)
\method{predict}{mmlt}(object, newdata, margins = 1:J, 
        type = c("trafo", "distribution", "survivor", "density", "hazard"), 
                 log = FALSE, args = object$args, ...)
\method{simulate}{mmlt}(object, nsim = 1L, seed = NULL, newdata, K = 50, ...)
}
\arguments{
  \item{\dots}{marginal transformation models, one for each response, for
\code{mmlt}. Additional arguments for the methods.
}
  \item{formula}{a model formula describing a model for the dependency
structure via the lambda parameters. The default is set to \code{~ 1} for constant lambdas.
}
  \item{data}{a data.frame.
}
  \item{conditional}{logical; parameters are defined conditionally (only
possible when all models are probit models). This is the default as
described by Klein et al. (2022). If \code{FALSE}, parameters can be
directly interpreted marginally, this is explained in Section 2.6 by Klein
et al. (2022). Using \code{conditional = FALSE} with probit-only models
gives the same likelihood but different parameter estimates.
}
  \item{theta}{an optional vector of starting values.
}
  \item{fixed}{an optional named numeric vector of predefined parameter values
or a logical (for \code{coef}) indicating to also return fixed parameters
(only when \code{type = "all"}).
}
  \item{scaleparm}{a logical indicating if (internal) scaling shall be applied
          to the model parameters.
}
  \item{optim}{a list of optimisers as returned by \code{\link{mltoptim}}
}
  \item{args}{a list of arguments for \code{\link[mvtnorm]{lpmvnorm}}.}
  \item{dofit}{logical; parameters are fitted by default, otherwise a list
               with log-likelihood and score function is returned.
}
  \item{domargins}{logical; all model parameters are fitted by default, 
                   including the parameters of marginal models.
}
  \item{object}{an object of class \code{mmlt}.
}
  \item{newdata}{an optional data.frame coefficients and predictions shall
                 be computed for.
}
  \item{type}{type of coefficient or prediction to be returned.
}
  \item{margins}{indices defining marginal models to be evaluated. Can be
single integers giving the marginal distribution of the corresponding
variable, or multiple integers (currently only \code{1:j} implemented).
}
  \item{log}{logical; return log-probabilities or log-densities if
\code{TRUE}.
}
  \item{nsim}{number of samples to generate.}
  \item{seed}{optional seed for the random number generator.}
  \item{K}{number of grid points to generate.}
}
\details{
  The function implements core functionality for fitting 
  multivariate conditional transformation models
  as described by \bibcitet{mlt::Klein:Hothorn:Barbanti:2020}.
}
\value{
An object of class \code{mmlt} with \code{coef} and \code{predict}
methods.

The model is always parameterised in terms of a unit lower-triangular matrix
\eqn{\Lambda}. If \code{conditional = FALSE} (the default), it is standardised to a
lower-triangular matrix \eqn{\Omega = \Lambda \text{diag}(\Lambda^{-1}
\Lambda^{-\top})} such that \eqn{\Sigma = \Omega^{-1} \Omega^{-\top}} is a
correlation matrix. Note that \eqn{Omega = \Lambda} when \code{conditional =
TRUE} was requested.

}
\references{\bibshow{*}}
\keyword{models}
