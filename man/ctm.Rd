\name{ctm}
\alias{ctm}
\alias{has_scale}
\title{
  Conditional Transformation Models
}
\description{
  Specification of conditional transformation models
}
\usage{
ctm(response, interacting = NULL, shifting = NULL, scaling = NULL, 
    scale_shift = FALSE, data = NULL, 
    todistr = c("Normal", "Logistic", "MinExtrVal", "MaxExtrVal", 
                "Exponential", "Laplace", "Cauchy"), 
    sumconstr = inherits(interacting, c("formula", "formula_basis")), ...)
has_scale(object)
}
\arguments{
  \item{response}{a basis function, ie, an object of class \code{basis}}
  \item{interacting}{a basis function, ie, an object of class \code{basis}}
  \item{shifting}{a basis function, ie, an object of class \code{basis}}
  \item{scaling}{a basis function, ie, an object of class \code{basis}}
  \item{scale_shift}{a logical choosing between two different model types
                     in the presence of a \code{scaling} term}
  \item{data}{either a \code{data.frame} containing the model variables
              or a formal description of these variables in an object of class \code{vars}}
  \item{todistr}{a character vector describing the distribution to be transformed}
  \item{sumconstr}{a logical indicating if sum constraints shall be applied}
  \item{object}{an object of class \code{ctm}}
  \item{\dots}{arguments to \code{as.basis} when \code{shifting} is a formula}
}
\details{

  This function only specifies the model which can then be fitted using
  \code{\link{mlt}}. The shift term is positive by default. All arguments except
  \code{response} can be missing (in this case an unconditional distribution
  is estimated). \bibcitet{mlt::Hothorn:Moest:Buehlmann:2017} explain the model class.

  Possible choices of the distributions the model transforms to (the inverse
  link functions \eqn{F_Z}) include the 
  standard normal (\code{"Normal"}), the standard logistic
  (\code{"Logistic"}), the standard minimum extreme value
  (\code{"MinExtrVal"}, also known as Gompertz distribution), and the
  standard maximum extreme value (\code{"MaxExtrVal"}, also known as Gumbel
  distribution) distributions. The exponential distribution
  (\code{"Exponential"}) can be used to fit Aalen additive hazard models.
  Laplace and Cauchy distributions are also available.

  Shift-scale models \bibcitep{mlt::Siegfried_Kook_Hothorn_2023} of the form
  \deqn{P(Y \le y \mid X = x) = F_Z(\sqrt{\exp(s(x)^\top \gamma)} [(a(y) \otimes b(x))^\top \vartheta] + d(x)^\top \beta)}
  (\code{scale_shift = FALSE}) or
  \deqn{P(Y \le y \mid X = x) = F_Z(\sqrt{\exp(s(x)^\top \gamma)} [(a(y) \otimes b(x))^\top \vartheta + d(x)^\top \beta])}
  (\code{scale_shift = TRUE})
  with bases \eqn{a(y)} (\code{response}), \eqn{b(x)} (\code{interacting}),
  \eqn{d(x)} (\code{shifting}), and \eqn{s(x)} (\code{scaling}) can be
  specified as well. 
}
\value{
  An object of class \code{ctm}; and a logical is returned by
  \code{has_scale} when a scale term is present in \code{object}.
}
\references{\bibshow{*}}
