\name{plot-predict-simulate}
\alias{predict.mlt}
\alias{simulate.mlt}
\alias{plot.mlt}
\title{
  Plots, Predictions and Samples from mlt Objects
}
\description{
  Plot, predict and sample from objects of class mlt
}
\usage{
\method{plot}{mlt}(x, newdata, type = c("distribution", "survivor", "density", 
     "logdensity", "hazard", "loghazard", "cumhazard", "quantile", "trafo"),
     q = NULL, p = 1:(K - 1) / K, K = 50, col = rgb(.1, .1, .1, .1), add = FALSE, 
     ...)
\method{predict}{mlt}(object, newdata = object$data, type = c("trafo", "distribution", 
         "survivor", "density", "logdensity", "hazard", "loghazard", "cumhazard", 
         "quantile"), terms = c("bresponse", "binteracting", "bshifting"), 
         q = NULL, p = NULL, K = 50, interpolate = TRUE, ...)
\method{simulate}{mlt}(object, nsim = 1, seed = NULL, newdata = object$data, K = 50, 
         interpolate = TRUE, ...)
}
\arguments{
  \item{object}{a fitted conditional transformation model as returned by \code{\link{mlt}}}
  \item{x}{a fitted conditional transformation model as returned by \code{\link{mlt}}}
  \item{newdata}{an optional data frame of observations}
  \item{type}{type of prediction or plot to generate}
  \item{q}{quantiles at which to evaluate the model}
  \item{p}{probabilities for the evaluation of the quantile function (\code{type = "quantile"})}
  \item{terms}{terms to evaluate for the predictions, corresponds to the argument
     \code{response}, \code{interacting} and \code{shifting} in \code{\link{ctm}}}
  \item{K}{number of grid points to generate (in the absence of \code{q})}
  \item{col}{color for the lines to plot}
  \item{add}{logical indicating if a new plot shall be generated (the default)}
  \item{interpolate}{logical indicating if quantiles shall be interpolated linearily}
  \item{nsim}{number of samples to generate}
  \item{seed}{optional seed for the random number generator}
  \item{\dots}{additional arguments}
}
\details{
  \code{plot} evaluates the transformation function over a grid of \code{q} values
  for all observations in \code{newdata} and plots these functions (according to 
  \code{type}). \code{predict} evaluates the transformation function over a grid
  of \code{q} values for all observations in \code{newdata} and returns the
  result as a matrix (where _columns_ correspond to _rows_ in \code{newdata}).
  \code{simulate} draws samples from \code{object} by numerical inversion of the
  quantile function.
}
