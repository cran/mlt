\name{plot-predict-simulate}
\alias{predict.ctm}
\alias{predict.mlt}
\alias{simulate.ctm}
\alias{simulate.mlt}
\alias{plot.ctm}
\alias{plot.mlt}
\title{
  Plots, Predictions and Samples from mlt Objects
}
\description{
  Plot, predict and sample from objects of class mlt
}
\usage{
\method{plot}{ctm}(x, newdata, type = c("distribution", "survivor", "density", 
     "logdensity", "hazard", "loghazard", "cumhazard", "logcumhazard", "odds",
     "logodds", "quantile", "trafo"),
     q = NULL, prob = 1:(K - 1) / K, K = 50, col = rgb(.1, .1, .1, .1), lty = 1, 
     add = FALSE, ...)
\method{plot}{mlt}(x, ...)
\method{predict}{ctm}(object, newdata, type = c("trafo", "distribution", 
         "survivor", "density", "logdensity", "hazard", "loghazard", "cumhazard", 
         "logcumhazard", "odds", "logodds", "quantile"), 
         terms = c("bresponse", "binteracting", "bshifting"), 
         q = NULL, prob = NULL, K = 50, interpolate = FALSE, ...)
\method{predict}{mlt}(object, newdata = object$data, ...)
\method{simulate}{ctm}(object, nsim = 1, seed = NULL, newdata, K = 50, q = NULL,
         interpolate = FALSE, bysim = TRUE, ...)
\method{simulate}{mlt}(object, nsim = 1, seed = NULL, newdata = object$data, bysim = TRUE, ...)
}
\arguments{
  \item{object}{a fitted conditional transformation model as returned by \code{\link{mlt}}
                or an unfitted conditional transformation model as returned by \code{\link{ctm}}}
  \item{x}{a fitted conditional transformation model as returned by \code{\link{mlt}}}
  \item{newdata}{an optional data frame of observations}
  \item{type}{type of prediction or plot to generate}
  \item{q}{quantiles at which to evaluate the model}
  \item{prob}{probabilities for the evaluation of the quantile function (\code{type = "quantile"})}
  \item{terms}{terms to evaluate for the predictions, corresponds to the argument
     \code{response}, \code{interacting} and \code{shifting} in \code{\link{ctm}}}
  \item{K}{number of grid points to generate (in the absence of \code{q})}
  \item{col}{color for the lines to plot}
  \item{lty}{line type for the lines to plot}
  \item{add}{logical indicating if a new plot shall be generated (the default)}
  \item{interpolate}{logical indicating if quantiles shall be interpolated
                     linearily. This unnecessary option is no longer implemented (starting with 1.2-1).}
  \item{nsim}{number of samples to generate}
  \item{seed}{optional seed for the random number generator}
  \item{bysim}{logical, if \code{TRUE} a list with \code{nsim} elements is returned,
               each element is of length \code{nrow(newdata)} and 
               contains one sample from the conditional distribution for each
               row of \code{newdata}. If \code{FALSE}, a list of length \code{nrow(newdata)}
               is returned, its ith element of length \code{nsim} contains \code{nsim} samples
               from the conditional distribution given \code{newdata[i,]}.} 
  \item{\dots}{additional arguments}
}
\details{
  \code{plot} evaluates the transformation function over a grid of \code{q} values
  for all observations in \code{newdata} and plots these functions (according to 
  \code{type}). \code{predict} evaluates the transformation function over a grid
  of \code{q} values for all observations in \code{newdata} and returns the
  result as a matrix (where _columns_ correspond to _rows_ in
  \code{newdata}, see examples). Lack of \code{type = "mean"} is a feature
  and not a bug.

  Argument \code{type} defines the scale of the plots or predictions:
  \code{type = "distribution"} means the cumulative distribution function, 
  \code{type = "survivor"} is the survivor function (one minus distribution
  function), \code{type = "density"} the absolute continuous or discrete
  density (depending on the response), \code{type = "hazard"}, \code{type =
  "cumhazard"}, and \code{type = "odds"} refers to the hazard (absolute
  continuous or discrete), cumulative hazard (defined as minus log-survivor
  function in both the absolute continuous and discrete cases), and odds
  (distribution divided by survivor) functions. The quantile function can be
  evaluated for probabilities \code{prob} by \code{type = "quantile"}.
  
  Note that the \code{predict} method for \code{ctm} objects requires all
  model coefficients to be specified in this unfitted model.
  \code{simulate} draws samples from \code{object} by numerical inversion of the
  quantile function.

  Note that offsets are ALWAYS IGNORED when computing predictions. If you
  want the methods to pay attention to offsets, specify them as a variable
  in the model with fixed regression coefficient using the \code{fixed}
  argument in \code{\link{mlt}}.

  More examples can be found in Hothorn (2018).

}
\references{

  Torsten Hothorn (2020), Most Likely Transformations: The mlt Package,
  \emph{Journal of Statistical Software}, \bold{92}(1), 1--68,
  \doi{10.18637/jss.v092.i01}

}
\examples{

  library("survival")
  op <- options(digits = 2)

  ### GBSG2 dataset
  data("GBSG2", package = "TH.data")

  ### right-censored response
  GBSG2$y <- with(GBSG2, Surv(time, cens))

  ### define Bernstein(log(time)) parameterisation
  ### of transformation function. The response
  ### is bounded (log(0) doesn't work, so we use log(1))
  ### support defines the support of the Bernstein polynomial
  ### and add can be used to make the grid wider (see below)
  rvar <- numeric_var("y", bounds = c(0, Inf), 
                      support = c(100, 2000))
  rb <- Bernstein_basis(rvar, order = 6, ui = "increasing")
  ### dummy coding of menopausal status
  hb <- as.basis(~ 0 + menostat, data = GBSG2)
  ### treatment contrast of hormonal treatment
  xb <- as.basis(~ horTh, data = GBSG2, remove_intercept = TRUE)

  ### set-up and fit Cox model, stratified by menopausal status
  m <- ctm(rb, interacting = hb, shifting = xb, todistr = "MinExtrVal")
  fm <- mlt(m, data = GBSG2)

  ### generate grid for all three variables
  ### note that the response grid ranges between 1 (bounds[1])
  ### and 2000 (support[2])
  (d <- mkgrid(m, n = 10))
  ### data.frame of menopausal status and treatment
  nd <- do.call("expand.grid", d[-1])

  ### plot model on different scales, for all four combinations
  ### of menopausal status and hormonal treatment
  typ <- c("distribution", "survivor", "density", "hazard", 
           "cumhazard", "odds")
  layout(matrix(1:6, nrow = 2))
  nl <- sapply(typ, function(tp) 
      ### K = 500 makes densities and hazards smooth
      plot(fm, newdata = nd, type = tp, col = 1:nrow(nd), K = 500))
  legend("topleft", lty = 1, col = 1:nrow(nd), 
         legend = do.call("paste", nd), bty = "n")

  ### plot calls predict, which generates a grid with K = 50
  ### response values
  ### note that a K x nrow(newdata) matrix is returned
  ### (for reasons explained in the next example)
  predict(fm, newdata = nd, type = "survivor")

  ### newdata can take a list, and evaluates the survivor
  ### function on the grid defined by newdata 
  ### using a linear array model formulation and is 
  ### extremely efficient (wrt computing time and memory)
  ### d[1] (the response grid) varies fastest
  ### => the first dimension of predict() is always the response,
  ### not the dimension of the predictor variables (like one 
  ### might expect)
  predict(fm, newdata = d, type = "survivor")

  ### owing to this structure, the result can be quickly stored in 
  ### a data frame as follows
  cd <- do.call("expand.grid", d)
  cd$surv <- c(S <- predict(fm, newdata = d, type = "survivor"))

  ### works for distribution functions
  all.equal(1 - S, predict(fm, newdata = d, type = "distribution"))
  ### cumulative hazard functions
  all.equal(-log(S), predict(fm, newdata = d, type = "cumhazard"))
  ### log-cumulative hazard functions (= trafo, for Cox models)
  all.equal(log(-log(S)), predict(fm, newdata = d, type = "logcumhazard"))
  all.equal(log(-log(S)), predict(fm, newdata = d, type = "trafo"))
  ### densities, hazards, or odds functions
  predict(fm, newdata = d, type = "density")
  predict(fm, newdata = d, type = "hazard")
  predict(fm, newdata = d, type = "odds")
  ### and quantiles (10 and 20%)
  predict(fm, newdata = d[-1], type = "quantile", prob = 1:2 / 10)

  ### note that some quantiles are only defined as intervals
  ### (> 2000, in this case). Intervals are returned as an "response" 
  ### object, see ?R. Unfortunately, these can't be stored as array, so
  ### a data.frame is returned where the quantile varies first
  p <- c(list(prob = 1:9/10), d[-1])
  np <- do.call("expand.grid", p)
  (Q <- predict(fm, newdata = d[-1], type = "quantile", prob = 1:9 / 10))
  np$Q <- Q
  np

  ### simulating from the model works by inverting the distribution 
  ### function; some obs are right-censored at 2000
  (s <- simulate(fm, newdata = nd, nsim = 3))
  ### convert to Surv
  sapply(s, as.Surv)

  ### generate 3 parametric bootstrap samples from the model
  tmp <- GBSG2[, c("menostat", "horTh")]
  s <- simulate(fm, newdata = tmp, nsim = 3)
  ### refit the model using the simulated response
  lapply(s, function(y) {
    tmp$y <- y
    coef(mlt(m, data = tmp))
  })

  options(op)

}
