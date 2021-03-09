\name{R}
\alias{R}
\alias{R.Surv}
\alias{R.factor}
\alias{R.ordered}
\alias{R.numeric}
\alias{R.integer}
\alias{R.list}
\alias{R.response}
\alias{as.Surv}
\alias{as.Surv.response}
\title{
    Response Variables
}
\description{
    Represent a possibly censored or truncated response variable
}
\usage{
R(object, ...)
\method{R}{numeric}(object = NA, cleft = NA, cright = NA,
   tleft = NA, tright = NA, tol = sqrt(.Machine$double.eps), 
   as.R.ordered = FALSE, ...)
\method{R}{ordered}(object, cleft = NA, cright = NA, ...)
\method{R}{integer}(object, cleft = NA, cright = NA, bounds = c(min(object), Inf), ...)
\method{R}{factor}(object, ...)
\method{R}{Surv}(object, as.R.ordered = FALSE, ...)
as.Surv(object)
\method{as.Surv}{response}(object)
}
\arguments{
  \item{object}{
  A vector of (conceptually) exact measurements or an object of class
  \code{response} (for \code{as.Surv}) or a list.}
  \item{cleft}{
  A vector of left borders of censored measurements}
  \item{cright}{
  A vector of right borders of censored measurements}
  \item{tleft}{
  A vector of left truncations}
  \item{tright}{
  A vector of right truncations}
  \item{tol}{
  Tolerance for checking if \code{cleft} < \code{cright}}
  \item{bounds}{
  Range of possible values for integers}
  \item{as.R.ordered}{logical, should numeric responses or right-censored
  (and possible left-truncated survival) times be coded as ordered factor?
  This leads to a parameterisation allowing to maximise the nonparametric
  maximum likelihood.}
  \item{\dots}{
  other arguments, ignored except for \code{tleft} and \code{tright} to
  \code{R.ordered} and \code{R.integer}}
}
\details{
  \code{R} is basically an extention of \code{\link[survival]{Surv}} for the
  representation of arbitrarily censored or truncated measurements at any scale.
  The \code{storage.mode} of \code{object} determines if models are fitted
  by the discrete likelihood (integers or factors) or the continuous
  likelihood (log-density for numeric \code{object}s). Interval-censoring
  is given by intervals (\code{cleft}, \code{cright}], also for integers and
  factors (see example below). Left- and right-truncation can be defined
  by the \code{tleft} and \code{tright} arguments. Existing \code{Surv}
  objects can be converted using \code{R(Surv(...))}$ and, in some cases, an
  \code{as.Surv()} method exists for representing \code{response} objects as
  \code{Surv} objects.

  \code{R} applied to a list calls \code{R} for each of the list elements
  and returns a joint object.

  More examples can be found in Hothorn (2018) and in 
  \code{vignette("tram", package = "tram")}.

}
\references{

  Torsten Hothorn (2020), Most Likely Transformations: The mlt Package,
  \emph{Journal of Statistical Software}, \bold{92}(1), 1--68,
  \doi{10.18637/jss.v092.i01}

}
\examples{

 library("survival")
 
 ### randomly right-censored continuous observations
 time <- as.double(1:9)
 event <- rep(c(FALSE, TRUE), length = length(time))

 Surv(time, event)
 R(Surv(time, event))

 ### right-censoring, left-truncation
 ltm <- 1:9 / 10
 Surv(ltm, time, event)
 R(Surv(ltm, time, event))

 ### interval-censoring
 Surv(ltm, time, type = "interval2")
 R(Surv(ltm, time, type = "interval2"))

 ### interval-censoring, left/right-truncation
 lc <- as.double(1:4)
 lt <- c(NA, NA, 7, 8)
 rt <- c(NA, 9, NA, 10)
 x <- c(3, NA, NA, NA)
 rc <- as.double(11:14)
 R(x, cleft = lt, cright = rt)
 as.Surv(R(x, cleft = lt, cright = rt))
 R(x, tleft = 1, cleft = lt, cright = rt)
 R(x, tleft = 1, cleft = lt, cright = rt, tright = 15)
 R(x, tleft = lc, cleft = lt, cright = rt, tright = rc)

 ### discrete observations: counts
 x <- 0:9
 R(x)
 ### partially interval-censored counts
 rx <- c(rep(NA, 6), rep(15L, 4))
 R(x, cright = rx)

 ### ordered factor
 x <- gl(5, 2, labels = LETTERS[1:5], ordered = TRUE)
 R(x)
 ### interval-censoring (ie, observations can have multiple levels)
 lx <- ordered(c("A", "A", "B", "C", "D", "E"), 
               levels = LETTERS[1:5], labels = LETTERS[1:5])
 rx <- ordered(c("B", "D", "E", "D", "D", "E"), 
               levels = LETTERS[1:5], labels = LETTERS[1:5])
 R(rx, cleft = lx, cright = rx)

 ### facilitate nonparametric maximum likelihood
 (y <- round(runif(10), 1))
 R(y, as.R.ordered = TRUE)

 R(Surv(time, event), as.R.ordered = TRUE)
 R(Surv(ltm, time, event), as.R.ordered = TRUE)

}
