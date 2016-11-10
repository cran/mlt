
\name{NEWS}
\title{NEWS file for the mlt package}

\section{Changes in Version 0.1-0 (2016-11-09)}{
  \itemize{
    \item{Add augmented lagrangian minimization (\code{auglag()} from package \pkg{alabama}).}
    \item{Make optimiation procedure more general and adaptive, allow users to change
          defaults and even add their own optimiser.}
 }
}


\section{Changes in Version 0.0-34 (2016-09-13)}{
  \itemize{
    \item{fix bug when calling \code{survfit} for computing initial probabilities.}
    \item{add \code{bysim} argument to \code{simulate}.}
    \item{make sure \code{checkGrad} is respected by \code{update}.}
    \item{\code{predict} computes \code{q} with \code{K} elements 
          if not given (as \code{plot} always did).}
 }
}

\section{Changes in Version 0.0-33 (2016-05-25)}{
  \itemize{
    \item{Make sure \code{times} are ordered before calling
          \code{survival::summary.survfit}.}
  }
}

\section{Changes in Version 0.0-32 (2016-05-24)}{
  \itemize{
    \item{Introduce \code{as.mlt} generic.}
    \item{Introduce a \code{coef} slot in \code{ctm} objects
          and a corresponding \code{coef<-} and \code{coef}
          method for setting and extracting coefficients
          to and from unfitted conditional transformation models.}
    \item{\code{predict}, \code{simulate} and \code{plot}
          work on \code{ctm} objects (with meaningful coefficients)
          now.}
  }
}

\section{Changes in Version 0.0-31 (2016-03-11)}{
  \itemize{
    \item{Some small improvements wrt run time and memory consumption.}
    \item{Use \code{theta = coef(object)} as default starting parameters in \code{update()}.}
    \item{\code{logLik} has a new \code{newdata} argument.}
    \item{\code{simulate} has a new \code{q} argument.}
  }
}

\section{Changes in Version 0.0-30 (2016-02-19)}{
  \itemize{
    \item{The \pkg{mlt} package was first 
          published on CRAN.}
  }
}

