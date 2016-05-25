
\name{NEWS}
\title{NEWS file for the mlt package}

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

