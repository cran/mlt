
\name{NEWS}
\title{NEWS file for the mlt package}

\section{Changes in Version 1.7-0 (2025-10-17)}{
  \subsection{Features}{
    \itemize{
      \item{\code{mltoptim()} features \code{nlminb} for unconstrained optimisation.}
      \item{\code{mltoptim()} features \code{optim} for unconstrained optimisation.}
      \item{\code{mltoptim()} defaults changed for all algorithms.}
      \item{\code{mltoptim()} refrains from generating random starting
            values now. Consequently, results are fully reproducible in the absence of
            user-defined seeds.}
      \item{The log-likelihood for infeasible parameter configurations
            is now always \code{-Inf}.}
      \item{Argument \code{scale} was renamed to \code{scaleparm} because
            of a conflict in \code{tram::Survreg}.}
      \item{Compute better starting values, especially for interval-censored
            outcomes using the fast Turnbull estimators provided 
            by \pkg{icenReg}.}
      \item{Implement (internal) pre-transformation based on nonparametric
            unconditional distribution to be used by \pkg{tram}.}
    }
  }
  \subsection{Bugfixes}{
    \itemize{
      \item{\code{summary.mmlt} works now.}
      \item{Improved error messages in \code{mmlt}.}
      \item{\code{R(, as.R.interval = TRUE)} produces correct empirical
            likelihood for numeric outcomes.}
      \item{\code{mmlt} only uses optimisers able to approximate the
            Hessian by default.}
      \item{\code{simulate.ctm} might have got internal dimensions wrong.}
    }
  }
}

\section{Changes in Version 1.6-6 (2025-06-17)}{
  \subsection{Features}{
    \itemize{
      \item{Simpler switch for computing numerically approximated Hessians in \code{mltoptim}.}
    }
  }
}

\section{Changes in Version 1.6-5 (2025-05-09)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{\code{rep} to \code{rep_len} when possible.}
      \item{Avoid partial argument matches.}
    }
  }
}


\section{Changes in Version 1.6-4 (2025-04-29)}{
  \subsection{Features}{
    \itemize{
      \item{Improve \code{R()} for interval censored data.}
      \item{Add \code{constrOptim} to \code{mltoptim}.}
    }
  }
  \subsection{Bugfixes}{
    \itemize{
      \item{Handle fixed starting values in \code{mmlt} correctly.}
    }
  }
}

\section{Changes in Version 1.6-3 (2025-03-20)}{
  \subsection{Features}{
    \itemize{
      \item{Add \code{as.Surv} method for numeric vectors.}
    }
  }
}


\section{Changes in Version 1.6-2 (2025-01-29)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Fall back to \pkg{coneproj} if \pkg{quadprog} fails to get starting values.}
      \item{Stabilise bugfixes and avoid ATLAS issues.}
    }
  }
}


\section{Changes in Version 1.6-1 (2024-11-18)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Try harder to obtain starting values meeting constraints, 
            issue a warning if this failed.}
    }
  }
}


\section{Changes in Version 1.6-0 (2024-09-23)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Allow interval censoring using \code{-Inf} and \code{Inf}.
            Starting values are affected but the corresponding observations
            contribute zero to the log-likelihood.}
     \item{Make sure \code{R()} does not drop any unused levels (unlike \code{factor()}).}
    }
  }
  \subsection{New features}{
    \itemize{
      \item{\code{mmlt} moved from \pkg{tram} allowing arbitrary mixes of continuous
          and discrete data, including exact and censored observations
          for the same variable and thus also allowing non-informative missings encoded as
          \code{-Inf, Inf}, see below. Corresponding \code{logLik} and \code{estfun}
          methods feature \code{newdata} and \code{weights} arguments, as the
          methods for \code{mlt} always did.}
      \item{\code{R()} translates missing response values to  \code{-Inf} and \code{Inf} interval censoring.
            This means that missing values are allowed in response variables
            and the corresponding observations do
            not contribute to the log-likelihood. In contrast to
            \code{na.omit}, the (technical) dimension of the training data is not
            changed.}
    }
  }
}

\section{Changes in Version 1.5-2 (2024-08-24)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Handle \code{Surv} objects with smaller than zero values in \code{R.Surv(..., as.R.interval = TRUE)}.}
    }
  }
  \subsection{New features}{
    \itemize{
      \item{Add \code{maxeval = 1000L} for \code{nloptr}.}
      \item{Report \code{nloptr} more detailed status messages.}
    }
  }
}

\section{Changes in Version 1.5-1 (2024-04-11)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Remove warning when dealing with binary factors.}
    }
  }
}

\section{Changes in Version 1.5-0 (2023-12-13)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Do no use zero as lower bound (when running \code{R(Surv(...), as.R.interval =
            TRUE)}) except when zero is first death time.}
    }
  }
  \subsection{New features}{
    \itemize{
      \item{Add slots to \code{mmlt} for externally computing transformations and derivatives 
            thereof, needed by \code{tram::mmlt}.}
    }
  }
}

\section{Changes in Version 1.4-9 (2023-08-21)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Make sure \code{R(Surv(...), as.R.ordered = TRUE)} matches 
            Kaplan-Meier for right-censored and possibly tied observations.
            Such responses give NPML estimates.}
    }
  }
  \subsection{New features}{
    \itemize{
      \item{\code{R(Surv(...), as.R.numeric = TRUE)} can be used to fit
            smooth models by maximising the nonparametric likelihood.}
    }
  }
}


\section{Changes in Version 1.4-8 (2023-06-27)}{
  \subsection{New features}{
    \itemize{
      \item{Better support for frailty parameter.}
    }
  }
}

\section{Changes in Version 1.4-7 (2023-05-22)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Starting values might end up being NA; set to zero.}
    }
  }
  \subsection{New features}{
    \itemize{
      \item{\code{mlt(..., theta = theta, dofit = TRUE)} returns a "fitted"
            \code{mlt} model with coefficients \code{theta}. 
            If \code{theta} is missing, an unfitted model is returned.}
    }
  }
}

\section{Changes in Version 1.4-6 (2023-04-14)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Register internal methods.}
    }
  }
  \subsection{New features}{
    \itemize{
      \item{Re-add support for \pkg{nloptr} in \code{mltoptim()}.}
    }
  }
}


\section{Changes in Version 1.4-5 (2023-03-09)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Fix S3 argument mismatches.}
      \item{Try to detect negative standard errors in case optimisation failed silently.}
    }
  }
}


\section{Changes in Version 1.4-4 (2022-11-02)}{
  \subsection{New features}{
    \itemize{
      \item{Improve numerical stability in \code{simulate.ctm}.}
    }
  }
}


\section{Changes in Version 1.4-3 (2022-10-31)}{
  \subsection{New features}{
    \itemize{
      \item{\code{predict} and \code{plot} have new types involving logs.}
    }
  }
  \subsection{Bugfixes}{
    \itemize{
      \item{\code{plot()} ignored  \code{lwd}.}
    }
  }
}

\section{Changes in Version 1.4-2 (2022-06-30)}{
  \subsection{New features}{
    \itemize{
      \item{\code{as.R.interval} allows evaluation of nonparametric
            likelihood for smoothly defined models.}
      \item{New Laplace link function, contributed by Ainesh Sewak.}
      \item{New cauchit link function.}
    }
  }
  \subsection{Bugfixes}{
    \itemize{
      \item{\code{drop = TRUE} in score for shift-scale models.}
      \item{Eliminate fallout of fix of PR#17616 in Rout.save files.}
    }
  }
}


\section{Changes in Version 1.4-1 (2022-03-31)}{
  \subsection{New features}{
    \itemize{
      \item{\code{update} allows fixing parameters.}
    }
  }
  \subsection{Bugfixes}{
    \itemize{
      \item{Fix discrete gradient for shift-scale models.}
    }
  }
}

\section{Changes in Version 1.4-0 (2022-01-14)}{
  \subsection{New features}{
    \itemize{
      \item{Add infrastructure for shift-scale transformation models.}
    }
  }
  \subsection{Documentation}{
    \itemize{
      \item{\code{Gradient} and \code{estfun} always returned 
            negative scores, this is documented now.}
    }
  }
  \subsection{Bugfixes}{
    \itemize{
      \item{Arguments passed to \code{mlt} via dots were silently ignored.}
      \item{Remove nloptr dependency for the time being.}
    }
  }
}


\section{Changes in Version 1.3-2 (2021-08-24)}{
  \subsection{New features}{
    \itemize{
      \item{Add \code{as.double} for \code{response} objects 
            (replaces \code{trtf:::.R2vec}).}
    }
  }
}

\section{Changes in Version 1.3-1 (2021-03-22)}{
  \subsection{New features}{
    \itemize{
      \item{Improved numerical stability for censored data.}
    }
  }
}

\section{Changes in Version 1.3-0 (2021-03-03)}{
  \subsection{New features}{
    \itemize{
      \item{New argument \code{as.R.ordered} allowing numeric and survival
            responses to be coded as ordered factors, for nonparametric maximum
            likelihood estimation.}
      \item{Allow sparse model matrices. This is useful for nonparametric
            maximum likelihood estimation with many distinct outcomes.}
      \item{Some speed-ups.}
    }
  }
  \subsection{Bugfixes}{
    \itemize{
      \item{Plotting of quantiles sometimes failed because inversion of cdf
            was not possible for certain quantiles. These are now removed before
            plotting.}
      \item{Fitting models to interval censored responses containing intervals
            \code{c(-Inf, Inf)} failed.}
      \item{Always return names score matrices and residuals.}
    }
  }
}


\section{Changes in Version 1.2-3 (2021-02-22)}{
  \subsection{New features}{
    \itemize{
      \item{Improve documentation.}
    }
  }
}


\section{Changes in Version 1.2-2 (2021-02-12)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Sampling from unconditional models did not pay attention to
            number of observations.}
    }
  }
}

\section{Changes in Version 1.2-1 (2021-02-03)}{
  \subsection{New features}{
    \itemize{
      \item{Quantiles and thus simulations are now computationally more exact and more robust.
            The unnecessary \code{interpolate} argument to \code{predict} and
            \code{simulate} is now ignored.}
      \item{Adjust contrasts a fixed parameter contributes to.}
      \item{Return numerically determined Hessians upon request.}
      \item{Implement frailty error distributions, experimentally and internal only.}
      \item{Implement cure mixture models, experimentally and internal only.}
      \item{Improve computations of log-probabilities.}
    }
  }
  \subsection{Bugfixes}{
    \itemize{
      \item{Discrete hazard functions were incorrect.}
    }
  }
}


\section{Changes in Version 1.2-0 (2020-05-11)}{
  \itemize{
    \item{Add exponential distribution (for Aalen additive hazards models).}
    \item{Pay attention to model class when computing cumulative hazards.}
    \item{Add log-cumulative hazards, log-odds, and odds for predictions and plots.}
  }
}


\section{Changes in Version 1.1-2 (2020-02-13)}{
  \itemize{
    \item{Allow permutations of single variables.}
    \item{Update citation info.}
  }
}

\section{Changes in Version 1.1-1 (2020-01-09)}{
  \itemize{
    \item{Try harder to invert Hessians.}
    \item{Update reference output.}
  }
}

\section{Changes in Version 1.1-0 (2019-11-06)}{
  \itemize{
    \item{Add support for \pkg{nloptr} (still experimental and thus switched off by default).}
    \item{Make sure \code{coef()} always returns named argument.}
    \item{Fix problem in \code{as.Surv} reported by Balint Tamasi.}
  }
}


\section{Changes in Version 1.0-7 (2019-10-07)}{
  \itemize{
    \item{Less paranoia in \file{bugfixes.R}.}
  }
}

\section{Changes in Version 1.0-6 (2019-10-04)}{
  \itemize{
    \item{Return Hessian for fixed parameters if requested.}
    \item{Fix subsetting problem in \code{R.numeric}.}
    \item{Allow to \code{update} offsets.}
  }
}


\section{Changes in Version 1.0-5 (2019-06-20)}{
  \itemize{
    \item{Add a \code{bread} method.}
    \item{Check response variable against observations in \code{data}.}
    \item{Make sure integers larger zero are handled correctly in \code{R}.}
    \item{Implement \code{resid} method, ie the score wrt a constant.}
    \item{Cox examples with Bernstein polynomials of log-time.}
  }
}

\section{Changes in Version 1.0-4 (2018-12-03)}{
  \itemize{
    \item{Arguments \code{K} and \code{cheat} where
          ignored by \code{confband} when \code{newdata}
          had multiple rows.}
    \item{Computation of starting values more robust now.}
    \item{Order of fixed parameters (\code{fixed} argument to \code{mlt})
          might have been wrong due to incomplete matching.}
  }
}

\section{Changes in Version 1.0-3 (2018-09-17)}{
  \itemize{
    \item{Add \code{lty} argument to \code{plot.ctm}.}
    \item{\code{update} needs free coefficients only.}
    \item{Internal interface changes.}
  }
}


\section{Changes in Version 1.0-2 (2018-06-22)}{
  \itemize{
    \item{Make sure transformation functions outside \code{bounds} are minus
          or plus \code{Inf}.}
    \item{Initial guestimates for ordered responses were incorrect
          and may potentially have led to nonsense results.}
    \item{Some smaller improvements in computation of log-likelihoods and
          scores with respect to accuracy and speed.}
    \item{\code{print} respects \code{options(digits)}.}
  }
}

\section{Changes in Version 1.0-1 (2018-04-20)}{
  \itemize{
    \item{\code{estfun, parm = coef(object, fixed = TRUE))} evaluates
          scores for all model parameters, including fixed ones.}
    \item{\code{logLik(..., newdata, w)} ignored weights \code{w} when
          \code{newdata} was given. Same problem was also fixed for
          \code{estfun}.}
  }
}


\section{Changes in Version 1.0-0 (2018-03-07)}{
  \itemize{
    \item{A paper describing version 1.0-0 of the \pkg{mlt}, \pkg{basefun}, and
          \pkg{variables} packages was accepted for publication in the Journal of
          Statistical Software 2018-03-05.}
    \item{Documentation updates.}
  }
}

\section{Changes in Version 0.2-2 (2018-02-23)}{
  \itemize{
    \item{Use \pkg{coneprog} for getting the starting values.}
    \item{\code{logLik} and \code{estfun} accept matrices as \code{parm}
          argument for the evalution of log-likelihoods and scores
          with subject-specific parameters (for example in transformation
          trees or forests and boosting procedures.}
  }
}

\section{Changes in Version 0.2-1 (2017-12-12)}{
  \itemize{
    \item{\code{q} is forwarded to \code{qmlt} by \code{predict.ctm} now.}
    \item{\code{p} is now \code{prob} in \code{qmlt} and thus \code{predict.ctm}.}
    \item{Update citation.}
  }
}


\section{Changes in Version 0.2-0 (2017-06-20)}{
  \itemize{
    \item{Most Likely Transformations will be published in
          the Scandinavian Journal of Statistics.}
    \item{Import package \pkg{alabama}.}
    \item{\code{as.Surv(R(Surv(...)))} returns \code{Surv(...)},
          useful for converting output by \code{simulate}
          to \code{Surv} objects.}
  }
}

\section{Changes in Version 0.1-3 (2017-04-04)}{
  \itemize{
    \item{Add \code{subset} argument to \code{update} (for faster transformation 
          trees and forests).}
    \item{Sum over score contributions with positive weight only when
          evaluating the gradient.}
    \item{Having all response observations being interval-censored is
          allowed again (too heavy checking was in place).}
    \item{Don't try to numerically check KKT conditions automatically.}
    \item{Check for unused arguments in dots where necessary.}
    \item{Make sure the score doesn't get too large (avoid division by near
          zero probabilities).}
    \item{Improve \code{survfit} to compute non-parametric unconditional
          probabilities for obtaining starting values in the presence of
          censoring and truncation.}
 }
}

\section{Changes in Version 0.1-2 (2017-02-14)}{
  \itemize{
    \item{\code{logLik} with \code{newdata} argument ignored \code{parm} 
          and \code{weights} arguments.}
    \item{\code{estfun} now also has a \code{newdata} argument.}
    \item{Correct axes labelling when plotting quantile functions.}
    \item{make sure names are correct in \code{coef(model, fixed = FALSE)}.}
    \item{check if any exact or interval-censored response with non-zero 
          weight exists before trying to fit the model.}
    \item{make checks a little more robust against huge diffs under Windows.}
 }
}

\section{Changes in Version 0.1-1 (2016-11-29)}{
  \itemize{
    \item{Fix two bugs in computation of log-likelihood for possibly
          left-truncated responses such as \code{Surv(start, time, status)}.}
 }
}

\section{Changes in Version 0.1-0 (2016-11-10)}{
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

