\name{soilhypfitNews}
\encoding{UTF-8}

\title{News for Package \pkg{soilhypfit}}

\section{Changes in version 0.1-7 (2022-08-31)}{
  \subsection{NEW FEATURES}{
  \itemize{

      \item new file \code{inst/NEWS.Rd}

      \item new internal functions for handling names of model parameters,
      fitting flags and parameter transformations

    }
  }

  \subsection{BUG FIXES}{
    \itemize{

      \item sanitized help pages complying with KaTeX

    }
  }
}


\section{Changes in version 0.1-5 (2022-01-07)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item issue of a warning if standard errors of nonlinear parameters
      are computed for a fit with global algorithm

      \item new functions \code{prfloglik_sample()} and
      \code{confint_prfloglik_sample()} for profile loglikelihood
      computations

      \item new default value \code{"ml"} for \code{methods} argument of
      \code{control_fit_wrc_hcc()}

      \item additional components in output of \code{fit_wrc_hcc()}

      \item new \code{confint} method for computing confidence intervals
      based on likelihood ratio test and the asymptotic normal distribution
      of ML estimates

    }

  }
  \subsection{BUG FIXES}{
    \itemize{

      \item adjust values of thetar and thetas if outside of allowed bounds

    }
  }
}


\section{Changes in version 0.1-4 (2021-12-07)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item new argument \code{head_saturation} for zero head values in
      functions \code{plot.fit_wrc_hcc()} and \code{lines.fit_wrc_hcc()}
      iss
      \item new defaults for local constrained and unconstrained algorithms

    }

  }
}
