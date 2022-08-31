% 2022-08-29 A. Papritz
% cd ~/R.user.home/soilhypfit/soilhypfit_01-7/pkg/man
% R CMD Rdconv -t html -o bla.html soilhypfit_internal_functions.Rd ; open bla.html; R CMD Rd2pdf --force soilhypfit_internal_functions.Rd
% hunspell -d en_GB-ise internal_functions.Rd

\encoding{utf8}

\name{internal_functions}
\alias{d1ineq_constraint_nlp}
\alias{d1l0}
\alias{d1lc}
\alias{estimate_lp}
\alias{estimate_nlp}
\alias{fit_wrc_hcc_fit}
\alias{gradient_nlp}
\alias{ineq_constraint_nlp}
\alias{model_names_nlp_lp}
\alias{model_fit_param_consistent}
\alias{model_param_tf_nlp_identity}
\alias{objective_nlp}
\alias{stop_cluster}

\title{Internal Functions of Package \pkg{soilhypfit}}

% description

\description{
The unexported internal functions
  \itemize{
    \item{\code{d1ineq_constraint_nlp}}
    \item{\code{d1l0}}
    \item{\code{d1lc}}
    \item{\code{estimate_lp}}
    \item{\code{estimate_nlp}}
    \item{\code{fit_wrc_hcc_fit}}
    \item{\code{gradient_nlp}}
    \item{\code{ineq_constraint_nlp}}
    \item{\code{model_names_nlp_lp}}
    \item{\code{model_fit_param_consistent}}
    \item{\code{model_param_tf_nlp_identity}}
    \item{\code{objective_nlp}}
    \item{\code{stop_cluster}}
}
are not intended for direct use.  However, as any unexported function, they
can be accessed by typing \code{soilhypfit:::\var{function-name}}.
}

% author

\author{
  Andreas Papritz \email{papritz@retired.ethz.ch}.
}

% see also

\seealso{
  \code{\link{soilhypfitIntro}} for a description of the models and a brief
  summary of the parameter estimation approach;

  \code{\link{fit_wrc_hcc}} for (constrained) estimation of parameters of
  models for soil water retention and hydraulic conductivity data;

  \code{\link{control_fit_wrc_hcc}} for options to control
  \code{fit_wrc_hcc};

  \code{\link{soilhypfitmethods}} for common S3 methods for class
  \code{fit_wrc_hcc};

  \code{\link{vcov}} for computing (co-)variances of the estimated
  nonlinear parameters;

  \code{\link{prfloglik_sample}} for profile loglikelihood
  computations;

  \code{\link{wc_model}} and \code{\link{hc_model}} for currently
  implemented models for soil water retention curves and hydraulic
  conductivity functions;

  \code{\link{evaporative-length}} for physically constraining parameter
  estimates of soil hydraulic material functions.
}

\keyword{internal}
