% 2021-12-20 A. Papritz
% cd ~/R.user.home/soilhypfit/soilhypfit_01-5/pkg/man
% R CMD Rdconv -t html -o bla.html soilhypfit_utility_functions.Rd ; open bla.html; R CMD Rd2pdf --force soilhypfit_utility_functions.Rd;
% hunspell -d en_GB-ise utility_functions.Rd

\encoding{utf8}

\name{utility-functions}
\alias{check_param_boundf}
\alias{convergence_message}
\alias{extract_error_messages}
\alias{select_failed_fits}

\concept{soil water}
\concept{soil water characteristic}
\concept{soil water retention curve}
\concept{hydraulic conductivity function}
\concept{saturated hydraulic conductivity}

\title{Utility functions for package soilhypfit}

% description

\description{
  This page documents the functions \code{convergence_message},
  \code{extract_error_messages}, \code{select_failed_fits} and
  \code{check_param_boundf}.
}

% usage

\usage{
convergence_message(x, sce = FALSE)

extract_error_messages(object, start = 1, stop = 80, prt = TRUE)

select_failed_fits(object)

check_param_boundf(y, compare_thetar_thetas = TRUE)
}

% arguments

\arguments{
  \item{x}{an integer scalar issued by the optimisers of
  \code{\link[nloptr]{nloptr}} or \code{\link[SoilHyP]{SCEoptim}} on\cr
  (non-)convergence}

  \item{sce}{a logical scalar to select the optimiser
  \code{\link[nloptr]{nloptr}} (\code{FALSE}, default) or
  \code{\link[SoilHyP]{SCEoptim}} (\code{TRUE}).}

  \item{object}{an object of class \code{fit_wrc_hcc}, see
  \code{\link{fit_wrc_hcc}}.}

  \item{prt}{a logical scalar controlling whether the error messages should
  be printed.}

  \item{start, stop}{integer scalar with the first and last character to
  print.}

  \item{y}{a named list of numeric vectors of length 2 that define the
  allowed lower and upper bounds (box constraints) for the parameters of
  the models, see argument \code{param_bound} of
  \code{\link{control_fit_wrc_hcc}}.}

  \item{compare_thetar_thetas}{logical scalar to control cross-comparison
  of valid ranges of parameters \code{thetar} and \code{thetas}.}

}

% details

\details{
The function \code{convergence_message} prints a message that explains
the convergence codes, for \code{\link[nloptr]{nloptr}}, see
\href{https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#return-values}{NLopt
return values}.  The function \code{extract_error_messages} extract the
error messages of estimations that failed and optionally prints sub-strings
of them.\cr  The function \code{select_failed_fits} returns the
\code{id}s of the soil samples for which parameter estimation failed.
The function \code{check_param_boundf} checks the validity and consistecy
of bounds of box constraints of model parameters.
}

% value

\value{
The function \code{convergence_message} and \code{extract_error_messages}
return invisibly the convergence code or the error messages.
}

% references

\references{
  Duan, Q., Sorooshian, S., and Gupta, V. K. (1994) Optimal use of the
  SCE-UA global optimisation method for calibrating watershed models,
  \emph{Journal of Hydrology} \bold{158}, 265--284,
  \doi{10.1016/0022-1694(94)90057-4}.

  Johnson, S.G. The NLopt nonlinear-optimisation package.
  \url{https://github.com/stevengj/nlopt}.
}

% author

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}.
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

% examples

\examples{
convergence_message(3)
convergence_message(2, sce = TRUE)
}
