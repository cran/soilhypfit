% 2021-10-13 A. Papritz
% cd ~/R.user.home/soilhypfit/soilhypfit_01-1/pkg/man
% R CMD Rdconv -t html -o bla.html soilhypfit_vcov.Rd ; open bla.html; R CMD Rd2pdf --force soilhypfit_vcov.Rd;
% hunspell -d en_GB-ise soilhypfitS3methods.Rd

\encoding{utf8}

\name{vcov}
\alias{vcov}
\alias{coef.vcov_fit_wrc_hcc}
\alias{vcov.fit_wrc_hcc}

\concept{soil water}
\concept{soil water characteristic}
\concept{soil water retention curve}
\concept{hydraulic conductivity function}
\concept{saturated hydraulic conductivity}

% \newcommand{\mbvec}{\eqn{\mbox{\boldmath$#1$\unboldmath}^\mathrm{#2} #4}{#1#3#4}}

\title{\code{vcov} Method for Class \code{fit_wrc_hcc}}

% description

\description{This page documents the method \code{vcov} for the class
\code{fit_wrc_hcc} and its \code{coef} method.  \code{vcov} extracts the
covariance matrices of the nonlinear parameters
\eqn{\widehat{\mbox{\boldmath$\nu$\unboldmath}}}{hat\nu} estimated by
maximum likelihood or maximum posterior density.
}

% usage

\usage{

\method{vcov}{fit_wrc_hcc}(object, subset = NULL, grad_eps,
    bound_eps = sqrt(.Machine$double.eps), ...)

\method{coef}{vcov_fit_wrc_hcc}(object, se = TRUE, correlation = se,
    status = FALSE, ...)

}

% arguments

\arguments{

  \item{object}{either an object of class \code{fit_wrc_hcc} for
  \code{vcov} or an object of class \code{vcov_fit_wrc_hcc} for
  \code{coef}.}

  \item{subset}{an integer, character or logical vector to the choose the
  soil samples for which covariance matrices should be extracted.
  Defaults to \code{NULL}, which extracts the covariances for all soil samples.}

  \item{grad_eps}{a numeric scalar defining a critical magnitude of the
  moduli of scaled gradient components so that they are considered to be
  approximately equal to zero, see \emph{Details}.}

  \item{bound_eps}{a numeric scalar defining the critical difference
  between parameter estimates and the boundaries of the parameter space so
  that the estimates are considered to be identical to the boundary values,
  see \emph{Details}.}

  \item{se}{a logical scalar to control whether standard errors of the
  estimated nonlinear parameters
  \eqn{\widehat{\mbox{\boldmath$\nu$\unboldmath}}}{hat\nu} should be
  returned (\code{TRUE}, default) or variances (\code{FALSE}).}

  \item{correlation}{a logical scalar to control whether correlations
  (\code{TRUE}, default) or covariances (\code{FALSE}) of the esitmated
  nonlinear parameters
  \eqn{\widehat{\mbox{\boldmath$\nu$\unboldmath}}}{hat\nu} should be
  returned.}

  \item{status}{a logical scalar to control whether diagnostics should be
  returned along with the results.}

  \item{...}{additional arguments passed to methods, currently not used.}

  }

% details

\details{

  The function \code{vcov} extracts (co-)variances of the nonlinear
  parameters from the inverse Hessian matrix of the objective function at
  the solution \eqn{\widehat{\mbox{\boldmath$\nu$\unboldmath}}}{hat\nu} for
  mpd and ml estimates, see \code{\link{soilhypfitIntro}} and \cite{Stewart
  and \enc{Sørensen}{Soerensen} (1981)}.

  \code{vcov} checks whether the gradient at the solution is approximately
  equal to zero and issues a warning if this is not the case.  This is
  controlled by the argument \code{grad_eps} which is the tolerable largest
  modulus of the scaled gradient (= gradient divided by the absolute value
  of objective function) at the solution.  The function
  \code{\link{control_fit_wrc_hcc}} selects a default value for
  \code{grad_eps} in the dependence of the chosen optimisation approach
  (argument \code{settings} of \code{\link{control_fit_wrc_hcc}}).

  \code{vcov} sets covariances equal to \code{NA} if the parameter
  estimates differ less than \code{bound_eps} from the boundaries of the
  parameter space as defined by \code{\link{param_boundf}}.

}

% value

\value{

  The method \code{vcov} returns an object of of class
  \code{vcov_fit_wrc_hcc}, which is a list of covariance matrices of the
  estimated nonlinear parameters for the soil samples.  The attribute
  \code{status} of the matrices qualifies the covariances.

  The \code{coef} method for class \code{vcov_fit_wrc_hcc} extracts the
  entries of the covariances matrices, optionally computes standard errors
  and correlation coefficients and returns the results in a dataframe.

}

% references

\references{
  Stewart, W.E. and \enc{Sørensen}{Soerensen}, J.P. (1981)
  Bayesian estimation of common
  parameters from multiresponse data with missing observations.
  \emph{Technometrics}, \bold{23}, 131--141, \cr
  \doi{10.1080/00401706.1981.10486255}.
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

%   \code{\link{vcov}} for computing (co-)variances of the estimated
%   nonlinear parameters;

  \code{\link{wc_model}} and \code{\link{hc_model}} for currently
  implemented models for soil water retention curves and hydraulic
  conductivity functions;

  \code{\link{evaporative-length}} for physically constraining parameter
  estimates of soil hydraulic material functions.
}

% examples

\examples{
\donttest{
# use of \donttest{} because execution time exceeds 5 seconds

data(sim_wrc_hcc)

# define number of cores for parallel computations
if(interactive()) ncpu <- parallel::detectCores() - 1L else ncpu <- 1L

# estimate parameters for 3 samples by unconstrained, global optimisation
# algorithm NLOPT_GN_MLSL
# sample 1: use only conductivity data
# sample 2: use only water content data
# sample 3: use both types of data
rfit_uglob <- fit_wrc_hcc(
  wrc_formula = wc ~ head | id,
  hcc_formula = hc ~ head | id,
  wrc_subset = id != 1,
  hcc_subset = id != 2,
  data = sim_wrc_hcc,
  control = control_fit_wrc_hcc(pcmp = control_pcmp(ncores = ncpu)))
print(rfit_uglob)
summary(rfit_uglob)
coef(rfit_uglob, what = "nonlinear")
coef(rfit_uglob, what = "linear", gof = TRUE)
coef(vcov(rfit_uglob), status = TRUE, se = FALSE)
op <- par(mfrow = c(3, 2))
plot(rfit_uglob)
on.exit(par(op))

}}

\keyword{methods}
\keyword{print}
