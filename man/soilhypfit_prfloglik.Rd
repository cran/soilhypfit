% 2022-12-17 A. Papritz
% cd ~/R.user.home/soilhypfit/soilhypfit_01-7/pkg/man
% R CMD Rdconv -t html -o bla.html soilhypfit_prfloglik.Rd ; open bla.html; R CMD Rd2pdf --force soilhypfit_prfloglik.Rd;
% hunspell -d en_GB-ise soilhypfit_prfloglik.Rd

\encoding{utf8}

\name{profile_loglikelihood}
\alias{prfloglik_sample}
\alias{confint_prfloglik_sample}
\alias{confint.fit_wrc_hcc}

\concept{soil water}
\concept{hydraulic conductivity function}
\concept{saturated hydraulic conductivity}

\newcommand{\mbvec}{\eqn{\boldsymbol{#1}^\mathrm{#2} #3}}

\title{
  Profile Loglikelihoods and Confidence Intervals for
  Parametric Modelling of Soil Hydraulic Properties}

% description

\description{
  The function \code{prfloglik_sample} computes for a single soil sample a
  loglikelihood profile as a function of the specified values for subsets
  of model parameters of soil water retention and/or soil hydraulic
  conductivity functions.  The function \code{confint_prfloglik_sample}
  computes a confidence interval of one model parameter based on the
  likelihood ratio test for a single soil sample, and the S3 method
  \code{confint} computes confidence intervals of \emph{nonlinear}
  model parameters for multiple soil samples.}

% usage

\usage{
prfloglik_sample(object, values, soil_sample,
    ncores = min(detectCores() - 1L, NROW(values)), verbose = 0)

confint_prfloglik_sample(object, parm = names(default_fit_param()),
    soil_sample, level = 0.95, test = c("F", "Chisq"),
    denominator_df = c("nonlinear", "all"), param_bound = NULL,
    root_tol = .Machine$double.eps^0.25, froot_tol = sqrt(root_tol),
    verbose = 0)

\method{confint}{fit_wrc_hcc}(object,
    parm = names(object[["control"]][["initial_param"]]), level = 0.95,
    subset = NULL, type = c("loglik", "normal"), test = c("F", "Chisq"),
    denominator_df = c("nonlinear", "all"),
    root_tol = .Machine$double.eps^0.25, froot_tol = sqrt(root_tol),
    ncores = detectCores() - 1L, verbose = 0, ...)

}

% arguments

\arguments{
  \item{object}{an object of class \code{fit_wrc_hcc}, see
  \code{\link{fit_wrc_hcc}}.}

  \item{values}{a \code{data.frame} or a \code{matrix} with the values of
  the conditionally linear (\mbvec{\mu}{}{}) and nonlinear parameters
  (\mbvec{\nu}{}{}) that should be kept fixed to compute the likelihood
  profile (mandatory argument, see \code{\link{soilhypfitIntro}},
  \code{\link{wc_model}} and \code{\link{hc_model}} for information about
  the parametrization of models for soil water retention curves and/or soil
  hydraulic conductivity functions.).  The names of the columns of
  \code{values} must match the names of model parameters.}

  \item{soil_sample}{a character scalar with the label of the soil sample
  for which the loglikelihood profile or the confidence interval should be
  computed.  If \code{object} contains parameter estimates for a single
  soil sample then \code{soil_sample} is ignored, otherwise
  \code{soil_sample} is a mandatory argument.}

  \item{ncores}{an integer defining the number of cores for parallel
  computations.  \code{ncores = 1} suppresses parallel computations.}

  \item{verbose}{positive integer controlling logging of diagnostic
  messages to the console and plotting of data and model curves
  during fitting, see \code{\link{fit_wrc_hcc}}.}

  \item{parm}{character scalar (\code{confint_prfloglik_sample}) or vector
  (\code{confint}) with name(s) of parameter(s) for which to compute the
  confidence interval(s).  Note that \code{confint_prfloglik_sample} allows
  to compute a confidence interval for \emph{all} parameters (including the
  linear ones), whereas the \code{confint} method computes confidence
  intervals for the \emph{nonlinear} parameters (\mbvec{\nu}{}{}) only,
  see \code{\link{soilhypfitIntro}} for information about the
  parametrization of models for soil water retention curves and/or soil
  hydraulic conductivity functions.}

  \item{level}{numeric scalar with the confidence level
  required to compute the confidence interval.}

  \item{test}{character keyword specifying whether to use the asymptotic
  \eqn{\chi^2}-distribution or the finite sample approximate F-distribution
  for the likelihood ratio test statistic when computing the confidence
  interval, see \emph{Details}.}

  \item{denominator_df}{character keyword specifying whether the
  denominator degrees of freedom for the F-distribution of the test
  statistic is equal to the number of estimated \emph{nonlinear} parameters
  (\code{"nonlinear"}, default) or equal to the total number of estimated
  parameters (\code{"all"}).}

  \item{param_bound}{a numeric vector of length 2 with the allowed range of
  the parameter for which the confidence interval is computed.  The limits
  of the confidence interval are searched within this range, see
  \emph{Details}.  When equal to \code{NULL} (default) \code{param_bound}
  is taken from to component \code{initial_objects} of the selected
  component \code{fit} of \code{object}, see \emph{Details}.}

  \item{root_tol}{a numeric scalar defining the desired accuracy (convergence
  tolerance) for root finding by \code{\link[stats]{uniroot}} when
  computing the confidence interval.}

  \item{froot_tol}{a numeric scalar defining the desired accuracy (function
  value tolerance) for deciding whether a root has been found.}

  \item{subset}{an integer, character or logical vector to the choose the
  soil samples for which confidence intervals should be computed.  Defaults
  to \code{NULL} which computes the intervals for all samples contained in
  \code{object}.}

  \item{type}{character keyword specifying whether to compute confidence
  intervals based on the likelihood ratio test (\code{"logik"}, default) by
  \code{confint_prfloglik_sample} or based on the asymptotic normal
  distribution of maximum likelihood estimates (\code{"normal"}), see
  \emph{Details}.}

  \item{...}{additional arguments passed to methods, currently not used.}

}


% details

\details{

  \subsection{Computing likelihood profiles}{

    The function \code{prfloglik_sample} computes the loglikelihood profile
    for a subset of the nonlinear (\mbvec{\nu}{}{}) (and linear
    \mbvec{\mu}{}{}) model parameters.  We denote the profiling
    parameters of interest by \eqn{\boldsymbol{\phi}}{\phi}
    and the nuisance parameters that are estimated when computing the
    profile by \eqn{\boldsymbol{\psi}}{\psi}, so that
    \eqn{
      (
        \boldsymbol{\phi}^\mathrm{T},
        \boldsymbol{\psi}^\mathrm{T}
      )^\mathrm{T}
    }{
      (\phi^T, \psi^T)^T
    }
    is a rearranged version of the parameter vector
      \eqn{
      (
        \boldsymbol{\mu}^\mathrm{T},
        \boldsymbol{\nu}^\mathrm{T}
      )^\mathrm{T}
    }{
      (\mu^T, \nu^T)^T
    }, see \code{\link{soilhypfitIntro}}.
    \code{prfloglik_sample} computes the estimates
    \eqn{
      \widehat{\boldsymbol{\psi}}(\boldsymbol{\phi})
    }{
       hat\psi(\phi)
    }
    of \eqn{\boldsymbol{\psi}}{\psi}
    and the profile loglikelihood
    \eqn{
      Q(\boldsymbol{\phi},
        \widehat{\boldsymbol{\psi}}(\boldsymbol{\phi});
        \boldsymbol{\theta}, \boldsymbol{K},
        \boldsymbol{h}
      )
    }{
      Q(\phi, hat\psi(\phi); \theta, K, h)
    }
    as a function of \eqn{\boldsymbol{\phi}}{\phi}.

    To compute the estimates,
    \eqn{\boldsymbol{\psi}}{\psi}
    is partitioned into the nonlinear and conditionally linear parameters, see
    \code{\link{soilhypfitIntro}} for details.
    \code{prfloglik_sample} uses the packages \pkg{parallel} and
    \pkg{snowfall} for parallelized computation of loglikelihood
    profiles.

  }

  \subsection{Computing likelihood based confidence intervals}{

    The function \code{confint_prfloglik_sample} computes the
    confidence interval of a single model parameter
    \eqn{
      \phi \in (
        \boldsymbol{\mu}^\mathrm{T},
        \boldsymbol{\nu}^\mathrm{T}
      )^\mathrm{T}
    }{
       \phi \in (\mu^T, \nu^T)^T
    }
    based on the likelihood ratio test.
    The interval is computed by assuming either

    \itemize{
      \item  that the test statistic
      \eqn{
        T(\phi) = \Delta Q(\phi) = 2(
        Q(\widehat{\phi},
        \widehat{\boldsymbol{\psi}};
        \boldsymbol{\theta}, \boldsymbol{K},
        \boldsymbol{h}
        ) -
        Q(\phi,
        \widehat{\boldsymbol{\psi}}(\phi);
        \boldsymbol{\theta}, \boldsymbol{K},
        \boldsymbol{h}
      ))
      }{
        T(\phi) = \Delta Q(\phi) = 2(Q(hat\phi, hat\psi; \theta, K, h) - Q(\phi, hat\psi(\phi); \theta, K, h))
      }
      follows a \eqn{\chi^2_1}-distribution with 1 degrees of freedom
      (possible choice when both water retention and/or
      hydraulic conductivity data were used to estimate the parameters), or

      \item that the transformed test statistic
      \eqn{
        T(\phi) = \left(\exp(\Delta Q(\phi) / n) - 1 \right) (n - p)
      }{
        T(\phi) = (exp(\Delta Q(\phi) / n) - 1) (n - p)
      }
      follows an \eqn{F(1,n-p)}-distribution where \eqn{n} is the number of
      water content or hydraulic conductivity measurements, and \eqn{p} is
      the number of estimated parameters (see Uusipaikka,
      2008, p. 115, for details).  \code{denominator_df} allows to control
      how \eqn{p} is chosen. Note that this test distribution can
      only be chosen (and is then the default) when only water retention or
      only hydraulic conductivty data were used to estimate the parameters.

    }

    \code{confint_prfloglik_sample} computes profile loglikelihoods
    \eqn{
      Q(\phi,
        \widehat{\boldsymbol{\psi}}(\phi);
        \boldsymbol{\theta}, \boldsymbol{K},
        \boldsymbol{h}
      )
    }{
      Q(\phi, hat\psi(\phi); \theta, K, h)
    }
    by\cr \code{prfloglik_sample}
    and then uses \code{\link[stats]{uniroot}} to search the roots
    of the equation
    \deqn{
      f(\phi) = q_T(\gamma)- T(\phi)
    }
    in the interval defined by \code{param_bound}. \eqn{q_T(\gamma)} is
    the \eqn{\gamma}-quantile of the chosen test distribution for \eqn{T}.
  }

  \subsection{Computing confidence intervals for several soil samples}{

    The \code{confint} method computes 2 types of confidence intervals (in
    dependence of \code{type}) of only the nonlinear parameters
    \mbvec{\nu}{}{}, possibly for
    multiple soil samples at the same time:

    \itemize{

      \item intervals based on the asymptotic normal distribution
      of maximum likelihood estimates with standard errors computed by the
      \code{vcov} method for class \code{fit_wrc_hcc} (see
      \code{\link{vcov.fit_wrc_hcc}},

      \item intervals based on the likelihood ratio test by
      \code{confint_prfloglik_sample}.  The intervals for several soil
      samples are computed in parallel.

    }

  }

  \subsection{Requirements for computing loglikelihood profiles}{

    The parameters contained in \code{object} must be estimated by maximum
    likelihood (\code{method = "ml"}, see \code{\link{soilhypfitIntro}} and
    \code{\link{control_fit_wrc_hcc}}.  Use of other estimation methods
    results an error.

    Preferably an unconstrained local algorithm (\code{settings =
    "ulocal"}, see \code{\link{soilhypfitIntro}} and
    \code{\link{control_fit_wrc_hcc}})) is used for minimizing the negative
    loglikelihood when generating \code{object}.  Use of other algorithms
    generates a warning message.

  }

}

% value

\value{

  The function \code{prfloglik_sample} returns a \code{data.frame}
  with the columns of \code{values} (parameters of interest
  \eqn{\boldsymbol{\phi}}{\phi}),
  a column \code{loglik} with the
  maximized profile loglikelihood
  \eqn{Q(\phi,
        \widehat{\boldsymbol{\psi}}(\phi);
        \boldsymbol{\theta}, \boldsymbol{K},
        \boldsymbol{h}
      )
  }{
    Q(\phi, hat\psi(\phi); \theta, K, h)
  }, columns with the estimated parameters
  \eqn{
      \widehat{\boldsymbol{\psi}}(\boldsymbol{\phi})
    }{
       hat\psi(\phi)^T
    },
  columns with the gradient of the loglikelihood with respect to the
  estimated nonlinear parameters (missing if all nonlinear parameters were
  fixed) and a column \code{converged}, indicating whether convergence has
  occurred (see \code{\link{convergence_message}}) when estimating the
  nonlinear parameters (equal to \code{NA} when all nonlinear parameters
  are fixed).

  The function \code{confint_prfloglik_sample} returns a numeric
  vector of length 2 with the lower and upper limits of the confidence
  interval.  If no roots were found then \code{NA} is returned.  The
  returned result further contains as attribute list \code{prfloglik}
  the parameter estimate \eqn{\widehat{\phi}}{hat\phi}
  (\code{param_estimate}), the maximized loglikelihood
  \eqn{
    Q(\widehat{\phi},
        \widehat{\boldsymbol{\psi}};
        \boldsymbol{\theta}, \boldsymbol{K},
        \boldsymbol{h}
    )
  }{
    Q(hat\phi, hat\psi; \theta, K, h)
  } (\code{loglik}), the quantile of the test distribution
  \eqn{q_{\mathrm{test}}(\gamma)}{_{test}(\gamma)} (\code{qtest}),
  the type of test distribution used (\code{test}),
  the significance \code{level}, the number
  of water content (\code{nobs_wrc}) and conductivity measurements
  (\code{nobs_hcc}) and the function values of \eqn{f(\phi)}
  evaluated at the roots returned by \code{\link[stats]{uniroot}}.

  The method \code{confint} returns a dataframe with the lower and upper
  limits of the confidence intervals for the estimated nonlinear
  parameters.


}

% references

\references{
  Uusipaikka, E. (2008) Confidence Intervals in Generalized Regression Models.
  Chapman & Hall/CRC Press, Boca Raton \doi{10.1201/9781420060386}.
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

%   \code{\link{prfloglik_sample}} for profile loglikelihood
%   computations;
%
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

if(!requireNamespace("lattice", quietly = TRUE)){
  install.packages("lattice", repos = "https://cloud.r-project.org/")
  requireNamespace("lattice", quietly = TRUE)
}

data(sim_wrc_hcc)

# define number of cores for parallel computations
if(interactive()) ncpu <- parallel::detectCores() - 1L else ncpu <- 1L

# estimate parameters for 3 samples simultaneously by ...

# ... unconstrained, global optimisation algorithm NLOPT_GN_MLSL (default)
rfit_uglob <- fit_wrc_hcc(
  wrc_formula = wc ~ head | id, hcc_formula = hc ~ head | id,
  data = sim_wrc_hcc,
  control = control_fit_wrc_hcc(param_bound = param_boundf(
      alpha = c(0.00001, 50), n = c(1.0001, 7), tau = c(-1, 5)
    ), pcmp = control_pcmp(ncores = ncpu)))

# ... unconstrained, local optimisation algorithm NLOPT_LN_BOBYQA,
rfit_uloc <- update(
  rfit_uglob,
  param = as.matrix(coef(rfit_uglob, what = "nonlinear")),
  control = control_fit_wrc_hcc(
    settings = "ulocal", param_tf = param_transf(
      alpha = "identity", n = "identity", tau = "identity"
    ), param_bound = param_boundf(
      alpha = c(0.00001, 50), n = c(1.0001, 7), tau = c(-1, 5)
    ), pcmp = control_pcmp(ncores = ncpu)))

# extract estimated parameters for sample id == "1"
coef_id_1 <- unlist(coef(rfit_uloc, gof = TRUE, se = TRUE, subset = "1"))

# compute loglikelihood profile of parameter alpha for sample id == "1"
rprfloglik_alpha <- prfloglik_sample(
  rfit_uloc, values = data.frame(alpha = seq(1.5, 3.0, length = 40L)),
  soil_sample = "1", ncores = ncpu)
# plot loglikelihood profile along with 95\% confidence intervals
plot(loglik ~ alpha, rprfloglik_alpha, type = "l")
abline(v = coef_id_1["alpha"])
# 95\% confidence intervall based on likelihood ratio test
abline(h = -coef_id_1["obj"] - 0.5 * qchisq(0.95, 1), lty = "dashed")
# 95\% confidence intervall based on asymptotic normal distribution
segments(
  x0 = coef_id_1["alpha"] + qnorm(0.025) * coef_id_1["se.alpha"],
  x1 = coef_id_1["alpha"] + qnorm(0.975) * coef_id_1["se.alpha"],
  y0 = min(rprfloglik_alpha$loglik)
)

# compute loglikelihood profile of parameter n for sample id == "1"
rprfloglik_n <- prfloglik_sample(
  rfit_uloc, values = data.frame(n = seq(1.7, 2.25, length = 40L)),
  soil_sample = "1", ncores = ncpu
)
# plot loglikelihood profile along with 95\% confidence intervals
plot(loglik ~ n, rprfloglik_n, type = "l")
abline(v = coef_id_1["n"])
# 95\% confidence intervall based on likelihood ratio test
abline(h = -coef_id_1["obj"] - 0.5 * qchisq(0.95, 1), lty = "dashed")
# 95\% confidence intervall based on asymptotic normal distribution
segments(
  x0 = coef_id_1["n"] + qnorm(0.025) * coef_id_1["se.n"],
  x1 = coef_id_1["n"] + qnorm(0.975) * coef_id_1["se.n"],
  y0 = min(rprfloglik_n$loglik)
)

# compute loglikelihood profile of parameters alpha and n for sample id == "1"
rprfloglik_alpha_n <- prfloglik_sample(
  rfit_uloc, values = expand.grid(
    alpha = seq(1.5, 3.0, length = 40L), n = seq(1.7, 2.25, length = 40L)),
  soil_sample = "1", ncores = ncpu
)

# joint confidence region of alpha and n based on likelihood ratio test
lattice::levelplot(loglik ~ alpha + n, rprfloglik_alpha_n,
  panel = function(x, y, z, subscripts, ...){
    lattice::panel.levelplot(x, y, z, subscripts, ...)
    lattice::panel.levelplot(x, y, z, subscripts, region = FALSE, contour = TRUE,
      at = -coef_id_1["obj"] - 0.5 * qchisq(0.95, 2),
      lty = "solid"
    )
    lattice::panel.levelplot(x, y, z, subscripts, region = FALSE, contour = TRUE,
      at = -coef_id_1["obj"] - 0.5 * qchisq(0.9, 2),
      lty = "dashed"
    )
    lattice::panel.levelplot(x, y, z, subscripts, region = FALSE, contour = TRUE,
      at = -coef_id_1["obj"] - 0.5 * qchisq(0.8, 2),
      lty = "dotted"
    )
    lattice::panel.points(coef_id_1["alpha"], coef_id_1["n"], pch = 16)
    lattice::panel.lines(
      x = rprfloglik_alpha[, c("alpha", "n")], col = "blue"
    )
    lattice::panel.lines(
      x = rprfloglik_n[, c("alpha", "n")], col = "magenta"
    )
  }, key = list(
    corner = c(1, 0.97),
    lines = list(
      lty = c(rep("solid", 3), "dashed", "dotted"),
      col = c("blue", "magenta", rep("black", 3))
    ),
    text = list(c(
        "estimate of n as a function of fixed alpha",
        "estimate of alpha as a function of fixed n",
        "95\% joint confidence region",
        "90\% joint confidence region",
        "80\% joint confidence region"
      ))
  ))

# compute 95\%-confidence interval
(ci.alpha <- confint_prfloglik_sample(
  rfit_uloc, parm = "alpha", soil_sample = "1"
))
# use limits to draw loglikelihood profile for alpha
rprfloglik_alpha <- prfloglik_sample(
  rfit_uloc, values = data.frame(
    alpha = seq(0.9 * ci.alpha[1], 1.1 * ci.alpha[2], length = 40L)),
  soil_sample = "1", ncores = ncpu)
plot(loglik ~ alpha, rprfloglik_alpha, type = "l")
lines(
  ci.alpha,
  rep(diff(unlist(attr(ci.alpha, "prfloglik")[c("qtest", "loglik")])) , 2)
)
abline(v = attr(ci.alpha, "prfloglik")["estimate"], lty = "dashed")

# 95\%-confidence intervals of all nonlinear parameters based for all
# soil samples asymptotic normal distribution of maximum likelihood estimates
confint(rfit_uloc, type = "normal")

# 95\%-confidence intervals of all nonlinear parameters based for all
# soil samples based on likelihood ratio test
confint(rfit_uloc, ncores = ncpu)

}}
