% 2022-08-29 A. Papritz
% cd ~/R.user.home/soilhypfit/soilhypfit_01-7/pkg/man
% R CMD Rdconv -t html -o bla.html soilhypfit_S3methods.Rd ; open bla.html; R CMD Rd2pdf --force soilhypfit_S3methods.Rd;
% hunspell -d en_GB-ise soilhypfitS3methods.Rd

\encoding{utf8}

\name{soilhypfitS3methods}
\alias{soilhypfitS3methods}
\alias{soilhypfitmethods}
\alias{coef.fit_wrc_hcc}
\alias{lines.fit_wrc_hcc}
\alias{plot.fit_wrc_hcc}
\alias{print.fit_wrc_hcc}
\alias{print.summary_fit_wrc_hcc}
\alias{summary.fit_wrc_hcc}
\concept{soil water}
\concept{soil water characteristic}
\concept{soil water retention curve}
\concept{hydraulic conductivity function}
\concept{saturated hydraulic conductivity}

\title{Common S3 Methods for Class \code{fit_wrc_hcc}}

% description

\description{This page documents the methods \code{coef},
\code{summary}, \code{print}, \code{plot} and \code{lines} for the class
\code{fit_wrc_hcc}.}

% usage

\usage{

\method{coef}{fit_wrc_hcc}(object, what = c("all", "nonlinear", "linear"),
    subset = NULL, residual_se = FALSE, se = FALSE, gof = FALSE, lc = FALSE,
    e0 = FALSE, bound = lc, ...)

\method{summary}{fit_wrc_hcc}(object, what = c("all", "nonlinear", "linear"),
    subset = NULL, gof = TRUE, lc = TRUE, ...)

\method{print}{fit_wrc_hcc}(x, ...)

\method{plot}{fit_wrc_hcc}(x, what = c("wrc", "hcc"), y = NULL,
    subset = NULL, ylim_wc = NULL, ylim_hc = NULL,
    head_saturation = 0.01,
    beside = identical(sum(par("mfrow")), 2L), pch = 1, col_points = "black",
    col_line_x = "blue", lty_x = "solid",
    col_line_y = "orange", lty_y = "dashed",
    xlab_wc = "head [m]", ylab_wc = "water content [-]",
    xlab_hc = "head [m]", ylab_hc = "hyd. conductivity [m/d]",
    draw_legend = TRUE, draw_parameter = FALSE, cex_legend = 0.7, ...)

\method{lines}{fit_wrc_hcc}(x, what = c("wrc", "hcc"), id = 1,
    head_saturation = 0.01, ...)

}

% arguments

\arguments{

  \item{object, x, y}{an object of class \code{fit_wrc_hcc}, see
  \code{\link{fit_wrc_hcc}}.}

  \item{what}{character keyword indicating the type of parameters to return
  (\code{coef}) or the type of data to \code{plot}.}

  \item{subset}{an integer, character or logical vector to the choose the
  soil samples for which data and model curves are displayed or extracted.
  Defaults to \code{NULL} which displays results for all soil samples.}

  \item{residual_se}{a logical scalar to control whether residual standard
  errors (= standard deviations of residuals) should be returned, see
  \emph{Details}.}

  \item{se}{a logical scalar to control whether standard errors of the
  nonlinear parameters \mbvec{\nu}{}{} should be returned, see
  \emph{Details} and \code{\link{vcov.fit_wrc_hcc}}.}

  \item{gof}{a logical scalar to control whether goodness-of-fit statistic
  should be returned.}

  \item{lc}{a logical scalar to control whether the characteristic evaporative
  length should be returned, see \code{\link{evaporative-length}}.}

  \item{e0}{a logical scalar to control whether the evaporation rate should
  be returned.  This is only effective for constrained estimation, see
  \code{\link{evaporative-length}}.}

  \item{bound}{a logical scalar to control whether the lower and upper
  bounds of the ratio \eqn{L_\mathrm{c}/L_\mathrm{t}}{L_c/L_t} should be
  returned. This is only effective for constrained estimation, see
  \code{\link{evaporative-length}}.}

  \item{ylim_wc}{optional numeric vector of length 2 to set the range of
  water content values displayed on the y-axis (default \code{NULL} for
  automatic axis scales).}

  \item{ylim_hc}{optional numeric vector of length 2 to set the range of
  hydraulic conductivity values displayed on the y-axis (default
  \code{NULL} for automatic axis scales).}

  \item{head_saturation}{head value (unit m) assigned to zero head values in
  plots with logarithmic head scale.}

  \item{beside}{a logical scalar controlling whether water retention curves
  and hydraulic conductivity functions of a sample should be plotted side
  by side.}

  \item{pch}{plotting `character', i.e., symbol to use for the
  measurements, see \code{\link[graphics]{points}}.}

  \item{col_points}{color code or name for symbol colors for the
  measurements, see \code{\link[graphics]{par}}.}

  \item{col_line_x}{color code or name for the line representing the
  fitted model \code{x}, see \code{\link[graphics]{par}}.}

  \item{lty_x}{type of line representing the fitted model \code{x}, see
  \code{\link[graphics]{par}}.}

  \item{col_line_y}{color code or name for the line representing the
  fitted model \code{y}, see \code{\link[graphics]{par}}.}

  \item{lty_y}{type of line representing the fitted model \code{y}, see
  \code{\link[graphics]{par}}.}

  \item{xlab_wc}{a character string with the annotation for the
  \var{x}-axis of a water retention curve.}

  \item{ylab_wc}{a character string with the annotation for the
  \var{y}-axis of a water retention curve.}

  \item{xlab_hc}{a character string with the annotation for the
  \var{x}-axis of a hydraulic conductivity function.}

  \item{ylab_hc}{a character string with the annotation for the
  \var{y}-axis of a hydraulic conductivity function.}

  \item{draw_legend}{a logical scalar controlling whether a legend with the
  values of the arguments \code{x} and \code{y} and the residual sums of
  squares is drawn if \code{y} is non-\code{NULL}.}

  \item{draw_parameter}{a logical scalar controlling whether the
  parameters are drawn (default \code{FALSE}).}

  \item{cex_legend}{a character expansion factor for annotations by
  \code{\link[graphics]{legend}}.}

  \item{id}{a character string or integer scalar to select the sample for
  which to plot the modelled water retention curve or hydraulic
  conductivity function.}

  \item{...}{additional arguments passed to methods.}
}

% details

\details{

  Residual standard errors, standard errors of the nonlinear parameters and
  confidence intervals based on the asymptotic normal distribution are
  computed only for mpd and ml estimates, see \cr
  \code{\link{soilhypfitIntro}}, \code{\link{control_fit_wrc_hcc}} and
  \code{\link{vcov}}.

  The \code{plot} method for class \code{fit_wrc_hcc} displays for each
  sample the measurements of the water retention curve and/or the hydraulic
  conductivity function, along with the fitted model curve(s).  Optionally,
  the curves of a second model fit (specified by the
  argument \code{y}) can be plotted for each sample.

  The \code{lines} method adds the curve of a fitted model to an existing
  plot.

}

% value

\value{

  The method \code{coef} returns a dataframe with the estimated parameters
  (and optionally standard errors), optionally the value of the
  objective function along with convergence code and/or information on the
  characteristic evaporative length.

  The method \code{summary} generates a list (of class
  \code{summary.fit_wrc_hcc}) with the following components:
  %
  \describe{

    \item{\code{data}}{a named integer vector with the total number of
    samples (\code{nsamp}) and the number of samples with only water
    retention (\code{nwrc}), only hydraulic conductivity (\code{nhcc}) and
    both type of measurements (\code{nwrchcc}).}

    \item{\code{control}}{a list with a subset (\code{settings},
    \code{nloptr}, \code{sce}, \code{approximation_alpha_k0},
    \code{param_bound}, \code{param_tf}) of the components of
    \code{object[["control"]]}, see \code{\link{fit_wrc_hcc}} and
    \code{\link{control_fit_wrc_hcc}}.}

    \item{\code{result}}{a dataframe with the estimated parameters and
    optionally the residual sum of squares along with convergence code
    and/or information on the characteristic evaporative length.}

    \item{\code{call}}{the \code{call} component of \code{object}.}

  }

  Note that only a \code{print} method is available for class
  \code{summary.fit_wrc_hcc}.

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

%   \code{\link{soilhypfitmethods}} for common S3 methods for class
%   \code{fit_wrc_hcc};
%
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
  wrc_formula = wc ~ head | id, hcc_formula = hc ~ head | id,
  wrc_subset = id != 1, hcc_subset = id != 2,
  data = sim_wrc_hcc, fit_param = default_fit_param(tau = TRUE),
  control = control_fit_wrc_hcc(param_bound = param_boundf(
      alpha = c(0.00001, 50), n = c(1.0001, 7), tau = c(-1, 5)
    ), pcmp = control_pcmp(ncores = ncpu)))
print(rfit_uglob)
summary(rfit_uglob)
coef(rfit_uglob, what = "nonlinear")
coef(rfit_uglob, what = "linear", gof = TRUE)
coef(vcov(rfit_uglob), status = TRUE, se = FALSE)
op <- par(mfrow = c(3, 2))
plot(rfit_uglob)
on.exit(par(op))

}
}

\keyword{methods}
\keyword{print}
