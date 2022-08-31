% 2022-08-29 A. Papritz
% cd ~/R.user.home/soilhypfit/soilhypfit_01-7/pkg/man
% R CMD Rdconv -t html -o bla.html soilhypfit_hcc_model.Rd ; open bla.html; R CMD Rd2pdf --force soilhypfit_hcc_model.Rd;
% hunspell -d en_GB-ise hcc_model.Rd

\encoding{utf8}

\name{hc_model}
\alias{hc_model}
\alias{hcrel_model}
\alias{hc_model}

\concept{soil water}
\concept{hydraulic conductivity function}
\concept{saturated hydraulic conductivity}

% \newcommand{\mbvec}{\eqn{\boldsymbol{#1}^\mathrm{#2} #3}}

\title{
Models for Soil Hydraulic Conductivity Functions
}

% description

\description{
  The functions \code{hc_model} and \code{hcrel_model} compute, for given
  capillary pressure head \eqn{h}, the \emph{hydraulic conductivity}
  \eqn{K(h)} and the \emph{relative hydraulic conductivity} \eqn{k(h)}
  respectively, of a soil by parametrical models.  }

% usage

\usage{
hcrel_model(h, nlp, precBits = NULL, hcc_model = "vgm")

hc_model(h, nlp, lp, precBits = NULL, hcc_model = "vgm")
}

% arguments

\arguments{
  \item{h}{a mandatory numeric vector with values of capillary pressure
  head for which to compute the hydraulic conductivity.  For consistency
  with other quantities, the unit of  head should be \bold{meter}
  [m].}

  \item{nlp}{a mandatory named numeric vector, currently with elements
  named \code{"alpha"}, \code{"n"} and \code{"tau"}, which are the
  \emph{nonlinear} parameters \mbvec{\nu}{T}{= (\alpha, n, \tau)}, where
  \eqn{\alpha}, \eqn{n} and \eqn{\tau} are the inverse air entry pressure,
  the shape and the tortuosity parameters, see \emph{Details}.  For
  consistency with other quantities, the unit of \eqn{\alpha} should be
  \bold{1/meter} [\eqn{\mathrm{m}^{-1}}{m^-1}].}

  \item{lp}{a mandatory named numeric vector, currently with a single
  element named \code{"k0"}, which is the saturated hydraulic conductivity
  \eqn{K_0}, the only \emph{linear} parameter of the model, see
  \emph{Details}.  For consistency with other quantities, the unit of
  \eqn{K_0} should be \bold{meter/day} [\eqn{\mathrm{m}\,\mathrm{d}^{-1}}{m
  d^-1}].}

  \item{precBits}{an optional integer scalar defining the maximal precision
  (in bits) to be used in high-precision computations by
  \code{\link[Rmpfr]{mpfr}}.  If equal to \code{NULL} (default) then
  \code{\link[Rmpfr]{mpfr}} is not used and the result of the function call
  is of storage mode \code{double}, see \code{\link{soilhypfitIntro}}.}

  \item{hcc_model}{a keyword denoting the parametrical model for the
  hydraulic conductivity function.  Currently, only the \emph{Van
  Genuchten-Mualem} model (\code{wrc_model = "vgm"}) is implemented, see
  \emph{Details}.}
}

% details

\details{
The functions \code{hcrel_model} and \code{hc_model} currently model soil
hydraulic conductivity functions only by the simplified form of the Van
Genuchten-Mualem model (\cite{Van Genuchten, 1980}) with the restriction
\eqn{m = 1 - \frac{1}{n}}{m=1 - 1/n}, i.e. by
%
\deqn{
  k_\mathrm{VGM}(h; \boldsymbol{\nu}) =
    S_\mathrm{VG}(h; \boldsymbol{\nu})^\tau \,
    \left[
      1 - \left(
        1 - S_\mathrm{VG}(h; \boldsymbol{\nu})^\frac{n}{n-1}
      \right)^\frac{n-1}{n}
    \right]^2,
}{k_VGM(h; \nu) =
  S_VG(h; \nu)^\tau * (1 - (1 - S_VG(h; \nu))^(n/(n-1)))^((n-1)/n),
}
%
\deqn{
  K_\mathrm{VGM}(h; \mbox{$\mu$}, \boldsymbol{\nu}) =
  K_0 \, k_\mathrm{VGM}(h; \boldsymbol{\nu}),
}{
  K_VGM(h; \mu, \nu) = K_0 * k_VGM(h; \nu),
}
%
where \eqn{\mu = K_0} is the saturated hydraulic conductivity (\eqn{K_0 >
0}),
\mbvec{\nu}{T}{= (\alpha, n, \tau)} are the inverse air entry pressure
(\eqn{\alpha > 0}), the shape (\eqn{n > 1}) and tortuosity parameter
(\eqn{\tau > -2}), and \eqn{S_\mathrm{VG}(h;
\boldsymbol{\nu})}{S_VG(h; \nu)} is the the volumetric
water saturation, see \code{sat_model} for an expression.

Note that \eqn{\mu} and \mbvec{\nu}{}{} are passed to the functions by
the arguments \code{lp} and \code{nlp}, respectively.}

% value

\value{
  A numeric vector with values of (relative) hydraulic conductivity.
}

% references

\references{
  Van Genuchten, M. Th.  (1980) A closed-form equation for predicting the
  hydraulic conductivity of unsaturated soils. \emph{Soil Science Society of
  America Journal}, \bold{44}, 892--898,
  \doi{10.2136/sssaj1980.03615995004400050002x}.
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

%   \code{\link{wc_model}} and \code{\link{hc_model}} for currently
%   implemented models for soil water retention curves and hydraulic
%   conductivity functions;
%
  \code{\link{evaporative-length}} for physically constraining parameter
  estimates of soil hydraulic material functions.
}

% examples

\examples{
## define capillary pressure head (unit meters)
h <- c(0.01, 0.1, 0.2, 0.3, 0.5, 1., 2., 5.,10.)

## compute (relative) hydraulic conductivity
hcrel <- hcrel_model(h, nlp = c(alpha = 1.5, n = 2, tau = 0.5))
hc <- hc_model(h, nlp = c(alpha = 1.5, n = 2, tau = 0.5), lp = c(k0 = 5))

## display hydraulic conductivity function
op <- par(mfrow = c(1, 2))
plot(hcrel ~ h, log = "xy", type = "l")
plot(hc ~ h, log = "xy", type = "l")
on.exit(par(op))
}
