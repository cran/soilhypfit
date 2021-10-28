% 2020-06-10 A. Papritz
% cd R.user.home/soilhypfit/soilhypfit_01-1/pkg/man
% R CMD Rdconv -t html -o bla.html soilhypfit_wc_model.Rd ; open bla.html; R CMD Rd2pdf --force soilhypfit_wc_model.Rd;
% hunspell -d en_GB-ise wc_model.Rd

\encoding{utf8}
\name{wc_model}
\alias{wc_model}
\alias{sat_model}
\alias{wc_model}

\concept{soil water}
\concept{soil water characteristic}
\concept{soil water retention curve}

% \newcommand{\mbvec}{\eqn{\mbox{\boldmath$#1$\unboldmath}^\mathrm{#2} #4}{#1#3#4}}

\title{
Models for Soil Water Retention Curves
}

% description

\description{
	The functions \code{sat_model} and \code{wc_model} compute, for given
	capillary pressure head \eqn{h}, the volumetric water saturation
	\eqn{S(h)} and the volumetric water content \eqn{\theta(h)},
	respectively, of a soil by parametrical models.  }

% usage

\usage{
sat_model(h, nlp, precBits = NULL, wrc_model = "vg")

wc_model(h, nlp, lp, precBits = NULL, wrc_model = "vg")

}

% arguments

\arguments{
	\item{h}{a mandatory numeric vector with values of capillary pressure
	head for which to compute the volumetric water saturation or content.
	For consistency with other quantities, the unit of pressure head
	should be \bold{meter} [m].}

	\item{nlp}{a mandatory named numeric vector, currently with elements
	named \code{"alpha"} and \code{"n"}, which are the \emph{nonlinear}
	parameters \mbvec{\nu}{T}{^T}{= (\alpha, n)}, where \eqn{\alpha} and \eqn{n}
	are the inverse air entry pressure and the shape parameter, see
	\emph{Details}.  For consistency with other quantities, the unit of
	\eqn{\alpha} should be \bold{1/meter} [\eqn{\mathrm{m}^{-1}}{m^-1}].}

	\item{lp}{a mandatory named numeric vector, currently with elements named
	\code{"thetar"} and \code{"thetas"}, which are the \emph{linear}
	parameters \mbvec{\mu}{T}{^T}{= (\theta_r, \theta_s)} where
	\eqn{\theta_r} and \eqn{\theta_s} are the residual and saturated water
	content, respectively, see \emph{Details}.}

	\item{precBits}{an optional integer scalar defining the maximal precision
	(in bits) to be used in high-precision computations by
	\code{\link[Rmpfr]{mpfr}}.  If equal to \code{NULL} (default) then
	\code{\link[Rmpfr]{mpfr}} is not used and the result of the function call
	is of storage mode \code{double}, see \code{\link{soilhypfitIntro}}.}

	\item{wrc_model}{a keyword denoting the parametrical model for the
	water retention curve.  Currently, only the \emph{Van Genuchten model}
	(\code{wrc_model = "vg"}) is implemented, see \emph{Details}.} }

% details

\details{
The functions \code{sat_model} and \code{wc_model} currently model soil
water retention curves only by the simplified form of the model by
\cite{Van Genuchten (1980)} with the restriction \eqn{m = 1 - \frac{1}{n}}{m=1
- 1/n}, i.e.
%
\deqn{
  S_\mathrm{VG}(h; \mbox{\boldmath$\nu$\unboldmath}) = (
  1 + (\alpha \ h)^n)^\frac{1-n}{n}
}{S_VG(h; \nu) = (1 + (\alpha * h)^n)^((1-n)/n)}
%
\deqn{
  \theta_\mathrm{VG}(h; \mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath}) =
	\theta_r + (\theta_s - \theta_r) \, S_\mathrm{VG}(h; \mbox{\boldmath$\nu$\unboldmath})
}{\theta_VG(h; \mu, \nu) = \theta_r + (\theta_s - \theta_r) * S_VG(h; \nu)}
%
where \mbvec{\mu}{T}{^T}{= (\theta_r, \theta_s)} are the residual and
saturated water content (\eqn{0 \le \theta_r \le \theta_s \le 1}),
respectively, and \mbvec{\nu}{T}{^T}{= (\alpha, n)} are the inverse air
entry pressure (\eqn{\alpha > 0}) and the shape parameter (\eqn{n > 1}).

Note that \mbvec{\mu}{}{}{} and
\mbvec{\nu}{}{}{} are passed to the functions by
the arguments \code{lp} and \code{nlp}, respectively.}

% value

\value{
	A numeric vector with values of volumetric water saturation
	(\code{sat_model}) or water content (\code{wc_model}).
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

## compute water saturation and water content
sat <- sat_model(h, nlp = c(alpha = 1.5, n = 2))
theta <- wc_model(
  h ,nlp = c(alpha = 1.5, n = 2), lp = c(thetar = 0.1, thetas = 0.5))

## display water retention curve
op <- par(mfrow = c(1, 2))
plot(sat ~ h, log = "x", type = "l")
plot(theta ~ h, log = "x", type = "l")
on.exit(par(op))
}
