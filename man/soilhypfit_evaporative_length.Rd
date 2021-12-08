% 2021-10-22 A. Papritz
% cd R.user.home/soilhypfit/soilhypfit_01-1/pkg/man
% R CMD Rdconv -t html -o bla.html soilhypfit_evaporative_length.Rd ; open bla.html; R CMD Rd2pdf --force soilhypfit_evaporative_length.Rd;
% hunspell -d en_GB-ise evaporative_length.Rd

\encoding{utf8}

\name{evaporative-length}
\alias{evaporative-length}
\alias{lc}
\alias{lt}

\concept{soil water}
\concept{soil water characteristic}
\concept{soil water retention curve}
\concept{hydraulic conductivity function}
\concept{saturated hydraulic conductivity}

% \newcommand{\mbvec}{\eqn{\mbox{\boldmath$#1$\unboldmath}^\mathrm{#2} #4}{#1#3#4}}

\title{Evaporative Characteristic Length}

% description

\description{
  The functions \code{lc} and \code{lt} compute the \emph{characteristic
  length} \eqn{L_\mathrm{c}}{L_c} of \emph{stage-I} evaporation from a soil
  and its ``target'' (expected) value \eqn{L_\mathrm{t}}{L_t}, used to
  constrain estimates of nonlinear parameters of the Van Genuchten-Mualem
  (VGM) model for water retention curves and hydraulic conductivity
  functions.}

% usage

\usage{
lc(alpha, n, tau, k0, e0, c0 = NULL, c3 = NULL, c4 = NULL)

lt(n, tau, e0, c0, c1, c2, c3, c4)
}

% arguments

\arguments{
  \item{alpha}{parameter \eqn{\alpha} (inverse air entry pressure) of the
  VGM model, see \code{\link{wc_model}} and \code{\link{hc_model}}.  For
  consistency with other quantities, the unit of \eqn{\alpha} should be
  \bold{1/meter} [\eqn{\mathrm{m}^{-1}}{m^-1}].  }

  \item{n}{parameter \eqn{n} (shape parameter) of the VGM model, see
  \code{\link{wc_model}} and \code{\link{hc_model}}.}

  \item{tau}{parameter \eqn{\tau} (tortuosity parameter) of the VGM model,
  see \code{\link{hc_model}}.}

  \item{k0}{saturated hydraulic conductivity \eqn{K_0}, see
  \code{\link{hc_model}}.  If \code{k0} is missing or equal to \code{NA}
  in calls of \code{lc} then \code{k0} is approximated by the same relation
  as used for \code{lt}, see \emph{Details}.  For consistency with other
  quantities, the unit of \eqn{K_0} should be \bold{meter/day}
  [\eqn{\mathrm{m}\,\mathrm{d}^{-1}}{m d^-1}].}

  \item{e0}{a numeric scalar with the \emph{stage-I} rate of evaporation
  \eqn{E_0} from a soil, see \emph{Details} and
  \code{\link{soilhypfitIntro}}.  For consistency with other quantities,
  the unit of \eqn{E_0} should be \bold{meter/day}
  [\eqn{\mathrm{m}\,\mathrm{d}^{-1}}{m d^-1}].}

  \item{c0, c1, c2, c3, c4}{numeric constants to approximate the parameter
  \eqn{\alpha} and the saturated hydraulic conductivity \eqn{K_0} when
  computing \eqn{L_\mathrm{t}}{L_t}, see \emph{Details} and
  \code{\link{control_fit_wrc_hcc}}.  For consistency with other quantities,
  the following units should be used for the constants: \itemize{ \item
  \code{c1}: \eqn{\mathrm{m}^{-1}}{m^-1}, \item \code{c3}:
  \eqn{\mathrm{m}\,\mathrm{d}^{-1}}{m d^-1}.  } The remaining constants are
  dimensionless.

  }

}

% details

\details{

The \emph{characteristic length} of \emph{stage-I evaporation}
\eqn{L_\mathrm{c}}{L_c} (\cite{Lehmann et al., 2008, 2020}) is defined by
%
\deqn{
  L_\mathrm{c}(\mbox{\boldmath$\nu$\unboldmath}, K_0, E_0) = \frac{
    \frac{1}{\alpha \,§ n} \left(\frac{2n-1}{n-1}\right)^\frac{2n-1}{n}
  }{1 + \frac{E_0}{K_\mathrm{eff}}}
}{
  L_c(\nu, K_0, E_0) =
  1 / (\alpha * n) * ((2n-1)/(n-1))^((2n-1)/n) / (1 / ((1 + E_0/K_eff)))
}
%
where \mbvec{\nu}{T}{^T}{= (\alpha, n, \tau)} are the nonlinear parameters of
the VGM model, \eqn{K_0} is the saturated hydraulic conductivity,
\eqn{E_0} the stage-I evaporation rate and
\eqn{K_\mathrm{eff} = 4\, K_\mathrm{VGM}(h_\mathrm{crit}; K_0,
\mbox{\boldmath$\nu$\unboldmath})}{K_eff = 4 * K_VGM(h_crit; \nu)} is the
effective hydraulic conductivity at the critical pressure %
\deqn{
  h_\mathrm{crit} = \frac{1}{\alpha} \,
  \left(\frac{n-1}{n}\right)^\frac{1-2n}{n},
  }{
  h_crit = 1/\alpha * ((n-1)/n)^((1-2n)/n),
}
%
see \code{hc_model} for the definition of
\eqn{K_\mathrm{VGM}(h; K_0,
\mbox{\boldmath$\nu$\unboldmath})}{K_VGM(h; K_0, \nu)}.

The quantity \eqn{L_\mathrm{t}}{L_t} is the expected value (``target'') of
\eqn{L_\mathrm{c}}{L_c} for given shape (\eqn{n}) and tortuosity
(\eqn{\tau}) parameters.  To evaluate \eqn{L_\mathrm{t}}{L_t}, the
parameters \eqn{\alpha} and \eqn{K_0} are approximated by the following
relations
%
\deqn{
 \widehat{\alpha} = g_\alpha(n; c_0, c_1, c_2) =
    c_1 \, \frac{n - c_0}{1 + c_2 \, (n - c_0)},
}{
  hat\alpha = g_\alpha(n; c_0, c_1, c_2) = c_1 * (n - c_0) / (1 + c_2 * (n - c_0)),
}
%
\deqn{
 \widehat{K}_0 = g_{K_0}(n; c_0, c_3, c_4) =
    c_3 \, (n - c_0)^{c_4}.
}{
  hat K_0 = g_{K_0}(n; c_0, c_3, c_4) = c_3 * (n - c_0)^c_4.
}
%
The default values for \eqn{c_0} to \eqn{c_4} (see argument
\code{approximation_alpha_k0} of
\code{\link{control_fit_wrc_hcc}}) were estimated with data on African
desert regions of the database \emph{ROSETTA3} (\cite{Zhang and Schaap,
2017}), see \cite{Lehmann et al. (2020)} for details.

}

% value

\value{
A numeric scalar with the characteristic evaporative length (\code{lc}) or
its expected value (\code{lt}).}

\references{

  Lehmann, P., Assouline, S., Or, D. (2008) Characteristic lengths
  affecting evaporative drying of porous media.  \emph{Physical Review E},
  \bold{77}, 056309, \doi{10.1103/PhysRevE.77.056309}.

  Lehmann, P., Bickel, S., Wei, Z., Or, D. (2020) Physical Constraints for
  Improved Soil Hydraulic Parameter Estimation by Pedotransfer Functions.
  \emph{Water Resources Research} \bold{56}, e2019WR025963,
  \doi{10.1029/2019WR025963}.

  Zhang, Y. , Schaap, M. G. 2017. Weighted recalibration of the Rosetta
  pedotransfer model with improved estimates of hydraulic parameter
  distributions and summary statistics (Rosetta3). \emph{Journal of Hydrology},
  \bold{547}, 39-53, \doi{10.1016/j.jhydrol.2017.01.004}.

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

  \code{\link{wc_model}} and \code{\link{hc_model}} for currently
  implemented models for soil water retention curves and hydraulic
  conductivity functions;
%
%   \code{\link{evaporative-length}} for physically constraining parameter
%   estimates of soil hydraulic material functions.
}

% examples

\examples{
\donttest{
# use of \donttest{} because execution time exceeds 5 seconds

# estimate parameters of 4 samples of the Swiss forest soil dataset
# that have water retention (theta, all samples), saturated hydraulic conductivity
# (ksat) and optionally unsaturated hydraulic conductivity data
# (ku, samples "CH2_4" and "CH3_1")

data(swissforestsoils)

# select subset of data
sfs_subset <- droplevels(
  subset(
    swissforestsoils,
    layer_id \%in\% c("CH2_3", "CH2_4", "CH2_6", "CH3_1")
  ))

# extract ksat measurements
ksat <- sfs_subset[!duplicated(sfs_subset$layer_id), "ksat", drop = FALSE]
rownames(ksat) <- levels(sfs_subset$layer_id)
colnames(ksat) <- "k0"

# define number of cores for parallel computations
if(interactive()) ncpu <- parallel::detectCores() - 1L else ncpu <- 1L


# unconstrained estimation (global optimisation algorithm NLOPT_GN_MLSL)
# k0 fixed at measured ksat values
rsfs_uglob <- fit_wrc_hcc(
  wrc_formula = theta ~ head | layer_id,
  hcc_formula = ku ~ head | layer_id,
  data = sfs_subset,
  param = ksat,
  fit_param = default_fit_param(k0 = FALSE),
  control = control_fit_wrc_hcc(
    settings = "uglobal", pcmp = control_pcmp(ncores = ncpu)))
summary(rsfs_uglob)
coef(rsfs_uglob, lc = TRUE, gof = TRUE)


# constrained estimation by restricting ratio Lc/Lt to [0.5, 2]
# (global constrained optimisation algorithm NLOPT_GN_MLSL)
# k0 fixed at measured ksat values
rsfs_cglob <- update(
  rsfs_uglob,
  control = control_fit_wrc_hcc(
    settings = "cglobal", nloptr = control_nloptr(ranseed = 1),
    pcmp = control_pcmp(ncores = ncpu)))
summary(rsfs_cglob)
coef(rsfs_cglob, lc = TRUE, gof = TRUE)

# get initial parameter values from rsfs_cglob
ini_param <- cbind(
  coef(rsfs_cglob)[, c("alpha", "n")],
  ksat
)
# constrained estimation by restricting ratio Lc/Lt to [0.5, 2]
# (local constrained optimisation algorithm NLOPT_LD_CCSAQ)
# k0 fixed at measured ksat values
rsfs_cloc <- update(
  rsfs_uglob,
  param = ini_param,
  control = control_fit_wrc_hcc(
    settings = "clocal", nloptr = control_nloptr(ranseed = 1),
    pcmp = control_pcmp(ncores = ncpu)))
summary(rsfs_cloc)
coef(rsfs_cloc, lc = TRUE, gof = TRUE)

op <- par(mfrow = c(4, 2))
plot(rsfs_uglob, y = rsfs_cglob)
on.exit(par(op))

op <- par(mfrow = c(4, 2))
plot(rsfs_uglob, y = rsfs_cloc)
on.exit(par(op))

}}
