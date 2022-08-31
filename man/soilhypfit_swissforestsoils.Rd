% 2022-08-29 A. Papritz
% cd ~/R.user.home/soilhypfit/soilhypfit_01-7/pkg/man
% R CMD Rdconv -t html -o bla.html soilhypfit_swissforestsoils.Rd ; open bla.html; R CMD Rd2pdf --force soilhypfit_swissforestsoils.Rd;
% hunspell -d en_GB-ise swissforestsoils.Rd

\encoding{utf8}

\name{swissforestsoils}
\alias{swissforestsoils}

\docType{data}

\title{
  Physical properties of Swiss forest soils}

\description{
  The data give basic physical properties, water content and hydraulic
  conductivity measurements (at given capillary pressure head) for 128 soil
  layers (horizons) measured at 23 forest sites in Switzerland.
}

\usage{data(swissforestsoils)}

\format{
  A data frame with 1373 observations on the following 21 variables.
  \describe{

    \item{\code{profile_id}}{a factor with short labels for the 23 sites.}

    \item{\code{profile}}{a factor with with long labels for the 23 sites.}

    \item{\code{longitude}}{a numeric vector with the latitude of the site
    in degree.}

    \item{\code{latitude}}{a numeric vector with the latitude of the site
    in degree.}

    \item{\code{layer_id}}{a factor with labels for the 128 soil layer.}

    \item{\code{layer_ub}, \code{layer_lb}}{numeric vectors with the upper
    and lower depth (unit cm) of the soil layer for the measurements of the
    basic physical properties (\code{particle_density}, ..., \code{ksat}).}

    \item{\code{particle_density}}{a numeric vector with the density of the
    solid soil material (unit \eqn{\mathrm{g} \,\mathrm{cm^{-3}}}{g
    cm^{-3}}).}

    \item{\code{bulk_density}}{a numeric vector with soil (bulk) density
    (unit \eqn{\mathrm{g} \,\mathrm{cm^{-3}}}{g cm^{-3}}).}

    \item{\code{porosity}}{a numeric vector with the soil porosity (unit
    volume percentage).}

    \item{\code{clay}}{a numeric vector with the clay content (unit mass
    percentage).}

    \item{\code{silt}}{a numeric vector with the silt content (unit mass
    percentage).}

    \item{\code{sand}}{a numeric vector with the sand content (unit mass
    percentage).}

    \item{\code{ksat}}{a numeric vector with the saturated hydraulic
    conductivity (unit \eqn{\mathrm{m} \,\mathrm{d^{-1}}}{m d^{-1}}).}

    \item{\code{head}}{a numeric vector with capillary pressure head at
    which \code{theta} (water retention curve) and/or \code{ku} (hydraulic
    conductivity function) were measured (unit m).}

    \item{\code{layer_ub_theta}, \code{layer_lb_theta}}{numeric vectors
    with the upper and lower depth (unit cm) of the soil layer for which
    the water retention curve was measured.}

    \item{\code{theta}}{a numeric vector with volumetric water content
    measurements (dimensionless) of the water retention curve.}

    \item{\code{layer_ub_ku}, \code{layer_lb_ku}}{a numeric vector with the
    upper and lower depth (unit cm) of the soil layer for which the water
    retention curve was measured.}

    \item{\code{ku}}{a numeric vector with hydraulic conductivity
    measurements (unit \eqn{\mathrm{m} \,\mathrm{d^{-1}}}{m d^{-1}}) of
    the hydraulic conductivity function.}
  }
}

\details{
  \code{clay}, \code{silt} and \code{sand} refer to soil particles with
  diameter less than 2, between 2 and 50 and larger than 50 \eqn{\mu}m.
}

\source{
  Richard, F. & Lüscher, P. 1978 -- 1987.  Physikalische Eigenschaften von
  Böden der Schweiz.  Lokalformen Bände 1 -- 4.  Eidgenössische Anstalt für
  das forstliche Versuchswesen, Birmensdorf.
}


\examples{
\donttest{
# use of \donttest{} because execution time exceeds 5 seconds

# estimate parameters using all samples (samples with water retention,
# hydraulic conductivity, or with both type of measurements)

data(swissforestsoils)

# define number of cores for parallel computations
if(interactive()) ncpu <- parallel::detectCores() - 1L else ncpu <- 1L

# unconstrained estimation (global optimisation algorithm NLOPT_GN_MLSL)
r_uglob <- fit_wrc_hcc(
  wrc_formula = theta ~ head | layer_id,
  hcc_formula = ku ~ head | layer_id,
  data = swissforestsoils,
  control = control_fit_wrc_hcc(
    settings = "uglobal", pcmp = control_pcmp(ncores = ncpu)))
summary(r_uglob)
coef(r_uglob)

}}
\keyword{datasets}
