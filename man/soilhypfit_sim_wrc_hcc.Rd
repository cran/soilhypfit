% 2022-12-17 A. Papritz
% cd ~/R.user.home/soilhypfit/soilhypfit_01-7/pkg/man
% R CMD Rdconv -t html -o bla.html soilhypfit_sim_wrc_hcc.Rd ; open bla.html; R CMD Rd2pdf --force soilhypfit_sim_wrc_hcc.Rd;
% hunspell -d en_GB-ise sim_wrc_hcc.Rd

\encoding{utf8}

\name{sim_wrc_hcc}
\alias{sim_wrc_hcc}

\docType{data}

\title{
  Simulated Soil Water Retention and Hydraulic Conductivity Data
}

\description{
  The data give simulated values of water content and hydraulic
  conductivity at given capillary pressure head for 3 soil samples.
}

\usage{data(sim_wrc_hcc)}

\format{
  A data frame with 28 observations on the following 4 variables.

    \describe{ \item{\code{id}}{a factor with levels \code{1}, \code{2},
    \code{3} coding the soil samples.}

    \item{\code{head}}{a numeric vector with values of capillary pressure
    head (unit \eqn{\mathrm{m}}{m}).}

    \item{\code{wc}}{a numeric vector with values of simulated volumetric
    water content (unit -)}

    \item{\code{hc}}{a numeric vector with values of simulated hydraulic
    conductivity (unit \eqn{\mathrm{m} \, \mathrm{d}^{-1}}{m d^-1}).} }
}

\details{

  The values of \code{wc} and \code{hc} were simulated by the model of
  \emph{Van Genuchten Mualem} (Van Genuchten, 1980, see
  \code{\link{wc_model}} and \code{\link{hc_model}}) with the following
  parameters:
%
  \tabular{rrrrrrr}{
    sample id \tab \eqn{\theta_r} \tab \eqn{\theta_s} \tab \eqn{K_0} [\eqn{\mathrm{m} \, \mathrm{d}^{-1}}{m d^-1})]
              \tab \eqn{\alpha} [\eqn{\mathrm{m}^{-1}}{m^-1}] \tab \eqn{n} \tab \eqn{\tau} \cr
    1 \tab 0.05 \tab 0.45 \tab 0.1 \tab 2   \tab 2   \tab 0.5 \cr
    2 \tab 0.1  \tab 0.5  \tab 5   \tab 1.5 \tab 1.5 \tab 0.5 \cr
    3 \tab 0.05 \tab 0.45 \tab 2   \tab 1.4 \tab 1.3 \tab 0.5 \cr
  }
%
  Normally distributed errors were added to the model values (\code{wc}:
  sd: 0.05; \code{log(hc)}: sd 0.5).

}

% \source{
% %%  ~~ reference to a publication or URL from which the data were obtained ~~
% }

\references{
  Van Genuchten, M. Th. (1980) A closed-form equation for predicting the
  hydraulic conductivity of unsaturated soils. \emph{Soil Science Society of
  America Journal}, \bold{44}, 892--898,
  \doi{10.2136/sssaj1980.03615995004400050002x}.
}

\examples{
data(sim_wrc_hcc)

if(!requireNamespace("lattice", quietly = TRUE)){
  install.packages("lattice", repos = "https://cloud.r-project.org/")
  requireNamespace("lattice", quietly = TRUE)
}

lattice::xyplot(wc ~ head|id, type = "l", sim_wrc_hcc, 
  scales = list( x = list(log=TRUE)))
lattice::xyplot(hc ~ head|id, type = "l", sim_wrc_hcc, 
  scales = list( log = TRUE))

}

\keyword{datasets}
