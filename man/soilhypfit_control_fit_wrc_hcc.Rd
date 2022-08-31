% 2022-08-29 A. Papritz
% cd ~/R.user.home/soilhypfit/soilhypfit_01-7/pkg/man
% R CMD Rdconv -t html -o bla.html soilhypfit_control_fit_wrc_hcc.Rd ; open bla.html; R CMD Rd2pdf --force soilhypfit_control_fit_wrc_hcc.Rd;
% hunspell -d en_GB-ise control_fit_wrc_hcc.Rd

\encoding{utf8}

\name{control_fit_wrc_hcc}
\alias{control_soilhypfit}
\alias{bwd_transf}
\alias{control_fit_wrc_hcc}
\alias{control_nloptr}
\alias{control_pcmp}
\alias{control_sce}
\alias{dfwd_transf}
\alias{fwd_transf}
\alias{param_boundf}
\alias{param_transf}

\concept{soil water}
\concept{soil water characteristic}
\concept{soil water retention curve}
\concept{hydraulic conductivity function}
\concept{saturated hydraulic conductivity}

% \newcommand{\mbvec}{\eqn{\boldsymbol{#1}^\mathrm{#2} #3}}

\title{Controlling fit_wrc_hcc}

% description

\description{
  This page documents options to control \code{\link{fit_wrc_hcc}}.  It
  describes the arguments of the functions \code{control_fit_wrc_hcc},
  \code{param_boundf}, \code{param_transf}, \code{fwd_transf},
  \code{dfwd_transf}, \code{bwd_transf}, \code{control_nloptr},
  \code{control_sce} and \code{control_pcmp}, which all serve to steer
  \code{\link{fit_wrc_hcc}}.}

% usage

\usage{
control_fit_wrc_hcc(
    settings = c("uglobal", "ulocal", "clocal", "cglobal", "sce"),
    method = c("ml", "mpd", "wls"), hessian,
    nloptr = control_nloptr(), sce = control_sce(),
    wrc_model = "vg", hcc_model = "vgm",
    initial_param = c(alpha = 2., n = 1.5, tau = 0.5),
    approximation_alpha_k0 =
        c(c0 = 1, c1 = 5.55, c2 = 1.204, c3 = 2.11, c4 = 1.71),
    variable_weight = c(wrc = 1, hcc = 1),
    gam_k = 6, gam_n_newdata = 101, precBits = 256,
    min_nobs_wc = 5, min_nobs_hc = 5,
    keep_empty_fits = FALSE,
    param_bound = param_boundf(), param_tf = param_transf(),
    fwd_tf = fwd_transf(), deriv_fwd_tfd = dfwd_transf(), bwd_tf = bwd_transf(),
    pcmp = control_pcmp())

param_boundf(alpha = c(0 + 10 * sqrt(.Machine$double.eps), 500),
    n = c(1 + 10 * sqrt(.Machine$double.eps), 20), tau = c(-2, 20),
    thetar = c(0, 1), thetas = c(0, 1), k0 = c(0, Inf))

param_transf(alpha = c("log", "identity"), n = c("log1", "identity"),
    tau = c("logitlu", "identity"), k0 = c("log", "identity"))

fwd_transf(...)

dfwd_transf(...)

bwd_transf(...)

control_nloptr(...)

control_sce(reltol = 1.e-8, maxtime = 20, ...)

control_pcmp(ncores = detectCores() - 1L,
    fork = !identical(.Platform[["OS.type"]], "windows"))

}

% arguments

\arguments{
  \item{settings}{a keyword with possible values \code{"uglobal"}
  (default), \code{"ulocal"}, \code{"clocal"}, \code{"cglobal"} or
  \code{"sce"} to choose the approach for the nonlinear optimisation, see
  \emph{Details} and \code{\link{soilhypfitIntro}}.}

  \item{method}{a keyword with possible values \code{"ml"} (maximum
  likelihood, default), \code{"mpd"} (maximum posterior density) or
  \code{"wls"} (weighted least squares) to choose the method for estimating
  the nonlinear parameters \mbvec{\nu}{}{}, see
  \code{\link{soilhypfitIntro}}.}

  \item{hessian}{a logical scalar controlling whether the Hessian matrix of
  the objective function should be computed with respect to the possibly
  transformed nonlinear parameters \mbvec{\nu}{}{} at the solution
  (default: \code{TRUE} if \code{settings \%in\% c("uglobal", "ulocal",
  "sce") && method \%in\% c("mpd", "ml")} and \code{FALSE} otherwise).}

  \item{nloptr}{a list of arguments passed to the optimiser
  \code{\link[nloptr]{nloptr}} by its argument \code{opts}, or a function
  such as \code{control_nloptr} that generates such a list.\cr Note that
  \code{control_fit_wrc_hcc} chooses sensible default values for the
  various components of the list in dependence of the value chosen for
  \code{settings}, but these defaults can be overridden by the arguments of
  \code{control_nloptr}, see \emph{Details}.}

  \item{sce}{a list of arguments passed to the optimiser
  \code{\link[SoilHyP]{SCEoptim}} by its argument \code{control}, or a
  function such as \code{control_sce} that generates such a list.\cr Note
  that \code{control_fit_wrc_hcc} chooses sensible default values for the
  various components of the list, but these defaults can be overridden by
  the arguments of \code{control_sce}, see \emph{Details}.}

  \item{wrc_model}{a keyword choosing the parametrical model for the
  water retention curve.  Currently, only the \emph{Van Genuchten} model
  (\code{"vg"}) is implemented, see \code{\link{wc_model}}
  and \code{\link{soilhypfitIntro}}.}

  \item{hcc_model}{a keyword choosing the parametrical model for the
  hydraulic conductivity function.  Currently, only the \emph{Van
  Genuchten-Mualem} model (\code{"vgm"}) is implemented, see
  \code{\link{hc_model}} and \code{\link{soilhypfitIntro}}.}

  \item{initial_param}{a named numeric vector with default initial values
  for the \emph{nonlinear} parameters \mbvec{\nu}{}{}
  of the models for the water retention curve and/or hydraulic conductivity
  function.  Currently, initial values must be defined for the parameters
  \eqn{\alpha} (default 1.5 [\eqn{\mathrm{m}^{-1}}{m^-1}]), \eqn{n}
  (default 2 [-]) and \eqn{\tau} (default 0.5 [-]), hence the elements of
  \code{initial_param} must be named \code{"alpha"}, \code{"n"} and
  \code{"tau"}.}

  \item{approximation_alpha_k0}{a named numeric vector with constants to
  approximate the parameter \eqn{\alpha} and the saturated hydraulic
  conductivity \eqn{K_0} when constraining the estimated \emph{nonlinear}
  parameters \mbvec{\nu}{}{} of the Van
  Genuchten-Mualem model by the \emph{characteristic evaporative length}
  (\cite{Lehmann et al., 2020}), see \code{\link{evaporative-length}} and
  \code{\link{soilhypfitIntro}}.  For consistency with other quantities,
  the following units should be used for the constants:
  \itemize{
    \item \code{c1}: \eqn{\mathrm{m}^{-1}}{m^-1},
    \item \code{c3}: \eqn{\mathrm{m}\,\mathrm{d}^{-1}}{m d^-1}.
  } The remaining constants are
  dimensionless.}

  \item{variable_weight}{a named numeric vector of length 2 that defines
  the weights of water content and hydraulic conductivity measurements in
  the objective function for \code{method = "wls"}.  If equal to 1
  (default) the weights of the variables are equal to the inverse variances
  of the water content and (possibly log-transformed) hydraulic
  conductivity measurements.  If different from 1, then the inverse
  variances are multiplied by \code{variable_weight} to get the variable
  weights \eqn{w_\theta} and \eqn{w_K}, see \code{\link{fit_wrc_hcc}} and
  \code{\link{soilhypfitIntro}}.  Note that for estimation methods mpd and
  ml the variable weights are equal to 1 but the case weights
  \eqn{w^\prime_{\theta,i}} and \eqn{w^\prime_{K,j}} may differ from 1.}

  \item{gam_k}{the dimension of the basis of the additive model for
  the water retention data when computing initial values of the
  parameters \eqn{\alpha} and \eqn{n}, see \code{\link[mgcv]{s}} and\cr
  \code{\link{soilhypfitIntro}}.}

  \item{gam_n_newdata}{the number of evaluation points of the additive
  model for the water retention data when computing initial values of the
  parameters \eqn{\alpha} and \eqn{n}, see \code{\link{soilhypfitIntro}}.}

  \item{precBits}{an integer scalar defining the default precision (in
  bits) to be used in high-precision computations by
  \code{\link[Rmpfr]{mpfr}}, see \emph{Details}.}

  \item{min_nobs_wc}{an integer scalar defining the minimum number of
  water content measurements per sample required for fitting a model to the
  water content data.}

  \item{min_nobs_hc}{an integer scalar defining the minimum number of
  hydraulic conductivity measurements per sample required for fitting a
  model to hydraulic conductivity data.}

  \item{keep_empty_fits}{a logical scalar controlling whether missing
  fitting results should be dropped for samples without any data or failed
  fits (\code{FALSE}, default) or kept (\code{TRUE}).}

  \item{param_bound}{a named list of numeric vectors of length 2 that
  define the allowed lower and upper bounds (box constraints) for the
  parameters of the models or a function such as \code{param_boundf} which
  generates this list, see \emph{Details}.}

  \item{param_tf}{a named vector of keywords that define the
  transformations to be applied to the model parameters before estimation
  or a function such as \code{param_transf}, which generates this vector,
  see \emph{Details}.}

  \item{fwd_tf}{a named list of invertible functions to be used to
  transform model parameters or a function such as \code{fwd_transf}, which
  generates this list, see \emph{Details}.}

  \item{deriv_fwd_tfd}{a named list of functions corresponding to the first
  derivatives of the functions defined in \code{fwd_tf} or a function such
  as \code{dfwd_transf}, which generates this list, see \emph{Details}.}

  \item{bwd_tf}{a named list of inverse functions corresponding the
  functions defined in \code{fwd_tf} or a function such as \code{bwd_transf},
  which generates this list, see \emph{Details}.}

  \item{pcmp}{a list to control parallel computations or a function such as
  \code{control_pcmp} that generates this list, see \code{control_pcmp} for
  allowed arguments.}

  \item{alpha}{either a numeric vector of length 2 defining the allowed
  lower and upper bounds for the nonlinear parameter \eqn{\alpha}
  (\code{param_boundf}), or a keyword defining the transformation
  to be applied to the parameter \eqn{\alpha} before estimation
  (\code{param_transf}), see \emph{Details}.}

  \item{n}{either a numeric vector of length 2 defining the allowed lower
  and upper bounds for the nonlinear parameter \eqn{n}
  (\code{param_boundf}), or a keyword defining the transformation
  to be applied to the parameter \eqn{n} before estimation
  (\code{param_transf}), see \emph{Details}.}

  \item{tau}{either a numeric vector of length 2 defining the allowed
  lower and upper bounds for the nonlinear parameter \eqn{\tau}
  (\code{param_boundf}), or a keyword defining the transformation
  to be applied to the parameter \eqn{\tau} before estimation
  (\code{param_transf}), see \emph{Details}.}

  \item{thetar}{a numeric vector of length 2 defining the allowed
  lower and upper bounds for the linear parameter \eqn{\theta_r}
  (\code{param_boundf}), see \emph{Details}.}

  \item{thetas}{a numeric vector of length 2 defining the allowed
  lower and upper bounds for the linear parameter \eqn{\theta_s}
  (\code{param_boundf}), see \emph{Details}.}

  \item{k0}{either a numeric vector of length 2 defining the allowed
  lower and upper bounds for the linear parameter \eqn{K_0}
  (\code{param_boundf}), or a keyword defining the transformation
  to be applied to the parameter \eqn{K_0} before estimation
  (\code{param_transf}), see \emph{Details}.}

  \item{reltol}{a numeric scalar defining (one possible) convergence
  criterion for the optimiser \code{SCEoptim}, see argument \code{control}
  of \code{\link[SoilHyP]{SCEoptim}} for details.}

  \item{maxtime}{a numeric scalar defining the maximum duration of
  optimisation in seconds by the optimiser \code{SCEoptim}, see see
  argument \code{control} of \code{\link[SoilHyP]{SCEoptim}} for details.}

  \item{ncores}{an integer defining the number of cores for
  parallel computations.  Defaults to the number of available cores
  minus one. \code{ncores = 1} suppresses parallel computations.}

  \item{fork}{a logical scalar controlling whether forking should be used
  for parallel computations (default: \code{TRUE} on Unix and MacOS and
  \code{FALSE} on Windows operating systems).  Note that stetting
  \code{fork = TRUE} on Windows suppresses parallel computations.}

  \item{...}{further arguments, such as details on parameter
  transformations (\code{fwd_transf}, \code{dfwd_transf},
  \code{bwd_transf}) or control options passed to the
  optimisers \code{\link[nloptr]{nloptr}} and
  \code{\link[SoilHyP]{SCEoptim}}, see \emph{Details}.}

%   \item{delta_sat_0}{a numeric scalar defining the maximum water saturation
%   that is considered to be different from zero.  Required for constrained
%   estimation of the residual (\eqn{\theta_r}) and saturated water content
%   (\eqn{\theta_s}).}
%
}

% details

\details{

  \subsection{Enforcing bounds on the estimated parameters}{

    Parameters of models for the water retention curve and the hydraulic
    conductivity function may vary only within certain bounds (see
    \code{\link{param_boundf}} for allowed ranges).  \code{\link{fit_wrc_hcc}}
    uses two mechanisms to constrain parameter estimates to permissible
    ranges:

    \enumerate{

        \item \emph{Parameter transformations}

        If a local algorithm is used for nonlinear optimisation
        (\code{settings = "ulocal"} or \code{settings = "clocal"}) and a
        transformation not equal to \code{"identity"} is specified in
        \code{param_tf} for any of the \emph{nonlinear} parameters
        \mbvec{\nu}{}{}, then the elements of \mbvec{\nu}{}{} are
        transformed by the functions given in \code{param_tf}.  The values
        of the transformed parameters vary then over the whole real line,
        and an unconstrained algorithm can be used for nonlinear
        optimisation.

        Note that the \emph{linear} parameters \eqn{\theta_r} (residual)
        and \eqn{\theta_s} (saturated water content) are never transformed
        and for the saturated hydraulic conductivity, \eqn{K_0}, only
        \code{"log"} (default) or \code{"identity"} can be chosen.
        Quadratic programming (see \code{\link[quadprog]{solve.QP}}) is
        employed to enforce the box constraints specified in the argument
        \code{param_bound} for \eqn{\theta_r} and \eqn{\theta_s}.
        Quadratic programming is also used to enforce the positivity
        constraint on \eqn{K_0} if \eqn{K_0} is not log-transformed
        (\code{"identity"}).  Otherwise, the logarithm of \eqn{K_0} is
        estimated unconstrainedly, see \code{\link{soilhypfitIntro}} for
        further details.

        \item \emph{Box constraints}

        If a global algorithm is used for the optimisation (\code{settings}
        equal to \code{"uglobal"} \code{"cglobal"} or \code{"sce"}) or if
        \code{"identity"} transformations are specified for all elements of
        \mbvec{\nu}{}{}, then an optimisation algorithm is deployed that
        respects the box constraints given in \code{param_bound}.  If
        parameters are estimated for several soil samples in a single call
        of \code{\link{fit_wrc_hcc}} and if sample-specific box constraints
        should be used then the lower and upper bounds of the
        box-constraints must be specified by the arguments
        \code{lower_param} and \code{upper_param} of the function
        \code{\link{fit_wrc_hcc}}, see explanations there.

        Further note that the transformations specified by \code{param_tf}
        for the nonlinear parameters \mbvec{\nu}{}{} are ignored when a
        global optimisation algorithm is used.

    }

  }

  \subsection{Parameter transformations}{

    The arguments \code{param_tf}, \code{fwd_tf}, \code{deriv_fwd_tfd},
    \code{bwd_tf} define how the model parameters are transformed for
    estimation by local optimisation algortihms (see above and
    \code{\link{soilhypfitIntro}}).  The following transformations are
    currently available:

    \describe{

      \item{\code{"log"}:}{\eqn{\log(x)}{log(x)},}

      \item{\code{"log1"}:}{\eqn{\log(x-1)}{log(x-1)},}

      \item{\code{"logitlu"}:}{\eqn{\log((x - l) / (u - x))} with \eqn{l}
      and \eqn{u} the allowed lower and upper bounds for a parameter, see
      \code{param_boundf},}

      \item{\code{"identity"}:}{no transformation.}

    }

    These are the possible values that the various arguments of the
    function \code{param_transf} accept (as quoted character strings), and
    these are the names of the list components returned by
    \code{fwd_transf}, \code{dfwd_transf} and \code{bwd_transf}.

    Additional transformations can be implemented by:

    \enumerate{

      \item Extending the function definitions by arguments like

      \code{fwd_tf = fwd_transf(my_fun = function(x) your transformation)},\cr
      \code{deriv_fwd_tfd = dfwd_transf(my_fun = function(x) your derivative)},\cr
      \code{bwd_tf = bwd_transf(my_fun = function(x) your back-transformation)},

      \item Assigning to a given argument of \code{param_transf} the name
      of the new function, e.g.\cr \code{alpha = "my_fun"}.  }

    Note that the values given to the arguments of \code{param_transf} must
    match the names of the functions returned by \code{fwd_transf},
    \code{dfwd_transf} and \code{bwd_transf}.

  }


  \subsection{High-precision numerical computations}{

    Estimation of \eqn{\log(K_0)}{log(K_0)} is somewhat delicate for large
    values of the shape parameter \eqn{n} and/or small values of
    \eqn{\alpha}.  The water saturation and the relative conductivity are
    then close to zero for capillary pressure head exceeding
    \eqn{1/\alpha}.  To avoid numerical problems caused by limited accuracy
    of double precision numbers, \code{\link{fit_wrc_hcc}} uses the
    function \code{\link[Rmpfr]{mpfr}} of the package \pkg{Rmpfr} for
    high-accuracy computations.  The argument \code{precBits} of
    \code{\link{control_fit_wrc_hcc}} controls the accuracy.  Increase its
    value if computation fails with a respective diagnostic message.

  }


  \subsection{Options to choose the approach for nonlinear optimisation}{

  The argument \code{settings} defines sets of default options to control
  the optimisers.  The following \code{settings} are currently defined:

  \describe{

     \item{\code{"uglobal"}:}{unconstrained optimisation by any of the
     global algorithms (named \code{"NLOPT_G..."}) of the NLopt library.}

     \item{\code{"cglobal"}:}{constrained optimisation by the global
     algorithm \code{"NLOPT_GN_ISRES"} of NLopt that allows for inequality
     constraints.}

     \item{\code{"ulocal"}:}{unconstrained optimisation by any of the
     local algorithms\cr (named \code{"NLOPT_L..."}) of
     NLopt.}

     \item{\code{"clocal"}:}{constrained optimisation by any of the local
     algorithms\cr (\code{"NLOPT_LN_COBYLA"}, \code{"NLOPT_LN_AUGLAG"},
     \code{"NLOPT_LD_AUGLAG"}, \code{"NLOPT_LD_SLSQP"},\cr
     \code{"NLOPT_LD_MMA"}), \code{"NLOPT_LD_CCSAQ"}) of NLopt that allow
     for inequality constraints.}

     \item{\code{"sce"}:}{unconstrained optimisation by the global
     algorithm implemented in \code{\link[SoilHyP]{SCEoptim}}.}

  }



  The functions \code{control_nloptr} and \code{control_sce} allow finer
  control of the optimisers.
  %
  \code{control_nloptr} and \code{control_sce} take any
  argument available to pass controlling options to the optimisers
  \code{\link[nloptr]{nloptr}} (by its argument \code{opts}) and
  \code{\link[SoilHyP]{SCEoptim}} (by its argument \code{control}),
  respectively.

    \subsection{Controlling nloptr}{

      The function \code{\link[nloptr]{nloptr.print.options}} prints all
      options available to control \code{\link[nloptr]{nloptr}} by its
      argument \code{opts}.  Detailed information on the options can be
      found in the
      \href{https://nlopt.readthedocs.io/en/latest/NLopt_manual/}{NLopt
      documentation}.\cr

      The function \code{control_fit_wrc_hcc} sets meaningful defaults for
      \code{opts} in dependence of the chosen optimisation approach as
      specified by the argument \code{settings}, and it checks the
      consistency of the arguments of \code{control_nloptr} if they are
      explicitly passed to \code{fit_wrc_hcc}.\cr

      The following defaults are set by \code{control_fit_wrc_hcc} for the
      argument \code{opts} of \code{\link{nloptr}} (:

      \enumerate{

        % unconstrained, global optimisation

        \item Unconstrained, global optimisation (\code{settings =
        "uglobal"}):

        \preformatted{
nloptr = control_nloptr(
  algorithm = "NLOPT_GN_MLSL_LDS",
  local_opts = list(
    algorithm = "NLOPT_LN_BOBYQA",
    xtol_rel = -1.,
    ftol_rel = 1.e-6
  ),
  xtol_rel = -1,
  ftol_rel = -1,
  maxeval = 125,
  maxtime = -1)
        }

        In addition, any parameter transformations specified by
        \code{param_tf} are overridden and the untransformed parameters
        (\code{"identity"}) are estimated when \code{settings = "uglobal"} is
        chosen.


        % constrained, global optimisation

        \item Constrained, global optimisation (\code{settings =
        "cglobal"}):

        \preformatted{
nloptr = control_nloptr(
  algorithm = "NLOPT_GN_ISRES",
  xtol_rel = -1,
  ftol_rel = -1,
  maxeval = 1000,
  maxtime = -1)
        }

        In addition, any parameter transformations specified by
        \code{param_tf} are overridden and the untransformed parameters
        (\code{"identity"}) are estimated when \code{settings = "cglobal"}
        is chosen.


        % unconstrained, local optimisation

        \item Unconstrained, local optimisation (\code{settings =
        "ulocal"}):

        \preformatted{
nloptr = control_nloptr(
  algorithm = "NLOPT_LN_BOBYQA",
  xtol_rel = -1,
  ftol_rel = 1.e-8,
  maxeval = 250,
  maxtime = -1)
        }

        % unconstrained, local optimisation

        \item Constrained, local optimisation (\code{settings = "clocal"}):

        \preformatted{
nloptr = control_nloptr(
  algorithm = "NLOPT_LD_CCSAQ",
  xtol_rel = -1,
  ftol_rel = 1.e-8,
  maxeval = 1000,
  maxtime = -1)
        }

        If the algorithm \code{"NLOPT_LD_AUGLAG"} is used for constrained,
        local optimisation then

        \preformatted{
nloptr = control_nloptr(
  algorithm = "NLOPT_LD_AUGLAG",
  local_opts = list(
    algorithm = "NLOPT_LD_LBFGS",
    xtol_rel = -1.,
    ftol_rel = 1.e-6
  ),
  xtol_rel = -1,
  ftol_rel = 1.e-8,
  maxeval = 1000,
  maxtime = -1)
        }


      }

      For other, unspecified elements of \code{opts} default values as
      listed by \code{\link[nloptr]{nloptr.print.options}} are used.

    }


    \subsection{Controlling SCEoptim}{

      The function \code{control_sce} sets meaningful defaults for the
      argument \code{control} of \code{\link[SoilHyP]{SCEoptim}}.
      Currently, the following defaults are defined:

      \preformatted{
sce = control_sce(
  reltol = 1e-08,
  maxtime = 20)
      }

      In addition, any parameter transformations specified by
      \code{param_tf} are overridden and the untransformed parameters
      (\code{"identity"}) are estimated when \code{settings = "sce"} is
      chosen.

    }




  }

}

% value

\value{

  \code{control_fit_wrc_hcc} creates a list with components
  \code{settings}, \code{hessian}, \code{method}, \code{nloptr},
  \code{sce}, \code{wrc_model}, \code{hcc_model}, \code{initial_param},
  \code{approximation_alpha_k0}, \code{variable_weight}, \code{gam_k},
  \code{gam_n_newdata}, \code{precBits}, \code{min_nobs_wc},
  \code{min_nobs_hc}, \code{keep_empty_fits}, \code{param_bound},
  \code{param_tf}, \code{fwd_tf}, \code{deriv_fwd_tfd}, \code{bwd_tf},
  \code{pcmp} corresponding to its arguments and some further components
  (\code{delta_sat_0}, \code{grad_eps}, \code{use_derivative}) that cannot
  be changed by the user.

  \code{control_nloptr} and \code{control_sce}
  create lists with control parameters passed to
  \code{\link[nloptr]{nloptr}} and \code{\link[SoilHyP]{SCEoptim}},
  respectively, see \emph{Details}.

  \code{param_boundf} generates a list with allowed lower and upper bounds of
  the model parameters.

  \code{param_transf} generates a list with keywords that
  define what transformations are used for estimating the model
  parameters, and \code{fwd_transf}, \code{bwd_transf} and
  \code{dfwd_transf} return lists of functions with forward and backward
  transformations and the first derivatives of the forward
  transformations, see \emph{Details}.

  \code{control_pcmp} generates a list with control parameters for parallel
  computations.

}

% references

\references{

  Johnson, S.G. The NLopt nonlinear-optimisation package.
  \url{https://github.com/stevengj/nlopt}.

  Lehmann, P., Assouline, S., Or, D. (2008) Characteristic lengths
  affecting evaporative drying of porous media.  \emph{Physical Review E},
  \bold{77}, 056309, \doi{10.1103/PhysRevE.77.056309}.

  Lehmann, P., Bickel, S., Wei, Z., Or, D. (2020) Physical Constraints for
  Improved Soil Hydraulic Parameter Estimation by Pedotransfer Functions.
  \emph{Water Resources Research} \bold{56}, e2019WR025963,
  \doi{10.1029/2019WR025963}.

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

%   \code{\link{control_fit_wrc_hcc}} for options to control
%   \code{fit_wrc_hcc};
%
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
\donttest{
# use of \donttest{} because execution time exceeds 5 seconds
data(sim_wrc_hcc)

# estimate parameters for a single soil sample by maximizing loglikelihood ...

# ... with unconstrained, global optimisation algorithm NLOPT_GN_MLSL
coef(
  fit1 <- fit_wrc_hcc(
    wrc_formula = wc ~ head, hcc_formula = hc ~ head,
    data = sim_wrc_hcc, subset = id == 2
  ), gof = TRUE)

# ... as fit1 but fitting parameter tau as well
coef(
  fit2 <- update(fit1,
    fit_param = default_fit_param(tau = TRUE)
  ), gof = TRUE)

plot(fit1, y = fit2)

# ... with unconstrained, local optimisation algorithm NLOPT_LN_BOBYQA,
#     initial values for alpha and n are computed from data and
#     transformed nonlinear parameters are estimated without box-constraints
coef(
  fit3 <- update(
    fit2,
    control = control_fit_wrc_hcc(settings = "ulocal"),
    verbose = 2), gof = TRUE)

# estimate parameters by unconstrained, weighted least squares minimisation with
#     algorithm NLOPT_LD_LBFGS, giving larger weight to conductivity data,
#     using specified initial values for alpha and n and
#     fitting untransformed nonlinear parameters with default box constraints
#     defined by param_boundf()
#     diagnostic output directly from nloptr
coef(
  fit4 <- update(
    fit2,
    param = c(alpha = 1.7, n = 2),
    control = control_fit_wrc_hcc(
      settings = "ulocal", method = "wls",
      variable_weight = c(wrc = 1, hcc = 2),
      nloptr = control_nloptr(algorithm = "NLOPT_LD_LBFGS", print_level = 3),
      param_tf = param_transf(alpha = "identity", n = "identity", tau = "identity")
    ), verbose = 0), gof = TRUE)

# ... as fit4 but giving even larger weight to conductivity data
coef(
  fit5 <- update(
    fit4,
    control = control_fit_wrc_hcc(
      settings = "ulocal", method = "wls",
      variable_weight = c(wrc = 1, hcc = 5),
      nloptr = control_nloptr(algorithm = "NLOPT_LD_LBFGS", print_level = 3),
      param_tf = param_transf(alpha = "identity", n = "identity", tau = "identity")
    ), verbose = 0), gof = TRUE)

plot(fit4, y = fit5)

}}
