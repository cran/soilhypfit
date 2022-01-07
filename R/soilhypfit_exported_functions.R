## ======================================================================
lc <- function(alpha, n, tau, k0, e0, c0 = NULL, c3 = NULL, c4 = NULL){

  ## function computes capillary length Lc based on VGM model and
  ## evaporation rate of a saturated soil; if k0 is not available then
  ## the approxmation k0 = c3 * (n - c0)^c4 is used

  ## 2019-11-27 A. Papritz

  ## approximation for missing k0

  if(missing(k0) || is.na(k0)){
    if(
      any(c(is.null(c0), is.null(c3), is.null(c4))) ||
      n < c0
    ){
      return(NA_real_)
    } else {
      t.k0 <- c3 * (n - c0)^c4
    }
  } else {
    t.k0 <- k0
  }

  ## auxiliary computations

  nm1 <- -1 + n
  m2p1on <- -2 + 1/n
  aux1 <- (1 + (((nm1)/n)^(m2p1on))^n)^(-1 + 1/n)

  ## Lc

  result <- ((nm1)/(-1 + 2*n))^(m2p1on) /
  ((1 + e0/(4.*t.k0*(-1 + (1 - (aux1)^(n/(nm1)))^((nm1)/n))^2*(aux1)^tau))*n*alpha)

  attr(result, "input") <- c(
    alpha = unname(alpha), n = unname(n), tau = unname(tau),
    k0 = unname(t.k0), e0 = unname(e0))

  unname(result)

}


## ======================================================================
lt <- function(n, tau, e0, c0, c1, c2, c3, c4){

  ## function computes constraining value for capillary length Lc based on
  ## VGM model and evaporation rate of a saturated soil and approximations
  ## for alpha and k0

  ## 2019-11-27 A. Papritz

  if(
    any(c(missing(c0), missing(c1), missing(c2), missing(c3), missing(c4))) ||
    n < c0
  ){
    return(NA_real_)
  } else {
    lc(
      alpha = c1*(n - c0)/(1 + c2*(n - c0)),
      n = n, tau = tau,
      k0 = c3 * (n - c0)^c4,
      e0 = e0
    )
  }

}


## ======================================================================
sat_model <- function(h, nlp, precBits = NULL, wrc_model = "vg"){

  ## function computes relative capillary saturation based on van
  ## Genuchten-Mulalem (vGM) model

  ## 2019-11-27 A. Papritz


  switch(
    wrc_model,
    vg = {
      alpha <- nlp["alpha"]
      n <- nlp["n"]
      if(!is.null(precBits)) h <- mpfr(h, precBits)
      (1 + (h * alpha)^n)^(-1 + 1/n)
    },
    stop("wrc model '", wrc_model, "' not implemented")
  )

}


## ======================================================================
wc_model <- function(h, nlp, lp, precBits = NULL, wrc_model = "vg"){

  ## function computes water content based on van Genuchten-Mulalem (vGM)
  ## model

  ## 2019-11-27 A. Papritz

  switch(
    wrc_model,
    vg = {
      lp["thetar"] +
      (lp["thetas"] - lp["thetar"]) * sat_model(h, nlp, precBits, wrc_model)
    }
  )

}


## ======================================================================
hcrel_model <- function(h, nlp, precBits = NULL, hcc_model = "vgm"){

  ## function computes relative hydraulic conductivity based on
  ## van Genuchten-Mulalem (vGM) model

  ## 2019-11-27 A. Papritz

  switch(
    hcc_model,
    vgm = {
      n <- nlp["n"]
      tau <- nlp["tau"]
      if(!is.null(precBits)){
        one.mpfr <- mpfr(1, precBits)
      } else {
        one.mpfr <- 1.
      }
      sat <- sat_model(h, nlp, precBits, wrc_model = "vg")
      (sat^tau) * (one.mpfr - (one.mpfr - sat^(n/(n-1)))^((n-1)/n))^2
    },
    stop("hcc model'", hcc_model, "' not implemented")
  )
}


## ======================================================================
hc_model <- function(h, nlp, lp, precBits = NULL, hcc_model = "vgm"){

  ## function computes hydraulic conductivity based on
  ## van Genuchten-Mulalem (vGM) model

  ## 2019-11-27 A. Papritz

  switch(
    hcc_model,
    vgm = {
      lp["k0"] * hcrel_model(h, nlp, precBits, hcc_model)
    }
  )

}


##  ##############################################################################

param_boundf <- function(
  alpha = c(0. + 10. * sqrt(.Machine$double.eps), 5.e2),
  n = c(1. + 10. * sqrt(.Machine$double.eps), 20.), tau = c(-2., 20.),
  thetar = c(0., 1.), thetas = c(0., 1.), k0 = c(0., Inf)
){

  ## function returns valid range of values for nonlinear parameters

  ## 2019-11-27 A. Papritz
  ## 2021-02-28 AP  checking of valid boundaries in new function
  ##                check_param_boundf

  param_bound <- list(
    alpha = alpha, n = n, tau = tau,
    thetar = thetar, thetas = thetas,
    k0 = k0
  )

  check_param_boundf(param_bound)

  param_bound

}


##  ##############################################################################

check_param_boundf <- function(
  y, compare_thetar_thetas = TRUE
){

  ## function checks whether boundaries for box constraints are valid

  ## 2021-02-27 A. Papritz

  ## 2021-12-26 AP possibility to check incomplete parameter
  ##               boundaries
  ## 2022-01-03 AP optional cross-comparison of boundaries of
  ##               thetar and thetas

  if(
    !is.null(y[["thetar"]]) && (
      y[["thetar"]][1] < 0. || y[["thetar"]][2] > 1. ||
      y[["thetar"]][1] > y[["thetar"]][2]
    )
  ) stop(
    "invalid range of possible values for parameter 'thetar'"
  )

  if(
    !is.null(y[["thetas"]]) && (
      y[["thetas"]][1] < 0. || y[["thetas"]][2] > 1. ||
      y[["thetas"]][1] > y[["thetas"]][2]
    )
  ) stop(
    "invalid range of possible values for parameter 'thetas'"
  )

  if(
    compare_thetar_thetas && !is.null(y[["thetar"]]) && !is.null(y[["thetas"]]) && (
      y[["thetas"]][1] < y[["thetar"]][1] ||
      y[["thetar"]][2] > y[["thetas"]][2]
    )
  ) stop(
    "invalid ranges of possible values for parameters 'thetar' and 'thetas'"
  )

  if(
    !is.null(y[["alpha"]]) && (
      y[["alpha"]][1] <= 0. || y[["alpha"]][1] > y[["alpha"]][2]
    )
  ) stop(
    "invalid range of possible values for parameter 'alpha'"
  )

  if(
    !is.null(y[["n"]]) && (
      y[["n"]][1] < 1. || y[["n"]][1] > y[["n"]][2]
    )
  ) stop(
    "invalid range of possible values for parameter 'n'"
  )

  if(
    !is.null(y[["k0"]]) && (
      y[["k0"]][1] < 0.|| y[["k0"]][1] > y[["k0"]][2]
    )
  ) stop(
    "invalid range of possible values for parameter 'k0'"
  )

  "no error"

}


## ======================================================================
param_transf <- function(
  #   thetar = "identity", thetas = "identity",
  alpha = c("log", "identity"),
  n   = c("log1", "identity"),
  tau = c("logitlu", "identity"),
  k0 = c("log", "identity")

){

  ## function sets meaningful defaults for transformation of model
  ## parameters

  ## 2019-11-27 A. Papritz

  alpha <- match.arg(alpha)
  n     <- match.arg(n)
  tau   <- match.arg(tau)
  k0  <- match.arg(k0)

  list(
    #     thetar = thetar, thetas = thetas,
    alpha = alpha, n = n,
    k0 = k0,
    tau = tau
  )

}


## ======================================================================
fwd_transf <- function(
  ...
){

  ## definition of forward transformation of model parameters

  ## 2019-11-27 A. Papritz

  list(
    log       = function(x, ...) log(x),
    log1      = function(x, ...) log(x - 1.),
    logit01   = function(x, ...) log(x / (1. - x)),
    logitlu   = function(x, lb, ub) log((x - lb) / (ub - x)),
    identity  = function(x, ...) x, ...
  )
}


## ======================================================================
dfwd_transf<- function(
  ...
){

  ## definition of first derivative of forward transformation of model
  ## parameters

  ## 2019-11-27 A. Papritz

  list(
    log       = function(x, ...) 1. / x,
    log1      = function(x, ...) 1. / (x - 1.),
    logit01   = function(x, ...) 1. / (x - x^2),
    logitlu   = function(x, lb, ub) 1. / (ub - x) + 1. / (x - lb),
    identity  = function(x, ...) rep(1., length(x)), ...
  )

}


## ======================================================================
bwd_transf <- function(
  ...
){

  ## definition of backward transformation of model parameters

  ## 2019-11-27 A. Papritz

  list(
    log       = function(x, ...) exp(x),
    log1      = function(x, ...) exp(x) + 1.,
    logit01   = function(x, ...){
      if(!is.finite(x)){
        if(sign(x) < 0.) 0. else 1.
      } else exp(x) / (1. + exp(x))
    },
    logitlu   = function(x, lb, ub){
      if(!is.finite(x)){
        if(sign(x) < 0.) 0. else 1.
      } else (ub * exp(x)  + lb) / (1. + exp(x))
    },
    identity = function(x, ...) x, ...
  )
}


## ======================================================================
control_nloptr <- function(
  ...
)
{

  ## function sets meaningful defaults for selected control arguments of
  ## function nloptr

  ## 2019-11-27 A. Papritz

  list(
    ...
  )
}


## ======================================================================
control_sce <- function(
  reltol = 1.e-8, maxtime = 20, ...
)
{

  ## function sets meaningful defaults for control arguments of function
  ## SCEoptim of package hydromad

  ## 2019-11-27 A. Papritz

  list(
    reltol = reltol,
    maxtime = 20,
    ...
  )
}


## ======================================================================
control_pcmp <-
  function(
    ncores = detectCores() - 1L,
    fork = !identical( .Platform[["OS.type"]], "windows" )
  )
{

  ## function sets meaningful defaults for parallelized computations

  ## 2014-07-29 AP
  ## 2015-06-30 AP function and arguments renamed
  ## 2015-07-29 AP changes for elimination of parallelized computation of gradient or estimating equations
  ## 2016-07-20 AP renamed function, separate ncores arguments various parallelized computations
  list(
    ncores = ncores,
    fork = fork
  )

}

## ======================================================================
### control_fit_wrc_hcc

control_fit_wrc_hcc <- function(
  settings = c("uglobal", "ulocal", "clocal", "cglobal", "sce"),
  method = c("ml", "mpd", "wls"),
  hessian,
  nloptr = control_nloptr(),
  sce = control_sce(),
  wrc_model = "vg", hcc_model = "vgm",
  initial_param = c(alpha = 2., n = 1.5, tau = 0.5),
  approximation_alpha_k0 = c(
    c0 = 1.0, c1 = 5.55, c2 = 1.204, c3 = 2.11, c4 = 1.71
  ),
  variable_weight = c(wrc = 1., hcc = 1.),
  gam_k = 6,
  gam_n_newdata = 101,
  precBits = 256,
  min_nobs_wc = 5,
  min_nobs_hc = 5,
  keep_empty_fits = FALSE,
  param_bound = param_boundf(),
  param_tf = param_transf(),
  fwd_tf = fwd_transf(),
  deriv_fwd_tfd = dfwd_transf(),
  bwd_tf = bwd_transf(),
  pcmp = control_pcmp()
){

  ## function set meaningful default for settings that control execution of
  ## function fit_wrc_hcc

  ## 2019-11-27 A. Papritz
  ## 2020-01-06 AP computation of Hessian matrix at solution
  ## 2020-01-23 AP new argument method
  ## 2020-01-29 AP new default value argument hessian
  ## 2021-05-21 AP new default values for maxeval for settings == "clocal"
  ## 2021-05-23 AP new list of allowed algorithms for settings == "cglobal"
  ## 2021-05-21 AP new default values for maxeval for settings == "cglobal|clocal"
  ## 2021-10-25 AP new default value "mpd" for methods argument
  ## 2021-12-06 AP new local constrained algorithm NLOPT_LD_CCSAQ
  ## 2021-12-07 AP new default local constrained algorithm NLOPT_LD_CCSAQ
  ## 2021-12-22 AP new default value "ml" for methods argument

#### -- check arguments

 ## match settings, method, wrc_model, hcc_model arguments

  settings  <- match.arg(settings)
  method    <- match.arg(method)
  wrc_model <- match.arg(wrc_model)
  hcc_model <- match.arg(hcc_model)

  ## check availability of functions for parameter transformations

  if(
    !(all(unlist(param_tf) %in% names(fwd_tf)) &&
      all(unlist(param_tf) %in% names(deriv_fwd_tfd)) &&
      all(unlist(param_tf) %in% names(bwd_tf))
    )
  ) stop(
    "undefined transformation of parameters; extend respective function definitions"
  )

  ## set value for hessian if missing

  if(missing(hessian)){
    hessian <- settings %in% c("uglobal", "ulocal", "sce") &&
    method %in% c("ml", "mpd")
  }

  ## check variable weights

  names(variable_weight) <- c("wrc", "hcc")
  sel <- variable_weight[!is.na(variable_weight)] < 0.
  if(length(sel) && any(sel)) stop("'variable_weight' must be positive")

#### -- prepare nloptr options

  ## match names of arguments of control_nloptr

  nloptr.dopt <- nloptr.get.default.options()

  names.opts <- NULL
  if(length(nloptr)) names(nloptr) <- names.opts <- match.arg(
    names(nloptr), c(nloptr.dopt[, "name"], "local_opts"), several.ok = TRUE
  )

  ## check whether specified algorithm is valid

  ## list of available algorithms

  nloptr.algorithms <- unlist(strsplit(
      nloptr.dopt[grep("algorithm", nloptr.dopt[, "name"]), "possible_values"],
      ", ", fixed = TRUE
    ))

  if("algorithm" %in% names.opts){
    nloptr[["algorithm"]] <- match.arg(
      toupper(nloptr[["algorithm"]]), choices = nloptr.algorithms
    )
  }

#### --- options for global NLopt optimization algorithms

  if(settings %in% c("uglobal", "cglobal")){

    ## specify algorithm

    if("algorithm" %in% names.opts){

      ## check whether specified algorithm does global optimization

      if(length(grep("^NLOPT_G", nloptr[["algorithm"]])) == 0L ) stop(
        "local optimization algorithm specified for global optimization"
      )

      if(settings == "cglobal" &&
        !nloptr[["algorithm"]] %in%
        c("NLOPT_GN_ISRES", "NLOPT_GN_ORIG_DIRECT", "NLOPT_GN_ORIG_DIRECT_L")
      ) stop(
        "wrong algorithm specified for constrained global optimization"
      )

    } else {

      ## choose default global optimization algorithm

      if(identical(settings, "uglobal")){
        nloptr[["algorithm"]] <- "NLOPT_GN_MLSL_LDS"
      } else {
        nloptr[["algorithm"]] <- "NLOPT_GN_ISRES"
      }

    }

#### ---- MLSL algorithm: check and specify local options

    if(length(grep("MLSL", nloptr[["algorithm"]]))){

      if("local_opts" %in% names.opts){

        ## local options present

        ## match names of local options

        names.local.opts <- names(nloptr[["local_opts"]])
        names.local.opts <- match.arg(
          names.local.opts, c("algorithm", "xtol_rel", "ftol_rel"),
          several.ok = TRUE
        )
        names(nloptr[["local_opts"]]) <- names.local.opts

        ## match name of algorithm for local solver

        if(length(grep("algorithm", names.local.opts))){

          ## match algorithm of local solver

          nloptr[["local_opts"]][["algorithm"]] <- match.arg(
            toupper(nloptr[["local_opts"]][["algorithm"]]),
            nloptr.algorithms[grep("NLOPT_L", nloptr.algorithms)]
          )

        }

        ## check and specify algorithm for local solver if missing

        if(length(grep("^NLOPT_GD", nloptr[["algorithm"]]))){

          ## outer algorithm uses derivatives

          if("algorithm" %in% names.local.opts){

            ## make sure that algorithm of local solver uses derivatives

            if(
              length(grep("NLOPT_LN", nloptr[["local_opts"]][["algorithm"]]))
            ) stop("wrong algorithm for local solver")

          } else {

            ## specify default algorithm of local solver

            nloptr[["local_opts"]][["algorithm"]] <- "NLOPT_LN_LBFGS"

          }

        } else {

          ## derivative-free outer algorithm

          if("algorithm" %in% names.local.opts){

            ## make sure that algorithm of local solver uses derivatives

            if(
              length(grep("NLOPT_LD", nloptr[["local_opts"]][["algorithm"]]))
            ) stop("wrong algoritm for local solver")

          } else {

            ## specify default algorithm of local solver

            nloptr[["local_opts"]][["algorithm"]] <- "NLOPT_LN_BOBYQA"

          }

        }

        ## default options controlling convergence of local solver

        if(!"xtol_rel" %in% names.local.opts) nloptr[["local_opts"]][["xtol_rel"]] <- -1.
        if(!"ftol_rel" %in% names.local.opts) nloptr[["local_opts"]][["ftol_rel"]] <- 1.e-6

      } else {

        ## local options missing: set default values

        nloptr[["local_opts"]] <- list(
          algorithm = if(length(grep("^NLOPT_GD", nloptr[["algorithm"]])))
          "NLOPT_LD_LBFGS" else "NLOPT_LN_BOBYQA",
          xtol_rel = -1.,
          ftol_rel = 1.e-6
        )

      }

    } else {

      ## omit local_opts for other algorithms

      nloptr <- nloptr[!names(nloptr) %in% "local_opts"]

    }

#### ---- overwrite default options

    ## overwrite options for parameter transformation for alpha and n and tau

    switch(
      wrc_model,
      vg = {
        param_tf[c("alpha", "n")] <- param_transf(
          alpha = "identity", n = "identity"
        )[c("alpha", "n")]
      },
      stop("wrc model '", wrc_model, "' not implemented")
    )

    switch(
      hcc_model,
      vgm = {
        param_tf[c("tau")] <- param_transf(tau = "identity")[c("tau")]
      },
      stop("hcc model '", hcc_model, "' not implemented")
    )

    ## overwrite default options controlling convergence

    if(!"xtol_rel" %in% names.opts) nloptr[["xtol_rel"]] <- -1.
    if(!"ftol_rel" %in% names.opts) nloptr[["ftol_rel"]] <- -1.
    if(!"maxtime"  %in% names.opts)  nloptr[["maxtime"]] <- -1.
    if(!"maxeval"  %in% names.opts){
      if(settings == "uglobal"){
        nloptr[["maxeval"]] <- 125
      } else {
        nloptr[["maxeval"]] <- 1000
      }
    }

  } else if(settings %in% c("ulocal", "clocal")){

#### --- options for local NLopt optimization algorithms

    ## specify algorithm

    if("algorithm" %in% names.opts){

      ## check whether specified algorithm does local optimization

      if(length(grep("^NLOPT_L", nloptr[["algorithm"]])) == 0L) stop(
        "global optimization algorithm specified for local optimization"
      )

      if(settings == "clocal" &&
        length(grep("COBYLA$|MMA$|CCSAQ$|SLSQP$|AUGLAG$", nloptr[["algorithm"]])) == 0L
      ) stop(
        "wrong algorithm specified for constrained local optimization"
      )

      if(
        all(param_tf %in% "identity") &&
        identical(nloptr[["algorithm"]], "NLOPT_LN_NEWUOA")
      ){
        warning(
          "Estimating untransformed parameters by algorithm 'NLOPT_LN_NEWUOA'",
          " may results in estimates that violate box constraints"
        )
      }

      ## AUGLAG algorithm: check and specify local options

      if(length(grep("AUGLAG$", nloptr[["algorithm"]]))){

        if("local_opts" %in% names.opts){

          ## local options present

          ## match names of local options

          names.local.opts <- names(nloptr[["local_opts"]])
          names.local.opts <- match.arg(
            names.local.opts, c("algorithm", "xtol_rel", "ftol_rel"),
            several.ok = TRUE
          )
          names(nloptr[["local_opts"]]) <- names.local.opts

          ## match name of algorithm for local solver

          if(length(grep("algorithm", names.local.opts))){

            ## match algorithm of local solver

            nloptr[["local_opts"]][["algorithm"]] <- match.arg(
              toupper(nloptr[["local_opts"]][["algorithm"]]), c(
                "NLOPT_LN_COBYLA", "NLOPT_LD_LBFGS",
                "NLOPT_LD_MMA", "NLOPT_LD_CCSAQ", "NLOPT_LD_SLSQP"
              )
            )

          }

          ## check and specify algorithm for local solver if missing

          if(length(grep("^NLOPT_LD", nloptr[["algorithm"]]))){

            ## outer algorithm uses derivatives

            if("algorithm" %in% names.local.opts){

              ## make sure that algorithm of local solver uses derivatives

              if(
                length(
                  grep("MMA$|CCSAQ$|SLSQP$|LBFGS", nloptr[["local_opts"]][["algorithm"]])
                ) == 0L
              ) stop("wrong algorithm for local solver")

            } else {

              ## specify default algorithm of local solver

              nloptr[["local_opts"]][["algorithm"]] <- "NLOPT_LD_LBFGS"

            }

          } else {

            ## derivative-free outer algorithm

            if("algorithm" %in% names.local.opts){

              ## make sure that algorithm of local solver uses derivatives

              if(
                length(grep("COBYLA$", nloptr[["local_opts"]][["algorithm"]])) == 0L
              ) stop("wrong algoritm for local solver")

            } else {

              ## specify default algorithm of local solver

              nloptr[["local_opts"]][["algorithm"]] <- "NLOPT_LN_COBYLA"

            }

          }

          ## default options controlling convergence of local solver

          if(!"xtol_rel" %in% names.local.opts) nloptr[["local_opts"]][["xtol_rel"]] <- -1.
          if(!"ftol_rel" %in% names.local.opts) nloptr[["local_opts"]][["ftol_rel"]] <- 1.e-6

        } else {

          ## local options missing: set default values

          nloptr[["local_opts"]] <- list(
            algorithm = if(length(grep("^NLOPT_LD", nloptr[["algorithm"]])))
              "NLOPT_LD_LBFGS" else "NLOPT_LN_COBYLA",
            xtol_rel = -1.,
            ftol_rel = 1.e-6
          )

        }

      } else {

        ## omit local_opts for other algorithms

        nloptr <- nloptr[!names(nloptr) %in% "local_opts"]

      }

    } else {

      ## choose default local optimization algorithm

      if(identical(settings, "ulocal")){
        nloptr[["algorithm"]] <- "NLOPT_LN_BOBYQA"
      } else {
        nloptr[["algorithm"]] <- "NLOPT_LD_CCSAQ"
      }

    }

    ## overwrite default options controlling convergence

    if(!"xtol_rel" %in% names.opts) nloptr[["xtol_rel"]] <- -1.
    if(!"ftol_rel" %in% names.opts) nloptr[["ftol_rel"]] <- 1.e-8
    if(!"maxtime"  %in% names.opts)  nloptr[["maxtime"]]  <- -1.
    if(!"maxeval"  %in% names.opts){
      if(settings == "ulocal"){
        nloptr[["maxeval"]] <- 250
      } else {
        nloptr[["maxeval"]] <- 1000
      }
    }

  }

#### -- prepare control arguments of SCEoptim

  if(identical(settings, "sce")){

    ## overwrite options for parameter transformation for alpha, n and tau:

    switch(
      wrc_model,
      vg = {
        param_tf[c("alpha", "n")] <- param_transf(
          alpha = "identity", n = "identity"
        )[c("alpha", "n")]
      },
      stop("wrc model '", wrc_model, "' not implemented")
    )

    switch(
      hcc_model,
      vgm = {
        param_tf[c("tau")] <- param_transf(tau = "identity")[c("tau")]
      },
      stop("hcc model '", hcc_model, "' not implemented")
    )

  }

  ## set flag for gradient method for NLopt optimizers

  use_derivative <- FALSE
  if(length(grep("_LD_|_GD_", nloptr[["algorithm"]]))) use_derivative <- TRUE

  ## set default values for grad_eps criterion

  grad_eps <- switch(
    settings,
    "sce" = ,
    "cglobal" = ,
    "uglobal" = 1.e-2,
    "clocal" = ,
    "ulocal" = if(
      length(grep("_LD_", nloptr[["algorithm"]]))
    ){
        1.e-5
    } else {
        1.e-3
    }
  )

#### -- collect all control arguments

  control <- list(
    settings = settings, hessian = hessian, method = method,
    nloptr = nloptr, sce = sce,
    wrc_model = wrc_model, hcc_model = hcc_model,
    initial_param = initial_param,
    approximation_alpha_k0 = approximation_alpha_k0,
    variable_weight = variable_weight,
    gam_k = gam_k, gam_n_newdata = gam_n_newdata,
    precBits = precBits, delta_sat_0 = 1.e-6,
    min_nobs_wc = min_nobs_wc, min_nobs_hc = min_nobs_hc,
    grad_eps = grad_eps,
    keep_empty_fits = keep_empty_fits,
    param_bound = param_bound,
    param_tf = param_tf, fwd_tf = fwd_tf, deriv_fwd_tfd = deriv_fwd_tfd, bwd_tf = bwd_tf,
    pcmp = pcmp,
    use_derivative = use_derivative
  )

  return(control)

}


# ## ======================================================================
# default.param <- function(
#   thetar = NA_real_, thetas = NA_real_,
#   alpha = NA_real_, n = NA_real_,
#   k0 = NA_real_,
#   tau = 2
# ){
#
#   ## function sets default flags for fitting parameters
#
#   ## 2019-03-22 A. Papritz
#
#   c(
#     thetar = thetar, thetas = thetas,
#     alpha = alpha, n = n,
#     k0 = k0,
#     tau = tau
#   )
#
# }


## ======================================================================
default_fit_param <- function(
  alpha = TRUE, n = TRUE, tau = FALSE,
  thetar = TRUE, thetas = TRUE, k0 = TRUE
){

  ## function sets default flags for fitting parameters

  ## 2019-11-27 A. Papritz
#
  c(
    alpha = alpha, n = n, tau = tau,
    thetar = thetar, thetas = thetas,
    k0 = k0
  )

}


## ======================================================================
### fit_wrc_hcc

fit_wrc_hcc <- function(
  wrc_formula, hcc_formula, data,
  subset = NULL, wrc_subset = subset, hcc_subset = subset,
  weights = NULL, wrc_weights = weights, hcc_weights = weights,
  na.action,
  param = NULL, lower_param = NULL, upper_param = NULL,
  fit_param = default_fit_param(),
  e0 = 2.5e-3,
  ratio_lc_lt_bound  = c(lower = 0.5, upper = 2.),
  control = control_fit_wrc_hcc(),
  verbose = 0
){

  ## user front-end function for estimating parameters of soil
  ## hydraulic functions

  ## 2019-11-27 A. Papritz
  ## 2021-02-27 AP new arguments lower_param, upper_param for specifying
  ##               sample-specific lower and upper parameter bounds
  ## 2021-04-11 AP correction of error when coding box-constraint for single parameter
  ## 2021-05-10 AP correction of error when coding box-constraints only for thetar
  ## 2021-05-22 AP check whether inequality constraints are satisfied
  ##               for constrained estimation
  ## 2021-06-04 AP correction of error when passing param and fit_param to input.data
  ## 2021-12-22 AP warning if lower_param or upper_param are vectors
  ## 2021-12-27 AP check that order of rows in param, fit_param, etc is the same
  ##               as order of samples in output
  ## 2021-12-29 AP sample_id_variable stored as additional component in output
  ##               new method to match input.data and sample-specific data for param,
  ##               fit_param, lower_param, upper_param, e0 ratio_lc_lt_bound

  if(!is.finite(verbose)) browser()

  ## flags for controlling estimation of water retention and hydraulic
  ## conductivity curves

  wrc <- !missing(wrc_formula)
  hcc <- !missing(hcc_formula)

  ## check whether all mandatory arguments have been provided

  if(missing(data) || (!wrc && !hcc)) stop(
    "some mandatory arguments are missing"
  )

  ## get call for later use and for building of model frames

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

#### -- prepare parameters and flags for fitting

  ## match names of param, fit_param, lower_param, upper_param
  ## !! CARE: code works only if fit_param has not been used before (lazy evaluation)

  all.param.name <- names(default_fit_param())

  ## param

  names.param <- NULL
  if(!is.null(param)){
    if(is.data.frame(param)){
      if(!all(sapply(param, is.numeric))) stop(
        "'param' must be a dataframe or a matrix with numeric variables"
      )
      param <- as.matrix(param)
    }
    if(is.vector(param)) param <- t(as.matrix(param))
    names.param <- colnames(param)
    names.param <- sapply(
      names.param,
      function(x, choices) match.arg(x, choices),
      choices = all.param.name
    )
    colnames(param) <- names.param
  }

  ## fit_param

  if(is.data.frame(fit_param)){
    if(!all(sapply(fit_param, is.logical))) stop(
      "'fit_param' must be a dataframe or a matrix with logical variables"
    )
    fit_param <- as.matrix(fit_param)
  }
  if(is.vector(fit_param)) fit_param <- t(as.matrix(fit_param))
  names.fit_param <- colnames(fit_param)
  if(!missing(fit_param)){
    names.fit_param <- sapply(
      names.fit_param,
      function(x, choices) match.arg(x, choices),
      choices = all.param.name
    )
    colnames(fit_param) <- names.fit_param
  }

  ## lower_param

  if(is.data.frame(lower_param)){
    if(!all(sapply(lower_param, is.numeric))) stop(
      "'lower_param' must be a dataframe or a matrix with numeric variables"
    )
    lower_param <- as.matrix(lower_param)
  }
  if(is.vector(lower_param)){
    warning(
      "'lower_param' is a vector; a better way to define parameter bounds ",
      "is to use the function 'param_boundf', see help(control_fit_wrc_hcc)"
    )
    lower_param <- t(as.matrix(lower_param))
  }

  names.lower_param <- colnames(lower_param)
  if(!missing(lower_param)){
    names.lower_param <- sapply(
      names.lower_param,
      function(x, choices) match.arg(x, choices),
      choices = all.param.name
    )
    colnames(lower_param) <- names.lower_param

    ## set lower bounds for thetas equal to lower bounds of thetar if
    ## only the latter have been specified

    if("thetar" %in% names.lower_param & !"thetas" %in% names.lower_param){
      lower_param <- cbind(
        lower_param,
        thetas = lower_param[, "thetar"]
      )
    }

  }

  ## upper_param

  if(is.data.frame(upper_param)){
    if(!all(sapply(upper_param, is.numeric))) stop(
      "'upper_param' must a dataframe or a matrix with numeric variables"
    )
    upper_param <- as.matrix(upper_param)
  }
  if(is.vector(upper_param)){
    warning(
      "'upper_param' is a vector; a better way to define parameter bounds ",
      "is to use the function 'param_boundf', see help(control_fit_wrc_hcc)"
    )
    upper_param <- t(as.matrix(upper_param))
  }

  names.upper_param <- colnames(upper_param)
  if(!missing(upper_param)){
    names.upper_param <- sapply(
      names.upper_param,
      function(x, choices) match.arg(x, choices),
      choices = all.param.name
    )
    colnames(upper_param) <- names.upper_param

    ## set upper bounds for thetar equal to upper bounds of thetas if
    ## only the latter have been specified

    if("thetas" %in% names.upper_param & !"thetar" %in% names.upper_param){
      upper_param <- cbind(
        upper_param,
        thetar = upper_param[, "thetas"]
      )
    }

  }

  #   print(param)
  #   print(fit_param)

  ## consistent coding of flags for fitting

  switch(
    control[["wrc_model"]],
    vg = {
      if(wrc){

        ## nonlinear parameters alpha and n

        if(!all(c("alpha", "n") %in% names.fit_param)){

          if(verbose >= 0.) warning(
            "no fitting control provided for parameter(s) '",
            paste(c("alpha", "n")[ !c("alpha", "n") %in% names.fit_param ], collapse= ", "),
            "': parameters will be fitted"
          )
          if(!"alpha" %in% names.fit_param) fit_param <- cbind(
            alpha = rep(TRUE, NROW(fit_param)), fit_param
          )
          if(!"n" %in% names.fit_param) fit_param <- cbind(
            n = rep(TRUE, NROW(fit_param)), fit_param
          )

        }

        if(!all(fit_param[, "alpha"]) && (is.null(param) || !("alpha" %in% names.param))) stop(
          "no value provided for fixed 'alpha'"
        )
        if(!all(fit_param[, "n"]) && ( !(is.null(param) || "n" %in% names.param))) stop(
          "no value provided for fixed 'n'"
        )

        ## linear parameters thetar and thetas

        if(!all(c("thetar", "thetas") %in% names.fit_param)){

          if(verbose >= 0.) warning(
            "no fitting control provided for parameter(s) '",
            paste(c("thetar", "thetas")[ !c("thetar", "thetas") %in% names.fit_param ], collapse= ", "),
            "': parameters will be fitted"
          )
          if(!"thetar" %in% names.fit_param) fit_param <- cbind(
            thetar = rep(TRUE, NROW(fit_param)), fit_param
          )
          if(!"thetas" %in% names.fit_param) fit_param <- cbind(
            thetas = rep(TRUE, NROW(fit_param)), fit_param
          )
        }

        if(!all(fit_param[, "thetar"]) && (is.null(param) || !("thetar" %in% names.param))) stop(
          "no value provided for fixed 'thetar'"
        )
        if(!is.null(param) && all(fit_param[, "thetar"])){
          param <- param[, !names.param %in% "thetar", drop = FALSE]
          names.param <- names.param[!names.param %in% "thetar"]
        }


        if(!all(fit_param[, "thetas"]) && (is.null(param) || !("thetas" %in% names.param))) stop(
          "no value provided for fixed 'thetas'"
        )
        if(!is.null(param) && all(fit_param[, "thetas"])){
          param <- param[, !names.param %in% "thetas", drop = FALSE]
          names.param <- names.param[!names.param %in% "thetas"]
        }

      }
    },
    stop("wrc model '", control[["wrc_model"]], "' not implemented")
  )


  #   print(param)
  #   print(fit_param)

  switch(
    control[["hcc_model"]],
    vgm = {
      if(hcc){

        ## nonlinear parameters alpha, n and tau

        if(!all(c("alpha", "n", "tau") %in% names.fit_param)){

          if(verbose >= 0.) warning(
            "no fitting control provided for parameter(s) '",
            paste(c("alpha", "n", "tau")[ !c("alpha", "n", "tau") %in% names.fit_param ], collapse= ", "),
            "': parameters will be fitted"
          )
          if(!"alpha" %in% names.fit_param) fit_param <- cbind(
            alpha = rep(TRUE, NROW(fit_param)), fit_param
          )
          if(!"n" %in% names.fit_param)     fit_param <- cbind(
            n = rep(TRUE, NROW(fit_param)), fit_param
          )
          if(!"tau" %in% names.fit_param)   fit_param <- cbind(
            tau = rep(TRUE, NROW(fit_param)), fit_param
          )

        }

        if(!all(fit_param[, "alpha"]) && (is.null(param) || !("alpha" %in% names.param))) stop(
          "no value provided for fixed 'alpha'"
        )
        if(!all(fit_param[, "n"]) && ( !(is.null(param) || "n" %in% names.param))) stop(
          "no value provided for fixed 'n'"
        )
        #     if(!all(fit_param[, "tau"]) && (is.null(param) || !("tau" %in% names.param))) stop(
        #       "no value provided for fixed 'tau'"
        #     )

        ## linear parameter k0

        if(!"k0" %in% names.fit_param){

          if(verbose >= 0.) warning(
            "no fitting control provided for parameter 'k0': parameter will be fitted"
          )
          fit_param <- cbind(k0 = rep(TRUE, NROW(fit_param)), fit_param)

        }

        if(!all(fit_param[, "k0"]) && (is.null(param) || !("k0" %in% names.param))) stop(
          "no value provided for fixed 'k0'"
        )
        if(!is.null(param) && all(fit_param[, "k0"])){
          param <- param[, !names.param %in% "k0", drop = FALSE]
          names.param <- names.param[!names.param %in% "k0"]
        }

      }
    },
    stop("hcc model '", control[["hcc_model"]], "' not implemented")
  )

#   print(param)
#   print(fit_param)

#### -- get variable names for water retention and hydraulic conductivity curves

  wrc_mf <- NULL
  if(wrc){

    ## variable names for wc, head and sample id

    wrc.variables <- as.character(wrc_formula)[-1]
    wrc.variables <- c(
      wrc.variables[1],
      gsub(
        "[\\(\\) ]",
        "",
        unlist(strsplit(wrc.variables[2], "|", fixed = TRUE))
      )
    )

    if(!all(wrc.variables %in% colnames(data))) stop(
      "missing variables for water retention curve(s)"
    )

    ## formula without sample id

    wrc_formula <- formula(
      paste(wrc.variables[1], wrc.variables[2], sep = "~"),
      env = environment(wrc_formula)
    )

    ## preparation for model frame

    m.wrc <- match(
      c(
        "wrc_formula", "subset", "wrc_subset",
        "weights", "wrc_weights", "na.action"
      ),
      names(mf), 0L
    )
    wrc_mf <- mf[c(1L, m.wrc)]

    ## get names of specified arguments

    nmes <- names(wrc_mf)

    ## handle subset and weights arguments

    if("wrc_subset" %in% nmes){
      if("subset" %in% nmes){
        warning(
          "both 'subset' and 'wrc_subset' arguments were specified:",
          "'subset' will be ignored"
        )
        ex <- match("subset", nmes)
        nmes <- nmes[!nmes == "wrc_subset"]
        wrc_mf <- wrc_mf[-ex]
      } else {
        nmes[nmes == "wrc_subset"] <- "subset"
      }
    }

    if("wrc_weights" %in% nmes){
      if("weights" %in% nmes){
        warning(
          "both 'weights' and 'wrc_weights' arguments were specified:",
          "'weights' will be ignored"
        )
        ex <- match("weights", nmes)
        nmes <- nmes[!nmes == "wrc_weights"]
        wrc_mf <- wrc_mf[-ex]
      } else {
        nmes[nmes == "wrc_weights"] <- "weights"
      }
    }

    ## rename wrc_formula and wrc_weights

    nmes[nmes == "wrc_formula"] <- "formula"

    ## set names of specified argument of function call

    names(wrc_mf) <- nmes

    ## set remaing argument of call and prepare evalution of model frame

    wrc_mf[["formula"]] <- wrc_formula
    wrc_mf[["drop.unused.levels"]] <- TRUE
    wrc_mf[["data"]] <- as.name("input.data")
    wrc_mf[[1L]] <- as.name("model.frame")

  } else {
    wrc_formula <- NULL
  }

  hcc_mf <- NULL
  if(hcc){

    ## variable names for hydraulic conductivity, head and sample id

    hcc.variables <- as.character(hcc_formula)[-1]
    hcc.variables <- c(
      hcc.variables[1],
      gsub(
        "[\\(\\) ]",
        "",
        unlist(strsplit(hcc.variables[2], "|", fixed = TRUE))
      )
    )

    if(!all(hcc.variables %in% colnames(data))) stop(
      "missing variables for hydraulic conductivity curve(s)"
    )

    if(wrc){

      if(!identical(length(wrc.variables), length(hcc.variables))) stop(
        "inconsistent formulae for water retention and hydraulic conductivity curves"
      )

      if(
        identical(length(wrc.variables), 3L) &&
        !identical(wrc.variables[3], hcc.variables[3])
      ) stop(
        "variable that defines samples differs for water retention and hydraulic conductivity curves"
      )

    }

    ## formula without sample id

    hcc_formula <- formula(
      paste(hcc.variables[1], hcc.variables[2], sep = "~"),
      env = environment(hcc_formula)
    )

    ## preparation for model frame

    m.hcc <- match(
      c(
        "hcc_formula", "subset", "hcc_subset",
        "weights", "hcc_weights", "na.action"
      ),
      names(mf), 0L
    )
    hcc_mf <- mf[c(1L, m.hcc)]

    ## get names of specified arguments

    nmes <- names(hcc_mf)

    ## handle subset and weights arguments

    if("hcc_subset" %in% nmes){
      if("subset" %in% nmes){
        warning(
          "both 'subset' and 'hcc_subset' arguments were specified:",
          "'subset' will be ignored"
        )
        ex <- match("subset", nmes)
        nmes <- nmes[!nmes == "hcc_subset"]
        hcc_mf <- hcc_mf[-ex]
      } else {
        nmes[nmes == "hcc_subset"] <- "subset"
      }
    }

    if("hcc_weights" %in% nmes){
      if("weights" %in% nmes){
        warning(
          "both 'weights' and 'hcc_weights' argument specified:",
          "'weights' will be ignored"
        )
        ex <- match("weights", nmes)
        nmes <- nmes[!nmes == "hcc_weights"]
        hcc_mf <- hcc_mf[-ex]
      } else {
        nmes[nmes == "hcc_weights"] <- "weights"
      }
    }

    ## rename hcc_formula and hcc_weights

    nmes[nmes == "hcc_formula"] <- "formula"

    ## set names of specified argument of function call

    names(hcc_mf) <- nmes

    ## set remaing argument of call and prepare evalution of model frame

    hcc_mf[["formula"]] <- hcc_formula
    hcc_mf[["drop.unused.levels"]] <- TRUE
    hcc_mf[["data"]] <- as.name("input.data")
    hcc_mf[[1L]] <- as.name("model.frame")

  } else {
    hcc_formula <- NULL
  }

#### -- prepare ratio_lc_lt_bound

  if(is.vector(ratio_lc_lt_bound)){
    tmp <- names(ratio_lc_lt_bound)
    ratio_lc_lt_bound <- matrix(ratio_lc_lt_bound, nrow = 1)
    colnames(ratio_lc_lt_bound) <- tmp
  }

  ## check consistency of bounds on ratio Lc/Lt

  if(any(apply(ratio_lc_lt_bound, 1, function(x) x[1] > x[2]))) stop(
    "lower bound of ratio Lc/Lt > upper bound"
  )

#### -- create list of dataframes where each component contains data for a
  ## single sample

  if(wrc){
    split.variables <- wrc.variables
  } else if(hcc){
    split.variables <- hcc.variables
  }

  if(identical(length(split.variables), 3L)){
    input.data <- split(data, data[, split.variables[3]])
  } else{
    input.data <- list(data)
  }

  if(identical(length(input.data), 1L)) names(input.data) <- ""

  nmes.samples <- names(input.data)

  ## check that param, fit_param, lower_param, upper_param, e0 and
  ## ratio_lc_lt_bound have valid (row)names if they contain
  ## sample-specific information, i.e. when they are matrices or a vector

  if(length(input.data) > 1L){

    if(!is.null(param) && NROW(param) > 1L){
      if(is.null(rownames(param))) stop(
        "'param' is a matrix with sample-specific data ",
        "but rownames are missing"
      )
      if(!all(nmes.samples %in% rownames(param))) stop(
        "'param' lacks data for some samples"
      )
    }

    if(NROW(fit_param) > 1L){
      if(is.null(rownames(fit_param))) stop(
        "'fit_param' is a matrix with sample-specific data ",
        "but rownames are missing"
      )
      if(!all(nmes.samples %in% rownames(fit_param))) stop(
        "'fit_param' lacks data for some samples"
      )
    }

    if(!is.null(lower_param) && NROW(lower_param) > 1L){
      if(is.null(rownames(lower_param))) stop(
        "'lower_param' is a matrix with sample-specific data ",
        "but rownames are missing"
      )
      if(!all(nmes.samples %in% rownames(lower_param))) stop(
        "'lower_param' lacks data for some samples"
      )
    }

    if(!is.null(upper_param) && NROW(upper_param) > 1L){
      if(is.null(rownames(upper_param))) stop(
        "'upper_param' is a matrix with sample-specific data ",
        "but rownames are missing"
      )
      if(!all(nmes.samples %in% rownames(upper_param))) stop(
        "'upper_param' lacks data for some samples"
      )
    }

    if(length(e0) > 1L){
      if(is.null(names(e0))) stop(
        "'e0' is a vector with sample-specific data ",
        "but names are missing"
      )
      if(!all(nmes.samples %in% names(e0))) stop(
        "'e0' lacks data for some samples"
      )
    }

    if(NROW(ratio_lc_lt_bound) > 1L){
      if(is.null(rownames(ratio_lc_lt_bound))) stop(
        "'ratio_lc_lt_bound' is a matrix with sample-specific data ",
        "but rownames are missing"
      )
      if(!all(nmes.samples %in% rownames(ratio_lc_lt_bound))) stop(
        "'ratio_lc_lt_bound' lacks data for some samples"
      )
    }

  }

  ## augment list of dataframes with param, fit_param, lower_param,
  ## upper_param, e0 and ratio_lc_lt_bound components

  input.data <- lapply(
    1:length(input.data),
    function(i){

      param.ii <- NULL
      if(!is.null(param)){
        if(NROW(param) > 1L){
          ii <- match(nmes.samples[i], rownames(param))
        } else {
          ii <- 1L
        }
        param.ii <- unname(param[ii, ])
        names(param.ii) <- colnames(param)
      }

      if(NROW(fit_param) > 1L){
        ii <- match(nmes.samples[i], rownames(fit_param))
      } else {
        ii <- 1L
      }
      fit_param.ii <- unname(fit_param[ii, ])
      names(fit_param.ii) <- colnames(fit_param)

      lower_param.ii <- NULL
      if(!is.null(lower_param)){
        if(NROW(lower_param) > 1L){
          ii <- match(nmes.samples[i], rownames(lower_param))
        } else {
          ii <- 1L
        }
        lower_param.ii <- unname(lower_param[ii, ])
        names(lower_param.ii) <- colnames(lower_param)
      }

      upper_param.ii <- NULL
      if(!is.null(upper_param)){
        if(NROW(upper_param) > 1L){
          ii <- match(nmes.samples[i], rownames(upper_param))
        } else {
          ii <- 1L
        }
        upper_param.ii <- unname(upper_param[ii, ])
        names(upper_param.ii) <- colnames(upper_param)
      }

      if(length(e0) > 1L){
        ii <- match(nmes.samples[i], names(e0))
      } else {
        ii <- 1L
      }
      e0.ii <- unname(e0[ii])

      if(NROW(ratio_lc_lt_bound) > 1L){
        ii <- match(nmes.samples[i], rownames(ratio_lc_lt_bound))
      } else {
        ii <- 1L
      }
      ratio_lc_lt_bound.ii <- unname(ratio_lc_lt_bound[ii, ])
      names(ratio_lc_lt_bound.ii) <- colnames(ratio_lc_lt_bound)

      list(
        data = input.data[[i]],
        param = param.ii,
        fit_param = fit_param.ii,
        lower_param = lower_param.ii,
        upper_param = upper_param.ii,
        e0 = e0.ii,
        ratio_lc_lt_bound = ratio_lc_lt_bound.ii
      )
    }
  )

  names(input.data) <- nmes.samples

#### -- preparations for parallel processing

  ## adjust number of cores for parallel processing

  control[["pcmp"]][["ncores"]] <- min(
    control[["pcmp"]][["ncores"]],
    length(input.data)
  )

  ## diagnostic output for parallel processing

  if(control[["pcmp"]][["ncores"]] > 1 && verbose >=1. ) verbose <- 0.5

  ## create environment for storing index of iteration, linear parameter
  ## and value of objective function

  common.env <- new.env()

  ## estimate VGM parameters for all samples in parallel

  ## auxiliary function to fit parameters of single sample

  f.aux <- function(i){

    ## note that input.data, wrc, wrc_formula, wrc_mf, hcc, hcc_formula, hcc_mf,
    ## all.param.name, control, verbose, split.variables
    ## are taken from parent environment

    ## extract data, param, fit_param, e0, ratio_lc_lt_bound

    input.data.i  <- input.data[[i]][["data"]]
    param.i       <- input.data[[i]][["param"]]
    fit_param.i   <- input.data[[i]][["fit_param"]]
    lower_param.i <- input.data[[i]][["lower_param"]]
    upper_param.i <- input.data[[i]][["upper_param"]]
    e0.i          <- input.data[[i]][["e0"]]
    ratio_lc_lt_bound.i <- input.data[[i]][["ratio_lc_lt_bound"]]

    if(
      verbose >= 2. && identical(length(split.variables), 3L)) cat(
      "\n\nprocessing data of sample",
      as.character(unique(input.data.i[, split.variables[3]])), "\n"
    )

    ## update default initial values of nonlinear parameter in
    ## control[["initial.param"]] if lie outside of valid range

    lupn <- unique(
      c(names(lower_param.i), names(upper_param.i))
    )

    if(length(lupn)){

      ## update parameter boundaries in control[["param_bound"]] and
      ## initial values in control[["initial_param"]]

      if(!is.null(lower_param.i)){
        tmp <- control[["param_bound"]][names(lower_param.i)]
        tmp <- lapply(
          names(tmp),
          function(i) c(lower_param.i[i], tmp[[i]][2])
        )
        control[["param_bound"]][names(lower_param.i)] <- unname(tmp)
      }

      if(!is.null(upper_param.i)){
        tmp <- control[["param_bound"]][names(upper_param.i)]
        tmp <- lapply(
          names(tmp),
          function(i) c(tmp[[i]][1],upper_param.i[i])
        )
        control[["param_bound"]][names(upper_param.i)] <- unname(tmp)
      }

      ## check validity of new parameter boundaries

      check_param_bound <- try(
        check_param_boundf(control[["param_bound"]]),
        silent = TRUE
      )

      if(identical(class(check_param_bound), "try-error")){

        ## some boundaries are invalid: return fit with an error

        message <- paste(
          "an error occurred when estimating parameters",
          if(identical(length(split.variables), 3L)) paste(
            "for sample", as.character(unique(input.data.i[, split.variables[3]]))
          ), "\n", as.character(check_param_bound)
        )
        warning(message)
        fit <- list(
          converged = NULL,
          message = message,
          objective = NULL,
          gradient = NULL,
          evaluation = NULL,
          lp = NULL, nlp = NULL,
          inequality_constraint = NULL,
          variable_weight = NULL,
          case.weight = NULL,
          model = NULL,
          initial_objects = NULL
        )

        return(fit)

      }

      ## update initial values of nonlinear parameters

      ipn <- names(control[["initial_param"]])

      bla <- switch(
        control[["hcc_model"]],
        vgm = {
          if(any(c("alpha", "n", "tau") %in% lupn)){
            sel.param <- ipn[ipn %in% lupn]
            if(length(sel.param)){
              result <- sapply(
                sel.param,
                function(xn){
                  x  <- control[["initial_param"]][xn]
                  xb <- control[["param_bound"]][[xn]]
                  unname(ifelse(
                    x < xb[1], xb[1] + 0.01 * abs(xb[1]),
                    ifelse(
                      x > xb[2], xb[2] - 0.01 * abs(xb[2]),
                      x
                    )
                  ))
                }
              )
              control[["initial_param"]][sel.param] <- result
            }
          }
        },
        stop("hcc model '", control[["hcc_model"]], "' not implemented")
      )

    }

    ## fit models

    fit <- try(
      fit_wrc_hcc_fit(
        i,
        input.data.i, param.i, fit_param.i,
        e0.i, ratio_lc_lt_bound.i,
        wrc, wrc_formula, wrc_mf,
        hcc, hcc_formula, hcc_mf,
        all.param.name,
        control,
        common.env,
        verbose
      )
      , silent = TRUE
    )

    if(identical(class(fit), "try-error")){

      ## model fit failed

      message <- paste(
        "an error occurred when estimating parameters",
        if(identical(length(split.variables), 3L)) paste(
          "for sample", as.character(unique(input.data.i[, split.variables[3]]))
        ), "\n", as.character(fit)
      )
      warning(message)
      fit <- list(
        converged = NULL,
        message = message,
        objective = NULL,
        gradient = NULL,
        evaluation = NULL,
        lp = NULL, nlp = NULL,
        inequality_constraint = NULL,
        variable_weight = NULL,
        case.weight = NULL,
        model = NULL,
        initial_objects = NULL
      )

    }

    fit

  }

#### -- loop over all samples

  if(
    control[["pcmp"]][["ncores"]] > 1L &&
    !control[["pcmp"]][["fork"]]
  ){

    ## create a SNOW cluster on windows OS

    options(error = stop_cluster)
    junk <- sfInit( parallel = TRUE, cpus = control[["pcmp"]][["ncores"]] )

    ## export required items to workers (debugging)

#     junk <- sfLibrary("quadprog", verbose = FALSE, character.only = TRUE)
#     junk <- sfLibrary("mgcv", verbose = FALSE, character.only = TRUE)
#     junk <- sfLibrary("nloptr", verbose = FALSE, character.only = TRUE)
#     junk <- sfLibrary("Rmpfr", verbose = FALSE, character.only = TRUE)
#     junk <- sfLibrary("SoilHyP", verbose = FALSE, character.only = TRUE)

    junk <- sfLibrary("soilhypfit", verbose = FALSE, character.only = TRUE)

    ## export required items to workers (if package is installed)

    #     junk <- sfExportAll()

    export.items <- c(
      "input.data", "split.variables",
      "wrc", "wrc_formula", "wrc_mf",
      "hcc", "hcc_formula", "hcc_mf",
      "all.param.name",
      "control",
      "common.env",
      "verbose"
    )

    junk <- sfExport(list = export.items)

    result <- sfLapply(1L:length(input.data), f.aux)

    junk <- stop_cluster()

  } else {

    ## fork child processes on non-windows OS

    result <- mclapply(
      1L:length(input.data),
      f.aux,
      mc.cores = control[["pcmp"]][["ncores"]],
      mc.allow.recursive = FALSE
    )

  }

#### -- prepare output object

  if(missing(na.action)) na.action <- NULL

  names(result) <- nmes.samples

  ## drop empty fits

  if(!control[["keep_empty_fits"]]){
    ex <- sapply(
      result,
      function(x) all(
        sapply(
          x[!names(x) %in% "message"],
          function(y) is.null(y)
        )
      )
    )
    if(identical(sum(ex), length(result))){
      cat(
        "\n\nonly 'empty' fits with following error messages:\n\n"
      )
      bla <- lapply(result, function(x) cat(x[["message"]]))

    }
    result <- result[!ex]
  }

  ## check wether inequality constraints are satified for
  #  constrained estimation algorithms

  if(control[["settings"]] %in% c("cglobal", "clocal")){

    not.ok <- sapply(
      result,
      function(x){
        any(
          sapply(
            x[["inequality_constraint"]],
            function(y) any(y > 0.00001)
          )
        )
      }
    )

    if(any(not.ok)) warning(
      "inequality constraints not satisfied for some samples"
    )

  }

  ## create list and set class attribute

  if(identical(length(split.variables), 3L)){
    tmp <- split.variables[3]

  } else {
    tmp <- NULL
  }

  result <- list(
    fit = result,
    sample_id_variable = tmp,
    wrc = wrc, wrc_formula = wrc_formula, wrc_mf = wrc_mf,
    hcc = hcc, hcc_formula = hcc_formula, hcc_mf = hcc_mf,
    control = control, call = cl, na.action = na.action
  )

  class(result) <- "fit_wrc_hcc"

  invisible(result)

}


## ======================================================================
convergence_message <- function(x, sce = FALSE){

  ## function prints messages corresponding to covergence codes of nloptr
  ## und SCEoptim

  ## 2019-11-27 A. Papritz

  xx <- as.character(x)

  if(sce){

    cat(switch(
        xx,
        "0" = "Change in solution over [tolsteps] less than specified tolerance (reltol).",
        "1" = "Maximum number of function evaluations or iterations reached.",
        "2" = "Exceeded maximum time.",
        "Unknown convergence code"
      ))

  } else {

    cat(switch(
        xx,
        "1" = "Generic success return value.",
        "2" = "Optimization stopped because stopval was reached.",
        "3" = "Optimization stopped because ftol_rel or ftol_abs was reached.",
        "4" = "Optimization stopped because xtol_rel or xtol_abs was reached.",
        "5" = "Optimization stopped because maxeval was reached.",
        "6" = "Optimization stopped because maxtime was reached.",
        "-1" = "Generic failure code.",
        "-2" = "Invalid arguments (e.g. lower bounds are bigger than upper bounds,\n  an unknown algorithm was specified, etc).",
        "-3" = "Ran out of memory.",
        "-4" = "Halted because roundoff errors limited progress. (In this case,\n  the optimization still typically returns a useful result).",
        "-5" = "Halted because of a forced termination.",
        "Unknown convergence code"
      ))

  }

  invisible(x)

}


## ======================================================================
select_failed_fits <- function(object){

  ## function finds soil samples for which the estimaton of parameter
  ## failed and returns an integer vector with the indices of those samples

  ## 2019-12-03 A. Papritz

  ## check consistency of object


  sel <- sapply(object[["fit"]], function(x) is.null(x[["objective"]]))

  names(object[["fit"]])[sel]

}


## ======================================================================
extract_error_messages <- function(object, start = 1, stop = 80, prt = TRUE){

  ## function extracts error messages of failed fits and prints a substring
  ## of them

  ## 2019-12-03 A. Papritz

  ## check consistency of object

  stopifnot(identical(class(object), "fit_wrc_hcc"))

  sel <- select_failed_fits(object)

  msg <- sapply(object[["fit"]][sel], function(x) x[["message"]])

  if(prt) cat(substr(msg, start, stop))

  invisible(msg)

}


## ======================================================================
### prfloglik_sample

prfloglik_sample <- function(
  object, values, soil_sample,
  ncores = min(detectCores() - 1L, NROW(values)),
  verbose = 0
){

  ## function to compute loglikelihood profile for a single sample of a
  ## fit_wrc_hcc object function is based on profilelogLik{georob}

  ## 2021-12-20 Andreas Papritz
  ## 2021-12-22 AP correction of errors for parallel computations on windows
  ## 2021-12-28 AP small changes
  ## 2021-12-31 AP new version calling directly fit_wrc_hcc_fit
  ## 2022-01-06 AP adjust values of thetar and thetas if outside 
  ##               of allowed bounds

  if(!is.finite(verbose)) browser()

#### -- auxiliary function

  ## auxiliary function to fit model and return maximized log-likelihood
  ## and conditionally estimated parameters

  f.aux <- function(i){

    ## values, soil_sample, input.data, param, fit_param, e0, ratio_lc_lt_bound,
    ## object, common.env, verbose are taken from parent environment

    ## set fixed initial values

    param[colnames(values)] <- unname(values[i, ])

    ## estimate parameters

    fit <- try(
      fit_wrc_hcc_fit(
        i = soil_sample,
        input.data = input.data,
        param = param, fit_param = fit_param,
        e0 = e0, ratio_lc_lt_bound = ratio_lc_lt_bound,
        wrc = object[["wrc"]], wrc_formula = object[["wrc_formula"]],
        wrc_mf = object[["wrc_mf"]],
        hcc = object[["hcc"]], hcc_formula = object[["hcc_formula"]],
        hcc_mf = object[["hcc_mf"]],
        all.param.name = names(default_fit_param()),
        control = object[["control"]],
        common.env = common.env,
        verbose = verbose
      )
      , silent = TRUE
    )

    ## result

    if(identical(class(fit), "try-error")){

      param[names(fit_param)[fit_param]] <- NA_real_
      gradient <- rep(NA_real_, length(object[["fit"]][[1]][["nlp"]]))
      names(gradient) <- names(object[["fit"]][[1]][["nlp"]])
      gradient <- gradient[names(gradient) %in% names(param)[is.na(param)]]
      names(gradient) <- paste0("gradient.", names(gradient))

      result <- c(
        loglik = NA_real_,
        param,
        gradient,
        converged = NA_real_
      )

    } else {

      result <- c(
        loglik = - fit[["objective"]],
        c(fit[["nlp"]], fit[["lp"]]),
        gradient = fit[["gradient"]],
        converged = fit[["converged"]]
      )

    }

    result

  }


#### -- check arguments and contents of object and values

  if(missing(object) || missing(values)) stop(
    "some mandatory arguments are missing"
  )

  stopifnot(identical(class(object)[1], "fit_wrc_hcc"))
  stopifnot(identical(length(verbose), 1L) && is.numeric(verbose) && verbose >= 0)
  stopifnot(identical(length(ncores), 1L) && is.numeric(ncores) && ncores >= 1)

  if(!(is.matrix(values) || is.data.frame(values))) stop(
    "'values' must be a dataframe or a matrix"
  )

  if(length(object[["fit"]]) > 1L){
    if(missing(soil_sample)) stop(
      "argument 'soil_sample' must be a character scalar because 'object' contains ",
      " parameter estimates for several samples"
    )
    if(length(soil_sample) > 1L || !is.character(soil_sample)) stop(
      "'soil_sample' must be a character scalar"
    )
  } else {
    if(!missing(soil_sample)){
      warning(
        "argument 'soil_sample' will be ignored because 'object' contains ",
        "parameter estimates for just one sample"
      )
    }
    soil_sample <- ""
  }

  if(!identical(object[["control"]][["method"]], "ml")) stop(
    "to compute profile loglikelihood, parameters must be estimated by ",
    "maximum likelihood method using control argument 'method = ml'"
  )

  if(object[["control"]][["settings"]] %in% c("uglobal", "cglobal", "sce")) warning(
    "parameters have been estimated by global algorithm"
  )

#### -- prepare objects required for call of fit_wrc_hcc_fit

  ## select component of object$fit for sample chosen by soil_sample

  if(length(object[["fit"]]) > 1L){
    if(soil_sample %in% names(object[["fit"]])){
      object[["fit"]] <- object[["fit"]][soil_sample]
    } else stop(
      "parameter estimates missing in 'object' missing for sample ",
      soil_sample
    )
  }

  ## check whether selected component of object$fit is valid

  if(is.null(object[["fit"]][[1]])) stop(
    "parameter estimation failed for sample ", soil_sample
  )

  ## update components wrc and hcc in object$fit[[1]] in dependence of
  ## available data

  ## no wrc data for selected sample

  if(object[["wrc"]] && is.null(object[["fit"]][[1]][["model"]][["wrc"]])){
    object[["wrc"]] <- FALSE
  }

  ## no hcc data for selected sample

  if(object[["hcc"]] && is.null(object[["fit"]][[1]][["model"]][["hcc"]])){
    object[["hcc"]] <- FALSE
  }

  ## get input.data, param, fit_param, e0, ratio_lc_bound for selected
  ## sample

  input.data <- NULL
  if(object[["wrc"]]){
    input.data <- object[["fit"]][[1]][["model"]][["wrc"]]
  }
  if(object[["hcc"]]){
    if(is.null(input.data)){
      input.data <- object[["fit"]][[1]][["model"]][["hcc"]]
    } else {
      input.data <- merge(
        input.data, object[["fit"]][[1]][["model"]][["hcc"]],
        all = TRUE
      )
    }
  }

  param <- c(
    object[["fit"]][[1]][["nlp"]], object[["fit"]][[1]][["lp"]]
  )

  fit_param <- object[["fit"]][[1]][["initial_objects"]][["fit_param"]]

  param_bound <- object[["fit"]][[1]][["initial_objects"]][["param_bound"]]

  e0 <- object[["fit"]][[1]][["initial_objects"]][["e0"]]

  ratio_lc_lt_bound <- object[["fit"]][[1]][["initial_objects"]][["ratio_lc_lt_bound"]]

  ## adjust values for thetar and thetas if estimates are outside of
  ## allowed bounds
  
  if("thetar" %in% names(param)){
    param["thetar"] <- min(
      max(param["thetar"], param_bound[["thetar"]][1]), 
      param["thetas"]
    )
  }
  
  if("thetas" %in% names(param)){
    param["thetas"] <- max(
      min(param["thetas"], param_bound[["thetas"]][2]), 
      param["thetar"]
    )
  }
  

  ## omit subset argument from modelframes

  if(!is.null(object[["wrc_mf"]])){
    x <- as.list(object[["wrc_mf"]])
    x <- x[!names(x) %in% "subset"]
    object[["wrc_mf"]] <- as.call(x)
  }

  if(!is.null(object[["hcc_mf"]])){
    x <- as.list(object[["hcc_mf"]])
    x <- x[!names(x) %in% "subset"]
    object[["hcc_mf"]] <- as.call(x)
  }

  ## check names of fixed parameters in values

  fixed.param <- colnames(values)
  fixed.param <- unname(sapply(
    fixed.param, match.arg,
    choices = c(names(default_fit_param()))
  ))
  if(any(!fixed.param %in% names(param))) stop(
    "column names of 'values' do not match names of parameters in 'object'"
  )
  colnames(values) <- fixed.param

  ## check whether values of fixed parameters are within allowed ranges

  minmax.values <- lapply(values, range)

  check_param_boundf(minmax.values, compare_thetar_thetas = FALSE)

  bla <- lapply(
    names(minmax.values),
    function(i){
      if(
        minmax.values[[i]][1] < param_bound[[i]][1] ||
        minmax.values[[i]][2] > param_bound[[i]][2]
      ) warning(
        "extrema of parameter '", i, "' in 'values' outside of allowed range"
      )
    }
  )

  ## change fit_param to FALSE for variables present in values

  fit_param[fixed.param] <- FALSE

  ## set ncores = 1L in object$control$pcmp

  object[["control"]][["pcmp"]][["ncores"]] <- 1L

  ## set hessian = FALSE in object$control

  object[["control"]][["hessian"]] <- FALSE

  ## create environment for storing index of iteration, linear parameter
  ## and value of objective function

  common.env <- new.env()


#### -- fit model for sets of parameter values

  ## loop over all elements of values

  values <- as.matrix(values)

  ## set default value for control of forking if missing (required for
  ## backward compatibility)

  if(
    ncores > 1L &&
    !object[["control"]][["pcmp"]][["fork"]]
  ){

    ## create a SNOW cluster on windows OS

    options(error = stop_cluster)
    junk <- sfInit( parallel = TRUE, cpus = ncores )

    ## export required items to workers

    junk <- sfLibrary("soilhypfit", verbose = FALSE, character.only = TRUE)

    ## export required items to workers (if package is installed)

    #     junk <- sfExportAll()
    #     export.items <- c("values", "object", "data", "change_function_call_set_onexxx_to_value")

    export.items <- c(
      "values",
      "soil_sample",
      "input.data", "param", "fit_param", "e0", "ratio_lc_lt_bound",
      "object", "common.env", "verbose"
    )

    junk <- sfExport(list = export.items)

    result <- sfLapply(1L:NROW(values), f.aux)

    junk <- stop_cluster()

  } else {

    ## fork child processes on non-windows OS

    result <- mclapply(
      1L:NROW(values),
      f.aux,
      mc.cores = ncores,
      mc.allow.recursive = FALSE
    )

  }


#### -- prepare output

  ## collect results

  result <- t(simplify2array(result))
  result <- result[, !colnames(result) %in% colnames(values), drop = FALSE]

  as.data.frame(cbind(values, result))

}


## ======================================================================
### confint_prfloglik_sample

confint_prfloglik_sample <- function(
  object, parm = names(default_fit_param()), soil_sample,
  level = 0.95, test = c("F", "Chisq"),
  denominator_df = c("nonlinear", "all"),
  param_bound = NULL,
  root_tol = .Machine$double.eps^0.25,
  froot_tol = sqrt(root_tol),
  verbose = 0
){

  ## function to compute lconfidence interval of a parameter for a single
  ## fit_wrc_hcc fit based on likelihood ratio test

  ## 2021-12-29 Andreas Papritz
  ## 2022-01-02 AP new version based on new prfloglik_sample calling
  ##               directly fit_wrc_hcc_fit
  ## 2022-01-03 AP new criterion to decide whether root has been found

  if(!is.finite(verbose)) browser()

#### -- auxiliary function with zero roots for confidence limits
  f.aux <- function(
    x, parm, object, qtest, maximized_loglik
  ){

    ## compute profile loglikelihood for parm = x

    values <- data.frame(x = x)
    colnames(values) <- parm

    tmp <- prfloglik_sample(
      object, values = values, ncores = 1L, verbose = verbose
    )

    ## compute and function value

    qtest - (maximized_loglik - tmp[1, "loglik"])

  }

#### -- check arguments and contents of object

  parm <- match.arg(parm)
  test  <- match.arg(test)
  denominator_df <- match.arg(denominator_df)

  if(missing(object)) stop(
    "some mandatory arguments are missing"
  )

  stopifnot(identical(class(object)[1], "fit_wrc_hcc"))
  stopifnot(identical(length(verbose), 1L) && is.numeric(verbose) && verbose >= 0)
  stopifnot(identical(length(level), 1L) && is.numeric(level) && level >= 0 & level <= 1)

  if(length(object[["fit"]]) > 1L){
    if(missing(soil_sample)) stop(
      "argument 'soil_sample' must be a character scalar because 'object' contains ",
      " parameter estimates for several samples"
    )
    if(length(soil_sample) > 1L || !is.character(soil_sample)) stop(
      "'soil_sample' must be a character scalar"
    )
  } else {
    if(!missing(soil_sample)){
      warning(
        "argument 'soil_sample' will be ignored because 'object' contains ",
        "parameter estimates for just one sample"
      )
    }
    soil_sample <- ""
  }

  if(!identical(object[["control"]][["method"]], "ml")) stop(
    "to compute profile loglikelihood, parameters must be estimated by ",
    "maximum likelihood method using control argument 'method = ml'"
  )

  if(object[["control"]][["settings"]] %in% c("uglobal", "cglobal", "sce")) warning(
    "parameters have been estimated by global algorithm"
  )

  if(!is.null(param_bound)){
    stopifnot(is.numeric(param_bound) && identical(length(param_bound), 2L))
    tmp <- list(param_bound)
    names(tmp) <- parm
    check_param_boundf(tmp)
  }


#### -- prepare objects for computing confidence interval

  ## select component of object$fit for sample chosen by soil_sample

  if(length(object[["fit"]]) > 1L){
    if(soil_sample %in% names(object[["fit"]])){
      object[["fit"]] <- object[["fit"]][soil_sample]
    } else stop(
      "parameter estimates missing in 'object' missing for sample ",
      soil_sample
    )
  }

  ## check whether selected component of object$fit is valid

  if(is.null(object[["fit"]][[1]])) stop(
    "parameter estimation failed for sample ", soil_sample
  )

  ## extract maximized loglik, parameter estimate and boundaries of parm

  tmp <- coef(object, gof = TRUE)
  param_estimate <- unname(tmp[, parm])
  maximized_loglik <- -unname(tmp[, "obj"])

  if(is.null(param_bound)){
    param_bound <- unname(
      object[["fit"]][[1]][["initial_objects"]][["param_bound"]][[parm]]
    )
  }

  ## determine number of measurements of wrc and/or hcc

  if(!is.null(object[["fit"]][[1]][["model"]][["wrc"]])){
    nobs_wrc <- NROW(object[["fit"]][[1]][["model"]][["wrc"]])
  } else {
    nobs_wrc <- 0L
  }

  if(!is.null(object[["fit"]][[1]][["model"]][["hcc"]])){
    nobs_hcc <- NROW(object[["fit"]][[1]][["model"]][["hcc"]])
  } else {
    nobs_hcc <- 0L
  }
  nobs <- nobs_wrc + nobs_hcc

  ## check type of test given type of data

  if(identical(test, "F") && nobs_wrc * nobs_hcc > 0){
    test <- "Chisq"
    warning(
      "Chisq-Test will be used because both water retention and", "
      hydraulic conductivity data were used to estimate the parameters"
    )
  }

  ## determine how many (nonlinear) parameters were fitted

  nme_param <- switch(
    denominator_df,
    nonlinear = names(object[["fit"]][[1]][["nlp"]]),
    all = c(
      names(object[["fit"]][[1]][["nlp"]]),
      names(object[["fit"]][[1]][["lp"]])
    ),
    stop("unknown value for 'denominator_df'")
  )
  n_nlp <- sum(object[["fit"]][[1]][["initial_objects"]][["fit_param"]][nme_param])

  ## compute (transformed) quantile of chi2- or f-distribution

  qtest <- switch(
    test,
    "Chisq" = 0.5 * qchisq(level, 1),
    "F" = 0.5 * nobs * log(1 + 1/(nobs - n_nlp) * qf(level, 1, nobs - n_nlp)),
    stop("unknown test")
  )


### -- compute limits of confidence interval

  ## lower limit

  flower <- f.aux(
    param_bound[1], parm, object, qtest, maximized_loglik
  )
  fupper <- f.aux(
    param_estimate, parm, object, qtest, maximized_loglik
  )

  if(identical(sign(flower * fupper), -1.)){

    t.lower <- uniroot(
      f.aux,
      lower = param_bound[1], upper = param_estimate,
      f.lower = flower, f.upper = fupper, tol = root_tol,
      parm = parm, object = object, qtest = qtest,
      maximized_loglik = maximized_loglik
    )

    if(abs(t.lower[["f.root"]]) < froot_tol ){
      lower_limit <- t.lower[["root"]]
    } else {
      warning("failure to compute lower confidence limit")
      lower_limit <- NA_real_
    }
    lower_froot <- t.lower[["f.root"]]

  } else {
    lower_limit <- lower_froot <- NA_real_
  }

  ## upper limit

  flower<- f.aux(
    param_estimate, parm, object, qtest, maximized_loglik
  )
  fupper <- f.aux(
    param_bound[2], parm, object, qtest, maximized_loglik
  )

  if(identical(sign(flower * fupper), -1.)){

    t.upper <- uniroot(
      f.aux,
      upper = param_bound[2], lower = param_estimate,
      f.lower = flower, f.upper = fupper, tol = root_tol,
      parm = parm, object = object, qtest = qtest,
      maximized_loglik = maximized_loglik
    )

    if(abs(t.upper[["f.root"]]) < froot_tol ){
      upper_limit <- t.upper[["root"]]
    } else {
      warning("failure to compute upper confidence limit")
      upper_limit <- NA_real_
    }
    upper_froot <- t.upper[["f.root"]]

  } else {
    upper_limit <- upper_froot <- NA_real_
  }

### -- prepare and return output

  result <- c(lower_limit, upper_limit)
  names(result) <- paste(parm, c("lower", "upper"), sep = ".")
  attr(result, "prfloglik") <- list(
    estimate = param_estimate,
    loglik = maximized_loglik,
    qtest = qtest,
    test = test, level = level, nobs_wrc = nobs_wrc, nobs_hcc = nobs_hcc,
    lower_froot = lower_froot, upper_froot = upper_froot
  )

  result

}
