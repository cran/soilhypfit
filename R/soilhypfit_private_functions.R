# ## ======================================================================
ineq_constraint_nlp <- function(
  adjustable.param, adjustable.param.name,
  fixed.param, fixed.param.name,
  param.name, method = NULL,
  param_tf, bwd_tf, deriv_fwd_tfd = NULL, param_bound,
  fit.linear.param = NULL, param_tf.k0 = NULL,
  precBits = NULL, delta_sat_0 = NULL,
  wrc = NULL, wc = NULL, head.wc = NULL, weights_wc = NULL, wrc_model,
  hcc = NULL, hc = NULL, head.hc = NULL, weights_hc = NULL, hcc_model,
  common.env, verbose = NULL,
  e0, param.lc, approximation_alpha_k0,
  ratio_lc_lt_bound,
  back.transform = TRUE
){

  ## function evaluates inequalities for bounds of Lc/Lt

  ## 2019-11-27 A. Papritz

  stopifnot(all(c(wrc_model == "vg", hcc_model == "vgm")))

  ## get linear parameters

  linear.param <- get("linear.param", common.env)

  ## restore names of parameters

  names(adjustable.param) <- adjustable.param.name
  if(length(fixed.param)) names(fixed.param) <- fixed.param.name
  param <- c(adjustable.param, fixed.param)[param.name]

  ## transform parameters back to original scale

  if(back.transform){
    param <- sapply(
      param.name,
      function(x, bwd_tf, param_tf, param, bounds){
        bwd_tf[[param_tf[[x]]]](param[x], bounds[[x]][1], bounds[[x]][2])
      },
      bwd_tf = bwd_tf,
      param_tf = param_tf,
      param = param,
      bounds = param_bound
    )
  }
  names(param) <- param.name

  ## check validity of c0

  if(param["n"] < approximation_alpha_k0["c0"]) warning(
    "estimated 'n' < 'c0': Lt and/or Lc cannot be computed\ndecrease value of 'c0'"
  )

  ## provide value for k0

  t.k0 <- param.lc["k0"]
  if("k0" %in% names(linear.param)) t.k0 <- linear.param["k0"]
  t.k0 <- unname(t.k0)

  ## provide value for tau

  t.tau <- param.lc["tau"]
  if("tau" %in% names(param)) t.tau <- param["tau"]
  t.tau <- unname(t.tau)

  ## values of Lc and Lt

  lc <- lc(
    alpha = param["alpha"], n = param["n"], tau = t.tau,
    k0 = t.k0, e0 = e0,
    c0 = approximation_alpha_k0["c0"],
    c3 = approximation_alpha_k0["c3"],
    c4 = approximation_alpha_k0["c4"]
  )

  lt <- lt(
    n = param["n"], tau = t.tau, e0 = e0,
    c0 = approximation_alpha_k0["c0"],
    c1 = approximation_alpha_k0["c1"],
    c2 = approximation_alpha_k0["c2"],
    c3 = approximation_alpha_k0["c3"],
    c4 = approximation_alpha_k0["c4"]
  )

  ## compute value of inequality constraint

  t.lower <- -lc / lt + ratio_lc_lt_bound[1]

  t.upper <-  lc / lt - ratio_lc_lt_bound[2]

  ## return result

  result <- c(lower = unname(t.lower), upper = unname(t.upper))

  attr(result, "lc") <- lc
  attr(result, "lt") <- lt
  attr(result, "ratio_lc_lt_bound") <- ratio_lc_lt_bound
  attr(result, "e0") <- c(e0 = e0)

  result

}


## ======================================================================
d1ineq_constraint_nlp <- function(
  adjustable.param, adjustable.param.name,
  fixed.param, fixed.param.name,
  param.name, method = NULL,
  param_tf, bwd_tf, deriv_fwd_tfd, param_bound,
  fit.linear.param = NULL, param_tf.k0 = NULL,
  precBits = NULL, delta_sat_0 = NULL,
  wrc = NULL, wc = NULL, head.wc = NULL, weights_wc = NULL, wrc_model,
  hcc = NULL, hc = NULL, head.hc = NULL, weights_hc = NULL, hcc_model,
  common.env = NULL, verbose = NULL,
  e0, param.lc, approximation_alpha_k0,
  ratio_lc_lt_bound,
  back.transform = TRUE
){

  # function evaluates jacobian of inequalities for bounds of Lc/Lt

  # 2019-11-27 A. Papritz

  stopifnot(all(c(wrc_model == "vg", hcc_model == "vgm")))

  ## get linear parameters

  linear.param <- get("linear.param", common.env)

  ## restore names of parameters

  names(adjustable.param) <- adjustable.param.name
  if(length(fixed.param)) names(fixed.param) <- fixed.param.name
  param <- c(adjustable.param, fixed.param)[param.name]

  ## transform parameters back to original scale

  if(back.transform){
    param <- sapply(
      param.name,
      function(x, bwd_tf, param_tf, param, bounds){
        bwd_tf[[param_tf[[x]]]](param[x], bounds[[x]][1], bounds[[x]][2])
      },
      bwd_tf = bwd_tf,
      param_tf = param_tf,
      param = param,
      bounds = param_bound
    )
  }

  names(param) <- param.name

  ## check validity of c0

  if(param["n"] < approximation_alpha_k0["c0"]) warning(
    "estimated 'n' < 'c0': Lt and/or Lc cannot be computed\ndecrease value of 'c0'"
  )

  ## provide value for k0

  t.k0 <- param.lc["k0"]
  if("k0" %in% names(linear.param)) t.k0 <- linear.param["k0"]
  t.k0 <- unname(t.k0)

  ## provide value for tau

  t.tau <- param.lc["tau"]
  if("tau" %in% names(param)) t.tau <- param["tau"]
  t.tau <- unname(t.tau)

  ## compute jacobian of inequality constraint

  target <- adjustable.param.name

  ## values of Lc and Lt

  lc <- lc(
    alpha = param["alpha"], n = param["n"], tau = t.tau,
    k0 = t.k0, e0 = e0,
    c0 = approximation_alpha_k0["c0"],
    c3 = approximation_alpha_k0["c3"],
    c4 = approximation_alpha_k0["c4"]
  )

  lt <- lt(
    n = param["n"], tau = t.tau, e0 = e0,
    c0 = approximation_alpha_k0["c0"],
    c1 = approximation_alpha_k0["c1"],
    c2 = approximation_alpha_k0["c2"],
    c3 = approximation_alpha_k0["c3"],
    c4 = approximation_alpha_k0["c4"]
  )

  ## derivative of Lc with respect to alpha, n, tau

  d1lc <- d1lc(
    alpha = param["alpha"], n = param["n"], tau = t.tau,
    k0 = t.k0, e0 = e0, target = adjustable.param.name,
    c0 = approximation_alpha_k0["c0"],
    c3 = approximation_alpha_k0["c3"],
    c4 = approximation_alpha_k0["c4"]
  )

  ## derivative of Lc with respect to transformed parameters

  if(back.transform){
    d1lc <- sapply(
      names(d1lc),
      function(x, pd, param, param_tf, deriv_fwd_tfd, bounds){
        pd[x] <- unname(
          pd[x] / deriv_fwd_tfd[[param_tf[[x]]]](param[x], bounds[[x]][1], bounds[[x]][2])
        )
      },
      pd = d1lc,
      param = param,
      param_tf = param_tf,
      deriv_fwd_tfd = deriv_fwd_tfd,
      bounds = param_bound
    )
  }


  ## derivative of Lt with respect to n, tau

  d1l0 <- d1l0(
    n = param["n"], tau = t.tau, e0 = e0,
    target = adjustable.param.name,
    c0 = approximation_alpha_k0["c0"],
    c1 = approximation_alpha_k0["c1"],
    c2 = approximation_alpha_k0["c2"],
    c3 = approximation_alpha_k0["c3"],
    c4 = approximation_alpha_k0["c4"]
  )

  ## derivative of Lc with respect to transformed parameters

  if(back.transform){
    d1l0 <- sapply(
      names(d1l0),
      function(x, pd, param, param_tf, deriv_fwd_tfd, bounds){
        pd[x] <- unname(
          pd[x] / deriv_fwd_tfd[[param_tf[[x]]]](param[x], bounds[[x]][1], bounds[[x]][2])
        )
      },
      pd = d1l0,
      param = param,
      param_tf = param_tf,
      deriv_fwd_tfd = deriv_fwd_tfd,
      bounds = param_bound
    )
  }


  ## compute partial derivatives

  t.upper <- (d1lc * lt - lc * d1l0) / lt^2
  t.lower <- -t.upper

#   ## compute partial derivatives
#
#   t.lower <- -d1lc / lc + d1l0 / lt
#   t.upper <-  d1lc / lc - d1l0 / lt
#
  ## return jacobian as matrix

  result <- rbind(lower = t.lower, upper = t.upper)

  attr(result, "lc") <- lc
  attr(result, "d1lc") <- d1lc
  attr(result, "lt") <- lt
  attr(result, "d1l0") <- d1l0

  result

}


## ======================================================================
d1l0 <- function(n, tau, e0, target, c0, c1, c2, c3, c4){

  ## function computes first partial derivatives of constraining value for
  ## capillary length Lc based on VGM model and evaporation rate of a
  ## saturated soil  and approximations
  ## for alpha and k0 with respect to VGM parameters alpha, n and tau

  ## 2019-11-27 A. Papritz

  stopifnot(all(target %in% c("alpha", "n", "tau")))

  nm1 <- (-1 + n)
  nonm1 <- n/nm1
  nm1on <- nm1/n
  m2p1on <- (-2 + 1/n)

  aux00 <- (nm1on)^m2p1on
  aux0 <- 1 + (aux00)^n
  aux1 <- (aux0)^(-1 + 1/n)
  aux2 <- (aux1)^tau
  aux3 <- (1 - (aux1)^(nonm1))^(nm1on)
  aux4 <- aux2*(-c0 + n)^c4
  aux5 <- (-1 + aux3)^2*aux4
  aux6 <- 4.*c3*aux5

  t.n <- t.tau <-  NULL

  if("n" %in% target){

    t.n <- (
      (2 - 1/n)^(2 - 1/n)*(nm1on)^(1/n)*n*(
        c2*nm1*(-c0 + n)*(1 + e0/(aux6)) -
        nm1*(1 - c0*c2 + c2*n)*(1 + e0/(aux6)) -
        (-c0 + n)*(1 - c0*c2 + c2*n)*(1 + e0/(aux6)) +
        (nm1*(-c0 + n)*(1 - c0*c2 + c2*n)*(1 + e0/(aux6))*(1 + log(2 - 1/n)))/n^2 -
        (nm1*(-c0 + n)*(1 - c0*c2 + c2*n)*(1 + e0/(aux6))*(1 + log(nm1on)))/n^2 -
        (
          e0*nm1*(-c0 + n)^(1 - c4)*(1 - c0*c2 + c2*n)*(
            (c4*(1 - aux3))/(c0 - n) +
            ((1 - aux3)*tau*((aux0)*log(aux0) + (aux00)^n*(1 - 2*n + nm1*n*log(aux00) - nm1*log(nm1on)))) / ((aux0)*n^2) +
            (
              2*(
                -((-1 + (aux1)^(nonm1))*log(1 - (aux1)^(nonm1))) +
                (
                  (aux1)^(nonm1)*(
                    (aux0)*n*log(aux1) - (1 - n)*((aux0)*log(aux0) + (aux00)^n*(1 - 2*n + nm1*n*log(aux00) - nm1*log(nm1on)))
                  )
                ) / ((aux0)*nm1))
            ) / ((1 - (aux1)^(nonm1))^(1/n)*n^2)
          )
        ) / (4.*c3*(1 - aux3)^3*aux2)
      )
    ) / (c1*(c0 - n)^2*nm1^3*(1 + e0/(aux6))^2)

  }

  if("tau" %in% target){

    t.tau <- (
      e0*(-c0 + n)^(-1 - c4)*(nm1/(-1 + 2*n))^m2p1on*(1 - c0*c2 + c2*n)*log(aux1)
    ) / (
      4.*c1*c3*(-1 + aux3)^2*aux2*n*(1 + e0/(aux6))^2
    )

  }

  c(alpha = 0., n = unname(t.n), tau = unname(t.tau))
#   c(n = unname(t.n), tau = unname(t.tau))

}



## ======================================================================
d1lc <- function(alpha, n, tau, k0, e0, target, c0 = NULL, c3 = NULL, c4 = NULL){

  ## function computes first partial derivatives of capillary length Lc
  ## based on VGM model and evaporation rate of a saturated soil with
  ## respect to VGM parameters alpha, n and tau; if k0 is not available then
  ## the approxmation k0 = c3 * (n - c0)^c4 is used

  ## 2019-11-27 A. Papritz

  stopifnot(all(target %in% c("alpha", "n", "tau")))

  ## auxiliary computations

  nm1 <- (-1 + n)
  nonm1 <- n/nm1
  nm1on <- nm1/n
  m2p1on <- (-2 + 1/n)

  aux00 <- (nm1on)^m2p1on
  aux0 <- 1 + (aux00)^n
  aux1 <- (aux0)^(-1 + 1/n)
  aux2 <- (aux1)^tau
  aux3 <- (1 - (aux1)^(nonm1))^(nm1on)


  t.alpha <- t.n <- t.tau <-  NULL

  if("alpha" %in% target){

    ## old approximation for k0 = c2 * n^c3
    ##    t.k0 <- k0
    ##    if(is.na(k0)) t.k0 <- c2 * n^c3

    t.k0 <- k0
    if(is.na(k0)) t.k0 <- c3 * (n - c0)^c4

    t.alpha <-  -(
      ((nm1)/(-1 + 2*n))^(m2p1on) / (
        (n + (e0*n)/(4.*t.k0*(-1 + (1 - (aux1)^(n/(nm1)))^((nm1)/n))^2*aux2))*alpha^2)
    )

  }

  if("n" %in% target){

    if(is.na(k0)){

      ## using approximation for k0

      aux4 <- aux2*(-c0 + n)^c4
      aux5 <- (-1 + aux3)^2*aux4
      aux6 <- 4.*c3*aux5

      t.n <- (
        (2 - 1/n)^(2 - 1/n)*(nm1on)^(1/n)*n*(
          -1 - e0/(aux6) + (nm1*(1 + e0/(aux6))*(1 + log(2 - 1/n)))/n^2 -
          (nm1*(1 + e0/(aux6))*(1 + log(nm1on)))/n^2 -
          (
            e0*nm1*(
              (c4*(1 - aux3))/(c0 - n) +
              ((1 - aux3)*tau*((aux0)*log(aux0) + (aux00)^n*(1 - 2*n + nm1*n*log(aux00) - nm1*log(nm1on)))) / ((aux0)*n^2) +
              (
                2*(
                  -((-1 + (aux1)^(nonm1))*log(1 - (aux1)^(nonm1))) +
                  (
                    (aux1)^(nonm1)*((aux0)*n*log(aux1)- (1 - n)*((aux0)*log(aux0) + (aux00)^n*(1 - 2*n + nm1*n*log(aux00) - nm1*log(nm1on))))
                  )/((aux0)*nm1))
              ) / ((1 - (aux1)^(nonm1))^(1/n)*n^2))
          ) / (4.*c3*(1 - aux3)^3*aux4)
        )
      ) / (
        nm1^3*(1 + e0/(aux6))^2*alpha
      )

    } else {

      ## using known k0

      aux4 <- 4.*k0*(-1 + aux3)^2*aux2
      aux5 <- e0/(aux4)
      aux6 <- (1 + aux5)*nm1

      t.n <- (
        (2 - 1/n)^(2 - 1/n)*(nm1on)^(1/n)*n*(
          -1 - aux5 + (aux6*(1 + log(2 - 1/n)))/n^2 - (aux6*(1 + log(nm1on)))/n^2 + (
            e0*(
              (1 - aux3)*(1 - n)*tau*((aux0)*log(aux0) + (aux00)^n*(1 - 2*n + nm1*n*log(aux00) - nm1*log(nm1on))) +
              (
                2*(
                  (-1 + (aux1)^(nonm1))*(aux0)*nm1*log(1 - (aux1)^(nonm1)) +
                  (aux1)^(nonm1)*(-((aux0)*n*log(aux1)) + (1 - n)*((aux0)*log(aux0) + (aux00)^n*(1 - 2*n + nm1*n*log(aux00) - nm1*log(nm1on))))
                )
              ) / (1 - (aux1)^(nonm1))^(1/n))
          ) / (
            4.*k0*(1 - aux3)^3*(aux0)*aux2*n^2
          )
        )
      ) / (
        (1 + aux5)^2*nm1^3*alpha
      )
    }

  }

  if("tau" %in% target){

    t.k0 <- k0
    if(is.na(k0)) t.k0 <- c3 * (n - c0)^c4

    t.tau <- (
      e0*(nm1/(-1 + 2*n))^m2p1on*log(aux1)
    ) / (
      4.*t.k0*(-1 + aux3)^2*(1 + e0/(4.*t.k0*(-1 + aux3)^2*aux2))^2*aux2*n*alpha
    )

  }

  c(alpha = unname(t.alpha), n = unname(t.n), tau = unname(t.tau))

}


## ======================================================================
### estimate_nlp

estimate_nlp <- function(
  nonlinear.param,
  fit.nonlinear.param,
  linear.param,
  fit.linear.param,
  param_tf.k0,
  e0, ratio_lc_lt_bound, param.lc,
  wrc, wc, head.wc, weights_wc, wrc_model,
  hcc, hc, head.hc, weights_hc, hcc_model,
  control,
  common.env,
  verbose
){

  ## function estimates the nonlinear parameters alpha, n and tau subject
  ## to physical constraints; the linear parameters thetar, thetas and k0
  ## are estimated as "by-product"

  ## 2019-11-27 A. Papritz
  ## 2020-01-27 AP ml and mpd estimation
  ## 2021-02-27 AP correction of error checking for checking validity of Hessian

#### -- transform parameters

  param.name <- names(nonlinear.param)

  transformed.param <- sapply(
    param.name,
    function(x, fwd_tf, param_tf, param, bounds){
      fwd_tf[[param_tf[[x]]]](param[x], bounds[[x]][1], bounds[[x]][2])
    },
    fwd_tf = control[["fwd_tf"]],
    param_tf = control[["param_tf"]],
    param = nonlinear.param,
    bounds = control[["param_bound"]]
  )

  names(transformed.param) <- param.name

  ## determine adjustable and fixed parameters

  sel <- names(fit.nonlinear.param)[fit.nonlinear.param]
  adjustable.param <- transformed.param[names(transformed.param) %in% sel]
  fixed.param <- transformed.param[!names(transformed.param) %in% sel]

#### -- initialize values for counters of iterations and gradient evaluations,
  ## of objective function and linear parameters

  assign("n.it",         1L,           envir = common.env)
  assign("n.grd",        0L,           envir = common.env)
  assign("linear.param", linear.param, envir = common.env)
  assign("obj",          Inf,          envir = common.env)

#### -- estimate parameters

  #   print("gradient")

  if(length(adjustable.param)){

    ## at least one parameter must be estimated estimated

    if(identical(control[["settings"]], "sce")){

#### --- sce

      ## Shuffled Complex Evolution (SCE) optimisation (unconstrained
      ## global algorithm, function SCEoptim of R package
      ## hydromad, cf. http://hydromad.catchment.org/#)

      result <- SCEoptim(
        FUN = objective_nlp,
        par = adjustable.param,
        lower = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[1]),
        upper = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[2]),
        control = control[["sce"]],
        adjustable.param.name = names(adjustable.param),
        fixed.param = fixed.param, fixed.param.name = names(fixed.param),
        param.name = param.name,
        method = control[["method"]],
        param_tf = control[["param_tf"]],
        bwd_tf = control[["bwd_tf"]],
        deriv_fwd_tfd = control[["deriv_fwd_tfd"]],
        param_bound = control[["param_bound"]],
        fit.linear.param = fit.linear.param,
        param_tf.k0 = param_tf.k0,
        precBits = control[["precBits"]],
        delta_sat_0 = control[["delta_sat_0"]],
        wrc = wrc, wc = wc, head.wc = head.wc,
        weights_wc = weights_wc, wrc_model = wrc_model,
        hcc = hcc, hc = hc, head.hc = head.hc,
        weights_hc = weights_hc, hcc_model = hcc_model,
        common.env = common.env,
        verbose = verbose,
        e0 = e0,
        param.lc = param.lc,
        approximation_alpha_k0 = control[["approximation_alpha_k0"]]
      )

      tmp <- result[["par"]]
      names(tmp) <- names(adjustable.param)
      transformed.param <- c(tmp, fixed.param)[param.name]

      result[["evaluation"]] = c(
        "function" = get("n.it", common.env),
        "gradient" = get("n.grd", common.env)
      )
      result[["objective"]] = result[["value"]]
      result[["message"]] <- paste("SCEoptim:", result[["message"]])

    } else {

      ## optimization by algorithms of NLOPT library (R package nloptr, cf.
      ## https://nlopt.readthedocs.io/en/latest/


      result <- switch(
        control[["settings"]],

        ## unconstrained local NLopt algorithm

        "ulocal" = {

          if(control[["use_derivative"]]){

#### --- ulocal LD algorithm

            if(all(control[["param_tf"]][names(adjustable.param)] %in% "identity")){

              ## specify bounds if all parameters are untransformed

              nloptr(
                x0 = adjustable.param,
                eval_f = objective_nlp,
                eval_grad_f = gradient_nlp,
                lb = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[1]),
                ub = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[2]),
                opts = control[["nloptr"]],
                adjustable.param.name = names(adjustable.param),
                fixed.param = fixed.param, fixed.param.name = names(fixed.param),
                param.name = param.name,
                method = control[["method"]],
                param_tf = control[["param_tf"]],
                bwd_tf = control[["bwd_tf"]],
                deriv_fwd_tfd = control[["deriv_fwd_tfd"]],
                param_bound = control[["param_bound"]],
                fit.linear.param = fit.linear.param,
                param_tf.k0 = param_tf.k0,
                precBits = control[["precBits"]],
                delta_sat_0 = control[["delta_sat_0"]],
                wrc = wrc, wc = wc, head.wc = head.wc,
                weights_wc = weights_wc, wrc_model = wrc_model,
                hcc = hcc, hc = hc, head.hc = head.hc,
                weights_hc = weights_hc, hcc_model = hcc_model,
                common.env = common.env,
                verbose = verbose,
                e0 = e0,
                param.lc = param.lc,
                approximation_alpha_k0 = control[["approximation_alpha_k0"]],
                ratio_lc_lt_bound = ratio_lc_lt_bound,
                back.transform = TRUE
              )

            } else {

              ## do not specify bounds if some parameters are transformed

              nloptr(
                x0 = adjustable.param,
                eval_f = objective_nlp,
                eval_grad_f = gradient_nlp,
                opts = control[["nloptr"]],
                adjustable.param.name = names(adjustable.param),
                fixed.param = fixed.param, fixed.param.name = names(fixed.param),
                param.name = param.name,
                method = control[["method"]],
                param_tf = control[["param_tf"]],
                bwd_tf = control[["bwd_tf"]],
                deriv_fwd_tfd = control[["deriv_fwd_tfd"]],
                param_bound = control[["param_bound"]],
                fit.linear.param = fit.linear.param,
                param_tf.k0 = param_tf.k0,
                precBits = control[["precBits"]],
                delta_sat_0 = control[["delta_sat_0"]],
                wrc = wrc, wc = wc, head.wc = head.wc,
                weights_wc = weights_wc, wrc_model = wrc_model,
                hcc = hcc, hc = hc, head.hc = head.hc,
                weights_hc = weights_hc, hcc_model = hcc_model,
                common.env = common.env,
                verbose = verbose,
                e0 = e0,
                param.lc = param.lc,
                approximation_alpha_k0 = control[["approximation_alpha_k0"]],
                ratio_lc_lt_bound = ratio_lc_lt_bound,
                back.transform = TRUE
              )

            }

          } else {

#### --- ulocal LN algorithm

            if(
              all(control[["param_tf"]][names(adjustable.param)] %in% "identity") &&
              !identical(control[["nloptr"]][["algorithm"]], "NLOPT_LN_NEWUOA")

            ){

              ## specify bounds if all parameters are untransformed

              nloptr(
                x0 = adjustable.param,
                eval_f = objective_nlp,
                lb = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[1]),
                ub = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[2]),
                opts = control[["nloptr"]],
                adjustable.param.name = names(adjustable.param),
                fixed.param = fixed.param, fixed.param.name = names(fixed.param),
                param.name = param.name,
                method = control[["method"]],
                param_tf = control[["param_tf"]],
                bwd_tf = control[["bwd_tf"]],
                deriv_fwd_tfd = control[["deriv_fwd_tfd"]],
                param_bound = control[["param_bound"]],
                fit.linear.param = fit.linear.param,
                param_tf.k0 = param_tf.k0,
                precBits = control[["precBits"]],
                delta_sat_0 = control[["delta_sat_0"]],
                wrc = wrc, wc = wc, head.wc = head.wc,
                weights_wc = weights_wc, wrc_model = wrc_model,
                hcc = hcc, hc = hc, head.hc = head.hc,
                weights_hc = weights_hc, hcc_model = hcc_model,
                common.env = common.env,
                verbose = verbose,
                e0 = e0,
                param.lc = param.lc,
                approximation_alpha_k0 = control[["approximation_alpha_k0"]],
                ratio_lc_lt_bound = ratio_lc_lt_bound,
                back.transform = TRUE
              )

            } else {

              ## do not specify bounds if some parameters are transformed

              nloptr(
                x0 = adjustable.param,
                eval_f = objective_nlp,
                opts = control[["nloptr"]],
                adjustable.param.name = names(adjustable.param),
                fixed.param = fixed.param, fixed.param.name = names(fixed.param),
                param.name = param.name,
                method = control[["method"]],
                param_tf = control[["param_tf"]],
                bwd_tf = control[["bwd_tf"]],
                deriv_fwd_tfd = control[["deriv_fwd_tfd"]],
                param_bound = control[["param_bound"]],
                fit.linear.param = fit.linear.param,
                param_tf.k0 = param_tf.k0,
                precBits = control[["precBits"]],
                delta_sat_0 = control[["delta_sat_0"]],
                wrc = wrc, wc = wc, head.wc = head.wc,
                weights_wc = weights_wc, wrc_model = wrc_model,
                hcc = hcc, hc = hc, head.hc = head.hc,
                weights_hc = weights_hc, hcc_model = hcc_model,
                common.env = common.env,
                verbose = verbose,
                e0 = e0,
                param.lc = param.lc,
                approximation_alpha_k0 = control[["approximation_alpha_k0"]],
                ratio_lc_lt_bound = ratio_lc_lt_bound,
                back.transform = TRUE
              )
            }

          }

        },

        ## constrained local NLopt algorithm

        "clocal" = {

          if(control[["use_derivative"]]){

#### --- clocal LD algorithm

            if(all(control[["param_tf"]][names(adjustable.param)] %in% "identity")){

              ## specify bounds if all parameters are untransformed

              nloptr(
                x0 = adjustable.param,
                eval_f = objective_nlp,
                eval_grad_f = gradient_nlp,
                eval_g_ineq = ineq_constraint_nlp,
                eval_jac_g_ineq = d1ineq_constraint_nlp,
                lb = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[1]),
                ub = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[2]),
                opts = control[["nloptr"]],
                adjustable.param.name = names(adjustable.param),
                fixed.param = fixed.param, fixed.param.name = names(fixed.param),
                param.name = param.name,
                method = control[["method"]],
                param_tf = control[["param_tf"]],
                bwd_tf = control[["bwd_tf"]], deriv_fwd_tfd = control[["deriv_fwd_tfd"]],
                param_bound = control[["param_bound"]],
                fit.linear.param = fit.linear.param,
                param_tf.k0 = param_tf.k0,
                precBits = control[["precBits"]],
                delta_sat_0 = control[["delta_sat_0"]],
                wrc = wrc, wc = wc, head.wc = head.wc,
                weights_wc = weights_wc, wrc_model = wrc_model,
                hcc = hcc, hc = hc, head.hc = head.hc,
                weights_hc = weights_hc, hcc_model = hcc_model,
                common.env = common.env,
                verbose = verbose,
                e0 = e0,
                param.lc = param.lc,
                approximation_alpha_k0 = control[["approximation_alpha_k0"]],
                ratio_lc_lt_bound = ratio_lc_lt_bound,
                back.transform = TRUE
              )

            } else {

              ## do not specify bounds if some parameters are transformed

              nloptr(
                x0 = adjustable.param,
                eval_f = objective_nlp,
                eval_grad_f = gradient_nlp,
                eval_g_ineq = ineq_constraint_nlp,
                eval_jac_g_ineq = d1ineq_constraint_nlp,
                opts = control[["nloptr"]],
                adjustable.param.name = names(adjustable.param),
                fixed.param = fixed.param, fixed.param.name = names(fixed.param),
                param.name = param.name,
                method = control[["method"]],
                param_tf = control[["param_tf"]],
                bwd_tf = control[["bwd_tf"]], deriv_fwd_tfd = control[["deriv_fwd_tfd"]],
                param_bound = control[["param_bound"]],
                fit.linear.param = fit.linear.param,
                param_tf.k0 = param_tf.k0,
                precBits = control[["precBits"]],
                delta_sat_0 = control[["delta_sat_0"]],
                wrc = wrc, wc = wc, head.wc = head.wc,
                weights_wc = weights_wc, wrc_model = wrc_model,
                hcc = hcc, hc = hc, head.hc = head.hc,
                weights_hc = weights_hc, hcc_model = hcc_model,
                common.env = common.env,
                verbose = verbose,
                e0 = e0,
                param.lc = param.lc,
                approximation_alpha_k0 = control[["approximation_alpha_k0"]],
                ratio_lc_lt_bound = ratio_lc_lt_bound,
                back.transform = TRUE
              )

            }

          } else {

#### --- clocal LN algorithm

            if(all(control[["param_tf"]][names(adjustable.param)] %in% "identity")){

              ## specify bounds if all parameters are untransformed

              nloptr(
                x0 = adjustable.param,
                eval_f = objective_nlp,
                eval_g_ineq = ineq_constraint_nlp,
                lb = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[1]),
                ub = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[2]),
                opts = control[["nloptr"]],
                adjustable.param.name = names(adjustable.param),
                fixed.param = fixed.param, fixed.param.name = names(fixed.param),
                param.name = param.name,
                method = control[["method"]],
                param_tf = control[["param_tf"]],
                bwd_tf = control[["bwd_tf"]], deriv_fwd_tfd = control[["deriv_fwd_tfd"]],
                param_bound = control[["param_bound"]],
                fit.linear.param = fit.linear.param,
                param_tf.k0 = param_tf.k0,
                precBits = control[["precBits"]],
                delta_sat_0 = control[["delta_sat_0"]],
                wrc = wrc, wc = wc, head.wc = head.wc,
                weights_wc = weights_wc, wrc_model = wrc_model,
                hcc = hcc, hc = hc, head.hc = head.hc,
                weights_hc = weights_hc, hcc_model = hcc_model,
                common.env = common.env,
                verbose = verbose,
                e0 = e0,
                param.lc = param.lc,
                approximation_alpha_k0 = control[["approximation_alpha_k0"]],
                ratio_lc_lt_bound = ratio_lc_lt_bound,
                back.transform = TRUE
              )

            } else {

              ## do not specify bounds if some parameters are transformed

              nloptr(
                x0 = adjustable.param,
                eval_f = objective_nlp,
                eval_g_ineq = ineq_constraint_nlp,
                opts = control[["nloptr"]],
                adjustable.param.name = names(adjustable.param),
                fixed.param = fixed.param, fixed.param.name = names(fixed.param),
                param.name = param.name,
                method = control[["method"]],
                param_tf = control[["param_tf"]],
                bwd_tf = control[["bwd_tf"]], deriv_fwd_tfd = control[["deriv_fwd_tfd"]],
                param_bound = control[["param_bound"]],
                fit.linear.param = fit.linear.param,
                param_tf.k0 = param_tf.k0,
                precBits = control[["precBits"]],
                delta_sat_0 = control[["delta_sat_0"]],
                wrc = wrc, wc = wc, head.wc = head.wc,
                weights_wc = weights_wc, wrc_model = wrc_model,
                hcc = hcc, hc = hc, head.hc = head.hc,
                weights_hc = weights_hc, hcc_model = hcc_model,
                common.env = common.env,
                verbose = verbose,
                e0 = e0,
                param.lc = param.lc,
                approximation_alpha_k0 = control[["approximation_alpha_k0"]],
                ratio_lc_lt_bound = ratio_lc_lt_bound,
                back.transform = TRUE
              )

            }

          }

        },

        ## global NLopt algorithms

        "uglobal"= {

          if(control[["use_derivative"]]){

#### --- uglobal GD algorithm

            nloptr(
              x0 = adjustable.param,
              eval_f = objective_nlp,
              eval_grad_f = gradient_nlp,
              lb = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[1]),
              ub = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[2]),
              opts = control[["nloptr"]],
              adjustable.param.name = names(adjustable.param),
              fixed.param = fixed.param, fixed.param.name = names(fixed.param),
              param.name = param.name,
              method = control[["method"]],
              param_tf = control[["param_tf"]],
              bwd_tf = control[["bwd_tf"]],
              deriv_fwd_tfd = control[["deriv_fwd_tfd"]],
              param_bound = control[["param_bound"]],
              fit.linear.param = fit.linear.param,
              param_tf.k0 = param_tf.k0,
              precBits = control[["precBits"]],
              delta_sat_0 = control[["delta_sat_0"]],
              wrc = wrc, wc = wc, head.wc = head.wc,
              weights_wc = weights_wc, wrc_model = wrc_model,
              hcc = hcc, hc = hc, head.hc = head.hc,
              weights_hc = weights_hc, hcc_model = hcc_model,
              common.env = common.env,
              verbose = verbose,
              e0 = e0,
              param.lc = param.lc,
              approximation_alpha_k0 = control[["approximation_alpha_k0"]],
              ratio_lc_lt_bound = ratio_lc_lt_bound,
              back.transform = TRUE
            )

          } else {

#### --- uglobal GN algorithm

            nloptr(
              x0 = adjustable.param,
              eval_f = objective_nlp,
              lb = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[1]),
              ub = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[2]),
              opts = control[["nloptr"]],
              adjustable.param.name = names(adjustable.param),
              fixed.param = fixed.param, fixed.param.name = names(fixed.param),
              param.name = param.name,
              method = control[["method"]],
              param_tf = control[["param_tf"]],
              bwd_tf = control[["bwd_tf"]],
              deriv_fwd_tfd = control[["deriv_fwd_tfd"]],
              param_bound = control[["param_bound"]],
              fit.linear.param = fit.linear.param,
              param_tf.k0 = param_tf.k0,
              precBits = control[["precBits"]],
              delta_sat_0 = control[["delta_sat_0"]],
              wrc = wrc, wc = wc, head.wc = head.wc,
              weights_wc = weights_wc, wrc_model = wrc_model,
              hcc = hcc, hc = hc, head.hc = head.hc,
              weights_hc = weights_hc, hcc_model = hcc_model,
              common.env = common.env,
              verbose = verbose,
              e0 = e0,
              param.lc = param.lc,
              approximation_alpha_k0 = control[["approximation_alpha_k0"]],
              ratio_lc_lt_bound = ratio_lc_lt_bound,
              back.transform = TRUE
            )

          }

        },

        "cglobal" = {

          if(control[["use_derivative"]]){

            ## GD algorithm

            stop("no constrained global NLOPT algorithm available")

          } else {

#### --- cglobal GN algorithm

            nloptr(
              x0 = adjustable.param,
              eval_f = objective_nlp,
              eval_g_ineq = ineq_constraint_nlp,
              lb = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[1]),
              ub = sapply(control[["param_bound"]][names(adjustable.param)], function(x) x[2]),
              opts = control[["nloptr"]],
              adjustable.param.name = names(adjustable.param),
              fixed.param = fixed.param, fixed.param.name = names(fixed.param),
              param.name = param.name,
              method = control[["method"]],
              param_tf = control[["param_tf"]],
              bwd_tf = control[["bwd_tf"]], deriv_fwd_tfd = control[["deriv_fwd_tfd"]],
              param_bound = control[["param_bound"]],
              fit.linear.param = fit.linear.param,
              param_tf.k0 = param_tf.k0,
              precBits = control[["precBits"]],
              delta_sat_0 = control[["delta_sat_0"]],
              wrc = wrc, wc = wc, head.wc = head.wc,
              weights_wc = weights_wc, wrc_model = wrc_model,
              hcc = hcc, hc = hc, head.hc = head.hc,
              weights_hc = weights_hc, hcc_model = hcc_model,
              common.env = common.env,
              verbose = verbose,
              e0 = e0,
              param.lc = param.lc,
              approximation_alpha_k0 = control[["approximation_alpha_k0"]],
              ratio_lc_lt_bound = ratio_lc_lt_bound,
              back.transform = TRUE
            )

          }

        },

        stop("unknown optimization approach")

      )

      tmp <- result[["solution"]]
      names(tmp) <- names(adjustable.param)
      transformed.param <- c(tmp, fixed.param)[param.name]

      result[["evaluation"]] = c(
        "function" = get("n.it", common.env),
        "gradient" = get("n.grd", common.env)
      )
      result[["convergence"]] = result[["status"]]

    }

  } else {

    ## all nonlinear parameters are fixed

    transformed.param <- fixed.param[param.name]

    result <- list(
      par = transformed.param,
      objective = NA_real_,
      convergence = NA_integer_,
      evaluation = c("function" = 1, "gradient" = 0),
      message = "all nonlinear parameters are fixed"
    )

  }

#### -- transform parameters back to original scale

  param <- sapply(
    param.name,
    function(x, bwd_tf, param_tf, param, bounds){
      bwd_tf[[param_tf[[x]]]](param[x], bounds[[x]][1], bounds[[x]][2])
    },
    bwd_tf = control[["bwd_tf"]],
    param_tf = control[["param_tf"]],
    param = transformed.param,
    bounds = control[["param_bound"]]
  )
  names(param) <- param.name

#### -- get value of objective function and linear parameters at the solution

  obj <- objective_nlp(
    adjustable.param = transformed.param[names(adjustable.param)],
    adjustable.param.name = names(adjustable.param),
    fixed.param = fixed.param, fixed.param.name = names(fixed.param),
    param.name = param.name,
    method = control[["method"]],
    param_tf = control[["param_tf"]],
    bwd_tf = control[["bwd_tf"]],
    param_bound = control[["param_bound"]],
    fit.linear.param = fit.linear.param,
    param_tf.k0 = param_tf.k0,
    precBits = control[["precBits"]],
    delta_sat_0 = control[["delta_sat_0"]],
    wrc = wrc, wc = wc, head.wc = head.wc,
    weights_wc = weights_wc, wrc_model = control[["wrc_model"]],
    hcc = hcc, hc = hc, head.hc = head.hc, weights_hc = weights_hc,
    hcc_model = control[["hcc_model"]],
    common.env = common.env,
    verbose = 0
  )

  attr(obj, "obj_wc") <- get("obj_wc", common.env)
  attr(obj, "obj_hc") <- get("obj_hc", common.env)
  attr(obj, "ssq_wc") <- get("ssq_wc", common.env)
  attr(obj, "ssq_hc") <- get("ssq_hc", common.env)

  linear.param <- get("linear.param", common.env)

#### -- compute gradient of objective with respect to nonlinear parameters at
  ## solution

  grad <- gradient_nlp(
    adjustable.param = transformed.param[names(adjustable.param)],
    adjustable.param.name = names(adjustable.param),
    fixed.param = fixed.param, fixed.param.name = names(fixed.param),
    param.name = param.name,
    method = control[["method"]],
    param_tf = control[["param_tf"]],
    bwd_tf = control[["bwd_tf"]],
    param_bound = control[["param_bound"]],
    fit.linear.param = fit.linear.param,
    param_tf.k0 = param_tf.k0,
    precBits = control[["precBits"]],
    delta_sat_0 = control[["delta_sat_0"]],
    wrc = wrc, wc = wc, head.wc = head.wc,
    weights_wc = weights_wc, wrc_model = control[["wrc_model"]],
    hcc = hcc, hc = hc, head.hc = head.hc, weights_hc = weights_hc,
    hcc_model = control[["hcc_model"]],
    common.env = common.env,
    verbose = 0
  )

#### -- compute values of Lc and inequality constraints

  ineq <- list(
    lc = ineq_constraint_nlp(
      adjustable.param = transformed.param[names(adjustable.param)],
      adjustable.param.name = names(adjustable.param),
      fixed.param = fixed.param, fixed.param.name = names(fixed.param),
      param.name = param.name,
      param_tf = control[["param_tf"]],
      bwd_tf = control[["bwd_tf"]],
      param_bound = control[["param_bound"]],
      wrc_model = control[["wrc_model"]], hcc_model = control[["hcc_model"]],
      common.env = common.env,
      e0 = e0,
      param.lc = param.lc,
      approximation_alpha_k0 = control[["approximation_alpha_k0"]],
      ratio_lc_lt_bound = ratio_lc_lt_bound,
      back.transform = TRUE
    )
  )

#### -- optionally compute hessian matrix
  ## with respect to transformed nonlinear parameters at solution if either
  ## only wrc or hcc data are available

  ## see Bates & Watts, 1988, eq. 4.6, p. 138

  if(control[["hessian"]]){

    hessian <- optimHess(
      par = transformed.param[names(adjustable.param)],
      fn = objective_nlp,
      gr = gradient_nlp,
      adjustable.param.name = names(adjustable.param),
      fixed.param = fixed.param, fixed.param.name = names(fixed.param),
      param.name = param.name,
      method = control[["method"]],
      param_tf = control[["param_tf"]],
      bwd_tf = control[["bwd_tf"]],
      param_bound = control[["param_bound"]],
      fit.linear.param = fit.linear.param,
      param_tf.k0 = param_tf.k0,
      precBits = control[["precBits"]],
      delta_sat_0 = control[["delta_sat_0"]],
      wrc = wrc, wc = wc, head.wc = head.wc,
      weights_wc = weights_wc, wrc_model = control[["wrc_model"]],
      hcc = hcc, hc = hc, head.hc = head.hc, weights_hc = weights_hc,
      hcc_model = control[["hcc_model"]],
      common.env = common.env,
      verbose = 0
    )

    if(
      length(hessian) > 0L && any( eigen(hessian)[["values"]] < 0. ) ) warning(
      "hessian not positive definite, check whether local maximum of sum of squares has been found"
    )

  } else {
    hessian <- NULL
  }


#### -- optionally draw data and fitted model curves

  if(verbose >= 1. && verbose < 3.){

    op <- par(mfrow = c(1, sum(c(wrc, hcc))))

    if(wrc){
      t.head <- exp(seq(
          min(log(head.wc)), max(log(head.wc)),
          length = control[["gam_n_newdata"]]
        ))
      plot(head.wc, wc, log="x")
      lines(
        t.head,
        as.double(wc_model(
            t.head, param, linear.param, control[["precBits"]],
            wrc_model
          ))
      )
    }

    if(hcc){
      t.head <- exp(seq(
          min(log(head.hc)), max(log(head.hc)),
          length = control[["gam_n_newdata"]]
        ))
      plot(head.hc, hc, log="xy")
      lines(
        t.head,
        as.double(hc_model(
            t.head, param, linear.param,
            control[["precBits"]], hcc_model
          ))
      )
    }

    par(op)
  }


#### -- return estimated  parameters, value of objective function, etc.

  result <- c(
    result[c("convergence", "message", "evaluation")],
    list(
      objective = obj,
      gradient = grad,
      linear.param = linear.param,
      param = param,
      transformed.param = transformed.param,
      inequality_constraint = ineq,
      hessian = hessian
    )
  )

  return(result)

}



## ======================================================================
### estimate_lp

estimate_lp <- function(
  wrc, wc, head.wc, weights_wc, wrc_model,
  hcc, hc, head.hc, weights_hc, hcc_model,
  nlp.est, lp.fixed, fit.linear.param,
  param_bound, param_tf.k0, precBits, delta_sat_0,
  verbose
){

  ## function estimates the linear parameters thetar, thetas and k0
  ## either unconstrainededly by least squared or constrainedly
  ## (hounouring box constraints for thetar, thetas and k0 and the
  ## inequality thetar < thetas) by quadratic programming.

  ## 2019-11-27 A. Papritz
  ## 2022-01-06 AP adjust values of constrained estimates if outside
  ##               of allowed bounds


#### -- water retention function

  if(wrc){

    ## wrc data available, estimation of thetar and thetas - thetar

    ## prepare data (water content, relative saturation, weights)

    d <- data.frame(
      wc = wc,
      sat = as.double(
        sat_model(head.wc, nlp.est, precBits,
          wrc_model)
      ),
      w = weights_wc
    )

    if(any(fit.linear.param[c("thetar", "thetas")])){

      ## estimate either thetar, thetas or both parameters

      ## design matrix

      XX <- model.matrix(~ 1 + sat, d)

      ## matrix and vector coding constraints for thetar and thetas

      Amat <- t(
        rbind(
          thetar.l       = c( 1,  0),  # thetar >= lower limit
          thetar.u       = c(-1,  0),  # thetar <= upper limit
          thetas.l       = c( 1,  1),  # thetas >= lower limit
          thetas.u       = c(-1, -1),  # thetas <= upper limit
          thetas.minus.thetar   = c( 0,  1)   # thetas - thetar >= 0
        )
      )
      rownames(Amat) <- c("thetar", "thetas")

      bvec <- c(
        thetar.l =  param_bound[["thetar"]][1],  # thetar >= lower limit
        thetar.u = -param_bound[["thetar"]][2],  # thetar <= upper limit
        thetas.l =  param_bound[["thetas"]][1],  # thetas >= lower limit
        thetas.u = -param_bound[["thetas"]][2],  # thetas <= upper limit
        thetas.minus.thetar = 0                  # theta.s - thetar >= 0
      )


#### ---  estimate thetar and thetas

      if(fit.linear.param["thetar"] && fit.linear.param["thetas"]){

        ## both parameters

        X <- XX
        y <- d[, "wc"]

        if(all(is.finite(unlist(param_bound[c("thetar", "thetas")])))){

          ## constrained estimation

          WX <- d[, "w"] * X

          Dmat <- crossprod(X, WX)
          dvec <- drop(t(WX) %*% y)

          (fit <- try(
              solve.QP(Dmat, dvec, Amat, bvec),
              silent = TRUE
            ))

          if(identical(class(fit), "try-error")){

            mean.y <- mean(y)
            if(all(abs(d$sat - 0.) <= delta_sat_0)){

              ## estimate thetar and thetas for zero saturation

              thetar <- mean.y
              thetas <- 1.

            } else if(all(abs(d$sat - 1.) <= delta_sat_0)){

              ## estimate thetar and thetas for full saturation

              thetar <- 0.
              thetas <- mean.y

            } else {
              message <- paste(
                "an error occurred when estimating 'thetar' and 'thetas': \n",
                as.character(fit),
                "\nvalues of nonlinear parameters:\n",
                paste(paste(names(nlp.est), nlp.est), collapse = ", "), "\n"
              )
              cat(message)
              cat("\ndata:\n")
              print(d)
              cat("\nabs(saturation - 1):\n")
              print(abs(d$sat - 1.))
              cat("\ndelta.sat.0:", delta_sat_0, "\n")
              stop(message)
            }

          } else {

            ## estimate estimate thetar and thetas for  0 < saturation < 1

            thetar <- fit[["solution"]][1]
            thetas <- sum(fit[["solution"]])

          }

          se.thetar <- NA_real_
          se.thetas <- NA_real_

        } else {

          tmp <- try(chol(t(X) %*% X), silent = TRUE)

          if(identical(class(tmp), "try-error")){

            ## rank-deficient design matrix

            message <- paste(
              "an error occurred when estimating 'thetar' and 'thetas': \n",
              as.character(tmp),
              "\nvalues of nonlinear parameters:\n",
              paste(paste(names(nlp.est), nlp.est), collapse = ", "), "\n"
            )
            cat(message)
            cat("\ndata:\n")
            print(d)
            cat("design matrix:\n")
            print(X)
            stop(message)

          }

          ## unconstrained estimation

          fit <- lm(y ~ X - 1, weights = d[, "w"])
          (fit.coef <- coef(fit))
          fit.vcov <- vcov(fit)

          thetar <- fit.coef[1]
          thetas <- sum(fit.coef)
          se.thetar <- sqrt(fit.vcov[1, 1])
          se.thetas <- sqrt(sum(fit.vcov))

        }


      } else if(fit.linear.param["thetar"] && !fit.linear.param["thetas"]){

#### ---  estimate only thetar

        X <- 1. - XX[, 2, drop = FALSE]
        y <- with(d, wc - lp.fixed["thetas"] * sat)

        if(all(is.finite(unlist(param_bound[c("thetar")])))){

          ## constrained estimation

          sel <- c("thetar.l", "thetar.u")
          WX <- d[, "w"] * X

          Dmat <- crossprod(X, WX)
          dvec <- drop(t(WX) %*% y)

          (fit <- try(
              solve.QP(
                Dmat, dvec,
                Amat["thetar", sel, drop = FALSE],
                bvec[sel]
              ), silent = TRUE
            )
          )

          if(identical(class(fit), "try-error")){

            mean.y <- mean(y)

            if(
              all(abs(d$sat - 0.) <= delta_sat_0) &&
              mean.y < lp.fixed["thetas"]
            ){

              ## estimate thetar for zero saturation

              thetar <- mean.y

            } else if(all(abs(d$sat - 1.) <= delta_sat_0)){

              ## estimate thetar for full saturation

              thetar <- 0.

            } else {
              message <- paste(
                "an error occurred when estimating 'thetar' and 'thetas': \n",
                as.character(fit),
                "\nvalues of nonlinear parameters:\n",
                paste(paste(names(nlp.est), nlp.est), collapse = ", "), "\n"
              )
              cat(message)
              cat("\ndata:\n")
              print(d)
              cat("\nabs(saturation - 1):\n")
              print(abs(d$sat - 1.))
              cat("\ndelta.sat.0:", delta_sat_0, "\n")
              stop(message)
            }

          } else {

            ## estimate estimate thetar for  0 < saturation < 1

            thetar <- fit[["solution"]]

          }

          se.thetar <- NA_real_

        } else {

          ## unconstrained estimation

          tmp <- try(chol(t(X) %*% X))

          if(identical(class(tmp), "try-error")){

            ## rank-deficient design matrix

            message <- paste(
              "an error occurred when estimating 'thetar' and 'thetas': \n",
              as.character(tmp),
              "\nvalues of nonlinear parameters:\n",
              paste(paste(names(nlp.est), nlp.est), collapse = ", "), "\n"
            )
            cat(message)
            cat("\ndata:\n")
            print(d)
            cat("design matrix:\n")
            print(X)
            stop(message)

          }

          (fit <- lm(y ~ X - 1, weights = d[, "w"]))
          thetar <- coef(fit)
          se.thetar <- c(sqrt(vcov(fit)))

        }

        thetas <- lp.fixed["thetas"]
        se.thetas <- 0.


      } else if(!fit.linear.param["thetar"] && fit.linear.param["thetas"]){

#### ---  estimate only thetas

        X <- XX[, 2, drop = FALSE]
        y <- with(d, wc - lp.fixed["thetar"] * (1. - sat))

        if(all(is.finite(unlist(param_bound[c("thetas")])))){

          ## constrained estimation

          sel <- c("thetas.l", "thetas.u")
          WX <- d[, "w"] * X

          Dmat <- crossprod(X, WX)
          dvec <- drop(t(WX) %*% y)

          (fit <- try(
              solve.QP(
                Dmat, dvec,
                Amat["thetas", sel, drop = FALSE],
                bvec[sel]
              )
            )
          )

          if(identical(class(fit), "try-error")){

            mean.y <- mean(y)

            if(all(abs(d$sat - 0.) <= delta_sat_0)){

              ## estimate thetar for zero saturation

              thetas <- 1.

            } else if(
              all(abs(d$sat - 1.) <= delta_sat_0) &&
              mean.y > lp.fixed["thetar"]
            ){

              ## estimate thetar for full saturation

              thetas <- mean.y

            } else {
              message <- paste(
                "an error occurred when estimating 'thetar' and 'thetas': \n",
                as.character(fit),
                "\nvalues of nonlinear parameters:\n",
                paste(paste(names(nlp.est), nlp.est), collapse = ", "), "\n"
              )
              cat(message)
              cat("\ndata:\n")
              print(d)
              cat("\nabs(saturation - 1):\n")
              print(abs(d$sat - 1.))
              cat("\ndelta.sat.0:", delta_sat_0, "\n")
              stop(message)
            }

          } else {

            thetas <- fit[["solution"]]
          }
          se.thetas <- NA_real_

        } else {

          ## unconstrained estimation

          tmp <- try(chol(t(X) %*% X))

          if(identical(class(tmp), "try-error")){

            ## rank-deficient design matrix

            message <- paste(
              "an error occurred when estimating 'thetar' and 'thetas': \n",
              as.character(tmp),
              "\nvalues of nonlinear parameters:\n",
              paste(paste(names(nlp.est), nlp.est), collapse = ", "), "\n"
            )
            cat(message)
            cat("\ndata:\n")
            print(d)
            cat("design matrix:\n")
            print(X)
            stop(message)
          }

          (fit <- lm(y ~ X - 1, weights = d[, "w"]))
          thetas <- coef(fit)
          se.thetas <- c(sqrt(vcov(fit)))

        }

        thetar <- lp.fixed["thetar"]
        se.thetar <- 0.

      }

      ## adjust values of thetar and thetas if estimates are outside of
      ## allowed bounds

      if(thetar < param_bound[["thetar"]][1] || thetar > thetas){
        thetar <- min(max(thetar, param_bound[["thetar"]][1]), thetas)
      }

      if(thetas > param_bound[["thetas"]][2] || thetas < thetar){
        thetas <- max(min(thetas, param_bound[["thetas"]][2]), thetar)
      }

    } else {

      ## thetar and thetas both fixed

      thetar <- lp.fixed["thetar"]
      thetas <- lp.fixed["thetas"]
      se.thetar <- 0.
      se.thetas <- 0.

    }

  } else {

    ## no wrc data

    thetar    <- NULL
    se.thetar <- NULL
    thetas    <- NULL
    se.thetas <- NULL

  }


#### -- hydraulic conductivity function

  if(hcc){

    ## hc data available

    ## prepare data (hydraulic conductivity, relative conductivity, weights)

    krel <- hcrel_model(head.hc, nlp.est, precBits, hcc_model)


    if(any(!is.finite(log(krel)))) stop(
      "log(relative conductivity) (partly) non-finite:\n",
      "values of nonlinear parameters:\n",
      paste(paste(names(nlp.est), nlp.est), collapse = ", "), "\n"
    )

    d <- data.frame(
      hc = hc,
      w = weights_hc
    )

    if(fit.linear.param["k0"]){

      if(identical(param_tf.k0, "log")){

#### --- estimation of log(k0)

        tmp <- as.double(with(
            d,
            sum(w * log(hc / krel)) / sum(w)
          ))
        var.tmp <- with(
          d,
          var(as.double(log(hc / krel))) * sum(w^2) / (sum(w))^2
        )

        ## back transformation

        k0 <- exp(tmp)
        se.k0 <- k0 * sqrt(exp(var.tmp) - 1.)

        ## handle non-finite estimates of k0

        if(!is.finite(k0)){
          k0 <- sqrt(ifelse(tmp > 0., .Machine$double.xmax, .Machine$double.xmin))
          se.k0 <- NA_real_
        }

      } else {

#### --- estimation of k0

        ## design matrix

        X <- model.matrix(~ -1 + as.double(krel), d)
        y <- d[, "hc"]

        if(any(is.finite(param_bound[["k0"]]))){

          ## constrained estimation

          ## matrix and vector coding constraints for k0

          ## lower bound

          Amat <- t(
            rbind(
              k0.l = c( 1)
            )
          )
          rownames(Amat) <- c("k0")
          bvec <- c(
            k0.l =  param_bound[["k0"]][1]
          )

          if(is.finite(param_bound[["k0"]][2])){

            ## finite upper bound

            Amat <- cbind(Amat, k0.u = -1)
            bvec <- c(bvec, -param_bound[["k0"]][2])

          }

          WX <- d[, "w"] * X

          Dmat <- crossprod(X, WX)
          dvec <- drop(t(WX) %*% y)

          (fit <- try(
              solve.QP(Dmat, dvec, Amat, bvec)
            )
          )

          if(identical(class(fit), "try-error")){
            message <- paste(
              "an error occurred when estimating 'thetar' and 'thetas': \n",
              as.character(fit),
              "\nvalues of nonlinear parameters:\n",
              paste(paste(names(nlp.est), nlp.est), collapse = ", "), "\n"
            )
            cat(message)
            print(nlp.est)
            cat("\ndata:\n")
            print(d)
            stop(message)

          } else {

            k0 <- fit[["solution"]]

          }

          se.k0 <- NA_real_

          ## adjust values of thetar and thetas if estimates are outside of
          ## allowed bounds

          if(k0 < param_bound[["k0"]][1] || k0 > param_bound[["k0"]][2]){
            k0 <- min(max(k0, param_bound[["k0"]][1]), param_bound[["k0"]][2])
          }

        } else {

          ## unconstrained estimation

          (fit <- lm(y ~ X - 1, weights = d[, "w"]))

          k0 <- coef(fit)
          se.k0 <- sqrt(vcov(fit))

        }

      }

    } else {

      ## k0 fixed

      k0 = lp.fixed["k0"]
      se.k0 <- 0.

    }

  } else {

    ## no hc data

    k0     <- NULL
    se.k0  <- NULL

  }

#### -- collect and return results

  list(
    param = c(
      thetar = unname(thetar), thetas = unname(thetas),
      k0 = unname(k0)
    ),
    se = c(
      thetar = unname(se.thetar), thetas = unname(se.thetas),
      k0 = unname(se.k0)
    )
  )

}



## ======================================================================
### fit_wrc_hcc_fit

fit_wrc_hcc_fit <- function(
  i,
  input.data, param, fit_param,
  e0, ratio_lc_lt_bound,
  wrc, wrc_formula, wrc_mf,
  hcc, hcc_formula, hcc_mf,
  all.param.name,
  control,
  common.env,
  verbose
){

  # function estimate parameters of van Genuchten-Mulalem (vGM) model from
  # data on water content and hydraulic conductivity for a single soil
  # sample

  ## 2019-11-27 A. Papritz
  ## 2019-11-30 A. Papritz check of ranges of head, water content and
  ##               hydraulic conductivity data
  ## 2020-01-06 AP computation of Hessian matrix at solution
  ## 2020-01-17 AP fitted values and residuals in output
  ## 2020-01-27 AP ml and mpd estimation
  ## 2021-02-22 AP correction of error when processing empty Hessian matrix
  ## 2021-02-26 AP correction of error for insufficient number of measurements
  ## 2022-01-17 AP changes for processing values of nonlinear parameters

#### -- get model. frame and required data items

  warn.message <- ""

  #   ... for water retention curve

  wc <- head.wc <- weights_wc <- NULL

  if(wrc){

    wrc_mf <- eval(wrc_mf, environment())

    ## check whether the numer of data is sufficient

    if(NROW(wrc_mf) < control[["min_nobs_wc"]]){
      warn.message <- "not enough data for fitting model to water retention curve"
      wrc <- FALSE
    }

    ## setting-up terms objects

    mt.wc <- terms(wrc_formula)

    ## check whether 'empty' models have been entered

    if(is.empty.model(mt.wc)) stop(
      "an 'empty' water retention curve has been provided"
    )

    if(wrc){

      ## extract response variable and weight

      wc <- model.response(wrc_mf, "numeric")
      if(any(range(wc) < 0. | range(wc) > 1.)) stop(
        "volumetric water content data not restricted to [0, 1]"
      )
      head.wc <- as.vector(model.matrix(update(wrc_formula, . ~ . -1), wrc_mf))
      if(any(head.wc < 0.)) stop(
        "negative head for water retention curve"
      )
      weights_wc <- as.vector(model.weights(wrc_mf))
      if(any(weights_wc < 0.)) stop(
        "negative weights for water retention curve"
      )

      if(is.null(weights_wc)) weights_wc = rep(1., length(wc))

    }

  }

  ##   ... for conductivity function

  hc <- head.hc <- weights_hc <- NULL

  if(hcc){

    hcc_mf <- eval(hcc_mf, environment())

    ## check whether the number of data is sufficient

    if(NROW(hcc_mf) < control[["min_nobs_hc"]]){
      tmp <- "not enough data for fitting model to hydraulic conductivity curve"
      warn.message <- if(nchar(warn.message)){
        paste0(warn.message, "; ", tmp)
      } else tmp
      hcc <- FALSE
    }

    ## setting-up terms objects

    mt.hc <- terms(hcc_formula)

    ## check whether 'empty' models have been entered

    if(is.empty.model(mt.hc))
    stop("an 'empty' hydraulic conductivity curve has been provided")

    if(hcc){

      ## extract response variable and weight

      hc <- model.response(hcc_mf, "numeric")
      if(any(hc < 0.)) stop(
        "negative hydraulic conductivity data"
      )
      head.hc <- as.vector(model.matrix(update(hcc_formula, . ~ . -1), hcc_mf))
      if(any(head.hc < 0.)) stop(
        "negative head for hydraulic conductivity data"
      )
      weights_hc <- as.vector(model.weights(hcc_mf))
      if(any(weights_wc < 0.)) stop(
        "negative weights for ydraulic conductivity data"
      )

      if(is.null(weights_hc)) weights_hc = rep(1., length(hc))

    }

  }

  ## check whether neither wrc nor hcc data have been provided

  if(!wrc && !hcc){
    return(
      list(
        converged = NULL,
        message = paste0(
          "an error occurred when estimating parameters for sample ", i,
          "\n ", warn.message, "\n"
        ),
        objective = NULL,
        gradient = NULL,
        evaluation = NULL,
        lp = NULL, nlp = NULL,
        inequality_constraint = NULL,
        variable_weight = NULL,
        weight = NULL,
        model = NULL,
        initial_objects = NULL
      )
    )
  }

#### -- prepare initial parameter values

#### --- wrc

  ## get names of nonlinear and linear parameters

  tmp <- model_names_nlp_lp(control[["wrc_model"]])

  names.nlp.wrc <- tmp[["names.nlp"]]
  names.lp.wrc  <- tmp[["names.lp"]]

  if(wrc){

    ## water retention data available

    sel.missing.nlp.wrc <- !names.nlp.wrc %in% names(param)

    if(any(sel.missing.nlp.wrc)){

      ## initial values are missing for some nonlinear parameters

      names.missing.nlp.wrc <- names.nlp.wrc[sel.missing.nlp.wrc]

      if(
        identical(control[["wrc_model"]], "vg") &
        length(grep("local", control[["settings"]]))
      ){

        ## initial values of nonlinear parameters alpha and n for Van
        ## Genuchten model and local algorithm

        ## compute missing initial values for alpha | n for local
        ## optimization algorithm

        if(verbose >= 0.) warning(
          "no initial value(s) provided for parameter(s) '",
          paste(names.missing.nlp.wrc, collapse= "', '"),
          "': initial values will be computed"
        )

        ## fit GAM sat ~ s(log(h))

        sel <- head.wc > 0
        log.head <- log(head.wc[sel])

        sat <- ((wc - min(wc)) / diff(range(wc)))[sel]

        fit.gam <- try(
          gam(sat ~ s(log.head, k = min(control[["gam_k"]], length(wc)))),
          silent = TRUE
        )

        if(identical(class(fit.gam), "try-error")){

          ## fitting gam failed, use default initial values
          if(verbose >= 0.) warning(
            "failed to compute initial values for '",
            paste0(names.missing.nlp.wrc, collapse = "', '"),
            "': default initial values will be used"
          )
          for(i in names.missing.nlp.wrc){
            param[i] <- unname(control[["initial_param"]][i])
          }

        } else {

          ## compute initial values from inflection point of fitted gam
          ## curve

          ## predict saturation from gam

          new.data <- data.frame(
            log.head = seq(
              min(log.head), max(log.head), length = control[["gam_n_newdata"]]
            )
          )
          new.data$sat <- predict(fit.gam, newdata = new.data)

          ## find head and saturation where wrc is steepest (inflection
          ## point)

          i.steep <- which.min(diff(new.data$sat))
          if(i.steep < NROW(new.data) - 1L) i.steep <- c(i.steep, i.steep + 2)

          sat.steep  <- mean(new.data$sat[i.steep])
          head.steep <- exp(mean(new.data$log.head[i.steep]))

          if(!"n" %in% names(param)){

            ## use sat.steep to compute missing initial value for n

            eps <- sqrt(.Machine$double.eps)
            m.initial <- try(
              uniroot(
                function(m, s) s - (1 + 1/m)^(-m),
                interval = c(eps, 1), s = sat.steep
              ), silent = TRUE
            )

            if(identical(class(m.initial), "try-error")){
              if(verbose >= 0.) warning(
                "failed to compute initial value for 'n': default initial value will be used"
              )
              param["n"] <- unname(control[["initial_param"]]["n"])
            } else {
              param["n"] <- 1. / (1 - m.initial$root)
            }
          }

          m.initial <- 1. - (1. / param["n"])

          ## use head.steep, and n to compute initial value for alpha

          if(!"alpha" %in% names(param)){
            param["alpha"] <- 1. / head.steep * (1. / m.initial)^(1. - m.initial)
          }

        }

      } else {

        ## initial values of nonlinear parameters for other models or
        ## global algorithms

        if(verbose >= 0.) warning(
          "no initial value(s) provided for parameter(s) '",
          paste0(names.missing.nlp.wrc, collapse= "', '"),
          "': default initial values will be used"
        )

        for(i in names.missing.nlp.wrc){
          param[i] <- unname(control[["initial_param"]][i])
        }

      }

    }

  } else {

    ## no water retention data available, drop unnecessary linear
    ## parameters

    param <- param[!names(param) %in% names.lp.wrc]
    fit_param <- fit_param[!names(fit_param) %in% names.lp.wrc]

  }


#### --- hcc

  ## get names of nonlinear parameters

  tmp <- model_names_nlp_lp(control[["hcc_model"]])

  names.nlp.hcc <- tmp[["names.nlp"]]
  names.lp.hcc  <- tmp[["names.lp"]]

  if(hcc){

    ## conductivity data available

    sel.missing.nlp.hcc <- !names.nlp.hcc %in% names(param)

    if(any(sel.missing.nlp.hcc)){

      ## initial values are missing for some nonlinear parameters

      names.missing.nlp.hcc <- names.nlp.hcc[sel.missing.nlp.hcc]

      ## set missing initial values of nonlinear parameters to defaults

      if(verbose >= 0.) warning(
        "no initial value(s) provided for parameter(s) '",
        paste(names.missing.nlp.hcc, collapse= "', '"),
        "': default initial values will be used"
      )
      for(i in names.missing.nlp.hcc){
        param[i] <- unname(control[["initial_param"]][i])
      }

    }

    param.lc <- NULL

  } else {

    if(identical(control[["hcc_model"]], "vgm")){

      ## specifiy values of k0 and tau for computing capillary length Lc

      param.lc <- c(
        k0 = NA_real_,
        tau = unname(control[["initial_param"]]["tau"])
      )

      if("k0" %in% names(param)){
        param.lc["k0"] <- param["k0"]
      }

      if("tau" %in% names(param)){
        param.lc["tau"] <- param["tau"]
      }

    } else {

      param.lc <- NULL

    }

    ## drop unnecessary linear and nonlinear extra parameters for hcc

    param <- param[!names(param) %in% c(
      names.lp.hcc,
      names.nlp.hcc[!names.nlp.hcc %in% names.nlp.wrc]
    )]
    fit_param <- fit_param[!names(fit_param) %in% c(
      names.lp.hcc, names.nlp.hcc[!names.nlp.hcc %in% names.nlp.wrc]
    )]

  }

  #     print(param)
  #     print(fit_param)

#### -- further checks and preparations of initial values

  ## check mode

  if(!all(is.numeric(param))) stop(
    "initial values of parameters must be of mode 'numeric'"
  )
  if(!all(is.logical(fit_param))) stop(
    "fitting control flags of parameters must be of mode 'logical'"
  )

  ##  rearrange initial parameters

  param <- param[all.param.name[all.param.name %in% names(param)]]
  fit_param <- fit_param[all.param.name[all.param.name %in% names(fit_param)]]

  param.name <- names(param)
  fit_param.name <- names(fit_param)

  ## check whether initial values of parameters are valid

  if(wrc){

    if(!all(names.nlp.wrc %in% param.name)) stop(
      "some initial values of parameters of water retention function are missing"
    )
    if(!all(c(names.lp.wrc, names.nlp.wrc) %in% fit_param.name)) stop(
      "some fitting flags of parameters of water retention function are missing"
    )

    lapply(
      c(names.nlp.wrc, names(param)[names(param) %in% names.lp.wrc]),
      function(x, param, bounds){
        if(param[x] < bounds[[x]][1] || param[x] > bounds[[x]][2]) stop(
          "\n  initial value of parameter '", x, "' not valid\n"
        )
      }, param = param, bounds = control[["param_bound"]]
    )

  }

  if(hcc){

    if(!all(names.nlp.hcc %in% param.name)) stop(
      "some initial values of parameters of hydraulic conductivity function are missing"
    )

    if(!all(c(names.lp.hcc, names.nlp.hcc) %in% fit_param.name)) stop(
      "some fitting flags of parameters of hydraulic conductivity function are missing"
    )

    lapply(
      c(names.nlp.hcc, names(param)[names(param) %in% names.lp.hcc]) ,
      function(x, param, bounds){
        if(param[x] < bounds[[x]][1] || param[x] > bounds[[x]][2]) stop(
          "\n  initial value of parameter '", x, "' not valid\n"
        )
      }, param = param, bounds = control[["param_bound"]]
    )

  }

#### -- parameter transformations

  ## check parameter transformation for k0

  if(
    !is.null(control[["param_tf"]][["k0"]]) &&
    !control[["param_tf"]][["k0"]] %in% c("identity", "log")
  ) stop(
    "only log-transformation implemented for 'k0' parameter"
  )

  ##  preparation for parameter transformations

  nlp.name <- param.name[!param.name %in% c(names.lp.wrc, names.lp.hcc)]

  all.param_tf <- control[["param_tf"]]

  t.sel <- match(nlp.name, names(all.param_tf))

  if(any(is.na(t.sel))){
    stop("transformation undefined for some parameters")
  } else {
    control[["param_tf"]] <- all.param_tf[t.sel]
  }
  names(control[["param_tf"]]) <- nlp.name

  ## final checks

  stopifnot(all(names(fit_param)[!fit_param] %in% names(param)))

  if(wrc){
    stopifnot(all(names.nlp.wrc %in% names(param)))
    stopifnot(all(c(names.lp.wrc, names.nlp.wrc) %in% names(fit_param)))
  }

  if(hcc){
    stopifnot(all(names.nlp.hcc %in% names(param)))
    stopifnot(all(c(names.lp.hcc, names.nlp.hcc) %in% names(fit_param)))
  }

  ## store initial parameter values

  initial_objects <- list(
    param = param,
    fit_param = fit_param,
    variable_weight = control[["variable_weight"]],
    param_bound = control[["param_bound"]],
    e0 = e0,
    ratio_lc_lt_bound = ratio_lc_lt_bound
  )

#### -- prepare case weights

  ## adjust case weights by variable weights if both wrc and hcc data are
  ## available

  if(wrc && hcc){

    bla <- switch(
      control[["method"]],

      ## multi-response maximum likelihood or maximum posterior density
      ## estimation

      "ml" =,
      "mpd" = control[["variable_weight"]] <- c("wrc" = 1., "hcc" = 1.),

      ## (weighted) least squares estimation

      "wls" = {

        ## scale variable weights by variances

        control[["variable_weight"]]["wrc"] <-
          control[["variable_weight"]]["wrc"] / var(wc)

        if(identical(all.param_tf[["k0"]], "log")){
          var.mean.k0 <- var(log(hc))
        } else {
          var.mean.k0 <- var(hc)
        }

        control[["variable_weight"]]["hcc"] <-
          control[["variable_weight"]]["hcc"] / var.mean.k0

      },
      "unknown estimation method"

    )

    ## adjust case weights weights

    weights_wc <- weights_wc * control[["variable_weight"]]["wrc"]
    weights_hc <- weights_hc * control[["variable_weight"]]["hcc"]

  }


#### -- estimate parameters

  ## initialization

  ## linear parameters

  lp.name <- NULL
  if(wrc) lp.name <- c("thetar", "thetas")
  if(hcc) lp.name <- c(lp.name, "k0")

  linear.param <- param[lp.name]
  linear.param[is.na(linear.param)] <- Inf
  names(linear.param) <- lp.name

  fit.linear.param <- fit_param[lp.name]

  ## nonlinear parameters

  nlp.est <- param[nlp.name]

  ## prepare logging output

  if(verbose >= 2. && any(fit_param[nlp.name]) && (wrc || hcc)){
    cat( "\n     ",
      format("Iteration", width = 10L, justify = "right"),
      format(
        c(
          names(nlp.est[names(nlp.est) %in% nlp.name]), lp.name,
          "Obj", "relDeltaObj", "lc/lt"
        ) , width = 12L, justify = "right"
      ),
      "\n", sep = ""
    )

  }

  nlp <- estimate_nlp(
    nonlinear.param = nlp.est,
    fit.nonlinear.param = fit_param[nlp.name],
    linear.param = linear.param,
    fit.linear.param = fit.linear.param,
    param_tf.k0 = all.param_tf[["k0"]],
    e0 = e0, ratio_lc_lt_bound = ratio_lc_lt_bound, param.lc = param.lc,
    wrc = wrc, wc = wc, head.wc = head.wc, weights_wc = weights_wc,
    wrc_model = control[["wrc_model"]],
    hcc = hcc, hc = hc, head.hc = head.hc, weights_hc = weights_hc,
    hcc_model = control[["hcc_model"]],
    control = control,
    common.env = common.env,
    verbose = verbose
  )

#### -- check whether inequality constraints are satisfied
#       in constrained estimation

  if(control[["settings"]] %in% c("cglobal", "clocal")){

  }



#### -- collect and return all results

  result <- list(
    converged             = nlp[["convergence"]],
    message               = nlp[["message"]],
    evaluation            = nlp[["evaluation"]],
    objective             = nlp[["objective"]],
    gradient              = nlp[["gradient"]],
    lp                    = nlp[["linear.param"]],
    nlp                   = nlp[["param"]],
    inequality_constraint = nlp[["inequality_constraint"]],
    hessian               = nlp[["hessian"]],
    variable_weight = control[["variable_weight"]],
    weight = list(
      weights_wc = weights_wc,
      weights_hc = weights_hc
    ),
    fitted = list(
      wrc = get("fitted.wc", common.env),
      hcc = get("fitted.hc", common.env)
    ),
    residuals = list(
      wrc = get("residuals.wc", common.env),
      hcc = get("residuals.hc", common.env)
    ),
    model = list(
      wrc = if(wrc) wrc_mf else NULL,
      hcc = if(hcc) hcc_mf else NULL
    ),
    initial_objects = initial_objects
  )

  #   print(str(result))

  result

}

## ======================================================================
gradient_nlp <- function(
  adjustable.param, adjustable.param.name,
  fixed.param, fixed.param.name,
  param.name, method,
  param_tf, bwd_tf, deriv_fwd_tfd = NULL, param_bound,
  fit.linear.param, param_tf.k0, precBits, delta_sat_0,
  wrc, wc, head.wc, weights_wc, wrc_model,
  hcc, hc, head.hc, weights_hc, hcc_model,
  common.env, verbose,
  back.transform = NULL,
  e0 = NULL, param.lc = NULL, approximation_alpha_k0 = NULL,
  ratio_lc_lt_bound = NULL
){

  # function computes gradient of objective function objective_nlp for
  # estimation of the parameters param^T=c(alpha, n, tau) of the van
  # Genuchten-Mulalem (vGM) model from data on water content and
  # hydraulic conductivity by finite differences

  ## 2019-11-27 A. Papritz
  ## 2020-01-27 AP ml and mpd estimation

  ## get number of gradient function calls

  n.grd        <- get("n.grd",        common.env)

  ## return NULL if no parameters are estimated

  if(!length(adjustable.param)) return(NULL)

  ## restore names of parameters

  names(adjustable.param) <- adjustable.param.name
  if(length(fixed.param)) names(fixed.param) <- fixed.param.name
  param <- c(adjustable.param, fixed.param)[param.name]

  ## store items required for numerical evaluation of gradient of
  ## objective_nlp with respect to untransformed nonlinear parameters in
  ## temporary enviroment

  t.env <- new.env()

  t.tau <- NULL
  if(hcc) t.tau <- param["tau"]

  assign("alpha",                 param["alpha"], envir = t.env)
  assign("n",                     param["n"], envir = t.env)
  assign("tau",                   t.tau, envir = t.env)

  assign("adjustable.param.name", adjustable.param.name, envir = t.env)
  assign("fixed.param",           fixed.param, envir = t.env)
  assign("fixed.param.name",      fixed.param.name, envir = t.env)
  assign("param.name",            param.name, envir = t.env)
  assign("method",                method, envir = t.env)

  assign("param_tf",              param_tf, envir = t.env)
  assign("bwd_tf",                bwd_tf, envir = t.env)
  assign("param_bound",           param_bound, envir = t.env)
  assign("fit.linear.param",      fit.linear.param, envir = t.env)
  assign("param_tf.k0",           param_tf.k0, envir = t.env)
  assign("precBits",              precBits, envir = t.env)
  assign("delta_sat_0",           delta_sat_0, envir = t.env)

  assign("wrc",                   wrc, envir = t.env)
  assign("wc",                    wc, envir = t.env)
  assign("head.wc",               head.wc, envir = t.env)
  assign("weights_wc",             weights_wc, envir = t.env)
  assign("wrc_model",             wrc_model, envir = t.env)

  assign("hcc",                   hcc, envir = t.env)
  assign("hc",                    hc, envir = t.env)
  assign("head.hc",               head.hc, envir = t.env)
  assign("weights_hc",             weights_hc, envir = t.env)
  assign("hcc_model",             hcc_model, envir = t.env)

  assign("common.env",            common.env, envir = t.env)

  ## numerically evaluate gradient of objective_nlp with respect to
  ## untransformed nonlinear parameters

  f.aux <- function(
    alpha, n, tau,
    adjustable.param.name,
    fixed.param, fixed.param.name,
    param.name, method,
    param_tf, bwd_tf, param_bound,
    fit.linear.param, param_tf.k0, precBits, delta_sat_0,
    wrc, wc, head.wc, weights_wc, wrc_model,
    hcc, hc, head.hc, weights_hc, hcc_model,
    common.env
  ){

    objective_nlp(
      adjustable.param = c(alpha, n, tau)[adjustable.param.name],
      adjustable.param.name = adjustable.param.name,
      fixed.param = c(alpha, n, tau)[fixed.param.name],
      fixed.param.name = fixed.param.name,
      param.name = param.name, method = method,
      param_tf = param_tf, bwd_tf = bwd_tf, param_bound = param_bound,
      fit.linear.param =  fit.linear.param, param_tf.k0 = param_tf.k0,
      precBits = precBits, delta_sat_0 = delta_sat_0,
      wrc = wrc, wc = wc, head.wc = head.wc,
      weights_wc = weights_wc, wrc_model = wrc_model,
      hcc = hcc, hc = hc, head.hc = head.hc,
      weights_hc = weights_hc, hcc_model = hcc_model,
      common.env = common.env, verbose = 0
    )
  }


  result <- numericDeriv(
    expr = quote(
      f.aux(
        alpha = alpha, n = n, tau = tau,
        adjustable.param.name = adjustable.param.name,
        fixed.param = fixed.param, fixed.param.name = fixed.param.name,
        param.name = param.name, method = method,
        param_tf = param_tf, bwd_tf = bwd_tf, param_bound = param_bound,
        fit.linear.param =  fit.linear.param, param_tf.k0 = param_tf.k0,
        precBits = precBits, delta_sat_0 = delta_sat_0,
        wrc = wrc, wc = wc, head.wc = head.wc,
        weights_wc = weights_wc, wrc_model = wrc_model,
        hcc = hcc, hc = hc, head.hc = head.hc,
        weights_hc = weights_hc, hcc_model = hcc_model,
        common.env = common.env
      )
    ),
    theta = adjustable.param.name,
    rho = t.env
  )

  grad <- as.vector(attr(result, "gradient"))
  names(grad) <- adjustable.param.name

  #   print(grad, digits = 14)

  ## update elements of common.env

  assign("n.grd", n.grd + 1L, envir = common.env)

  return(grad)

}


## ======================================================================
objective_nlp <- function(
  adjustable.param, adjustable.param.name,
  fixed.param, fixed.param.name,
  param.name, method,
  param_tf, bwd_tf, deriv_fwd_tfd = NULL, param_bound,
  fit.linear.param, param_tf.k0, precBits, delta_sat_0,
  wrc, wc, head.wc, weights_wc, wrc_model,
  hcc, hc, head.hc, weights_hc, hcc_model,
  common.env, verbose,
  back.transform = TRUE,
  e0 = NULL, param.lc = NULL, approximation_alpha_k0 = NULL,
  ratio_lc_lt_bound = NULL
){

  # function computes objective function for estimation of the parameters
  # param^T=c(alpha, n, tau) of the van Genuchten-Mulalem (vGM) model from
  # data on water content and hydraulic conductivity

  ## 2019-11-27 A. Papritz
  ## 2020-01-27 AP ml and mpd estimation
  ## 2021-10-13 AP correction of degrees of freedom for ml and mpd method

  ## get number iteration, linear parameters and value of objective of
  ## previous iteration

  n.it         <- get("n.it",         common.env)
  linear.param <- get("linear.param", common.env)
  obj.old      <- get("obj",          common.env)

  ## restore names of parameters

  names(adjustable.param) <- adjustable.param.name
  if(length(fixed.param)) names(fixed.param) <- fixed.param.name
  param <- c(adjustable.param, fixed.param)[param.name]

  ##  transform parameters back to original scale

  param <- sapply(
    param.name,
    function(x, bwd_tf, param_tf, param, bounds){
      bwd_tf[[param_tf[[x]]]](param[x], bounds[[x]][1], bounds[[x]][2])
    },
    bwd_tf = bwd_tf,
    param_tf = param_tf,
    param = param,
    bounds = param_bound
  )
  names(param) <- param.name

  ## optionally output parameters to console

  if(verbose >= 2. && length(adjustable.param) >= 1L){

    cat( "     ",
      format(n.it, width = 10L, justify = "right"),
      format(
        signif(
          param,
          digits = 5L
        ),
        scientific = TRUE, width = 12L
      ), sep = ""
    )

  }

  ## estimate linear parameters

  lp <- estimate_lp(
    wrc = wrc, wc = wc, head.wc = head.wc,
    weights_wc = weights_wc, wrc_model = wrc_model,
    hcc = hcc, hc = hc, head.hc = head.hc,
    weights_hc = weights_hc, hcc_model = hcc_model,
    nlp.est = param,
    lp.fixed = linear.param[!fit.linear.param],
    fit.linear.param = fit.linear.param,
    param_bound = param_bound, param_tf.k0 = param_tf.k0,
    precBits = precBits, delta_sat_0 = delta_sat_0,
    verbose = verbose
  )

  ## compute value of objective function

  if(verbose >= 3.) op <- par(mfrow = c(1, sum(c(wrc, hcc))))

  ## sum of squares for wrc

  if(wrc){

    n.wc <- length(wc)

    thetavgm <- as.double(wc_model(
      head.wc, param, lp[["param"]], precBits, wrc_model
    ))
    ssq_wc <- sum(
      weights_wc * (wc - thetavgm )^2
    )

    obj_wc <- switch(
      method,

      ## contribution to negative loglikelihood (marginal posterior
      ## density (\ propto negative loglikelihood), Stewart et al., 1992, eqs 6 & 11)
      "ml" = 0.5 * n.wc * log(ssq_wc),

      ## contribution to negative logarithm of modified posterior density
      ## (Stewart et al., 1992, eqs 7 & 12)
      "mpd" = 0.5 * (n.wc + 2L) * log(ssq_wc),

      ## contribution to weighted sum of squares
      "wls" = ssq_wc,

      "unknown estimation method"
    )

    if(verbose >= 3.){
      t.head <- exp(seq(min(log(head.wc)), max(log(head.wc)),
          length = 101))
      plot(head.wc, wc, log="x")
      lines(
        t.head,
        as.double(wc_model(
          t.head, param, lp[["param"]], precBits, wrc_model
        )), col = "grey"
      )

    }

  } else {
    ssq_wc <- NULL
    obj_wc <- 0.
  }

  ## sum of squares for hc

  if(hcc){

    n.hc <- length(hc)

    hcvgm <- hc_model(
      head.hc, param, lp[["param"]], precBits, hcc_model
    )
    if(any(!is.finite(log(hcvgm)))) stop(
      "conductivity (partly) non-finite: increase precBits for more accurate numerics"
    )

    if(identical(param_tf.k0, "log")){
      ssq_hc <- as.double(sum(weights_hc * log( hc / hcvgm )^2))
    } else {
      ssq_hc <- as.double(sum(weights_hc * (hc - hcvgm)^2))
    }

    obj_hc <- switch(
      method,

      ## contribution to negative loglikelihood (marginal posterior
      ## density (\ propto negative loglikelihood), Stewart et al., 1992, eqs 6 & 11)
      "ml" = 0.5 * n.hc * log(ssq_hc),

      ## contribution to negative logarithm of modified posterior density
      ## (Stewart et al., 1992, eqs 7 & 12)
      "mpd" = 0.5 * (n.hc + 2L) * log(ssq_hc),

      ## contribution to weighted sum of squares
      "wls" = ssq_hc,

      "unknown estimation method"
    )


    if(verbose >= 3.){

      t.head <- exp(seq(min(log(head.hc)), max(log(head.hc)),
          length = 101))
      plot(head.hc, hc, log="xy")
      lines(
        t.head,
        as.double(hc_model(
          t.head, param, lp[["param"]], precBits, hcc_model
        )), col = "grey"
      )

    }

  } else {
    ssq_hc <- NULL
    obj_hc <- 0.
  }

  if(verbose >= 3.) par(op)

  ## overall sum of squares

  obj <- obj_wc + obj_hc

  ## optionally output objective and ratio Lc/Lt etc to console

  if(verbose >= 2. && length(adjustable.param) >= 1L){

    ## compute values of Lc and inequality constraints

    tmp <- ineq_constraint_nlp(
      adjustable.param = adjustable.param,
      adjustable.param.name = adjustable.param.name,
      fixed.param = fixed.param, fixed.param.name = fixed.param.name,
      param.name = param.name,
      param_tf = param_tf,
      bwd_tf = bwd_tf,
      param_bound = param_bound,
      wrc_model = wrc_model, hcc_model = hcc_model,
      common.env = common.env,
      e0 = e0,
      param.lc = param.lc,
      approximation_alpha_k0 = approximation_alpha_k0,
      ratio_lc_lt_bound = ratio_lc_lt_bound,
      back.transform = back.transform
    )
    t.ratio <- attr(tmp, "lc") / attr(tmp, "lt")
    attributes(t.ratio) <- NULL

    cat(
      format(
        signif(
          c( linear.param, obj, (obj - obj.old) / abs(obj.old), t.ratio),
          digits = 5L
        ),
        scientific = TRUE, width = 12L
      ), "\n" , sep = ""
    )

  }

  ## computed fitted values and residuals

  if(wrc){
    fitted.wc <- thetavgm
    residuals.wc <- wc - thetavgm
  } else {
    obj_wc <- NULL
    fitted.wc <- NULL
    residuals.wc <- NULL
  }

  if(hcc){
    if(identical(param_tf.k0, "log")){
      fitted.hc <- as.double(log(hcvgm))
      residuals.hc <- as.double(log( hc / hcvgm ))
    } else {
      fitted.hc <- hcvgm
      residuals.hc <- hc - hcvgm
    }
  } else {
    obj_hc <- NULL
    fitted.hc <- NULL
    residuals.hc <- NULL
  }

  ## update elements of common.env

  n.it <- n.it + 1L

  assign("n.it",            n.it,             envir = common.env)
  assign("obj",             obj,              envir = common.env)
  assign("obj_wc",          obj_wc,           envir = common.env)
  assign("obj_hc",          obj_hc,           envir = common.env)
  assign("ssq_wc",          ssq_wc,           envir = common.env)
  assign("ssq_hc",          ssq_hc,           envir = common.env)
  assign("linear.param",    lp[["param"]],    envir = common.env)
  assign("fitted.wc",       fitted.wc,        envir = common.env)
  assign("residuals.wc",    residuals.wc,     envir = common.env)
  assign("fitted.hc",       fitted.hc,        envir = common.env)
  assign("residuals.hc",    residuals.hc,     envir = common.env)

  ## return obj

  return(obj)

}


################################################################################

stop_cluster <- function(){

  ## function to stop snowfall clusters

  ## 2019-11-27 A. Papritz

  if(sfIsRunning()){
    sfStop()
  }

  options(error = NULL)

}


################################################################################

### model-dependent functions

#### -- model_names_nlp_lp

model_names_nlp_lp <- function(
  model
){

  ## function to define names of nonlinear and linear model parameters
  ## called in model_fit_param_consistent

  ## 2022-01-09 A. Papritz

  switch(
    model,

    ## wrc_models

    vg = {
      ## Van Genuchten model
      names.nlp <- c("alpha", "n")
      names.lp  <- c("thetar", "thetas")
    },

    ## hcc_models

    vgm = {
      ## Van Genuchten-Mualem model
      names.nlp <- c("alpha", "n", "tau")
      names.lp  <- c("k0")
    },

    stop("model '", model, "' not implemented")

  )

  list(names.nlp = names.nlp, names.lp = names.lp)

}

#### -- model_fit_param_consistent

model_fit_param_consistent <- function(
  model,
  fit_param, names.fit_param,
  param, names.param,
  verbose
){

  ## function for consistent coding of flags for fitting for wrc_models
  ## called in fit_wrc_hcc

  ## 2022-01-13 A. Papritz

  ## define names of nonlinear and linear model parameters

  tmp <- model_names_nlp_lp(model)

  names.nlp <- tmp[["names.nlp"]]
  names.lp <- tmp[["names.lp"]]

  ## check and adjust fit_param for nonlinear parameters

  sel.missing.nlp <- !names.nlp %in% names.fit_param

  if(any(sel.missing.nlp)){

    names.missing.nlp <- names.nlp[sel.missing.nlp]

    if(verbose >= 0.) warning(
      "no fitting control provided for parameter(s) '",
      paste(names.missing.nlp, collapse= "', '"),
      "': parameters will be fitted"
    )
    for(i in names.missing.nlp){
      fit_param <- cbind(rep(TRUE, NROW(fit_param)), fit_param)
      colnames(fit_param) <- c(i, colnames(fit_param)[-1])
    }
    names.fit_param <- colnames(fit_param)

  }

  #   bla <- sapply(
  #     names.nlp,
  #     function(i){
  #       if(!all(fit_param[, i]) && (is.null(param) || !(i %in% names.param))) stop(
  #         "no value(s) provided for fixed '", i, "'"
  #       )
  #     }
  #   )

  ## check and adjust fit_param for linear parameters

  sel.missing.lp <- !names.lp %in% names.fit_param

  if(any(sel.missing.lp)){

    names.missing.lp <- names.lp[sel.missing.lp]

    if(verbose >= 0.) warning(
      "no fitting control provided for parameter(s) '",
      paste(names.missing.lp, collapse= "', '"),
      "': parameters will be fitted"
    )
    for(i in names.missing.lp){
      fit_param <- cbind(rep(TRUE, NROW(fit_param)), fit_param)
      colnames(fit_param) <- c(i, colnames(fit_param)[-1])
    }
    names.fit_param <- colnames(fit_param)

  }

  for(i in names.lp){

    if(!all(fit_param[, i]) && (is.null(param) || !(i %in% names.param))) stop(
      "no value provided for fixed '", i, "'"
    )

    ## exclude initial value of linear parameter if they are estimated

    if(!is.null(param) && all(fit_param[, i])){
      param <- param[, !names.param %in% i, drop = FALSE]
    }

  }

  ## return result

  list(fit_param = fit_param, param = param)

}

#### -- model_param_tf_nlp_identity

model_param_tf_nlp_identity <- function(
  model
){

  ## function to choose identity transformation for nonlinear model
  ## parameters
  ## called in control_fit_wrc_hcc

  ## 2022-01-08 A. Papritz

  switch(
    model,

    ## wrc_models

    ## Van Genuchten model
    vg = param_transf(
      alpha = "identity", n = "identity"
    )[c("alpha", "n")],

    ## hcc_models

    ## Van Genuchten-Mualem model
    vgm = param_transf(
      alpha = "identity", n = "identity", tau = "identity"
    )[c("alpha", "n", "tau")],

    stop("model '", model, "' not implemented")

  )

}
