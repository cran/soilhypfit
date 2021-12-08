# ## ======================================================================
print.fit_wrc_hcc <- function(x, ...){

  ## function extracts fitted coefficients from an object of class fit_wrc_hcc

  ## 2019-11-27 A. Papritz

  ## check consistency of object

  stopifnot(identical(class(x), "fit_wrc_hcc"))

  ## extract number of samples with wrc and/or hcc curves

  nsamp <- length(x[["fit"]])
  type.of.data <- sapply(
    x[["fit"]],
    function(xx) c(
      wrc= !is.null(xx[["model"]][["wrc"]]),
      hcc= !is.null(xx[["model"]][["hcc"]])
    )
  )

  nwrc <- sum(type.of.data[1, ] & !type.of.data[2, ])
  nhcc <- sum(!type.of.data[1, ] & type.of.data[2, ])
  nwrchcc <- sum(type.of.data[1, ] & type.of.data[2, ])

  ## print number of samples
  cat(
    "\nSoil hydraulic parameters estimated for", nsamp, "soil samples, of which",
    "\n  have only water content data                           :", nwrc,
    "\n  have only hydraulic conductivity data                  :", nhcc,
    "\n  have both water content and hydraulic conductivity data:", nwrchcc
  )

  invisible(x)

}


# ## ======================================================================
summary.fit_wrc_hcc <- function(object,
  what = c("all", "nonlinear", "linear"),
  subset = NULL, gof = TRUE, lc = TRUE, ...
){

  ## function extracts fitted coefficients from an object of class fit_wrc_hcc

  ## 2019-11-27 A. Papritz

  ## check consistency of object

  stopifnot(identical(class(object), "fit_wrc_hcc"))

  ## select subsets

  if(!is.null(subset)){
    stopifnot(is.character(subset) || is.numeric(subset) || is.logical(subset))
    stopifnot(length(subset) <= length(object[["fit"]]))
    object[["fit"]] <- object[["fit"]][subset]
  }

  ## extract number of samples with wrc and/or hcc curves

  nsamp <- length(object[["fit"]])
  type.of.data <- sapply(
    object[["fit"]],
    function(xx) c(
      wrc= !is.null(xx[["model"]][["wrc"]]),
      hcc= !is.null(xx[["model"]][["hcc"]])
    )
  )

  nwrc <- sum(type.of.data[1, ] & !type.of.data[2, ])
  nhcc <- sum(!type.of.data[1, ] & type.of.data[2, ])
  nwrchcc <- sum(type.of.data[1, ] & type.of.data[2, ])

  ## extract algorithmic details

  control <- object[["control"]][c(
    "settings", "nloptr", "sce",
    "approximation_alpha_k0",
    "param_bound", "param_tf"
  )]


  ## extract fitted parameters, goodness-of-fit and Lc

  result <- coef(object, subset = NULL, gof = gof, lc = lc, ...)

  ## collect all results

  ans <- list(
    data = c(
      nsamp = nsamp, nwrc = nwrc, nhcc = nhcc, nwrchcc = nwrchcc
    ),
    control = control,
    result = result,
    call = object[["call"]]
  )

  class(ans) <- "summary_fit_wrc_hcc"

  ans

}


# ## ======================================================================
print.summary_fit_wrc_hcc <- function(
  x, digits = max(5, getOption("digits") - 3), ...
){

  ## function extracts fitted coefficients from an object of class fit_wrc_hcc

  ## 2019-11-27 A. Papritz
  ## 2019-12-20 A. Papritz  print only a summary of result
  ## 2020-06-15 A. Papritz  output parameter space also for local
  ##                        algorithm and untransformed parameters

  ## call

  cat("\nCall:\n  ")
  cat( paste( deparse(x[["call"]]), sep = "\n", collapse = "\n"),  "\n", sep = "" )

  ## details on algorithm

  cat("\nEstimation method:\n")

  ## algorithm and convergence criteria

  switch(
    x[["control"]][["settings"]],
    cglobal = cat("  global algorithm, constrained"),
    clocal  = cat("  local algorithm, constrained"),
    uglobal = cat("  global algorithm, unconstrained"),
    sce     = cat("  global algorithm, unconstrained"),
    ulocal  = cat("  local algorithm, unconstrained")
  )
  cat(" optimization\n")

  switch(
    x[["control"]][["settings"]],
    sce = {
      cat("  algorithm shuffled Complex Evolution (SCE) optimisation\n")

      cat("  convergence criteria\n")
      bla <- lapply(
        c("reltol", "tolsteps", "maxeval", "maxtime"),
        function(y){
          if(!is.null(x[["control"]][["sce"]][[y]]) &&
            x[["control"]][["sce"]][[y]] > 0.) cat(
            "    ", y, x[["control"]][["sce"]][[y]], "\n"
          )
        }
      )
    },
    {
      cat("  algorithm", x[["control"]][["nloptr"]][["algorithm"]], "\n")

      cat("  convergence criteria\n")
      if(
        !is.null(x[["control"]][["nloptr"]][["stopval"]]) &&
        x[["control"]][["nloptr"]][["stopval"]] > -Inf) cat(
        "    stopval", x[["control"]][["nloptr"]][["stopval"]], "\n"
      )
      bla <- lapply(
        c("ftol_rel", "ftol_abs", "xtol_rel", "xtol_abs", "maxtime", "maxeval"),
        function(y){
          if(!is.null(x[["control"]][["nloptr"]][[y]]) &&
            x[["control"]][["nloptr"]][[y]] > 0.) cat(
            "    ", y , x[["control"]][["nloptr"]][[y]], "\n"
          )
        }
      )

    }

  )

  ## parameter transformations

  cat("  parameter transformations")
  sel <- names(x[["control"]][["param_tf"]]) %in% colnames(x[["result"]])
  bla <- lapply(
    names(x[["control"]][["param_tf"]][sel]),
    function(y) cat(
      "\n    ",
      format(y, width = 8L, justify = "right"),
      format(x[["control"]][["param_tf"]][[y]], width = 10L)
    )
  )

  #   ## constraints for constrained optimization
  #
  #   if(length(grep("^c", x[["control"]][["settings"]]))) cat(
  #     "\n  bounds of ratio lc/lt:",
  #     "\n    lower", x[["control"]][["ratio_lc_lt_bound"]]["lower"],
  #     "\n    upper", x[["control"]][["ratio_lc_lt_bound"]]["upper"]
  #   )

  ## parameter space for global optimization or local
  ## optimization with all nonlinear parameters untransformed

  if(
    length(grep("global$", x[["control"]][["settings"]])) || (
      length(grep("local$", x[["control"]][["settings"]])) &&
      all(
        unlist(x[["control"]][["param_tf"]])[c("alpha", "n", "tau")] ==
        "identity"
      )
    )
  ){
    cat("\n  parameter space")
    sel <- names(x[["control"]][["param_bound"]]) %in% colnames(x[["result"]])
    bla <- lapply(
      names(x[["control"]][["param_bound"]][sel]),
      function(y) cat(
        "\n    ",
        format(y, width = 8L, justify = "right"),
        format(signif(x[["control"]][["param_bound"]][[y]], digits = 3L), width = 10L)
      )
    )
  }

  ## details on data

  cat(
    "\n\nSoil data:",
    "\n  # of samples with water content data                           :", x[["data"]]["nwrc"],
    "\n  # of samples with hydraulic conductivity data                  :", x[["data"]]["nhcc"],
    "\n  # of samples with water content and hydraulic conductivity data:", x[["data"]]["nwrchcc"]  )

  ## estimated parameters and goodness-of-fit

  cat("\n\nSummary of estimated parameters:\n")
  rownames(x[["result"]]) <- paste0("  ", rownames(x[["result"]]))
  print(summary(x[["result"]]), digits = digits)
  invisible(x)

}


# ## ======================================================================
vcov.fit_wrc_hcc <- function(
  object, subset = NULL,
  grad_eps, bound_eps = sqrt(.Machine$double.eps), ...
){

  ## function computes covariance matrix of estimated nonlinear parameters
  ## from a fit_wrc_hcc object

  ## 2020-02-03 A. Papritz
  ## 2020-02-26 A. Papritz correction error when computing standard errors
  ##                       for nonlinear parameters when some parameters
  ##                       are fixed
  ## 2021-10-13 A. Papritz correction of error in checking whether estimates
  ##                       are on boundary of parameter space or whether the
  ##                       gradient is not approximately equal to zero

  #   ## check consistency of object
  #
  #   stopifnot(identical(class(object), "fit_wrc_hcc"))

  ## check arguments

  if(!missing(grad_eps)) stopifnot(
    identical(length(grad_eps), 1L) && is.numeric(grad_eps) && grad_eps > 0
  )
  stopifnot(identical(length(bound_eps), 1L) && is.numeric(bound_eps) && bound_eps > 0)

  ## select subsets

  if(!is.null(subset)){
    stopifnot(is.character(subset) || is.numeric(subset) || is.logical(subset))
    stopifnot(length(subset) <= length(object[["fit"]]))
    object[["fit"]] <- object[["fit"]][subset]
  }

  ## check settings and method

  if(object[["control"]][["settings"]] %in% c("cglobal", "clocal")) warning(
    "computing covariance matrix of nonlinear parameters for constrained estimates"
  )

  if(!object[["control"]][["method"]] %in% c("ml", "mpd")) stop(
    "covariance matrix of nonlinear parameters only available for 'mpd' or 'ml' estimation"
  )

  if(is.null(object[["fit"]][[1]][["hessian"]])) stop(
    "covariance matrix cannot be computed because hessian matrix is missing"
  )

  ## set value for grad_eps if missing

  if(missing(grad_eps)) grad_eps <-
    object[["control"]][["grad_eps"]]

  ## extract boundaries of parameter space

  param_bound <- object[["control"]][["param_bound"]]

  ## compute covariance matrices, cf Stewart et al., 1981, section 3.4

  result <- lapply(
    1:length(object[["fit"]]),
    function(i){

      ## select current sample

      x <- object[["fit"]][[i]]

      if(length(x[["hessian"]]) > 0L){

        ## cholesky decomposition of hessian matrix

        tmp <- try(chol(x[["hessian"]]), silent = TRUE)

        ## compute and return covariance matrix

        if(!identical(class(tmp), "try-error")){

          status <- "hessian positive definite"

          ## invert positive definite hessian matrix

          res <- chol2inv(tmp)
          dimnames(res) <- dimnames(x[["hessian"]])

          ## set (co-)variances to missing for parameters with values on the
          ## boundary of the parameter space

          sel <- which(
            sapply(
              colnames(res),
              function(j){
                any(abs(param_bound[[j]] - x[["nlp"]][j]) < bound_eps)
              }
            )
          )

          if(length(sel)){

            status <- paste(
              status, "some estimates are on boundary of parameter space", sep = ", "
            )
            warning("sample ", names(object[["fit"]][i]), ": ", status)

            res[sel, ] <- NA_real_
            res[, sel] <- NA_real_

          }

          ## check whether scaled gradient is approximately equal to zero

          if(length(colnames(res))){

            ## exclude parameters on boundary of parameter space
            if(length(sel)){
              tmp <- x[["gradient"]][-sel]
            } else {
              tmp <- x[["gradient"]]
            }

            if(max(abs(tmp / x[["objective"]])) > grad_eps){
              status <- paste(
                status, paste("max(abs(scaled gradient)) >", grad_eps), sep = ", "
              )
              warning("sample ", names(object[["fit"]][i]), ": ", status)

            }
          }


        } else {

          status <- "hessian not positive definite"
          warning("sample ", names(object[["fit"]][i]), ": ", status)

          res <- matrix(
            NA_real_, nrow = NROW(x[["hessian"]]), ncol = NCOL(x[["hessian"]])
          )
          dimnames(res) <- dimnames(x[["hessian"]])

        }

      } else {

        status <- "hessian empty (all nonlinear parameters fixed)"
        warning("sample ", names(object[["fit"]][i]), ": ", status)

        res <- matrix(
          NA_real_, nrow = NROW(x[["hessian"]]), ncol = NCOL(x[["hessian"]])
        )
        dimnames(res) <- dimnames(x[["hessian"]])

      }

      attr(res, "status") <- status

      res

    }
  )

  names(result) <- names(object[["fit"]])
  class(result) <- "vcov_fit_wrc_hcc"

  result

}


## ======================================================================
coef.vcov_fit_wrc_hcc <- function(object, se = TRUE,
  correlation = se, status = FALSE, ...
){

  ## function extracts variances (or standard errors) and covariances (or
  ## correlations) of nonlinear parameters from a fitted vcov_fit_wrc_hcc
  ## object

  ## 2020-01-31 A. Papritz
  ## 2020-06-10 A. Papritz change of default value for correlation
  ## 2020-06-17 A. Papritz correction of error
  ## 2020-02-26 A. Papritz correction error when computing standard errors
  ##                       for nonlinear parameters when some parameters
  ##                       are fixed

  ## check arguments

  stopifnot(identical(length(se), 1L) && is.logical(se))
  stopifnot(identical(length(correlation), 1L) && is.logical(correlation))
  stopifnot(identical(length(status), 1L) && is.logical(status))

  ## get names of parameters from largest hessian matrix

  ## names for variances

  dims <- sapply(object, function(x) prod(dim(x)))
  names.var <- colnames(object[[which.max(dims)]])

  ## names for covariances/correlations

  if(length(names.var) > 1L){
    tmp <- outer(names.var, names.var, function(x,y) paste(x, y, sep = "."))
    names.cov <- tmp[upper.tri(tmp)]
  } else {
    names.cov <- NULL
  }


  ## prefixes for output

  prefix.se <- if(se) "se." else "var."
  prefix.corr <- if(correlation) "corr." else "cov."

  tmp <- lapply(
    object,
    function(x, correlation){

      diag.elements <- rep(NA_real_, length(names.var))
      names(diag.elements) <- names.var

      offdiag.elements <- rep(NA_real_, length(names.cov))
      names(offdiag.elements) <- names.cov

      ## variances or standard errors

      if(prod(dim(x))){

        ## variances/standard errors

        tmp <- diag(x)
        if(se) tmp <- sqrt(tmp)
        diag.elements[names(tmp)] <- tmp

        ## covariances/correlations

        if(NCOL(x) > 1L){

          ## temporarily suppress warnings for computing correlations

          if(correlation) oo <- options(warn = -1)

          ## covariances or correlations

          sel <- upper.tri(x)
          tmp <- if(correlation){
            cov2cor(x)[sel]
          } else {
            x[sel]
          }
          names.tmp <- outer(
            colnames(x),
            colnames(x),
            function(x, y) paste(x, y, sep = ".")
          )
          names(tmp) <- names.tmp[upper.tri(names.tmp)]

          offdiag.elements[names(tmp)] <- tmp

          ## restore warnings

          if(correlation) options(oo)

        }

      }

      ## set names

      if(length(diag.elements) > 0L){
        names(diag.elements)    <- paste0(prefix.se,   names(diag.elements))
      }

      if(length(offdiag.elements) > 0L){
        names(offdiag.elements) <- paste0(prefix.corr, names(offdiag.elements))
      }


      ## return result

      res <- list(
        diag.elements,
        offdiag.elements
      )

      if(status){
        res <- c(res, attr(x, "status"))
      }

      res

    }, correlation = correlation
  )

  ## convert to dataframe

  result <- as.data.frame(t(sapply(
        tmp,
        function(x){
          if(length(x[[1]]) > 0L) {
            unlist(x[1:2])
          } else NULL
        }
      )))


  if(identical(length(tmp[[1]]), 3L)){
    result[, "status"] <- sapply(
      tmp,
      function(x) x[[3]]
    )
  }

  result

}


# # ## ======================================================================
# print.ci_fit_wrc_hcc <- function(x,  ...
# ){
#
#   ## function prints confidence intervals for nonlinear parameters from a
#   ## fitted fit_wrc_hcc object
#
#   ## 2020-01-27 A. Papritz
#
#   #   ## check consistency of object
#   #
#   #   stopifnot(identical(class(object), "fit_wrc_hcc"))
#
#   stop("still to be written")
#
#   invisible(x)
#
# }



# ## ======================================================================
coef.fit_wrc_hcc <- function(
  object, what = c("all", "nonlinear", "linear"),
  subset = NULL, residual_se = FALSE, se = FALSE,
  gof = FALSE, lc = FALSE,
  e0 = FALSE, bound = lc, ...
){

  ## function extracts fitted coefficients from an object of class
  ## fit_wrc_hcc

  ## 2019-11-27 A. Papritz
  ## 2020-01-27 A. Papritz computation of residual standard for
  ##               ml and mpd estimates
  ## 2020-02-03 A. Papritz optional computation of standard errors
  ##               of nonlinear parameters
  ## 2020-02-26 A. Papritz correction error when computing standard errors
  ##                       for nonlinear parameters when some parameters
  ##                       are fixed
  ## 2021-10-13 AP correction of degrees of freedom for ml and mpd method


  #   ## check consistency of object
  #
  #   stopifnot(identical(class(object), "fit_wrc_hcc"))

  ## check arguments

  stopifnot(identical(length(residual_se), 1L) && is.logical(residual_se))
  stopifnot(identical(length(se), 1L) && is.logical(se))
  stopifnot(identical(length(gof), 1L) && is.logical(gof))
  stopifnot(identical(length(lc), 1L) && is.logical(lc))
  stopifnot(identical(length(e0), 1L) && is.logical(e0))
  stopifnot(identical(length(bound), 1L) && is.logical(bound))

  what <- match.arg(what)

  if(residual_se && identical(object[["control"]][["method"]], "wls")){
    stop(
      "computing residual standard is only meaningful for 'mpd' or 'ml' estimates"
    )
    residual_se = FALSE
  }

  if(se && is.null(object[["fit"]][[1]][["hessian"]])){
    se <- FALSE
    warning(
      "standard errors of nonlinear parameters cannot be computed\n",
      "because Hessian matrix is missing"
    )
  }

  if(se && identical(what, "linear")){
    se <- FALSE
    warning(
      "standard errors of linear parameters cannot be computed"
    )
  }

  ## select subsets

  if(!is.null(subset)){
    stopifnot(is.character(subset) || is.numeric(subset) || is.logical(subset))
    stopifnot(length(subset) <= length(object[["fit"]]))
    object[["fit"]] <- object[["fit"]][subset]
  }

  ## extract (fitted) parameters

  result <- lapply(
    object[["fit"]],
    function(x, what, settings){

      ## parameter estimates
      res <- switch(
        what,
        nonlinear = x[["nlp"]],
        linear = x[["lp"]],
        all = c(x[["nlp"]], x[["lp"]]),
        stop("unknown value of argument 'what'")
      )

      ## residual standard errors (Stewart et al., 1992, eqs 8 & 16)
      if(residual_se){
        tmp <- NULL
        if(!is.null(x[["residuals"]][["wrc"]])){
          if(identical(object[["control"]][["method"]], "ml")){
            n.wc <- length(x[["residuals"]][["wrc"]])
          } else {
            n.wc <- length(x[["residuals"]][["wrc"]]) + 2L
          }
          tmp <- c(tmp, se.eps.wc = sqrt(
              attr(x[["objective"]], "ssq_wc") / n.wc
            )
          )
        }
        if(!is.null(x[["residuals"]][["hcc"]])){
          if(identical(object[["control"]][["method"]], "ml")){
            n.hc <- length(x[["residuals"]][["hcc"]])
          } else {
            n.hc <- length(x[["residuals"]][["hcc"]]) + 2L
          }
          tmp <- c(tmp, se.eps.hc = sqrt(
              attr(x[["objective"]], "ssq_hc") / n.hc
            )
          )
        }
        res <- c(res, tmp)
      }

      ## Lc etc
      if(lc){

        ineq.lc <- x[["inequality_constraint"]][["lc"]]
        res <- c(
          res,
          lc = attr(ineq.lc, "lc"),
          lt = attr(ineq.lc, "lt"),
          ratio.lc.lt = attr(ineq.lc, "lc") / attr(ineq.lc, "lt")
        )

        if(e0){
          res <- c(res, attr(ineq.lc, "e0"))
        }

        if(
          bound && length(grep("^c", settings))){
          res <- c(res, attr(ineq.lc, "ratio_lc_lt_bound"))
        }

      }

      ## goodness-of-fit statistics
      if(gof) res <- c(
        res,
        obj = x[["objective"]], convergence = x[["converged"]]
      )

      res

    },
    what = what, settings = object[["control"]][["settings"]]
  )

  result <- as.data.frame(t(sapply(
    result,
    function(x, all.nmes){
      res <- rep(NA_real_, length(all.nmes))
      names(res) <- all.nmes
      sel <- match(names(x), all.nmes)
      res[sel] <- x
      res
    }, all.nmes = names(result[[which.max(sapply(result, length))]])
  )))

  ## add standard errors of estimated nonlinear parameters

  ncol.hessian <- sapply(object[["fit"]], function(x) NCOL(x[["hessian"]]))

  if(se && max(ncol.hessian) > 0L){

    sel <- which.max(ncol.hessian)
    nmes.nlp <- colnames(object[["fit"]][[sel]][["hessian"]])

    result.se <- coef(vcov(object, subset = subset), se = se, status = FALSE)
    result[, paste0("se.", nmes.nlp)] <-
      result.se[, paste0("se.", nmes.nlp), drop = FALSE]

    sel <- match(
      c(matrix(
          c(nmes.nlp, paste0("se.", nmes.nlp)),
          ncol = length(nmes.nlp), byrow = TRUE
        )),
      colnames(result)
    )

    result <- cbind(
      result[, sel, drop = FALSE],
      result[, -sel, drop = FALSE]
    )

  }

  ## return result

  result

}


# ## ======================================================================
plot.fit_wrc_hcc <- function(
  x, what = c("wrc", "hcc"), y = NULL, 
  subset = NULL, ylim_wc = NULL, ylim_hc = NULL,
  head_saturation = 0.01,
  beside = identical(sum(par("mfrow")), 2L),
  pch = 1, col_points = "black",
  col_line_x = "blue", lty_x = "solid",
  col_line_y = "orange", lty_y = "dashed",
  xlab_wc = "head [m]", ylab_wc = "water content [-]",
  xlab_hc = "head [m]", ylab_hc= "hyd. conductivity [m/d]",
  draw_legend = TRUE, draw_parameter = FALSE, cex_legend = 0.7,
  ...

){

  ## function plots measurements and fitted VGM models curves of water
  ## retention and hydraulic conductivity data

  ## 2019-11-27 A. Papritz
  ## 2019-12-01 A. Papritz, changes in legend text annotation
  ## 2020-06-15 A. Papritz, changes in legend text annotation
  ## 2020-08-11 A. Papritz, optional output of parameter values
  ## 2021-05-26 A. Papritz, function specific ellipsis arguments
  ## 2021-06-12 A. Papritz, correction of error in generation of arg.value.x, arg.value.y
  ## 2021-10-20 A. Papritz, optional user-defined ylim_wc and ylim_hc argument
  ## 2021-11-22 A. Papritz, argument head_saturation for zero head values

  ## get names of objects assigned to x and y argument

  cl <- match.call(expand.dots = FALSE)
  nmes <- names(cl)
  sel <- match("x", nmes)
  if(length(cl[[sel]]) > 1){
    arg.value.x <- as.character(deparse(cl[[sel]][[2]]))
  } else {
    arg.value.x <- as.character(deparse(cl[[sel]]))
  }


  arg.value.y <- NULL
  if(!is.null(y)){
    sel <- match("y", nmes)
    if(length(cl[[sel]]) > 1){
      arg.value.y <- as.character(deparse(cl[[sel]][[2]]))
    } else {
      arg.value.y <- as.character(deparse(cl[[sel]]))
    }
  }

  ## select subsets

  if(!is.null(subset)){
    stopifnot(is.character(subset) || is.numeric(subset) || is.logical(subset))
    stopifnot(length(subset) <= length(x[["fit"]]))
    x[["fit"]] <- x[["fit"]][subset]
    if(!is.null(y)) y[["fit"]] <- y[["fit"]][subset]
  }

  ## check consistency of x and y

  #   stopifnot(identical(class(x), "fit_wrc_hcc"))

  if(!is.null(y)){
    #     stopifnot(identical(class(y), "fit_wrc_hcc"))
    stopifnot(identical(length(y[["fit"]]), length(x[["fit"]])))
  }

  ## get control component

  control <- x[["control"]]

  ## drop control component

  x <- x[["fit"]]
  if(!is.null(y)) y <- y[["fit"]]

  ## get wrc_model and hcc_model

  wrc_model <- control[["wrc_model"]]
  hcc_model <- control[["hcc_model"]]

  ## eliminate samples in x without data

  sel <- sapply(
    x,
    function(x) !is.null(x[["model"]])
  )
  x <- x[sel]
  if(!is.null(y)) y <- y[sel]

  ## check availability of wrc and/or hcc data

  wrc.hcc <- sapply(
    x,
    function(xx) c(
      wrc = !is.null(xx[["model"]][["wrc"]]),
      hcc = !is.null(xx[["model"]][["hcc"]])
    )
  )

  wrc.hcc <- apply(wrc.hcc, 1, any)

  ## determine what to plot

  what <- match.arg(what, several.ok = TRUE)

  what <- what[wrc.hcc[what]]
  if(identical(length(what), 2L)) what <- what[match(what, c("wrc", "hcc"))]

  gam_n_newdata <- control[["gam_n_newdata"]]
  precBits <- control[["precBits"]]

  ## match ellipsis arguements

  ellipsis <- list(...)

  ellipsis.plot   <- ellipsis[names(ellipsis) %in% names(formals(plot.default))]
  ellipsis.lines  <- ellipsis[names(ellipsis) %in% names(formals(lines.default))]
  ellipsis.legend <- ellipsis[names(ellipsis) %in% names(formals(legend))]

 ## loop over all samples

  bla <- lapply(
    1:length(x),
    function(i, arg.value.x, arg.value.y, ylim_wc, ylim_hc){

      ## current sample

      id <- names(x[i])
      xx <- x[[i]]
      if(!is.null(y)) yy <- y[[i]]

      ## prepare data for current sample

      ## parameters

      if(!is.null(xx[["nlp"]])){
        nlp.x <- xx[["nlp"]]
        lp.x  <- xx[["lp"]]
      } else {
        return(NULL)
      }

      if(!is.null(y) && !is.null(yy[["nlp"]])){
        nlp.y <- yy[["nlp"]]
        lp.y  <- yy[["lp"]]
      } else {
        nlp.y <- NULL
        lp.y  <- NULL
      }

      ## measurements and modelled values

      ## wrc

      wrc <- "wrc" %in% what && !is.null(xx[["model"]][["wrc"]])
      if(wrc){
        t.terms <- attr(xx[["model"]][["wrc"]], "terms")
        head.wc  <- model.matrix(t.terms, xx[["model"]][["wrc"]])[, 2]
        zero.head <- head.wc <= 0.
        if(any(zero.head)){
          warning("zero head values replaced by 'head_saturation'")
          head.wc[zero.head] <- head_saturation
        }
        wc <- model.response(model.frame(t.terms, xx[["model"]][["wrc"]]))
        head.wc.fit <- exp(seq(
            min(log(head.wc)), max(log(head.wc)),
            length = gam_n_newdata
          ))
        wc.fit.x <- as.numeric(wc_model(
            head.wc.fit, nlp.x, lp.x, precBits, wrc_model
          ))
        if(!is.null(nlp.y)){
          wc.fit.y <- as.numeric(wc_model(
              head.wc.fit, nlp.y, lp.y, precBits, wrc_model
            ))
        } else wc.fit.y <- NULL
      }

      ## hcc

      hcc <- ("hcc" %in% what && !is.null(xx[["model"]][["hcc"]]))
      if(hcc){
        t.terms <- attr(xx[["model"]][["hcc"]], "terms")
        head.hc  <- model.matrix(t.terms, xx[["model"]][["hcc"]])[, 2]
        zero.head <- head.hc <= 0.
        if(any(zero.head)){
          warning("zero head values replaced by 'head_saturation'")
          head.hc[zero.head] <- head_saturation
        }
        hc <- model.response(model.frame(t.terms, xx[["model"]][["hcc"]]))
        head.hc.fit <- exp(seq(
            min(log(head.hc)), max(log(head.hc)),
            length = gam_n_newdata
          ))
        hc.fit.x <- as.numeric(hc_model(
            head.hc.fit, nlp.x, lp.x, precBits, hcc_model
          ))
        if(!is.null(nlp.y)){
          hc.fit.y <- as.numeric(hc_model(
              head.hc.fit, nlp.y, lp.y, precBits, hcc_model
            ))
        } else hc.fit.y <- NULL
      }

      ## title

      t.main <- NULL
      if(length(x) > 1) t.main <- paste("sample", id)

      ## plot measurements and fit

      if(identical(length(what), 2L)){

        if(beside) op <- par(mfrow = c(1, 2))

        if(wrc){

          ## plot wrc data

          if(is.null(ylim_wc)) ylim_wc <- range(wc, wc.fit.x, wc.fit.y) * c(0.96, 1.04)

          do.call(
            plot,
            c(
              list(
                x = head.wc, y = wc, ylim = ylim_wc, pch = pch, col = col_points,
                log = "x", xlab = xlab_wc, ylab = ylab_wc,
                main = t.main
              ),
              ellipsis.plot
            )
          )
          do.call(
            lines,
            c(
              list(x = head.wc.fit, y = wc.fit.x, col = col_line_x, lty = lty_x),
              ellipsis.lines
            )
          )
          if(draw_parameter) do.call(
            legend,
            c(
              list(
                x = "topright", legend = paste(names(c(nlp.x, lp.x)), signif(c(nlp.x, lp.x), digits = 4)),
                horiz = FALSE, text.col = col_line_x, bty = "n", cex = cex_legend
              ),
              ellipsis.legend
            )
          )
          if(!is.null(wc.fit.y)){
            do.call(
              lines,
              c(
                list(
                  x = head.wc.fit, y = wc.fit.y, col = col_line_y, lty = lty_y
                ),
                ellipsis.lines
              )
            )
            if(draw_parameter) do.call(
              legend,
              c(
                list(
                  x = "right", legend = paste(names(c(nlp.y, lp.y)), signif(c(nlp.y, lp.y), digits = 4)),
              horiz = FALSE, text.col = col_line_y, bty = "n", cex = cex_legend
                ),
                ellipsis.legend
              )
            )
            if(draw_legend) do.call(
              legend,
              c(
                list(
                  x = "bottomleft", lty = c(lty_x, lty_y, NA), col = c(col_line_x, col_line_y),
                  legend = c(
                    arg.value.x,
                    paste0(
                      arg.value.y, " (relative SSE y/x: ", signif(
                        attr(yy[["objective"]], "ssq_wc") / attr(xx[["objective"]], "ssq_wc"), 3
                      ), ")"
                    )
                  ), bty = "n", cex = cex_legend
                ),
                ellipsis.legend
              )
            )
          }


        } else {

          ## skip one plot for n_r * 2 layout

          mfg <- par("mfg")
          if((mfg[4]%%2) == 0) plot.new()

        }


        ## plot hcc data

        if(hcc){

          if(is.null(ylim_hc)) ylim_hc <- range(hc, hc.fit.x, hc.fit.y) * c(0.96, 1.04)

          do.call(
            plot,
            c(
              list(x = head.hc, y = hc, ylim = ylim_hc, pch = pch, col = col_points,
                log = "xy", xlab = xlab_hc, ylab = ylab_hc,
                main = t.main
              ),
              ellipsis.plot
            )
          )
          do.call(
            lines,
            c(
              list(x = head.hc.fit, y = hc.fit.x, col = col_line_x, lty = lty_x),
              ellipsis.lines)
          )
          if(draw_parameter) do.call(
            legend,
            c(
              list(
                x = "topright", legend = paste(names(c(nlp.x, lp.x)), signif(c(nlp.x, lp.x), digits = 4)),
                horiz = FALSE, text.col = col_line_x, bty = "n", cex = cex_legend
              ),
              ellipsis.legend
            )
          )
          if(!is.null(hc.fit.y)){
            do.call(
              lines,
              c(
                list(x = head.hc.fit, y = hc.fit.y, col = col_line_y, lty = lty_y),
                ellipsis.lines
              )
            )
            if(draw_parameter) do.call(
              legend,
              c(
                list(
                  x = "right", legend = paste(names(c(nlp.y, lp.y)), signif(c(nlp.y, lp.y), digits = 4)),
                  horiz = FALSE, text.col = col_line_y, bty = "n", cex = cex_legend
                ),
                ellipsis.legend
              )
            )
            if(draw_legend) do.call(
              legend,
              c(
                list(
                  x = "bottomleft", lty = c(lty_x, lty_y), col = c(col_line_x, col_line_y),
                  legend = c(
                    arg.value.x,
                    paste0(
                      arg.value.y, " (relative SSE y/x: ", signif(
                        attr(yy[["objective"]], "ssq_hc") / attr(xx[["objective"]], "ssq_hc"), 3
                      ), ")"
                    )
                  ), bty = "n", cex = cex_legend
                ),
                ellipsis.legend
              )
            )
          }

        } else {

          ## skip one plot for n_r * 2 layout

          mfg <- par("mfg")
          if((mfg[4]%%2) == 0) plot.new()
        }

        if(beside) par(op)

      } else {

        ## plot wrc data

        if(wrc){
          if(is.null(ylim_wc)) ylim_wc <- range(wc, wc.fit.x, wc.fit.y) * c(0.96, 1.04)

          do.call(
            plot,
            c(
              list(x = head.wc, y = wc, ylim = ylim_wc, pch = pch, col = col_points,
                log = "x", xlab = xlab_wc, ylab = ylab_wc,
                main = t.main),
              ellipsis.plot
            )
          )
          do.call(
            lines,
            c(
              list(
                x = head.wc.fit, y = wc.fit.x, col = col_line_x, lty = lty_x
              ),
              ellipsis.lines
            )
          )
          if(draw_parameter) do.call(
            legend,
            c(
              list(x = "topright", legend = paste(names(c(nlp.x, lp.x)), signif(c(nlp.x, lp.x), digits = 4)),
                horiz = FALSE, text.col = col_line_x, bty = "n", cex = cex_legend),
              ellipsis.legend
            )
          )
          if(!is.null(wc.fit.y)){
            do.call(
              lines,
              c(
                list(x = head.wc.fit, y = wc.fit.y, col = col_line_y, lty = lty_y),
                ellipsis.lines
              )
            )
            if(draw_parameter) do.call(
              legend,
              c(
                list(
                  x = "right", legend = paste(names(c(nlp.y, lp.y)), signif(c(nlp.y, lp.y), digits = 4)),
                  horiz = FALSE, text.col = col_line_y, bty = "n", cex = cex_legend
                ),
                ellipsis.legend
              )
            )
            if(draw_legend) do.call(
              legend,
              c(
                list(
                  x = "bottomleft", lty = c(lty_x, lty_y), col = c(col_line_x, col_line_y),
                  legend = c(
                    arg.value.x,
                    paste0(
                      arg.value.y, " (relative SSE y/x: ", signif(
                        attr(yy[["objective"]], "ssq_wc") / attr(xx[["objective"]], "ssq_wc"), 3
                      ), ")"
                    )
                  ), bty = "n", cex = cex_legend
                )
              )

            )
          }
        }

        ## plot hcc data

        if(hcc){
          if(is.null(ylim_hc)) ylim_hc <- range(hc, hc.fit.x, hc.fit.y) * c(0.96, 1.04)

          do.call(
            plot,
            c(
              list(x = head.hc, y = hc, ylim = ylim_hc, pch = pch, col = col_points,
                log = "xy", xlab = xlab_hc, ylab = ylab_hc,
                main = t.main),
              ellipsis.plot
            )
          )
          do.call(
            lines,
            c(
              list(
                x = head.hc.fit, hc.fit.x, col = col_line_x, lty = lty_x
              ),
              ellipsis.lines
            )
          )
          if(draw_parameter) do.call(
            legend,
            c(
              list(
                x = "topright", legend = paste(names(c(nlp.x, lp.x)), signif(c(nlp.x, lp.x), digits = 4)),
                horiz = FALSE, text.col = col_line_x, bty = "n", cex = cex_legend
              ),
              ellipsis.legend
            )
          )
          if(!is.null(hc.fit.y)){
            do.call(
              lines,
              c(
                list(
                  x = head.hc.fit, hc.fit.y, col = col_line_y, lty = lty_y
                ),
                ellipsis.lines
              )
            )
            if(draw_parameter) do.call(
              legend,
              c(
                list(
                  x = "right", legend = paste(names(c(nlp.y, lp.y)), signif(c(nlp.y, lp.y), digits = 4)),
                  horiz = FALSE, text.col = col_line_y, bty = "n", cex = cex_legend
                ),
                ellipsis.legend
              )
            )
            if(draw_parameter) do.call(
              legend,
              c(
                list(
                  x = "bottomleft", lty = c(lty_x, lty_y), col = c(col_line_x, col_line_y),
                  legend = c(
                    arg.value.x,
                    paste0(
                      arg.value.y, " (relative SSE y/x: ", signif(
                        attr(yy[["objective"]], "ssq_hc") / attr(xx[["objective"]], "ssq_hc"), 3
                      ), ")"
                    )
                  ), bty = "n", cex = cex_legend
                ),
                ellipsis.legend
              )
            )
          }
        }

      }

    }, arg.value.x = arg.value.x, arg.value.y = arg.value.y, ylim_wc = ylim_wc, ylim_hc = ylim_hc

  )

  invisible()

}


# ## ======================================================================

lines.fit_wrc_hcc <- function(
  x, what = c("wrc", "hcc"), id = 1, head_saturation = 0.01, ...
){

  ## function adds a modelled water retention curve or hydraulic
  ## conductivity function to a respective plot

  ## 2019-11-27 A. Papritz
  ## 2021-11-22 A. Papritz, argument head_saturation for zero head values

  ## get wrc_model and hcc_model

  wrc_model <- x[["control"]][["wrc_model"]]
  hcc_model <- x[["control"]][["hcc_model"]]

  ## drop control component

  x <- x[["fit"]]

  ## eliminate samples in x without data

  sel <- sapply(
    x,
    function(x) !is.null(x[["model"]])
  )
  x <- x[sel]

  ## determine what to plot

  what <- match.arg(what)

  control <- control_fit_wrc_hcc()
  gam_n_newdata <- control[["gam_n_newdata"]]
  precBits <- control[["precBits"]]

  ## current sample

  if(is.character(id)){
    if(!id %in% names(x)) stop("sample id not found in 'x'")
  } else if(is.numeric(id)){
    if(id < 1 || id > length(x)) stop("sample id not found in 'x'")
  } else stop("'id' must be of mode numeric or character")

  xx <- x[[id]]

  ## prepare data for current sample

  ## parameters

  if(!is.null(xx[["nlp"]])){
    nlp.x <- xx[["nlp"]]
    lp.x  <- xx[["lp"]]
  } else {
    return(NULL)
  }

  ## wrc

  wrc <- "wrc" %in% what && !is.null(xx[["model"]][["wrc"]])
  if(wrc){
    t.terms <- attr(xx[["model"]][["wrc"]], "terms")
    head.wc  <- model.matrix(t.terms, xx[["model"]][["wrc"]])[, 2]
    zero.head <- head.wc <= 0.
    if(any(zero.head)){
      warning("zero head values replaced by 'head_saturation'")
      head.wc[zero.head] <- head_saturation
    }
    wc <- model.response(model.frame(t.terms, xx[["model"]][["wrc"]]))
    head.wc.fit <- exp(seq(
        min(log(head.wc)), max(log(head.wc)),
        length = gam_n_newdata
      ))
    wc.fit.x <- as.numeric(wc_model(
        head.wc.fit, nlp.x, lp.x, precBits, wrc_model
      ))
    lines(head.wc.fit, wc.fit.x, ...)
  }

  ## hcc

  hcc <- ("hcc" %in% what && !is.null(xx[["model"]][["hcc"]]))
  if(hcc){
    t.terms <- attr(xx[["model"]][["hcc"]], "terms")
    head.hc  <- model.matrix(t.terms, xx[["model"]][["hcc"]])[, 2]
    zero.head <- head.hc <= 0.
    if(any(zero.head)){
      warning("zero head values replaced by 'head_saturation'")
      head.hc[zero.head] <- head_saturation
    }
    hc <- model.response(model.frame(t.terms, xx[["model"]][["hcc"]]))
    head.hc.fit <- exp(seq(
        min(log(head.hc)), max(log(head.hc)),
        length = gam_n_newdata
      ))
    hc.fit.x <- as.numeric(hc_model(
        head.hc.fit, nlp.x, lp.x, precBits, hcc_model
      ))
    lines(head.hc.fit, hc.fit.x, ...)
  }

  invisible()

}
