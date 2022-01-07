% 2021-12-22 A. Papritz
% cd ~/R.user.home/soilhypfit/soilhypfit_01-5/pkg/man
% R CMD Rdconv -t html -o bla.html soilhypfit_package.Rd;  open bla.html; R CMD Rd2pdf --force soilhypfit_package.Rd;
% hunspell -d en_GB-ize soilhypfit_package.Rd

\encoding{utf8}

\name{soilhypfit-package}
\alias{soilhypfit-package}
\alias{soilhypfit}
\alias{soilhypfitIntro}

\concept{soil water}
\concept{soil water characteristic}
\concept{soil water retention curve}
\concept{hydraulic conductivity function}
\concept{saturated hydraulic conductivity}

\docType{package}

% \newcommand{\mbvec}{\eqn{\mbox{\boldmath$#1$\unboldmath}^\mathrm{#2} #4}{#1#3#4}}

\title{The soilhypfit Package}

% description

\description{This is a summary of the features and functionality of
  \pkg{soilhypfit}, a package in \R for parameteric modelling of soil water
  retention and hydraulic conductivity data.}

% details

\details{

  % estimation approach

  \subsection{Estimation approach}{

    The main function, \code{\link{fit_wrc_hcc}}, estimates parameters of
    models for soil water retention and hydraulic conductivity by
    \emph{maximum likelihood} (ml, default), \emph{maximum posterior
    density} (mpd) \emph{estimation} (\cite{Stewart et al., 1992}) or
    nonlinear \emph{weighted least squares} (wls) from data on volumetric
    \emph{soil water content}, \mbvec{\theta}{T}{^T}{=(\theta_1, \theta_2,
    ..., \theta_{n_\theta})}, and/or \emph{soil hydraulic conductivity},
    \mbvec{K}{T}{^T}{=(K_1, K_2, ..., K_{n_K})}, both measured at given
    capillary pressure head, \mbvec{h}{T}{^T}{=(h_1, h_2, ...)}.

    For mpd and ml estimation, the models for the measurements are
    \deqn{
      \theta_i = \theta(h_i;\mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath})
      + \varepsilon_{\theta,i}, \ \ i=1, 2, ..., n_\theta,
    }{
      \theta_i = \theta(h_i; \mu, \nu) + \epsilon_{\theta,i}, i = 1, 2, ..., n_\theta,
    }
    \deqn{
      \log(K_j) = \log(K(h_j;\mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath}))
      + \varepsilon_{K,j}, \ \ j=1, 2, ..., n_K,
    }{
      log(K_j) = log(K(h_j; \mu, \nu)) + \epsilon_{K,j}, j = 1, 2, ..., n_K,
    }
    where
    \eqn{
      \theta(h_i; \mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath})
    }{\theta(h_i; \mu, \nu)} and
    \eqn{
      K(h_j; \mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath})
    }{K(h_j; \mu, \nu)}
    denote modelled water content and conductivity,
    \mbvec{\mu}{}{}{}{} and \mbvec{\nu}{}{}{}{} are the
    \emph{conditionally linear} and \emph{nonlinear}
    model parameters (see below and \cite{Bates and Watts, 1988, sec. 3.3.5}),
    and
    \eqn{\varepsilon_{\theta,i}}{\epsilon_{\theta,i}} and
    \eqn{\varepsilon_{K,j}}{\epsilon_{K,j}} are independent,
    normally distributed errors
    with zero means and variances \eqn{\sigma^2_\theta}
    and \eqn{\sigma^2_K}, respectively.


    % objective function

    Let
    \deqn{
      Q_{\theta}(\mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath};
      \mbox{\boldmath$\theta$\unboldmath},
      \mbox{\boldmath$h$\unboldmath}) =\left(
        \mbox{\boldmath$\theta$\unboldmath} -
        \mbox{\boldmath$\theta$\unboldmath}(
          \mbox{\boldmath$h$\unboldmath};
          \mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath}
        )
      \right)^\mathrm{T}
      \mbox{\boldmath$W$\unboldmath}_\theta
      \left(
        \mbox{\boldmath$\theta$\unboldmath} -
        \mbox{\boldmath$\theta$\unboldmath}(
          \mbox{\boldmath$h$\unboldmath};
          \mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath}
        )
      \right),
    }{
      Q_\theta(\mu, \nu; \theta, h)
      = (\theta - \theta(h; \mu, \nu))^T W_\theta (\theta - \theta(h; \mu, \nu)),
    }
    and
    \deqn{
      Q_{K}(\mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath};
       \mbox{\boldmath$K$\unboldmath}, \mbox{\boldmath$h$\unboldmath}) = \left(
        \log(\mbox{\boldmath$K$\unboldmath}) -
        \log(\mbox{\boldmath$K$\unboldmath}(
          \mbox{\boldmath$h$\unboldmath};
          \mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath}
        ))
      \right)^\mathrm{T}
      \mbox{\boldmath$W$\unboldmath}_K
      \left(
        \log(\mbox{\boldmath$K$\unboldmath}) -
        \log(\mbox{\boldmath$K$\unboldmath}(
          \mbox{\boldmath$h$\unboldmath};
          \mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath}
        ))
      \right),
    }{
      Q_K(\mu, \nu; K, h) = (log(K) - log(K(h; \mu, \nu)))^T W_K (log(K) - log(K(h; \mu, \nu))),
    }
    denote the (possibly weighted) residual sums of squares between
    measurements and modelled values.
    \eqn{\mbox{\boldmath$\theta$\unboldmath}(
        \mbox{\boldmath$h$\unboldmath};
        \mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath}
      ) )
    }{\theta(h; \mu, \nu))} and
    \eqn{\mbox{\boldmath$K$\unboldmath}(
        \mbox{\boldmath$h$\unboldmath};
        \mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath}
      )
    }{K(h; \mu, \nu)}
    are vectors with modelled values of water content and hydraulic
    conductivity, and
    \eqn{\mbox{\boldmath$W$\unboldmath}_\theta}{W_\theta} and
    \eqn{\mbox{\boldmath$W$\unboldmath}_K}{W_K} are optional diagonal matrices
    with weights \eqn{w_{\theta,i}} and \eqn{w_{K,j}}, respectively.
    The weights are products of \emph{case weights} \eqn{w^\prime_{\theta,i}}
    and \eqn{w^\prime_{K,j}} and \emph{variable weights} \eqn{w_\theta}, \eqn{w_K}, hence
    \eqn{
      w_{\theta,i} = w_\theta\, w^\prime_{\theta,i}
    }{w_{\theta,i} = w_\theta * w^\prime_{\theta,i}} and
    \eqn{
      w_{K,j} = w_K\, w^\prime_{K,j}
    }{w_{K,j} = w_K * w^\prime_{K,j}}.

    The objective function for ml and mpd estimation
    is equal to (\cite{Stewart et al., 1992},
    eqs 14 and 15, respectively)
    \deqn{
      Q(\mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath};
        \mbox{\boldmath$\theta$\unboldmath}, \mbox{\boldmath$K$\unboldmath},
        \mbox{\boldmath$h$\unboldmath}
        ) = \frac{\kappa_\theta}{2}
        \log( Q_{\theta}(\mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath};
        \mbox{\boldmath$\theta$\unboldmath}, \mbox{\boldmath$h$\unboldmath})) +
        \frac{\kappa_K}{2}
        \log( Q_{K}(\mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath};
        \mbox{\boldmath$K$\unboldmath}, \mbox{\boldmath$h$\unboldmath})),
    }{
       Q(\mu, \nu; \theta, K, h) = \kappa_\theta / 2 Q_\theta(\mu, \nu; \theta, h) +
         \kappa_K / 2 Q_K(\mu, \nu; K, h),
    }
    where \eqn{\kappa_v = n_v + 2} for mpd and \eqn{\kappa_v = n_v} for ml
    estimation, \eqn{v \in (\theta, K)},
    and weights \eqn{w_{\theta,i} = w_{K,j} = 1}.
    Note that
    \eqn{
      Q(\mbox{\boldmath$\mu$\unboldmath},
      \mbox{\boldmath$\nu$\unboldmath}; \mbox{\boldmath$\theta$\unboldmath},
      \mbox{\boldmath$K$\unboldmath}, \mbox{\boldmath$h$\unboldmath})
    }{
      Q(\mu, \nu; \theta, K, h)
    } is equal to the negative logarithm of the
    concentrated loglikelihood or the concentrated posterior density,
    obtained by replacing
    \eqn{\sigma^2_\theta} and \eqn{\sigma^2_K} by their conditional maximum
    likelihood or maximum density estimates
    \eqn{
      \widehat{\sigma}^2_\theta(\mbox{\boldmath$\mu$\unboldmath}),
      \mbox{\boldmath$\nu$\unboldmath})
    }{hatsigma^2_\theta(\mu, \nu)} and
    \eqn{
      \widehat{\sigma}^2_K(\mbox{\boldmath$\mu$\unboldmath})
    }{hat\sigma^2_K(\mu, \nu)}
    respectively (\cite{Stewart et al., 1992}, p. 642).

    For wls the objective function is equal to
    \deqn{
      Q(\mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath};
        \mbox{\boldmath$\theta$\unboldmath}, \mbox{\boldmath$K$\unboldmath},
        \mbox{\boldmath$h$\unboldmath}
        ) = Q_{\theta}(\mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath};
        \mbox{\boldmath$\theta$\unboldmath}, \mbox{\boldmath$h$\unboldmath}) +
        Q_{K}(\mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath};
        \mbox{\boldmath$K$\unboldmath}, \mbox{\boldmath$h$\unboldmath}).
    }{
       Q(\mu, \nu; \theta, K, h) = Q_\theta(\mu, \nu; \theta, h) +
         Q_K(\mu, \nu; K, h).
    }
    If either water content or conductivity data are not available, then the respective
    terms are omitted from
    \eqn{Q(\mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath};
        \mbox{\boldmath$\theta$\unboldmath}, \mbox{\boldmath$K$\unboldmath},
        \mbox{\boldmath$h$\unboldmath}
        )
    }{Q(\mu, \nu; \theta, K, h)}.


    % estimating linear parameters

    The function \code{\link{fit_wrc_hcc}} does not attempt to estimate the parameters
    by minimising \cr
    \eqn{
      Q(\mbox{\boldmath$\mu$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath};
      \mbox{\boldmath$\theta$\unboldmath}, \mbox{\boldmath$K$\unboldmath},
      \mbox{\boldmath$h$\unboldmath})
    }{Q(\mu, \nu; \theta, K, h)} directly with respect to
    \mbvec{\mu}{}{}{}{} and \mbvec{\nu}{}{}{}{}.
    Rather, it exploits the fact that for \emph{given nonlinear parameters}
    \mbvec{\nu}{}{}{}{}, the \emph{conditionally linear parameters}
    \eqn{
      \mbox{\boldmath$\mu$\unboldmath}^T = (\theta_r, \theta_s, \log(K_0))
    }{
      \mu^T = (\theta_r, \theta_s, log(K_0))
    }
    can be estimated straightforwardly
    by minimising the \emph{conditional residual sums of squares}

    \deqn{
      Q_\theta^*(\theta_r, \theta_s;
      \mbox{\boldmath$\theta$\unboldmath},
      \mbox{\boldmath$h$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath})
      =\left(
        \mbox{\boldmath$\theta$\unboldmath} - [
          \mbox{\boldmath$1$\unboldmath},
          \mbox{\boldmath$S$\unboldmath}(
            \mbox{\boldmath$h$\unboldmath};
            \mbox{\boldmath$\nu$\unboldmath}
          )] \, \left[
            \begin{array}{c} \theta_r \\ \theta_s - \theta_r \end{array}
          \right]
      \right)^\mathrm{T}
      \mbox{\boldmath$W$\unboldmath}_\theta
      \left(
        \mbox{\boldmath$\theta$\unboldmath} - [
          \mbox{\boldmath$1$\unboldmath},
          \mbox{\boldmath$S$\unboldmath}(
            \mbox{\boldmath$h$\unboldmath};
            \mbox{\boldmath$\nu$\unboldmath}
          )] \, \left[
            \begin{array}{c} \theta_r \\ \theta_s - \theta_r \end{array}
          \right]
      \right)
    }{
      Q*_\theta(\theta_r, \theta_s; \theta,  h, \nu) =
      (\theta - [1, S(h;\nu)] [\theta_r, \theta_s - \theta_r]^T])^T W_\theta
      (\theta - [1, S(h;\nu)] [\theta_r, \theta_s - \theta_r]^T])
    }
    with respect to \eqn{\theta_r} and  \eqn{\theta_s-\theta_r} and/or
    \deqn{
      Q_K^*(K_0; \mbox{\boldmath$K$\unboldmath},
        \mbox{\boldmath$h$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath})=
      \left(
        \log(\mbox{\boldmath$K$\unboldmath}) - \log(K_0  \,
          \mbox{\boldmath$k$\unboldmath}(
            \mbox{\boldmath$h$\unboldmath};
            \mbox{\boldmath$\nu$\unboldmath}
          ))
      \right)^\mathrm{T}
      \mbox{\boldmath$W$\unboldmath}_K
      \left(
        \log(\mbox{\boldmath$K$\unboldmath}) - \log(K_0  \,
          \mbox{\boldmath$k$\unboldmath}(
            \mbox{\boldmath$h$\unboldmath};
            \mbox{\boldmath$\nu$\unboldmath}
          ))
      \right),
    }{
      Q*_K(K_0; K,  h, \nu) =
      (log(K) - log(K_0  k(h; \nu)))^T W_K
      (log(K) - log(K_0  k(h; \nu))),
    }
    with respect to \eqn{\log(K_0)}{log(K_0)}, where
    \mbvec{1}{}{}{}{} is a vector of ones,
    \eqn{\mbox{\boldmath$S$\unboldmath}(
            \mbox{\boldmath$h$\unboldmath};
            \mbox{\boldmath$\nu$\unboldmath}
          )^\mathrm{T} = (
            S(h_1; \mbox{\boldmath$\nu$\unboldmath}),
            ...,
            S(h_{n_\theta}; \mbox{\boldmath$\nu$\unboldmath})
          )
    }{
      S(h;\nu)^T = (S(h_1;\nu), S(h_2;\nu), ..., S(h_{n_\theta};\nu))
    } and
    \eqn{\mbox{\boldmath$k$\unboldmath}(
            \mbox{\boldmath$h$\unboldmath};
            \mbox{\boldmath$\nu$\unboldmath}
          )^\mathrm{T} = (
            k(h_1; \mbox{\boldmath$\nu$\unboldmath}),
            ...,
            k(h_{n_K}; \mbox{\boldmath$\nu$\unboldmath})
          )
    }{
      k(h;\nu)^T = (k(h_1;\nu), k(h_2;\nu), ..., k(h_{n_K};\nu))
    }
    are vectors of modelled \emph{water saturation}
    and modelled \emph{relative conductivity} values,
    \eqn{\theta_r} and \eqn{\theta_s} are the \emph{residual}
    and \emph{saturated water content}, and
    \eqn{K_0} is the \emph{saturated hydraulic conductivity}.

    Unconstrained conditional estimates, say \eqn{
    \widehat{\theta}_r(\mbox{\boldmath$\nu$\unboldmath})
    }{hat\theta_r(\nu)}, \eqn{
    \widehat{\theta}_s(\mbox{\boldmath$\nu$\unboldmath}) -
    \widehat{\theta}_r(\mbox{\boldmath$\nu$\unboldmath}) }{hat\theta_s(\nu)
    - hat\theta_r(\nu)} and \eqn{
    \widehat{\log(K_0)}(\mbox{\boldmath$\nu$\unboldmath}) }{hatlog(K_0)(\nu)} can be
    easily obtained from the normal equations of the respective (weighted)
    least squares problems, and quadratic programming yields conditional (weighted)
    least squares estimates that honour the inequality constraints \eqn{0
    \le \theta_r \le \theta_s \le 1}.


    % estimating nonlinear parameters

    Let \eqn{
     \widehat{\mbox{\boldmath$\mu$\unboldmath}}(\mbox{\boldmath$\nu$\unboldmath})^\mathrm{T} =
     (
     \widehat{\theta}_r(\mbox{\boldmath$\nu$\unboldmath}),
     \widehat{\theta}_s(\mbox{\boldmath$\nu$\unboldmath}),
     \widehat{\log(K_0)}(\mbox{\boldmath$\nu$\unboldmath})
     )
    }{
      hat\mu(\nu)^T = (hat\theta_r(\nu), hat\theta_s(\nu), hatlog(K_0)(\nu))
    }
    be the conditional estimates of the linear parameters
    obtained by minimising
    \eqn{
      Q_\theta^*(\theta_r, \theta_s;
      \mbox{\boldmath$\theta$\unboldmath},
      \mbox{\boldmath$h$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath})
    }{Q*_\theta(\theta_r, \theta_s; \theta, h, \nu)},
    and
    \eqn{
      Q_K^*(K_0;
      \mbox{\boldmath$K$\unboldmath},
      \mbox{\boldmath$h$\unboldmath}, \mbox{\boldmath$\nu$\unboldmath})
    }{Q*_K(K_0; K, h, \nu)}, respectively.
    \code{\link{fit_wrc_hcc}} then estimates the
    nonlinear parameters by minimising
    \eqn{
      Q(
      \widehat{\mbox{\boldmath$\mu$\unboldmath}}(\mbox{\boldmath$\nu$\unboldmath}),
      \mbox{\boldmath$\nu$\unboldmath};
      \mbox{\boldmath$\theta$\unboldmath}, \mbox{\boldmath$K$\unboldmath},
      \mbox{\boldmath$h$\unboldmath})
    }{Q(hat\mu(\nu), \nu; \theta, K, h)} with respect to \mbvec{\nu}{}{}{}{}
    by a nonlinear optimisation algorithm.

    % estimating the error variances

    For mpd and ml estimation the variances of the model errors are estimated by
    (\cite{Stewart et al., 1992}, eq. 16)
    \deqn{
      \widehat{\sigma}_\theta^2 = \frac{Q_{\theta}(
        \widehat{\mbox{\boldmath$\mu$\unboldmath}}(\widehat{\mbox{\boldmath$\nu$\unboldmath}}),
        \widehat{\mbox{\boldmath$\nu$\unboldmath}};
      \mbox{\boldmath$\theta$\unboldmath},
      \mbox{\boldmath$h$\unboldmath})}{\kappa_\theta},
    }{
      hat\sigma^2_\theta = Q_\theta(hat\mu(hat\nu), hat\nu; \theta, h) / \kappa_\theta,
    }
    and
    \deqn{
      \widehat{\sigma}_K^2 = \frac{Q_{K}(
        \widehat{\mbox{\boldmath$\mu$\unboldmath}}(\widehat{\mbox{\boldmath$\nu$\unboldmath}}),
        \widehat{\mbox{\boldmath$\nu$\unboldmath}};
      \mbox{\boldmath$K$\unboldmath},
      \mbox{\boldmath$h$\unboldmath})}{\kappa_K}.
    }{
      hat\sigma^2_K= Q_K(hat\mu(hat\nu), hat\nu; K, h) / \kappa_K.
    }

    % covariances of nonlinear parameters

    Furthermore, for mpd and ml estimation, the covariance matrix of the estimated
    nonlinear parameters may be approximated by the inverse Hessian matrix
    of \eqn{
      Q(\widehat{\mbox{\boldmath$\mu$\unboldmath}}(\mbox{\boldmath$\nu$\unboldmath}),
      \mbox{\boldmath$\nu$\unboldmath};
      \mbox{\boldmath$\theta$\unboldmath}, \mbox{\boldmath$K$\unboldmath},
      \mbox{\boldmath$h$\unboldmath})
    }{Q(hat\mu(\nu), \nu; \theta, K, h)} at the solution
    \eqn{
      \widehat{\mbox{\boldmath$\nu$\unboldmath}}
    }{hat\nu}
    (\cite{Stewart and \enc{Sørensen}{Soerensen}, 1981}), i.e.

    \deqn{
      \mathrm{Cov}[
        \widehat{\mbox{\boldmath$\nu$\unboldmath}},
        \widehat{\mbox{\boldmath$\nu$\unboldmath}}^T] \approx
        \mbox{\boldmath$A$\unboldmath}^{-1},
    }{
      Cov[hat\nu, hat\nu^T] ≈ A^{-1},
    }
    where
    \deqn{
      [\mbox{\boldmath$A$\unboldmath}]_{kl} =
      \frac{\partial^2}{\partial\nu_k\, \partial\nu_l} \left.
      Q(\widehat{\mbox{\boldmath$\mu$\unboldmath}}(\mbox{\boldmath$\nu$\unboldmath}),
      \mbox{\boldmath$\nu$\unboldmath};
      \mbox{\boldmath$\theta$\unboldmath}, \mbox{\boldmath$K$\unboldmath},
      \mbox{\boldmath$h$\unboldmath})\right|_{\nu=\hat{\nu}}.
    }{
      [A]_{k,l} = \partial^2 / (\partial\nu_k, \partial\nu_l)
      Q(hat\mu(\nu), \nu; \theta, K, h)|_(\nu = hat\nu).
    }

  }

  % estimation details

  \subsection{Details on parameter estimation}{

    % models

    \subsection{Models for water retention curves and hydraulic conductivity functions}{

      Currently, \code{\link{fit_wrc_hcc}} allows to estimate the parameters of the
      simplified form of the \emph{Van Genuchten-Mualem} (VGM) model (Van
      Genuchten, 1980) with the restriction \eqn{m = 1 - \frac{1}{n}}{m = 1
      - 1/n}, see \code{\link{wc_model}} and \code{\link{hc_model}}.
      This model has the following parameters:

      \itemize{
        \item \mbvec{\mu}{T}{^T}{= (\theta_r, \theta_s, K_0)} (see above) and

        \item \mbvec{\nu}{T}{^T}{= (\alpha, n, \tau)} where \eqn{\alpha} is
        the inverse air entry pressure, \eqn{n} the shape and \eqn{\tau}
        the tortuosity parameter.

      }

      Any of these parameters can be either estimated from data or kept
      fixed at the specified initial values (see arguments \code{param} and
      \code{fit_param} of \code{\link{fit_wrc_hcc}}).


    }

    % constraining parameter estimates

    \subsection{Imposing physical constraints on the estimated parameters}{

      Parameters of models for the water retention curve and the hydraulic
      conductivity function may vary only within certain bounds (see
      \code{\link{wc_model}}, \code{\link{hc_model}} and
      \code{\link{param_boundf}} for allowed ranges).
      \code{\link{fit_wrc_hcc}} either estimates \emph{transformed
      parameters} that vary over the whole real line and can therefore be
      estimated without constraints (see \code{\link{param_transf}}), or it
      uses algorithms (quadratic programming for estimating
      \mbvec{\mu}{}{}{}, nonlinear optimisation algorithms with box
      constraints for estimating \mbvec{\nu}{}{}{}) that restrict estimates
      to permissible ranges, see \emph{Details} section of
      \code{\link{control_fit_wrc_hcc}}.

      In addition, for natural soils, the parameters of the VGM model
      cannot vary independently from each other over the allowed ranges.
      Sets of fitted parameters should always result in soil hydraulic
      quantities that are physically meaningful.  One of these quantities
      is the \emph{characteristic length} \eqn{L_\mathrm{c}}{L_c} of
      \emph{stage-I} evaporation from a soil (Lehmann et al., 2008).
      \eqn{L_\mathrm{c}}{L_c} can be related to the parameters of the VGM
      model, see Lehmann et al.  (2008, 2020) and
      \code{\link{evaporative-length}}.

      Using several soil hydrological databases, Lehmann et al.  (2020)
      analysed the mutual dependence of VGM parameters and proposed
      regression equations to relate the inverse air entry pressure
      \eqn{\alpha} and the saturated hydraulic \eqn{K_0} to the shape
      parameter \eqn{n}, which characterises the width of the pore size
      distribution of a soil.  Using these relations, Lehmann et al.
      (2020) then computed the expected value (``target'')
      \eqn{L_\mathrm{t}}{L_t} of \eqn{L_\mathrm{c}}{L_c} for given \eqn{n}
      and tortuosity parameter \eqn{\tau}, see
      \code{\link{evaporative-length}}.  \code{\link{fit_wrc_hcc}} allows
      to constrain estimates of the nonlinear parameters \mbvec{\nu}{}{}{}
      by defining permissible lower and upper bounds for the ratio
      \eqn{L_\mathrm{c}/L_\mathrm{t}}{L_c/L_t}, see arguments
      \code{ratio_lc_lt_bound} of \code{\link{fit_wrc_hcc}} and
      \code{settings} of \code{\link{control_fit_wrc_hcc}}.

    }

    % choice of optimisation algorithm

    \subsection{Choice of optimisation algorithm for estimating the nonlinear parameters}{

      To estimate \mbvec{\nu}{}{}{}, \code{\link{fit_wrc_hcc}} minimises
      \eqn{ Q(
      \widehat{\mbox{\boldmath$\mu$\unboldmath}}(\mbox{\boldmath$\nu$\unboldmath}),
      \mbox{\boldmath$\nu$\unboldmath};
      \mbox{\boldmath$\theta$\unboldmath}, \mbox{\boldmath$K$\unboldmath},
      \mbox{\boldmath$h$\unboldmath}) }{Q(hat\mu(\nu), \nu; \theta, K, h)}
      either by a nonlinear optimisation algorithm available in the library
      \emph{NLopt} (Johnson, see \code{\link[nloptr]{nloptr}}) or by the Shuffled
      Complex Evolution (SCE) optimisation algorithm (Duan et al., 1994,
      see \code{\link[SoilHyP]{SCEoptim}}).  The choice of the algorithm
      is controlled by the argument \code{settings} of the function
      \code{\link{control_fit_wrc_hcc}}:

      \enumerate{

        \item global optimisation without constraints for the
        ratio \eqn{L_\mathrm{c}/L_\mathrm{t}}{L_c/L_t} \cr(\code{settings =
        "uglobal"} or \code{settings = "sce"}),

        \item global optimisation with inequality constraints
        for the ratio \eqn{L_\mathrm{c}/L_\mathrm{t}}{L_c/L_t}\cr
        (\code{settings = "cglobal"}),

        \item local optimisation without constraints for the
        ratio \eqn{L_\mathrm{c}/L_\mathrm{t}}{L_c/L_t}\cr (\code{settings =
        "ulocal"}),

        \item local optimisation with inequality constraints
        for the ratio \eqn{L_\mathrm{c}/L_\mathrm{t}}{L_c/L_t}\cr
        (\code{settings = "clocal"}).

      }

      The \code{settings} argument also sets reasonable default values
      for the termination (= convergence) criteria for the various
      algorithms, see
      \href{https://nlopt.readthedocs.io/en/latest/NLopt_Introduction/#termination-conditions}{NLopt
      documentation, section Termination conditions}.  The NLopt
      documentation contains a very useful discussion of
      \href{https://nlopt.readthedocs.io/en/latest/NLopt_Introduction/#optimisation-problems}{(constrained)
      optimisation problems} in general,
      \href{https://nlopt.readthedocs.io/en/latest/NLopt_Introduction/#global-versus-local-optimisation}{global
      vs.  local optimisation} and
      \href{https://nlopt.readthedocs.io/en/latest/NLopt_Introduction/#gradient-based-versus-derivative-free-algorithms}{gradient-based
      vs.  derivative-free algorithms}.

      Note that besides the \code{settings} argument of
      \code{\link{control_fit_wrc_hcc}}, the arguments \code{nloptr} and
      \code{sce} along with the functions \code{\link{control_nloptr}} and
      \code{\link{control_sce}} allow to fully control the nonlinear
      optimisation algorithms, see \code{\link{control_fit_wrc_hcc}} for
      details.




    }

    % initial values

    \subsection{Computing initial values of parameters}{

      For local optimisation algorithms ``good'' initial values of
      \mbvec{\nu}{}{}{} are indispensable for successful estimation.
      \code{\link{fit_wrc_hcc}} allows to compute initial values of
      \eqn{\alpha} and \eqn{n} from water retention data by the following
      procedure:

      \enumerate{

         \item Smooth the water retention data,
         \eqn{
           (\theta_i, y_i= \log(h_i)), i=1,2, ... n_\theta,
          }{
           (\theta_i, y_i = log(h_i)), i=1,2, ... n_\theta,
          } by an additive model.

         \item Determine the saturation, \eqn{S^*}{S*}, and the logarithm of
         capillary pressure head, \eqn{y^* = \log(h^*)}{y* = log(h*)}, at
         the inflection point of the additive model fit.

         \item Find the root, say \eqn{\widehat{m}}{hatm}, of
         \eqn{S^* = (1 + 1/m)^{-m}}{S* = (1 + 1/m)^{-m}}.
         One obtains the right-hand side of this equation by solving
         \eqn{\frac{\partial^2}{\partial y^2}
           \left[
             S_\mathrm{VG}(\exp(y); \mbox{\boldmath$\nu$\unboldmath})
           \right] = 0
         }{d^2/dy^2[S_VG(exp(y); \nu)] = 0}
         for \eqn{y} and plugging the result into the expression for
         \eqn{
           S_\mathrm{VG}(\exp(y); \mbox{\boldmath$\nu$\unboldmath}),
         }{S_VG(exp(y); \nu),}
         see \code{\link{wc_model}}.

         \item Compute
         \eqn{
           \widehat{n} = 1 / (1 - \widehat{m})
         }{hatn = 1 /(1-hatm)}
         and
         \eqn{
           \widehat{\alpha} = 1 / \exp(y^*) \,(1/\widehat{m})^{1-\widehat{m}}.
         }{
           hat\alpha = 1 / exp(y*) * (1/hatm)^{1-hatm}.
         }
         The second expression is again a result of solving
         \eqn{\frac{\partial^2}{\partial y^2}
           \left[
             S_\mathrm{VG}(\exp(y); \mbox{\boldmath$\nu$\unboldmath})
           \right] = 0.
         }{d^2/dy^2[S_VG(exp(y); \nu)] = 0.}

      }

      Initial values for local optimisation algorithms can of course also
      be obtained by first estimating the parameters by a global
      algorithm.  These estimates can be ``refined'' in a second
      step by a local unconstrained algorithm, possibly
      followed by a third call of \code{\link{fit_wrc_hcc}} to constrain
      the estimated parameters by the ratio
      \eqn{L_\mathrm{c}/L_\mathrm{t}}{L_c/L_t}.  The method
      \code{\link{coef.fit_wrc_hcc}} can be used to extract the estimated
      parameters from an object of class \code{fit_wrc_hcc} and to pass
      them as initial values to \code{\link{fit_wrc_hcc}}, see
      \code{\link{fit_wrc_hcc}} for examples.

    }

  }

}

% references

\references{
  Bates, D. M., and Watts, D. G. (1988) Nonlinear Regression Analysis and its
  Applications. John Wiley & Sons, New York \doi{10.1002/9780470316757}.

  Duan, Q., Sorooshian, S., and Gupta, V. K. (1994) Optimal use of the
  SCE-UA global optimisation method for calibrating watershed models,
  \emph{Journal of Hydrology} \bold{158}, 265--284,
  \doi{10.1016/0022-1694(94)90057-4}.

  Johnson, S.G. The NLopt nonlinear-optimisation package.
  \url{https://github.com/stevengj/nlopt}.

  Lehmann, P., Assouline, S., Or, D. (2008) Characteristic lengths
  affecting evaporative drying of porous media.  \emph{Physical Review E},
  \bold{77}, 056309, \doi{10.1103/PhysRevE.77.056309}.

  Lehmann, P., Bickel, S., Wei, Z., Or, D. (2020) Physical Constraints for
  Improved Soil Hydraulic Parameter Estimation by Pedotransfer Functions.
  \emph{Water Resources Research} \bold{56}, e2019WR025963,
  \doi{10.1029/2019WR025963}.

  Stewart, W.E., Caracotsios, M. \enc{Sørensen}{Soerensen}, J.P. 1992.  Parameter
  estimation from multiresponse data.  \emph{AIChE Journal}, \bold{38},
  641--650, \doi{10.1002/aic.690380502}.

  Stewart, W.E. and \enc{Sørensen}{Soerensen}, J.P. (1981)
  Bayesian estimation of common
  parameters from multiresponse data with missing observations.
  \emph{Technometrics}, \bold{23}, 131--141, \cr
  \doi{10.1080/00401706.1981.10486255}.

  Van Genuchten, M. Th. (1980) A closed-form equation for predicting the
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
%   \code{\link{soilhypfitIntro}} for a description of the models and a brief
%   summary of the parameter estimation approach;
%
  \code{\link{fit_wrc_hcc}} for (constrained) estimation of parameters of
  models for soil water retention and hydraulic conductivity data;

  \code{\link{control_fit_wrc_hcc}} for options to control
  \code{fit_wrc_hcc};

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

% \examples{
% }
