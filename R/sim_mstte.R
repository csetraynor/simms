#' Simulate multi-state survival data
#' @export
sim_idm <- function(dist01 = c("weibull", "exponential", "gompertz")  ,
                   dist02 = c("weibull", "exponential", "gompertz")  ,
                   dist12 = c("weibull", "exponential", "gompertz")  ,
                   lambdas01, gammas01, betas01,
                   lambdas02, gammas02, betas02,
                   lambdas12, gammas12, betas12,
                   hazard01, hazard02, hazard12,
                   loghazard01, loghazard02, loghazard12,
                   tde, tdefunction = NULL,
                   mixture = FALSE, pmix = 0.5, hazard, loghazard,
                   cumhazard, logcumhazard,
                   x, cens , idvar = NULL, ids = NULL, nodes = 15,
                   interval = c(1E-8, 500),
                   seed = sample.int(.Machine$integer.max, 1), ...){

  if (missing(lambdas01))
    lambdas01 <- NULL
  if (missing(lambdas02))
    lambdas02 <- NULL
  if (missing(lambdas12)){
    lambdas12 <- NULL
  }

  if (missing(tde)){
    tde <- NULL
  }

  if (missing(gammas01))
    gammas01 <- NULL
  if (missing(gammas02))
    gammas02 <- NULL
  if (missing(gammas12))
    gammas12 <- NULL

  if (missing(betas01))
    betas01 <- NULL
  if (missing(betas02))
    betas02 <- NULL
  if (missing(betas12))
    betas12 <- NULL

  if (missing(hazard01))
    hazard01 <- NULL
  if (missing(hazard02))
    hazard02 <- NULL
  if (missing(hazard12))
    hazard12 <- NULL

  if (missing(loghazard01))
    loghazard01 <- NULL
  if (missing(loghazard02))
    loghazard02 <- NULL
  if (missing(loghazard12))
    loghazard12 <- NULL

  if (!is.null(ids) == is.null(idvar))
    stop("Both 'idvar' and 'ids' must be supplied together.")
  if (!is.null(ids)) {
    N <- length(ids) # number of individuals
    if (any(duplicated(ids)))
      stop("The 'ids' vector must specify unique ID values.")
  } else {
    N <- nrow(x) # number of individuals
    ids <- seq(N)
  }

  if(is.null(hazard01) && is.null(loghazard01)){
    inverted_surv01 <- get_inverted_surv(dist = dist01,
                                         lambdas = lambdas01,
                                         gammas = gammas01)

    R <- sapply(ids, function(i) {
      x_i <- subset_df(x, i, idvar = idvar)
      betas_i <- subset_df(betas01, i, idvar = idvar)
      u_i <- stats::runif(1)
      t_i <- inverted_surv01(u = u_i, x = x_i, betas = betas_i)
      return(t_i)
    })

    inverted_surv02 <- get_inverted_surv(dist = dist02,
                                         lambdas = lambdas02, gammas = gammas02)

    D <- sapply(ids, function(i) {
      x_i <- subset_df(x, i, idvar = idvar)
      betas_i <- subset_df(betas02, i, idvar = idvar)
      u_i <- stats::runif(1)
      t_i <- inverted_surv02(u = u_i, x = x_i, betas = betas_i)
      return(t_i)
    })

    yesR <- R < D

    inverted_surv12 <- get_inverted_surv(dist = dist12,
                                         lambdas = lambdas12, gammas = gammas12)

    tt12 <- sapply(seq_along(ids[yesR]), function(i) {
      x_i <- subset_df(x[yesR, ], i, idvar = idvar)
      betas_i <- subset_df(betas12, i, idvar = idvar)
      u_i <- stats::runif(1)
      t_i <- inverted_surv12(u = u_i, x = x_i, betas = betas_i)
      return(t_i)
    })
  } else { # user-defined hazard or log hazard
    if (!is.null(tde))
      stop("'tde' cannot be specified with a user-defined [log] hazard ",
           "function; please just incorporate the time dependent effects ",
           "into the [log] hazard function you are defining.")
    if (!is.null(loghazard01)) {
      hazard01 <- function(t, x, betas, ...) exp(loghazard01(t, x, betas, ...))
    }
    if (!is.null(loghazard02)) {
      hazard02 <- function(t, x, betas, ...) exp(loghazard02(t, x, betas, ...))
    }
    if (!is.null(loghazard12)) {
      hazard12 <- function(t, x, betas, ...) exp(loghazard12(t, x, betas, ...))
    }
    qq <- get_quadpoints(nodes)
    R <- sapply(ids, function(i) {
      x_i <- subset_df(x, i, idvar = idvar)
      betas_i <- subset_df(betas01, i, idvar = idvar)
      log_u_i <- log(stats::runif(1))
      # check whether S(t) is still greater than random uniform variable u_i at the
      # upper limit of uniroot's interval (otherwise uniroot will return an error)
      at_limit <- rootfn_hazard(interval[2], hazard = hazard01, x = x_i,
                                betas = betas_i, log_u = log_u_i, qq = qq)
      if (is.nan(at_limit)) {
        STOP_nan_at_limit()
      } else if (at_limit > 0) { # individual will be censored anyway, so just return interval[2]
          return(interval[2])
      } else {
        t_i <- stats::uniroot(
          rootfn_hazard, hazard = hazard01, x = x_i, betas = betas_i,
          log_u = log_u_i, qq = qq,  interval = interval)$root
      }
      return(t_i)
    })

    D <- sapply(ids, function(i) {
      x_i <- subset_df(x, i, idvar = idvar)
      betas_i <- subset_df(betas02, i, idvar = idvar)
      log_u_i <- log(stats::runif(1))
      # check whether S(t) is still greater than random uniform variable u_i at the
      # upper limit of uniroot's interval (otherwise uniroot will return an error)
      at_limit <- rootfn_hazard(interval[2], hazard = hazard02, x = x_i,
                                betas = betas_i, log_u = log_u_i, qq = qq)
      if (is.nan(at_limit)) {
        STOP_nan_at_limit()
      } else if (at_limit > 0) { # individual will be censored anyway, so just return interval[2]
        return(interval[2])
      } else {
        t_i <- stats::uniroot(
          rootfn_hazard, hazard = hazard02, x = x_i, betas = betas_i,
          log_u = log_u_i, qq = qq,  interval = interval)$root
      }
      return(t_i)
    })

    yesR <- R < D

    tt12 <- sapply(seq_along(ids[yesR]), function(i) {
      x_i <- subset_df(x, i, idvar = idvar)
      betas_i <- subset_df(betas12, i, idvar = idvar)
      log_u_i <- log(stats::runif(1))
      # check whether S(t) is still greater than random uniform variable u_i at the
      # upper limit of uniroot's interval (otherwise uniroot will return an error)
      at_limit <- rootfn_hazard(interval[2], hazard = hazard12, x = x_i,
                                betas = betas_i, log_u = log_u_i, qq = qq)
      if (is.nan(at_limit)) {
        STOP_nan_at_limit()
      } else if (at_limit > 0) { # individual will be censored anyway, so just return interval[2]
        return(interval[2])
      } else {
        t_i <- stats::uniroot(
          rootfn_hazard, hazard = hazard12, x = x_i, betas = betas_i,
          log_u = log_u_i, qq = qq,  interval = interval)$root
      }
      return(t_i)
    })
  }

  D[yesR] <- R[yesR] + tt12

  delta1 <- rep(NA, N)
  delta2 <- rep(NA, N)
  y1 <- R
  y2 <- D


  if(missing(cens)){
    Cen <- c(max(c(y1, y2)) + 1, max(c(y1, y2)) + 2)
  } else {
    Cen <- runif(N, cens[1], cens[2])
  }

  ind01 <- which(D < R & D < Cen)
  y1[ind01] <- D[ind01]
  delta1[ind01] <- 0
  delta2[ind01] <- 1
  ind10 <- which(R < D & R < Cen & D >= Cen)
  y2[ind10] <- Cen[ind10]
  delta1[ind10] <- 1
  delta2[ind10] <- 0
  ind00 <- which(R >= Cen & D >= Cen)
  y1[ind00] <- Cen[ind00]
  y2[ind00] <- Cen[ind00]
  delta1[ind00] <- 0
  delta2[ind00] <- 0
  ind11 <- which(R < Cen & D < Cen & R < D)
  delta1[ind11] <- 1
  delta2[ind11] <- 1
  ret <- data.frame(cbind(id = ids, df_time = y1, df_event = delta1, os_time = y2, os_event = delta2))
  ret <- suppressMessages(dplyr::left_join(ret, x))
  #ret$time_diff <- ret$os_time - ret$df_time
  return(ret)
}





#
#
# rsimms_fn <- function(gamma01 = 0, shape01 = 0, beta01 = 0, covs = 0, maxt = 50, lambdas02 = 0, shape02 = 0, beta02 = 0, lambdas12 = 0, shape12 = 0, beta12 = 0, x = 0){
#
#   x = am(x)
#
#   fn01 <- function(t, x, betas,...){
#     (-1 + 0.02 * t - 0.03 * t ^ 2 + 0.005 * t ^ 3 + exp(beta01 * x))
#   }
#   s1 <- simsurv(loghazard = fn01, x = covs, maxt = 1.5)
#   colnames(s1)[2:3] <- c("yr", "dr")
#
#   # simulate non-terminal event
#   fn02 <- function(t, x, betas, ...){
#     (-1 + 0.02 * t - 0.03 * t ^ 2 + 0.005 * t ^ 3 + exp(beta02 * x))
#   }
#   s2 <- simsurv(loghazard = fn02, x = covs, maxt = 1.5)
#   head(s2)
#
#   colnames(s2)[2:3] <- c("yt", "dt")
#   # s2 and s1
#   s <- dplyr::left_join(s1, s2)
#
#   s$ostime <- with(s, pmin(yt, yr))
#   s$dftime <- with(s, pmin(yt, yr))
#   s$osevent <- s$dt
#   s$osevent[s$ostime == s$yr] <- 0
#   s$dfevent <- s$dr
#   s$dfevent[s$dftime == s$yt] <- 0
#
#   covs3 <- covs[s$dfevent == 1, ]
#   x3 <- x[s$dfevent == 1, ]
#
#   fn12 <- function(t, x, betas ,...){
#     (-1 + 0.02 * t - 0.03 * t ^ 2 + 0.005 * t ^ 3 + exp(beta12 * x3))
#   }
#
#   s3 <- simsurv(loghazard = fn12, x = covs3, maxt = 0.02)
#   head(s3)
#   colnames(s2)[2:3] <- c("yt", "dt")
#   s$ostime[s$dfevent == 1] <- s3$eventtime + s$ostime[s$dfevent == 1]
#   s$osevent[s$dfevent == 1] <- s3$status
#
#   s$trt <- covs$trt
#   s$timediff <- s$ostime - s$dftime
#   return(s)
# }
