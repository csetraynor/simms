#' Simulate a illness-death set with longitudinal data
#' @export
sim_idm_jm <- function(n = 200, M = 1, n_trans = 3,
                  fixed_trajectory = "cubic",
                  random_trajectory = "linear",
                  assoc = c("etavalue"),
                  basehaz = c("weibull"),
                  betaLong_intercept = 10,
                  betaLong_binary = -1,
                  betaLong_continuous = 1,
                  betaLong_linear = -0.25,
                  betaLong_quadratic = 0.03,
                  betaLong_cubic = -0.0015,
                  betaLong_aux = 0.5,
                  betaEvent_intercept01 = -7.5,
                  betaEvent_intercept02 = -6.5,
                  betaEvent_intercept12 = -5.5,
                  betaEvent_binary01 = -0.5,
                  betaEvent_binary02 = -0.4,
                  betaEvent_binary12 = -0.3,
                  betaEvent_continuous01 = 0.5,
                  betaEvent_continuous02 = 0.4,
                  betaEvent_continuous12 = 0.3,
                  betaEvent_assoc01 = 0.55,
                  betaEvent_assoc02 = 0.57,
                  betaEvent_assoc12 = 0.5,
                  betaEvent_aux01 = 1.2,
                  betaEvent_aux02 = 0.8,
                  betaEvent_aux12 = 0.9,
                  b_sd = c(1.5, 0.07), b_rho = -0.2,
                  prob_Z1 = 0.5,
                  mean_Z2 = 0, sd_Z2 = 1,
                  max_yobs = 10,
                  cens = c(199, 199),
                  max_fuptime = 20,
                  balanced = FALSE,
                  family = gaussian,
                  clust_control = NULL,
                  return_eta = FALSE,
                  seed = sample.int(.Machine$integer.max, 1),
                  interval = c(1E-8, 200)) {
  
  #----- Preliminaries
  
  set.seed(seed)
  
  basehaz <- match.arg(basehaz)
  
  if (max_yobs < 1)
    stop("'max_yobs' must be at least 1.")
  
  # Check input assoc is valid
  ok_assocs <- c("etavalue", "etaslope", "etaauc", "muvalue",
                 "null", "shared_b(1)", "shared_coef(1)",
                 "shared_b(2)", "shared_coef(2)")
  assoc <- maybe_broadcast(assoc, M)
  if (!all(assoc %in% ok_assocs))
    stop("'assoc' must be one of: ", paste(ok_assocs, collapse = ", "))
  
  # Check input to trajectory type is valid
  ok_trajs  <- c("none", "linear", "quadratic", "cubic")
  fixed_trajectory  <- maybe_broadcast(fixed_trajectory,  M)
  random_trajectory <- maybe_broadcast(random_trajectory, M)
  assoc <- maybe_broadcast(assoc, n_trans)
  basehaz <- maybe_broadcast(basehaz, n_trans)
  if (!all(fixed_trajectory %in% ok_trajs))
    stop("'fixed_trajectory' must be one of: ", paste(ok_trajs, collapse = ", "))
  if (!all(random_trajectory %in% ok_trajs))
    stop("'random_trajectory' must be one of: ", paste(ok_trajs, collapse = ", "))
  if (!length(fixed_trajectory) == M)
    stop("'fixed_trajectory' is the wrong length.")
  if (!length(random_trajectory) == M)
    stop("'random_trajectory' is the wrong length.")
  
  # Check family is valid
  if (!is(family, "list"))
    family <- list(family)
  family <- maybe_broadcast(family, M)
  family <- lapply(family, validate_family)
  ok_families <- c("gaussian", "binomial")
  lapply(family, function(x) if (!x$family %in% ok_families)
    stop("'family' must be one of: ", paste(ok_families, collapse = ", ")))
  
  # Error check to ensure that the random effects structure isn't
  # more complex than the corresponding fixed effects structure
  fixed_traj_index  <- match(fixed_trajectory,  ok_trajs)
  random_traj_index <- match(random_trajectory, ok_trajs)
  sel <- which(random_traj_index > fixed_traj_index)
  if (length(sel))
    stop("The 'random_trajectory' cannot be more complex than the ",
         "corresponding 'fixed_trajectory'. This problem was encountered ",
         "for longitudinal submodel(s): ", paste(sel, collapse = ", "))
  
  # Error checks for assoc type
  for (m in 1:M) {
    shared_slope <- (assoc[m] %in% c("shared_b(2)", "shared_coef(2)"))
    if (shared_slope && (!random_trajectory[m] == "linear"))
      stop("Can only use 'shared_b(2)' or 'shared_coef(2)' when ",
           "'random_trajectory' is linear.")
  }
  
  # Check specified structure for any lower-level clustering
  has_clust <- !missing(clust_control)
  if (has_clust) { # has lower-level clustering within the individual
    if (!is(clust_control, "list"))
      stop("'clust_control' should be a named list.")
    clust_nms <- names(clust_control)
    ok_clust_args <- c("L", "assoc", "random_trajectory", "u_sd", "u_rho")
    if (!all(clust_nms %in% ok_clust_args))
      stop("'clust_control' should only include the following named arguments: ",
           paste(ok_clust_args, collapse = ", "))
    req_clust_args <- c("L", "assoc", "random_trajectory", "u_sd")
    if (!all(req_clust_args %in% clust_nms))
      stop("'clust_control' must include the following named arguments: ",
           paste(req_clust_args, collapse = ", "))
    if (clust_control$L < 1)
      stop("In clust_control, 'L' should be a positive integer.")
    if (!is.numeric(clust_control$u_sd) || any(clust_control$u_sd < 0))
      stop("In clust_control, 'u_sd' must be a positive scalar.")
    ok_clust_assocs <- c("sum", "mean", "max", "min")
    if (!clust_control$assoc %in% ok_clust_assocs)
      stop("In clust_control, 'assoc' must be one of: ",
           paste(ok_clust_assocs, collapse = ", "))
    ok_clust_trajs  <- c("none", "linear")
    if (!(clust_control$random_trajectory %in% ok_clust_trajs))
      stop("'clust_control$random_trajectory' must be one of: ",
           paste(ok_clust_trajs, collapse = ", "))
    marker1_traj_types <- c(random_trajectory[1L], clust_control$random_trajectory)
    if (!any(marker1_traj_types == "none")) {
      stop("If lower-level clustering within the individual is specified, ",
           "then the random effects must be limited to a random intercept ",
           "(ie. 'random_trajectory = none') at either the individual-level ",
           "or the clustering-level within individuals.")
    }
    grp_assoc <- clust_control$assoc
  } else {
    grp_assoc <- NULL
  }
  
  #----- Parameters
  
  # Broadcast user-specified parameters
  betaLong_intercept  <- maybe_broadcast(betaLong_intercept,  M)
  betaLong_binary     <- maybe_broadcast(betaLong_binary,     M)
  betaLong_continuous <- maybe_broadcast(betaLong_continuous, M)
  betaLong_linear     <- maybe_broadcast(betaLong_linear,     M)
  betaLong_quadratic  <- maybe_broadcast(betaLong_quadratic,  M)
  betaLong_cubic      <- maybe_broadcast(betaLong_cubic,      M)
  betaLong_aux        <- maybe_broadcast(betaLong_aux,        M)
  betaEvent_assoc01     <- maybe_broadcast(betaEvent_assoc01,     M)
  betaEvent_assoc02     <- maybe_broadcast(betaEvent_assoc02,     M)
  betaEvent_assoc12     <- maybe_broadcast(betaEvent_assoc12,     M)
  
  # Draw individual-level REs
  b_dim <- sapply(random_trajectory, function(x)
    switch(x,
           none      = 1L, # random intercepts model
           linear    = 2L, # random slopes model
           quadratic = 3L, # random quadratic terms
           cubic     = 4L) # random cubic terms
  )
  b_dim_total <- sum(b_dim) # total num of individual-level REs
  
  # Validate b_sd
  if (length(unique(b_dim)) == 1L) { # same num of REs in each submodel
    if (length(b_sd) == b_dim[1]) {
      b_sd <- rep(b_sd, times = M)
    }
  }
  if (!length(b_sd) == b_dim_total)
    stop("'b_sd' appears to be the wrong length.")
  
  # Validate b_rho
  if (b_dim_total == 1) { # only one RE, no corr matrix needed
    b_corr_mat = matrix(1,1,1)
  } else { # >1 RE, requires corr matrix
    if (is.scalar(b_rho) && b_dim_total > 1L) {
      # user supplied a constant correlation term for REs
      b_corr_mat <- matrix(rep(b_rho, b_dim_total ^ 2), ncol = b_dim_total)
      diag(b_corr_mat) <- 1
    } else {
      # user supplied a correlation matrix for REs
      b_corr_mat <- validate_corr_matrix(b_rho)
    }
  }
  
  # Draw standardised REs and scale them by b_sd
  b_dd <- MASS::mvrnorm(n = n, mu = rep(0, b_dim_total), Sigma = b_corr_mat)
  b <- sapply(1:length(b_sd), function(x) b_sd[x] * b_dd[,x])
  colnames(b) <- paste0("b", 1:b_dim_total)
  
  # Draw random effects if lower level clustering within individuals
  # NB lower level clustering only applies to the first longitudinal submodel
  if (has_clust) {
    L <- clust_control$L
    u_sd <- clust_control$u_sd
    u_dim <- switch(clust_control$random_trajectory,
                    none   = 1L, # random intercepts model
                    linear = 2L) # random slopes model
    if (!length(u_sd) == u_dim)
      stop("In clust_control, 'u_sd' appears to be the wrong length. ",
           "Should be length ", u_dim, ".")
    Li <- as.integer(stats::runif(n, 1, L+1)) # num units within each individual
    if (u_dim == 1) { # only one RE, no corr matrix needed
      u_corr_mat = matrix(1,1,1)
    } else { # >1 RE, requires corr matrix
      u_rho <- clust_control$u_rho
      if (is.null(u_rho))
        stop("'u_rho' must be provided in the clust_control list when there is ",
             "more than one random effect for the lower level clustering factor.")
      if (is.scalar(u_rho) && u_dim > 1L) {
        # user supplied a constant correlation term for REs
        u_corr_mat <- matrix(rep(u_rho, u_dim ^ 2), ncol = u_dim)
        diag(u_corr_mat) <- 1
      } else {
        # user supplied a correlation matrix for REs
        u_corr_mat <- validate_corr_matrix(u_rho)
      }
    }
    u_dd <- MASS::mvrnorm(n = sum(Li), mu = rep(0, u_dim), Sigma = u_corr_mat)
    u <- sapply(1:length(clust_control$u_sd), function(x) u_sd[x] * u_dd[,x])
    colnames(u) <- paste0("u", 1:u_dim)
  }
  
  # Construct data frame of parameters
  betas <- list()
  
  # Longitudinal submodel parameters
  for (m in 1:M) {
    nm <- paste0("Long", m)
    betas[[nm]] <- data.frame(id = 1:n)
    
    # fixed effect parameters
    betas[[nm]][[paste0("betaLong_intercept", m)]] <- rep(betaLong_intercept[m], n)
    betas[[nm]][[paste0("betaLong_binary", m)]] <- rep(betaLong_binary[m], n)
    betas[[nm]][[paste0("betaLong_continuous", m)]] <- rep(betaLong_continuous[m], n)
    if (fixed_trajectory[m] %in% c("linear", "quadratic", "cubic")) {
      # fixed effect linear
      betas[[nm]][[paste0("betaLong_linear", m)]] <- rep(betaLong_linear[m], n)
    }
    if (fixed_trajectory[m] %in% c("quadratic", "cubic")) {
      # fixed effect quadratic
      betas[[nm]][[paste0("betaLong_quadratic", m)]] <- rep(betaLong_quadratic[m], n)
    }
    if (fixed_trajectory[m] %in% c("cubic")) {
      # fixed effect cubic
      betas[[nm]][[paste0("betaLong_cubic", m)]] <- rep(betaLong_cubic[m], n)
    }
    
    # add on subject-specific intercept
    shift <- if (m == 1) 0 else sum(b_dim[1:(m - 1)])
    b_idx <- shift + 1
    betas[[nm]][[paste0("betaLong_intercept",  m)]] <-
      betas[[nm]][[paste0("betaLong_intercept",  m)]] + b[, b_idx]
    
    # add on subject-specific linear term
    if (random_trajectory[m] %in% c("linear", "quadratic", "cubic")) {
      b_idx <- shift + 2
      betas[[nm]][[paste0("betaLong_linear",  m)]] <-
        betas[[nm]][[paste0("betaLong_linear",  m)]] + b[, b_idx]
    }
    
    # add on subject-specific quadratic term
    if (random_trajectory[m] %in% c("quadratic", "cubic")) {
      b_idx <- shift + 3
      betas[[nm]][[paste0("betaLong_quadratic",  m)]] <-
        betas[[nm]][[paste0("betaLong_quadratic",  m)]] + b[, b_idx]
    }
    
    # add on subject-specific cubic term
    if (random_trajectory[m] %in% c("cubic")) {
      b_idx <- shift + 4
      betas[[nm]][[paste0("betaLong_cubic",  m)]] <-
        betas[[nm]][[paste0("betaLong_cubic",  m)]] + b[, b_idx]
    }
  }
  
  # Additional clustering factors (only allowed for longitudinal submodel 1)
  if (has_clust) { # expand rows and add on cluster-specific random effects
    betas[["Long1"]] <-
      betas[["Long1"]][rep(row.names(betas[["Long1"]]), Li), , drop = FALSE]
    # cluster id within each individual
    betas[["Long1"]]$clust <- sequence(Li)
    # unique cluster id
    betas[["Long1"]]$clust_id <-
      paste(betas[["Long1"]][["id"]], betas[["Long1"]]$clust, sep = "_")
    # add on cluster-specific intercept
    betas[["Long1"]][["betaLong_intercept1"]] <-
      betas[["Long1"]][["betaLong_intercept1"]] + u[, 1]
    # add on cluster-specific linear slope
    if (clust_control$random_trajectory %in% "linear") {
      betas[["Long1"]][["betaLong_linear1"]] <-
        betas[["Long1"]][["betaLong_linear1"]] + u[, 2]
    }
  }
  
  # Event submodel parameters
  betas[["Event01"]] <- data.frame(
    id = 1:n,
    
    betaEvent_intercept  = rep(betaEvent_intercept01,  n),
    betaEvent_binary    = rep(betaEvent_binary01,     n),
    betaEvent_continuous = rep(betaEvent_continuous01, n),
    betaEvent_aux        = rep(betaEvent_aux01,        n)
  )
  betas[["Event02"]] <- data.frame(
    id = 1:n,
    
    betaEvent_intercept  = rep(betaEvent_intercept02,  n),
    betaEvent_binary     = rep(betaEvent_binary02,     n),
    betaEvent_continuous = rep(betaEvent_continuous02, n),
    betaEvent_aux        = rep(betaEvent_aux02,        n)
  )
  betas[["Event12"]] <- data.frame(
    id = 1:n,
    
    betaEvent_intercept  = rep(betaEvent_intercept12,  n),
    betaEvent_binary     = rep(betaEvent_binary12,     n),
    betaEvent_continuous = rep(betaEvent_continuous12, n),
    betaEvent_aux        = rep(betaEvent_aux12,        n)
  )
  
  for (m in 1:M) {
    betas[["Event01"]][[paste0("betaEvent_assoc", m)]] <- rep(betaEvent_assoc01[m], n)
    betas[["Event02"]][[paste0("betaEvent_assoc", m)]] <- rep(betaEvent_assoc02[m], n)
    betas[["Event12"]][[paste0("betaEvent_assoc", m)]] <- rep(betaEvent_assoc12[m], n)
  }
  
  #----- Data
  
  # Generate baseline covariates - binary
  Z1 <- stats::rbinom(n, 1, prob_Z1)      # covariate value for each subject
  
  # Generate baseline covariates - continuous
  Z2 <- stats::rnorm(n, mean_Z2, sd_Z2)   # covariate value for each subject
  
  # Construct data frame of baseline covariates
  covs <- data.frame(id = 1:n, Z1, Z2)
  
  # Generate survival times
  transitions <- c("Event01", "Event02", "Event12")
  
  ss01 <- simsurv(# the following arguments apply to 'simsurv'
    hazard = jmmstte_hazard, x = covs, betas = betas, 
    idvar = "id", ids = covs$id,
    maxt = max_fuptime, interval = interval,
    # the following arguments apply to 'jm_hazard'
    basehaz = basehaz[1], trans = 1, transitions = transitions,
    M = M, trajectory = fixed_trajectory,
    assoc = assoc, family = family, grp_assoc = grp_assoc)
  
  ss02 <- simsurv(# the following arguments apply to 'simsurv'
    hazard = jmmstte_hazard, x = covs, betas = betas, 
    idvar = "id", ids = covs$id,
    maxt = max_fuptime, interval = interval,
    # the following arguments apply to 'jmmstte_hazard'
    basehaz = basehaz[2], trans = 2, transitions = transitions,
    M = M, trajectory = fixed_trajectory,
    assoc = assoc, family = family, grp_assoc = grp_assoc)
  
  R <- ss01$eventtime
  D <- ss02$eventtime
  yesR <- R < D
  
  covs12 <- covs[yesR, ]
  betas12 <- lapply(betas, function(d)  d[yesR, ])
  tstart <- ss01$eventtime[yesR]
  
  ss12 <- simsurv(# the following arguments apply to 'simsurv'
    hazard = jmmstte_hazard, x = covs12, betas = betas12, 
    idvar = "id", ids = covs12$id, tstart = tstart,
    maxt = max_fuptime, interval = interval,
    # the following arguments apply to 'jmmstte_hazard'
    basehaz = basehaz[3], trans = 3, transitions = transitions,
    M = M, trajectory = fixed_trajectory,
    assoc = assoc, family = family, grp_assoc = grp_assoc)

  D[yesR] <- R[yesR] + ss12$eventtime
  
  delta1 <- rep(NA, n)
  delta2 <- rep(NA, n)
  y1 <- R
  y2 <- D
  # if(missing(cens)){
  #   Cen <- c(max(c(y1, y2)) + 1, max(c(y1, y2)) + 2)
  # } else {
  #   Cen <- runif(n, cens[1], cens[2])
  # }
  
  Cen <- runif(n, cens[1], cens[2])
  
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
  ss <- data.frame(cbind(id = covs$id, df_time = y1, df_event = delta1, os_time = y2, os_event = delta2))
  
  # Construct data frame of event data
  dat <- list(
    Event = dplyr::left_join( covs, ss, by = "id") # single row per subject
  )
  
  # Construct data frame of longitudinal data
  for (m in 1:M) {
    nm <- paste0("Long", m)
    dat[[nm]] <- merge(betas[[nm]], dat[["Event"]])
    dat[[nm]] <- merge(dat[[nm]], covs)
    dat[[nm]] <- dat[[nm]][rep(row.names(dat[[nm]]), max_yobs), ] # multiple row per subject
    # create observation times
    if (balanced) { # longitudinal observation times balanced across individuals
      tij_seq  <- max_fuptime * 0:(max_yobs - 1) / max_yobs # baseline and post-baseline
      dat[[nm]]$tij <- rep(tij_seq, each = nrow(betas[[nm]]))
    } else { # longitudinal observation times unbalanced across individuals
      tij_seq1 <- rep(0, nrow(betas[[nm]])) # baseline
      tij_seq2 <- stats::runif(nrow(dat[[nm]]) - nrow(betas[[nm]]), 0, max_fuptime) # post-baseline
      dat[[nm]]$tij <- c(tij_seq1, tij_seq2)
    }
    if (has_clust && m == 1) { # sort on ID, cluster ID and time
      dat[[nm]] <- dat[[nm]][order(dat[[nm]]$id, dat[[nm]]$clust_id, dat[[nm]]$tij), ]
    } else { # sort on ID and time
      dat[[nm]] <- dat[[nm]][order(dat[[nm]]$id, dat[[nm]]$tij), ]
    }
    dat[[nm]][[paste0("Xij_", m)]] <-
      dat[[nm]][[paste0("betaLong_intercept", m)]] +
      dat[[nm]][[paste0("betaLong_binary", m)]] * dat[[nm]][["Z1"]] +
      dat[[nm]][[paste0("betaLong_continuous", m)]] * dat[[nm]][["Z2"]]
    if (fixed_trajectory[m] %in% c("linear", "quadratic", "cubic")) {
      dat[[nm]][[paste0("Xij_", m)]] <- dat[[nm]][[paste0("Xij_", m)]] +
        dat[[nm]][[paste0("betaLong_linear", m)]] *
        dat[[nm]][["tij"]]
    }
    if (fixed_trajectory[m] %in% c("quadratic", "cubic")) {
      dat[[nm]][[paste0("Xij_", m)]] <- dat[[nm]][[paste0("Xij_", m)]] +
        dat[[nm]][[paste0("betaLong_quadratic", m)]] *
        dat[[nm]][["tij"]] *
        dat[[nm]][["tij"]]
    }
    if (fixed_trajectory[m] %in% c("cubic")) {
      dat[[nm]][[paste0("Xij_", m)]] <- dat[[nm]][[paste0("Xij_", m)]] +
        dat[[nm]][[paste0("betaLong_cubic", m)]] *
        dat[[nm]][["tij"]] *
        dat[[nm]][["tij"]] *
        dat[[nm]][["tij"]]
    }
    fam <- family[[m]]$family
    invlink <- family[[m]]$linkinv
    mu <- invlink(dat[[nm]][[paste0("Xij_", m)]])
    if (fam == "gaussian") {
      sigma <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- stats::rnorm(length(mu), mu, sigma)
    } else if (fam == "binomial") {
      trials <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- stats::rbinom(length(mu), trials, mu)
    } else if (fam == "poisson") {
      dat[[nm]][[paste0("Yij_", m)]] <- stats::rpois(length(mu), mu)
    } else if (fam == "neg_binomial_2") {
      size <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- stats::rnbinom(length(mu), size, mu)
    } else if (fam == "Gamma") {
      shape <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- stats::rgamma(length(mu), shape, rate = shape / mu)
    } else if (fam == "inverse.gaussian") {
      lambda <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- .rinvGauss(length(mu), mu, lambda)
    }
  }
  
  # Prepare final datasets
  ret <- lapply(dat, function(x) {
    if ("tij" %in% colnames(x))
      x <- x[x$tij <= x$os_time, ] # only keep rows before absorving time
    if (return_eta) {
      sel <- grep("^id$|^clust|^Z|^tij|^Yij|^Xij|os_time|os_event|df_time|df_event", colnames(x))
    } else {
      sel <- grep("^id$|^clust|^Z|^tij|^Yij|os_time|os_event|df_time|df_event", colnames(x))
    }
    x <- x[, sel, drop = FALSE]
    rownames(x) <- NULL
    return(x)
  })
  
  # Only keep individuals who have at least one measurement for each biomarker
  # (NB Following the issues around delayed entry, this is now all individuals,
  #     since every individual is forced to have a baseline measurement. This
  #     may change in the future though, if there is a need to test whether
  #     rstanarm is correctly accommodating delayed entry, i.e. the situation
  #     in which we exclude individuals who do not have a baseline measurement.)
  commonids <- ret[[1L]]$id
  for (i in 1:length(ret))
    commonids <- intersect(commonids, ret[[i]]$id)
  ret <- lapply(ret, function(x) x[x$id %in% commonids, , drop = FALSE])
  
  # Store 'true' parameter values
  long_params <- nlist(
    betaLong_intercept,
    betaLong_binary,
    betaLong_continuous,
    betaLong_aux
  )
  if (any(fixed_trajectory %in% c("linear", "quadratic", "cubic")))
    long_params$betaLong_linear <- betaLong_linear
  if (any(fixed_trajectory %in% c("quadratic", "cubic")))
    long_params$betaLong_quadratic <- betaLong_quadratic
  if (any(fixed_trajectory %in% c("cubic")))
    long_params$betaLong_cubic <- betaLong_cubic
  event_params <- nlist(
    betaEvent_intercept01,
    betaEvent_binary01,
    betaEvent_continuous01,
    betaEvent_assoc01,
    betaEvent_aux01,
    
    betaEvent_intercept02,
    betaEvent_binary02,
    betaEvent_continuous02,
    betaEvent_assoc02,
    betaEvent_aux02,
    
    betaEvent_intercept12,
    betaEvent_binary12,
    betaEvent_continuous12,
    betaEvent_assoc12,
    betaEvent_aux12
  )
  re_params <- nlist(
    b_sd,
    b_corr = b_corr_mat
  )
  cov_params <- nlist(
    prob_Z1, mean_Z2, sd_Z2
  )
  
  # Return object
  structure(ret,
            params = c(long_params, event_params, re_params, cov_params),
            n = length(unique(ret$Event$id)),
            M = M,
            max_yobs = max_yobs,
            max_fuptime = max_fuptime,
            balanced = balanced,
            assoc = assoc,
            family = family,
            fixed_trajectory = fixed_trajectory,
            random_trajectory = random_trajectory,
            return_eta = return_eta,
            clust_control = clust_control,
            seed = seed,
            class = c("simjm", class(ret)))
}


#--------------------- internal

# Return the hazard at time t, for a shared parameter joint model
#
# @param t The time variable
# @param x A named vector of covariate values
# @param betas A named vector of parameter values
# @param basehaz The type of baseline hazard
# @param M The number of longitudinal submodels
# @param assoc A vector of character strings, indicating the association
#   structure to use for each longitudinal submodel
# @param A list of family objects, providing the family (and link function)
#   for each longitudinal submodel
# @return A scalar
jmmstte_hazard <- function(t, x, betas, trans, transitions, basehaz = "weibull", M = 1,
                      trajectory = "linear", assoc = "etavalue", tstart_i,
                      family = list(gaussian()), grp_assoc = NULL) {
  
  if (t == 0)
    return(0) # boundary condition
  
  event = transitions[trans]
  
  # Baseline hazard
  if (basehaz == "weibull") {
    h0 <- betas[[event]][["betaEvent_aux"]] *
      (t ^ (betas[[event]][["betaEvent_aux"]] - 1))
  }
  
  # Time-invariant part of event submodel eta
  etaevent <-
    betas[[event]][["betaEvent_intercept"]] +
    betas[[event]][["betaEvent_binary"]] * x[["Z1"]] +
    betas[[event]][["betaEvent_continuous"]] * x[["Z2"]]
  
  # Add start time
  t = t + tstart_i$tstart
  
  # Association structure
  for (m in 1:M) {
    
    nm <- paste0("Long", m)
    
    etabaseline_m <- etavalue_m <-
      betas[[nm]][[paste0("betaLong_intercept", m)]] +
      betas[[nm]][[paste0("betaLong_binary", m)]] * x[["Z1"]] +
      betas[[nm]][[paste0("betaLong_continuous", m)]] * x[["Z2"]]
    if (trajectory[m] %in% c("linear", "quadratic", "cubic")) {
      etavalue_m <- etavalue_m +
        betas[[nm]][[paste0("betaLong_linear", m)]] * t
    }
    if (trajectory[m] %in% c("quadratic", "cubic")) {
      etavalue_m <- etavalue_m +
        betas[[nm]][[paste0("betaLong_quadratic", m)]] * (t * t)
    }
    if (trajectory[m] %in% c("cubic")) {
      etavalue_m <- etavalue_m +
        betas[[nm]][[paste0("betaLong_cubic", m)]] * (t * t * t)
    }
    
    if (!is.null(grp_assoc)) {
      if (grp_assoc == "sum") {
        res_etavalue_m <- sum(etavalue_m)
      } else if (grp_assoc == "mean") {
        res_etavalue_m <- mean(etavalue_m)
      }
    } else res_etavalue_m <- etavalue_m
    
    if (assoc[m] == "etavalue") {
      # eta value
      res_etavalue_m <- collapse_across_clusters(etavalue_m, grp_assoc)
      etaevent <- etaevent +
        betas[[event]][[paste0("betaEvent_assoc", m)]] * res_etavalue_m
    } else if (assoc[m] == "etaslope") {
      # eta slope
      if (trajectory[m] == "none") {
        etaslope_m <- 0
      } else if (trajectory[m] == "linear") {
        etaslope_m <-
          betas[[nm]][[paste0("betaLong_linear", m)]]
      } else if (trajectory[m] == "quadratic") {
        etaslope_m <-
          betas[[nm]][[paste0("betaLong_linear",  m)]] +
          betas[[nm]][[paste0("betaLong_quadratic",  m)]] * 2 * t
      } else if (trajectory[m] == "cubic") {
        etaslope_m <-
          betas[[nm]][[paste0("betaLong_linear",  m)]] +
          betas[[nm]][[paste0("betaLong_quadratic",  m)]] * 2 * t +
          betas[[nm]][[paste0("betaLong_cubic",  m)]] * 3 * I(t ^ 2)
      }
      res_etaslope_m <- collapse_across_clusters(etaslope_m, grp_assoc)
      etaevent <- etaevent +
        betas[[event]][[paste0("betaEvent_assoc", m)]] * res_etaslope_m
    } else if (assoc[m] == "etaauc") {
      # eta auc
      if (trajectory[m] == "none") {
        etaauc_m <-
          etabaseline_m * t
      } else if (trajectory[m] == "linear") {
        etaauc_m <-
          etabaseline_m * t +
          I(1/2) * betas[[nm]][[paste0("betaLong_linear", m)]] * I(t ^ 2)
      } else if (trajectory[m] == "quadratic") {
        etaauc_m <-
          etabaseline_m * t +
          I(1/2) * betas[[nm]][[paste0("betaLong_linear", m)]] * I(t ^ 2) +
          I(1/3) * betas[[nm]][[paste0("betaLong_linear", m)]] * I(t ^ 3)
      } else if (trajectory[m] == "cubic") {
        etaauc_m <-
          etabaseline_m * t +
          I(1/2) * betas[[nm]][[paste0("betaLong_linear", m)]] * I(t ^ 2) +
          I(1/3) * betas[[nm]][[paste0("betaLong_linear", m)]] * I(t ^ 3) +
          I(1/4) * betas[[nm]][[paste0("betaLong_linear", m)]] * I(t ^ 4)
      }
      res_etaauc_m <- collapse_across_clusters(etaauc_m, grp_assoc)
      etaevent <- etaevent +
        betas[[event]][[paste0("betaEvent_assoc", m)]] * res_etaauc_m
    } else if (assoc[m] == "muvalue") {
      # mu value
      invlink <- family[[m]]$invlink
      muvalue_m <- invlink(etavalue_m)
      res_muvalue_m <- collapse_across_clusters(muvalue_m, grp_assoc)
      etaevent <- etaevent +
        betas[[event]][[paste0("betaEvent_assoc", m)]] * res_muvalue_m
    }
    
  }
  
  # Calculate hazard
  h <- h0 * exp(etaevent)
  
  # Validate and return hazard
  if (!length(h) == 1) {
    stop("Bug found: returned hazard should be a scalar.")
  }
  return(h)
}

# Apply summary function in association structure when there is lower-level
# clustering within patients
#
# @param x A scalar or numeric vector with the 'etavalue', 'etaslope', 'etaauc'
#   value (or whatever quantity is being used in the association structure) for
#   one patient.
# @param collapse_fun A function to match and then apply to 'x'. If NULL, then
#   'x' is just returned, in which case 'x' should just be a scalar (i.e. there
#   should only be a patient-level value in the association structure and no
#   lower-level clusters within a patient).
# @return A scalar, used as the association term in the event submodel.
collapse_across_clusters <- function(x, collapse_fun = NULL) {
  if (is.null(collapse_fun)) {
    return(x)
  } else {
    fn <- match.fun(collapse_fun)
    return(fn(x))
  }
}