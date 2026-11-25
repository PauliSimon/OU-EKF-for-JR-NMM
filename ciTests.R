robust_fit <- function(nll, start, n_extra = 2, sd_log = 0.05) {
  f <- function(phi) { v <- nll(exp(phi)); if (is.finite(v)) v else 1e10 }
  starts <- c(list(start), lapply(seq_len(n_extra), function(i) start * exp(rnorm(length(start), 0, sd_log))))
  fits <- lapply(starts, function(s) optim(log(s), f, method = "BFGS",
                                           control = list(reltol = 1e-13, maxit = 500, ndeps = rep(1e-5, length(s)))))
  vals <- sapply(fits, `[[`, "value")
  best <- fits[[which.min(vals)]]
  list(par = setNames(exp(best$par), names(start)), value = best$value, all_values = vals)
}


.ci_names <- function(level) {
  paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2), trim = TRUE), "%")
}

# Per-step terms function
make_terms_fun <- function(full_param_list, variable_indices, ...,
                           link_fun = identity, nSubsteps = 1L, burnin = NULL) {
  force(full_param_list); force(variable_indices)
  nSubsteps <- as.integer(nSubsteps)
  mats_vary <- any(c(9, 11, 13) %in% variable_indices)
  function(pars) {
    pars <- link_fun(pars)
    for (idx in seq_along(variable_indices)) {
      full_param_list[[variable_indices[idx]]] <- pars[idx]
    }
    args <- .ekf_args(full_param_list, nSubsteps, recompute_mats = mats_vary)  # same filter as the nll
    bi   <- if (!is.null(burnin)) as.integer(burnin) else .ekf_burnin(args)
    do.call(construct_ekf_terms, c(args[1:15], list(nSubsteps = nSubsteps, burnin = bi)))
  }
}

# Finite-diff Hessian
.fd_stencil <- function(fun, phi, h = 1e-3, want_scores = TRUE) {
  p <- length(phi)
  E <- diag(p) * h
  S0 <- sum(fun(phi))
  Tp <- Tm <- vector("list", p)
  Sp <- Sm <- numeric(p)
  for (j in seq_len(p)) {
    Tp[[j]] <- fun(phi + E[j, ]); Tm[[j]] <- fun(phi - E[j, ])
    Sp[j] <- sum(Tp[[j]]);        Sm[j] <- sum(Tm[[j]])
  }
  H <- matrix(0, p, p)
  for (j in seq_len(p)) H[j, j] <- (Sp[j] + Sm[j] - 2 * S0) / h^2
  if (p > 1) for (i in 1:(p - 1)) for (j in (i + 1):p) {
    Spp <- sum(fun(phi + E[i, ] + E[j, ])); Smm <- sum(fun(phi - E[i, ] - E[j, ]))
    Spm <- sum(fun(phi + E[i, ] - E[j, ])); Smp <- sum(fun(phi - E[i, ] + E[j, ]))
    H[i, j] <- H[j, i] <- (Spp - Spm - Smp + Smm) / (4 * h^2)
  }
  g <- NULL
  if (want_scores) {
    g <- matrix(0, length(Tp[[1]]), p)
    for (j in seq_len(p)) g[, j] <- (Tp[[j]] - Tm[[j]]) / (2 * h)
  }
  list(H = H, g = g, S0 = S0)
}

# Newey-West
.hac_J <- function(g, L) {
  g <- scale(g, center = TRUE, scale = FALSE)
  N <- nrow(g)
  J <- crossprod(g)
  L <- min(L, N - 2)
  for (l in seq_len(L)) {
    w <- 1 - l / (L + 1)
    Gl <- crossprod(g[(l + 1):N, , drop = FALSE], g[1:(N - l), , drop = FALSE])
    J <- J + w * (Gl + t(Gl))
  }
  J
}

.safe_inverse <- function(H) {
  H <- (H + t(H)) / 2
  ev <- eigen(H, symmetric = TRUE)
  if (min(ev$values) <= 1e-10 * max(abs(ev$values))) {
    warning("Hessian not positive definite at the optimum; eigenvalues were clipped. ",
            "Check convergence / identifiability.", call. = FALSE)
    ev$values <- pmax(ev$values, 1e-10 * max(abs(ev$values)))
  }
  ev$vectors %*% diag(1 / ev$values, nrow = length(ev$values)) %*% t(ev$vectors)
}


.adj_profile_ci <- function(nll_phi, phi_hat, V_model, V_sand, level = 0.95,
                            step_se = 0.75, max_steps = 8, reltol = 1e-10) {
  p <- length(phi_hat)
  z <- qnorm(1 - (1 - level) / 2)
  nll_hat <- nll_phi(phi_hat)
  cj <- diag(V_sand) / diag(V_model)
  se_s <- sqrt(diag(V_sand))
  lims <- matrix(NA_real_, p, 2)
  n_opt <- 0L
  
  one_side <- function(j, s) {
    pts_phi <- phi_hat[j]; pts_r <- 0
    prev <- phi_hat
    for (m in seq_len(max_steps)) {
      target <- phi_hat[j] + s * m * step_se * se_s[j]
      start <- prev
      if (p > 1) start[-j] <- prev[-j] + V_model[-j, j] / V_model[j, j] * (target - prev[j])
      start[j] <- target
      if (p > 1) {
        f <- function(q) { ph <- start; ph[-j] <- q; nll_phi(ph) }
        opt <- optim(start[-j], f, method = "BFGS",
                     control = list(reltol = reltol, maxit = 200, ndeps = rep(1e-4, p - 1)))
        n_opt <<- n_opt + 1L
        val <- opt$value; sol <- start; sol[-j] <- opt$par
      } else {
        val <- nll_phi(target); sol <- start
      }
      Lam <- max(2 * (val - nll_hat), 0)
      r <- s * sqrt(Lam / cj[j])
      pts_phi <- c(pts_phi, target); pts_r <- c(pts_r, r)
      prev <- sol
      if (abs(r) > z) break
    }
    if (abs(tail(pts_r, 1)) <= z) return(NA_real_)          # never crossed
    o <- order(pts_phi); xx <- pts_phi[o]; yy <- pts_r[o]
    fr <- tryCatch(splinefun(xx, yy, method = "monoH.FC"),
                   error = function(e) approxfun(xx, yy))
    n <- length(xx)
    lo <- if (s > 0) xx[n - 1] else xx[1]
    hi <- if (s > 0) xx[n]     else xx[2]
    tryCatch(uniroot(function(x) fr(x) - s * z, c(lo, hi), tol = 1e-8)$root,
             error = function(e) NA_real_)
  }
  
  for (j in seq_len(p)) {
    lims[j, 1] <- one_side(j, -1)
    lims[j, 2] <- one_side(j, +1)
  }
  list(phi_ci = lims, n_opt = n_opt, c_j = cj)
}

#everything after the point estimate
.ekf_ci_engine <- function(theta_hat, nll_fun, terms_fun = NULL,
                           method = c("sandwich", "adjprofile", "wald"),
                           level = 0.95, hac_lag = NULL, h = 1e-3,
                           step_se = 0.75, max_steps = 8) {
  method <- match.arg(method)
  p <- length(theta_hat)
  z <- qnorm(1 - (1 - level) / 2)
  phi_hat <- log(theta_hat)
  nms <- names(theta_hat)
  
  have_terms <- !is.null(terms_fun)
  if (!have_terms && method != "wald") {
    warning("No per-step terms function available: falling back to model-based Wald (log scale).",
            call. = FALSE)
    method <- "wald"
  }
  
  if (have_terms) {
    tf <- function(phi) as.numeric(terms_fun(exp(phi))[, 1])
    st <- .fd_stencil(tf, phi_hat, h = h, want_scores = TRUE)
    nll0 <- nll_fun(theta_hat)
    # if (is.finite(nll0) && abs(st$S0 - nll0) > 1e-6 * max(1, abs(nll0))) {
    #   warning(sprintf("sum(terms) = %.6f differs from nll = %.6f; terms/nll are not the same filter.",
    #                   st$S0, nll0), call. = FALSE)
    # }
    N <- nrow(st$g)
  } else {
    st <- .fd_stencil(function(phi) nll_fun(exp(phi)), phi_hat, h = h, want_scores = FALSE)
    N <- NA_integer_
  }
  
  Hm <- (st$H + t(st$H)) / 2
  Hinv <- .safe_inverse(Hm)                    
  V_model <- Hinv
  if (have_terms) {
    L <- if (is.null(hac_lag)) ceiling(N^(1/3)) else hac_lag
    J <- .hac_J(st$g, L)
    V_sand <- Hinv %*% J %*% Hinv
    V_sand <- (V_sand + t(V_sand)) / 2
  } else V_sand <- V_model
  
  mk <- function(V) {
    se <- sqrt(pmax(diag(V), 0))
    out <- cbind(exp(phi_hat - z * se), exp(phi_hat + z * se))
    dimnames(out) <- list(nms, .ci_names(level)); out
  }
  ci_wald <- mk(V_model)
  ci_sand <- if (have_terms) mk(V_sand) else NULL
  
  ci_adj <- NULL; adj_info <- NULL
  if (method == "adjprofile") {
    nll_phi <- function(phi) { v <- nll_fun(exp(phi)); if (is.finite(v)) v else 1e10 }
    ap <- .adj_profile_ci(nll_phi, phi_hat, V_model, V_sand, level, step_se, max_steps)
    ci_adj <- exp(ap$phi_ci)
    
    miss <- is.na(ci_adj)
    if (any(miss)) {
      warning("Adjusted profile did not cross the cutoff for some limits; sandwich limits used there.",
              call. = FALSE)
      ci_adj[miss] <- ci_sand[miss]
    }
    dimnames(ci_adj) <- list(nms, .ci_names(level))
    adj_info <- list(n_opt = ap$n_opt, c_j = ap$c_j)
  }
  
  ci <- switch(method, sandwich = ci_sand, adjprofile = ci_adj, wald = ci_wald)

  diagn <- NULL
  if (have_terms) {
    zi <- terms_fun(theta_hat)[, 2]
    diagn <- c(mean = mean(zi), var = var(zi),
               acf1 = acf(zi, lag.max = 1, plot = FALSE)$acf[2],
               acf1_sq = acf(zi^2, lag.max = 1, plot = FALSE)$acf[2])
  }
  
  list(ci = ci, ci_wald = ci_wald, ci_sandwich = ci_sand, ci_adjprofile = ci_adj,
       se_log_model = sqrt(diag(V_model)), se_log_sandwich = sqrt(diag(V_sand)),
       se_ratio = sqrt(diag(V_sand) / diag(V_model)),
       V_model = V_model, V_sand = V_sand, hac_lag = if (have_terms) L else NA,
       adj_info = adj_info, innovation_diagnostics = diagn, method = method)
}

.attach_ci_info <- function(res, eng, ci_method) {
  attr(res, "ci_method") <- ci_method
  attr(res, "ci_all") <- list(wald = eng$ci_wald, sandwich = eng$ci_sandwich,
                              adjprofile = eng$ci_adjprofile)
  attr(res, "se_log_model") <- eng$se_log_model
  attr(res, "se_log_sandwich") <- eng$se_log_sandwich
  attr(res, "se_ratio") <- eng$se_ratio
  attr(res, "innovation_diagnostics") <- eng$innovation_diagnostics
  attr(res, "adj_info") <- eng$adj_info
  res
}

# Legacy
run_profile_ci_timeout <- function(obj, full_param_list, variable_indices,
                                   timeout = 3600, optim_args, setup_file = NULL,
                                   optimizer_type = c("ekf", "em", "wiener"),
                                   nSubsteps = 1L) {
  if (!requireNamespace("R.utils", quietly = TRUE)) stop("Package 'R.utils' is required.")
  optimizer_type <- match.arg(optimizer_type)
  timed_out <- FALSE; failed <- FALSE; error_message <- NULL
  value <- tryCatch({
    optimizer_builder <- switch(optimizer_type,
                                ekf = function(fpl, vi) create_optimizer_function(fpl, vi, nSubsteps = nSubsteps),
                                em = em_optimizer, wiener = wiener_optimizer)
    nll_NULL <- optimizer_builder(full_param_list, variable_indices)
    loglik <- function(pars) -nll_NULL(pars)
    R.utils::withTimeout(
      expr = profileCI::profileCI(obj, loglik = loglik, faster = TRUE, mult = 32,
                                  epsilon = 0, optim_args = optim_args),
      timeout = timeout, cpu = Inf, elapsed = timeout, onTimeout = "error")
  },
  TimeoutException = function(e) { timed_out <<- TRUE; NULL },
  error = function(e) { failed <<- TRUE; error_message <<- conditionMessage(e); NULL })
  list(timed_out = timed_out, failed = failed, error = error_message, value = value)
}

.legacy_profile <- function(fit, V, full_param_list, variable_indices, profile_timeout,
                            setup_file, optimizer_type, level, nSubsteps = 1L) {
  dimnames(V) <- list(names(fit$par), names(fit$par))
  obj <- structure(list(coefficients = fit$par, vcov = V), class = "fastmle")
  pr <- run_profile_ci_timeout(obj, full_param_list, variable_indices, timeout = profile_timeout,
                               setup_file = setup_file, optimizer_type = optimizer_type,
                               nSubsteps = nSubsteps,
                               optim_args = list(method = "L-BFGS-B", lower = c(10, 10, 10),
                                                 control = opt_control))
  if (pr$timed_out || pr$failed) {
    warning(if (pr$timed_out) sprintf("Profile CI exceeded %.1f min. Using Wald CI instead.", profile_timeout / 60)
            else paste("Profile CI failed:", pr$error, "Using Wald CI instead."), call. = FALSE)
    z <- qnorm(1 - (1 - level) / 2); se <- sqrt(pmax(diag(V), 0))
    ci <- cbind(pmax(fit$par - z * se, 10), fit$par + z * se)
    dimnames(ci) <- list(names(fit$par), .ci_names(level)); ci
  } else pr$value
}

# final thing
estimation_profile_ci <- function(unknown_parameters, observation,
                                  initial_values = NULL,
                                  nObs = length(observation),
                                  full_param_list = NULL, ...,
                                  calc.ci = TRUE,
                                  ci_method = c("sandwich", "adjprofile", "wald", "profile"),
                                  level = 0.95, hac_lag = NULL,
                                  profile_timeout = 3600,
                                  setup_file = "variableSetup.R",
                                  nSubsteps = 1L, n_extra = 2) {
  ci_method <- match.arg(ci_method)
  
  cPlace <- which("C" == unknown_parameters)
  if (length(cPlace) > 0) unknown_parameters[cPlace] <- "C1"
  variable_indices <- unname(sapply(unknown_parameters, function(i) which(param_strings %in% i)))
  if (is.null(initial_values)) initial_values <- initialGuesses[unknown_parameters]
  opt_control <- c(opt_control, list(parscale = initial_values))
  if (is.null(full_param_list)) {
    full_param_list <- usual_param_value
    full_param_list[[1]] <- observation
    full_param_list[[5]] <- nObs
    stdString <- which(param_strings == "obsStd")
    full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  }
  
  nll_NULL <- create_optimizer_function(full_param_list, variable_indices, nSubsteps = nSubsteps)
  
  rf  <- robust_fit(nll_NULL, initial_values, n_extra = n_extra)
  fit <- list(par = rf$par, value = rf$value, all_values = rf$all_values)
  if (ci_method == "profile") fit$hessian <- optimHess(fit$par, nll_NULL)
  
  if (!calc.ci) return(cbind(fit$par))
  
  if (ci_method == "profile") {
    ci <- .legacy_profile(fit, solve(fit$hessian), full_param_list, variable_indices,
                          profile_timeout, setup_file, "ekf", level, nSubsteps = nSubsteps)
    return(cbind(fit$par, ci))
  }
  
  terms_fun <- make_terms_fun(full_param_list, variable_indices, nSubsteps = nSubsteps)
  eng <- .ekf_ci_engine(fit$par, nll_NULL, terms_fun, method = ci_method,
                        level = level, hac_lag = hac_lag)
  res <- cbind(fit$par, eng$ci)
  colnames(res)[1] <- ""
  # res <- list(res,fit$all_values)
  .attach_ci_info(res, eng, eng$method)
}

################################################################
#################### WIENER EM STUFF ##########################
################################################################




# wiener_inference_plus_ci <- function(observation, full_param_list = NULL, ...,
#                                      calc.ci = TRUE,
#                                      ci_method = c("wald", "sandwich", "adjprofile", "profile"),
#                                      terms_fun = NULL, level = 0.95, hac_lag = NULL,
#                                      profile_timeout = 3600,
#                                      setup_file = file.path(getwd(), "variableSetup.R")) {
#   ci_method <- match.arg(ci_method)
#   unknown_parameters <- c("C1", "mu", "sigma")
#   initial_values <- initialGuesses[unknown_parameters]
#   variable_indices <- unname(sapply(unknown_parameters, function(i) which(param_strings %in% i)))
#   
#   if (is.null(full_param_list)) {
#     full_param_list <- usual_param_value
#     full_param_list[[2]] <- stepMatHalf.
#     full_param_list[[1]] <- observation
#     full_param_list[[3]] <- WIEcovMat(full_param_list[[13]])
#     full_param_list[[5]] <- length(observation)
#   }
#   stdString <- which(param_strings == "obsStd")
#   full_param_list[[stdString]] <- full_param_list[[stdString]]^2
#   
#   nll_NULL <- wiener_optimizer(full_param_list, variable_indices)
#   fit <- optim(par = initial_values, fn = nll_NULL, method = "L-BFGS-B",
#                lower = c(10, 10, 10), control = opt_control,
#                hessian = (calc.ci && ci_method == "profile"))
#   
#   if (!calc.ci) return(cbind(fit$par))
#   
#   if (ci_method == "profile") {
#     ci <- .legacy_profile(fit, solve(fit$hessian), full_param_list, variable_indices,
#                           profile_timeout, setup_file, "wiener", level)
#     return(cbind(fit$par, ci))
#   }
#   
#   eng <- .ekf_ci_engine(fit$par, nll_NULL, terms_fun, method = ci_method,
#                         level = level, hac_lag = hac_lag)
#   res <- cbind(fit$par, eng$ci)
#   colnames(res)[1] <- ""
#   .attach_ci_info(res, eng, eng$method)
# }
.make_terms_fun_simple <- function(full_param_list, variable_indices, cpp_fun,
                                   cov_fun, burnin = NULL) {
  force(full_param_list); force(variable_indices); force(cpp_fun); force(cov_fun)
  function(pars) {
    for (idx in seq_along(variable_indices)) {
      full_param_list[[variable_indices[idx]]] <- pars[idx]
    }
    full_param_list[[3]] <- cov_fun(full_param_list[[13]])     # [[13]] = sigma
    args <- unname(full_param_list[-13])                        # same positional layout as the nll
    extra <- if (!is.null(burnin)) list(burnin = as.integer(burnin))
    else if (length(args) >= 16) list(burnin = as.integer(args[[16]]))
    else list()                                        # -> C++ default (2), same as the nll
    do.call(cpp_fun, c(args[1:15], extra))
  }
}

make_terms_fun_wiener <- function(full_param_list, variable_indices, burnin = NULL) {
  .make_terms_fun_simple(full_param_list, variable_indices, wi_ekf_terms,
                         function(sigma) WIEcovMat(sigma), burnin)
}

make_terms_fun_em <- function(full_param_list, variable_indices, burnin = NULL) {
  .make_terms_fun_simple(full_param_list, variable_indices, em_ekf_terms,
                         function(sigma) diag(c(rep(0, 3), tau., sigma, tau.)^2) * Delta., burnin)
}

.simple_filter_ci <- function(type = c("wiener", "em"), observation,
                              unknown_parameters = c("C", "mu", "sigma"),
                              initial_values = NULL,
                              nObs = length(observation),
                              full_param_list = NULL,
                              lower = NULL,
                              use_parscale = FALSE,
                              square_obs_std = is.null(full_param_list),
                              calc.ci = TRUE,
                              ci_method = c("sandwich", "adjprofile", "wald", "profile"),
                              level = 0.95, hac_lag = NULL,
                              profile_timeout = 3600,
                              setup_file = "variableSetup.R") {
  type <- match.arg(type)
  ci_method <- match.arg(ci_method)
  
  cPlace <- which("C" == unknown_parameters)
  if (length(cPlace) > 0) unknown_parameters[cPlace] <- "C1"
  variable_indices <- unname(sapply(unknown_parameters, function(i) which(param_strings %in% i)))
  if (is.null(initial_values)) initial_values <- initialGuesses[unknown_parameters]
  
  if (is.null(lower)) lower <- rep(10, length(initial_values))
  
  opt_control <- opt_control                                         
  if (use_parscale) opt_control$parscale <- initial_values
  
  if (is.null(full_param_list)) {
    full_param_list <- usual_param_value
    full_param_list[[1]] <- observation 
    full_param_list[[5]] <- nObs
    if (type == "wiener") {
      full_param_list[[2]] <- stepMatHalf.
      full_param_list[[3]] <- WIEcovMat(full_param_list[[13]])
    }
  }

  if (type == "wiener") {
    if (!is.null(full_param_list[[2]]) && !isTRUE(all.equal(unname(full_param_list[[2]]), unname(stepMatHalf.))))
      message("wiener: slot 2 of full_param_list was not stepMatHalf.; replacing it.")
    full_param_list[[2]] <- stepMatHalf.
    full_param_list[[3]] <- WIEcovMat(full_param_list[[13]])
  }
  
  if (square_obs_std) {
    stdString <- which(param_strings == "obsStd")
    full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  }
  
  nll_NULL <- if (type == "wiener") wiener_optimizer(full_param_list, variable_indices)
  else                  em_optimizer(full_param_list, variable_indices)
  
  fit <- optim(par = initial_values, fn = nll_NULL, method = "L-BFGS-B",
               lower = lower, control = opt_control,
               hessian = (calc.ci & ci_method == "profile"))
  
  if (!calc.ci) return(cbind(fit$par))
  
  if (ci_method == "profile") {
    ci <- .legacy_profile(fit, solve(fit$hessian), full_param_list, variable_indices,
                          profile_timeout, setup_file, type, level)
    return(cbind(fit$par, ci))
  }
  
  terms_fun <- if (type == "wiener") make_terms_fun_wiener(full_param_list, variable_indices)
  else                  make_terms_fun_em(full_param_list, variable_indices)
  eng <- tryCatch(
    .ekf_ci_engine(fit$par, nll_NULL, terms_fun, method = ci_method,
                   level = level, hac_lag = hac_lag),
    error = function(e) e)
  if (inherits(eng, "error")) {
    warning("CI engine failed: ", conditionMessage(eng), ". Returning estimate with NA CIs.",
            call. = FALSE)
    res <- cbind(fit$par, NA_real_, NA_real_)
    colnames(res) <- c("", .ci_names(level))
    attr(res, "ci_error") <- conditionMessage(eng)
    attr(res, "convergence") <- fit$convergence
    return(res)
  }
  res <- cbind(fit$par, eng$ci)
  colnames(res)[1] <- ""
  res <- .attach_ci_info(res, eng, eng$method)
  attr(res, "convergence") <- fit$convergence
  attr(res, "nll") <- fit$value
  res
}


wiener_inference_plus_ci <- function(observation, full_param_list = NULL, ...,
                                     unknown_parameters = c("C", "mu", "sigma"),
                                     calc.ci = TRUE,
                                     ci_method = c("sandwich", "adjprofile", "wald", "profile"),
                                     level = 0.95, hac_lag = NULL,
                                     profile_timeout = 3600,
                                     setup_file = "variableSetup.R") {
  .simple_filter_ci("wiener", observation, unknown_parameters = unknown_parameters,
                    full_param_list = full_param_list, calc.ci = calc.ci,
                    ci_method = match.arg(ci_method), level = level, hac_lag = hac_lag,
                    profile_timeout = profile_timeout, setup_file = setup_file, ...)
}

em_inference_plus_ci <- function(observation, full_param_list = NULL, ...,
                                 unknown_parameters = c("C", "mu", "sigma"),
                                 calc.ci = TRUE,
                                 ci_method = c("sandwich", "adjprofile", "wald", "profile"),
                                 level = 0.95, hac_lag = NULL,
                                 profile_timeout = 3600,
                                 setup_file = "variableSetup.R") {
  .simple_filter_ci("em", observation, unknown_parameters = unknown_parameters,
                    full_param_list = full_param_list, calc.ci = calc.ci,
                    ci_method = match.arg(ci_method), level = level, hac_lag = hac_lag,
                    profile_timeout = profile_timeout, setup_file = setup_file, ...)
}
