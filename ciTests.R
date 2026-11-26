
.get_link <- function(link = c("log", "identity")) {
  link <- match.arg(link)
  if (link == "log") list(name = "log", fwd = log, inv = exp)
  else               list(name = "identity", fwd = identity, inv = identity)
}


.default_lower <- function(par_names, link = "identity") {
  if (link == "identity") ifelse(par_names == "obsStd", 0, 1e-6)
  else rep(0, length(par_names))
}


.lrt_cut <- function(level = 0.95) qchisq(max(1 - 2 * (1 - level), 0.5), 1) / 2


.snap_to_bounds <- function(f, par, value, lower, scale, snap_tol = 0.05) {
  p <- length(par)
  snapped <- par <= lower
  for (j in seq_len(p)) {
    if (par[j] <= lower[j] || (par[j] - lower[j]) > 0.01 * scale[j]) next
    free <- setdiff(seq_len(p), j)
    cand <- par; cand[j] <- lower[j]
    if (length(free) == 0L) {
      v <- f(cand)
    } else {
      o <- optim(par[free], function(x) { th <- cand; th[free] <- x; f(th) },
                 method = "L-BFGS-B", lower = lower[free],
                 control = list(maxit = 500, factr = 1e5, pgtol = 1e-10,
                                parscale = scale[free], ndeps = rep(1e-5, length(free))))
      cand[free] <- o$par; v <- o$value
    }
    if (is.finite(v) && v - value <= snap_tol) { par <- cand; value <- v; snapped[j] <- TRUE }
  }
  list(par = par, value = value, snapped = snapped)
}


.boundary_upper <- function(nll_fun, theta_hat, j, lower, level = 0.95, scale = pmax(abs(theta_hat), 1e-2)) {
  p <- length(theta_hat); free <- setdiff(seq_len(p), j)
  f0 <- nll_fun(theta_hat)
  cut <- .lrt_cut(level)
  prof <- function(v) {
    th <- theta_hat; th[j] <- v
    if (length(free) == 0L) return(nll_fun(th))
    o <- tryCatch(optim(theta_hat[free],
                        function(x) { t2 <- th; t2[free] <- x; val <- nll_fun(t2); if (is.finite(val)) val else 1e10 },
                        method = "L-BFGS-B", lower = lower[free],
                        control = list(maxit = 100, factr = 1e7, pgtol = 1e-6,
                                       parscale = scale[free], ndeps = rep(1e-4, length(free)))),
                  error = function(e) NULL)
    if (is.null(o)) NA_real_ else o$value
  }
  g <- function(v) prof(v) - f0 - cut
  hi <- lower[j] + 1e-6 * max(scale[j], 1)
  for (k in 1:30) {
    gv <- g(hi)
    if (is.na(gv)) return(NA_real_)
    if (gv > 0) break
    hi <- lower[j] + 4 * (hi - lower[j])
  }
  if (!is.finite(gv) || gv <= 0) return(NA_real_)
  tryCatch(uniroot(g, c(lower[j], hi), tol = max(1e-9, 1e-3 * (hi - lower[j])))$root,
           error = function(e) NA_real_)
}

robust_fit <- function(nll, start, n_extra = 2, sd_log = 0.05,
                       link = c("log", "identity"), lower = NULL, snap_tol = 0.05) {
  link <- match.arg(link)
  p <- length(start); nms <- names(start)
  starts <- c(list(start),
              lapply(seq_len(n_extra), function(i) start * exp(rnorm(p, 0, sd_log))))
  if (link == "log") {
    f <- function(phi) { v <- nll(exp(phi)); if (is.finite(v)) v else 1e10 }
    fits <- lapply(starts, function(s) optim(log(s), f, method = "BFGS",
                                             control = list(reltol = 1e-13, maxit = 500, ndeps = rep(1e-5, length(s)))))
    pars <- lapply(fits, function(x) exp(x$par))
  } else {
    if (is.null(lower)) lower <- .default_lower(nms, "identity")
    lower <- rep_len(lower, p)
    f <- function(theta) { v <- nll(theta); if (is.finite(v)) v else 1e10 }
    scale <- pmax(abs(start), 1e-2)                       # parscale: parameters differ by orders of magnitude
    fits <- lapply(starts, function(s) optim(pmax(s, lower), f, method = "L-BFGS-B", lower = lower,
                                             control = list(maxit = 500, factr = 1e5, pgtol = 1e-10,
                                                            parscale = scale, ndeps = rep(1e-5, p))))
    pars <- lapply(fits, `[[`, "par")
  }
  vals <- sapply(fits, `[[`, "value")
  b <- which.min(vals)
  par_b <- pars[[b]]; val_b <- vals[b]
  if (link == "identity" && snap_tol > 0) {
    sn <- .snap_to_bounds(f, par_b, val_b, lower, pmax(abs(start), 1e-2), snap_tol)
    par_b <- sn$par; val_b <- sn$value
  }
  list(par = setNames(par_b, nms), value = val_b, all_values = vals,
       convergence = fits[[b]]$convergence)
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

# Finite-diff Hessian (h may be a vector: one step size per parameter)
.fd_stencil <- function(fun, phi, h = 1e-3, want_scores = TRUE) {
  p <- length(phi)
  h <- rep_len(h, p)
  E <- diag(h, nrow = p)
  S0 <- sum(fun(phi))
  Tp <- Tm <- vector("list", p)
  Sp <- Sm <- numeric(p)
  for (j in seq_len(p)) {
    Tp[[j]] <- fun(phi + E[j, ]); Tm[[j]] <- fun(phi - E[j, ])
    Sp[j] <- sum(Tp[[j]]);        Sm[j] <- sum(Tm[[j]])
  }
  H <- matrix(0, p, p)
  for (j in seq_len(p)) H[j, j] <- (Sp[j] + Sm[j] - 2 * S0) / h[j]^2
  if (p > 1) for (i in 1:(p - 1)) for (j in (i + 1):p) {
    Spp <- sum(fun(phi + E[i, ] + E[j, ])); Smm <- sum(fun(phi - E[i, ] - E[j, ]))
    Spm <- sum(fun(phi + E[i, ] - E[j, ])); Smp <- sum(fun(phi - E[i, ] + E[j, ]))
    H[i, j] <- H[j, i] <- (Spp - Spm - Smp + Smm) / (4 * h[i] * h[j])
  }
  g <- NULL
  if (want_scores) {
    g <- matrix(0, length(Tp[[1]]), p)
    for (j in seq_len(p)) g[, j] <- (Tp[[j]] - Tm[[j]]) / (2 * h[j])
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

.ekf_ci_engine <- function(theta_hat, nll_fun, terms_fun = NULL,
                           method = c("sandwich", "adjprofile", "wald"),
                           level = 0.95, hac_lag = NULL, h = 1e-3,
                           step_se = 0.75, max_steps = 8,
                           link = c("log", "identity"), lower = NULL, boundary_ci = TRUE) {
  method <- match.arg(method)
  lk <- .get_link(link)
  p <- length(theta_hat)
  z <- qnorm(1 - (1 - level) / 2)
  nms <- names(theta_hat)
  if (is.null(lower)) lower <- rep(0, p)
  lower <- rep_len(lower, p)
  
  
  bound_tol <- if (lk$name == "identity") 3 * h * pmax(abs(theta_hat), 1e-2) else 1e-8 * pmax(1, abs(lower))
  at_bound <- (theta_hat - lower) <= bound_tol
  
  if (lk$name == "identity" && boundary_ci && any(is.finite(lower) & !at_bound)) {
    fsafe <- function(th) { v <- nll_fun(th); if (is.finite(v)) v else 1e10 }
    sn <- tryCatch(.snap_to_bounds(fsafe, theta_hat, fsafe(theta_hat), lower,
                                   pmax(abs(theta_hat), 1e-2), snap_tol = .lrt_cut(level)),
                   error = function(e) NULL)
    if (!is.null(sn)) at_bound <- at_bound | sn$snapped
  }
  free <- which(!at_bound)
  if (length(free) == 0L) stop("All parameters are on their lower bound; no confidence intervals can be computed.")
  if (any(at_bound))
    warning("Parameter(s) ", paste(nms[at_bound], collapse = ", "),
            " on, or statistically indistinguishable from (one-sided LR test), the lower bound: treated as fixed in the curvature/sandwich computation; ",
            "reported as c(lower, upper) with a one-sided profile-likelihood upper limit.", call. = FALSE)
  pf <- length(free)
  bound_hi <- rep(NA_real_, p)
  if (lk$name == "identity" && boundary_ci && any(at_bound)) {
    bscale <- pmax(abs(theta_hat), 1e-2)
    for (j in which(at_bound))
      bound_hi[j] <- tryCatch(.boundary_upper(nll_fun, theta_hat, j, lower, level, bscale),
                              error = function(e) NA_real_)
  }
  
  to_theta <- function(phi_f) { th <- theta_hat; th[free] <- lk$inv(phi_f); th }   # free params -> full theta
  phi_hat <- lk$fwd(theta_hat[free])
  hstep <- if (lk$name == "log") rep(h, pf) else
    pmin(h * pmax(abs(theta_hat[free]), 1e-2), (theta_hat[free] - lower[free]) / 3)  # relative step, stay inside bounds
  
  have_terms <- !is.null(terms_fun)
  if (!have_terms && method != "wald") {
    warning("No per-step terms function available: falling back to model-based Wald.",
            call. = FALSE)
    method <- "wald"
  }
  
  if (have_terms) {
    tf <- function(phi_f) as.numeric(terms_fun(to_theta(phi_f))[, 1])
    st <- .fd_stencil(tf, phi_hat, h = hstep, want_scores = TRUE)
    N <- nrow(st$g)
  } else {
    st <- .fd_stencil(function(phi_f) nll_fun(to_theta(phi_f)), phi_hat, h = hstep, want_scores = FALSE)
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
  dimnames(V_model) <- dimnames(V_sand) <- list(nms[free], nms[free])
  
  mk <- function(V) {
    se <- sqrt(pmax(diag(V), 0))
    lo <- hi <- rep(NA_real_, p)
    if (lk$name == "log") {
      lo[free] <- exp(phi_hat - z * se); hi[free] <- exp(phi_hat + z * se)
    } else {
      lo[free] <- pmax(theta_hat[free] - z * se, lower[free]); hi[free] <- theta_hat[free] + z * se
    }
    lo[at_bound] <- lower[at_bound]; hi[at_bound] <- bound_hi[at_bound]   # boundary: one-sided profile upper limit (NA if not computed)
    out <- cbind(lo, hi); dimnames(out) <- list(nms, .ci_names(level)); out
  }
  full_se <- function(V) { s <- rep(NA_real_, p); s[free] <- sqrt(pmax(diag(V), 0)); names(s) <- nms; s }
  ci_wald <- mk(V_model)
  ci_sand <- if (have_terms) mk(V_sand) else NULL
  
  ci_adj <- NULL; adj_info <- NULL
  if (method == "adjprofile") {
    nll_phi <- function(phi_f) { v <- nll_fun(to_theta(phi_f)); if (is.finite(v)) v else 1e10 }
    ap <- .adj_profile_ci(nll_phi, phi_hat, V_model, V_sand, level, step_se, max_steps)
    lims <- ap$phi_ci
    if (lk$name == "log") lims <- exp(lims) else {
      lims[, 1] <- pmax(lims[, 1], lower[free]); lims[, 2] <- pmax(lims[, 2], lower[free])
    }
    ci_adj <- matrix(NA_real_, p, 2); ci_adj[free, ] <- lims
    ci_adj[at_bound, 1] <- lower[at_bound]; ci_adj[at_bound, 2] <- bound_hi[at_bound]
    
    miss <- is.na(ci_adj); miss[at_bound, ] <- FALSE
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
       se_log_model = full_se(V_model), se_log_sandwich = full_se(V_sand),
       se_ratio = sqrt(full_se(V_sand)^2 / full_se(V_model)^2),
       V_model = V_model, V_sand = V_sand, hac_lag = if (have_terms) L else NA,
       adj_info = adj_info, innovation_diagnostics = diagn, method = method,
       link = lk$name, at_bound = setNames(at_bound, nms))
}

.attach_ci_info <- function(res, eng, ci_method) {
  attr(res, "ci_method") <- ci_method
  attr(res, "link") <- eng$link
  attr(res, "at_bound") <- eng$at_bound
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
                            setup_file, optimizer_type, level, nSubsteps = 1L,
                            lower = rep(10, length(fit$par))) {
  dimnames(V) <- list(names(fit$par), names(fit$par))
  obj <- structure(list(coefficients = fit$par, vcov = V), class = "fastmle")
  pr <- run_profile_ci_timeout(obj, full_param_list, variable_indices, timeout = profile_timeout,
                               setup_file = setup_file, optimizer_type = optimizer_type,
                               nSubsteps = nSubsteps,
                               optim_args = list(method = "L-BFGS-B", lower = lower,
                                                 control = opt_control))
  if (pr$timed_out || pr$failed) {
    warning(if (pr$timed_out) sprintf("Profile CI exceeded %.1f min. Using Wald CI instead.", profile_timeout / 60)
            else paste("Profile CI failed:", pr$error, "Using Wald CI instead."), call. = FALSE)
    z <- qnorm(1 - (1 - level) / 2); se <- sqrt(pmax(diag(V), 0))
    ci <- cbind(pmax(fit$par - z * se, lower), fit$par + z * se)
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
                                  nSubsteps = 1L, n_extra = 2,
                                  link = c("log", "identity"),
                                  lower = NULL, snap_tol = "lrt", boundary_ci = TRUE) {
  ci_method <- match.arg(ci_method)
  link <- match.arg(link)
  if (identical(snap_tol, "lrt")) snap_tol <- .lrt_cut(level)
  
  cPlace <- which("C" == unknown_parameters)
  if (length(cPlace) > 0) unknown_parameters[cPlace] <- "C1"
  variable_indices <- unname(sapply(unknown_parameters, function(i) which(param_strings %in% i)))
  if (is.null(initial_values)) initial_values <- initialGuesses[unknown_parameters]
  names(initial_values) <- unknown_parameters
  if (is.null(lower)) lower <- .default_lower(unknown_parameters, link)
  lower <- rep_len(lower, length(unknown_parameters))
  if (is.null(full_param_list)) {
    full_param_list <- usual_param_value
    full_param_list[[1]] <- observation
    full_param_list[[5]] <- nObs
    stdString <- which(param_strings == "obsStd")
    full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  }
  
  nll_NULL <- create_optimizer_function(full_param_list, variable_indices, nSubsteps = nSubsteps)
  
  rf  <- robust_fit(nll_NULL, initial_values, n_extra = n_extra, link = link, lower = lower, snap_tol = snap_tol)
  fit <- list(par = rf$par, value = rf$value, all_values = rf$all_values)
  if (link == "log" && any(unknown_parameters == "obsStd") && fit$par["obsStd"] < 1e-6)
    warning("obsStd (observation-noise variance) is ~0 with link = \"log\", where it can never reach 0 ",
            "and its confidence interval is not meaningful. Consider link = \"identity\".", call. = FALSE)
  if (ci_method == "profile") fit$hessian <- optimHess(fit$par, nll_NULL)
  
  if (!calc.ci) return(cbind(fit$par))
  
  if (ci_method == "profile") {
    ci <- .legacy_profile(fit, solve(fit$hessian), full_param_list, variable_indices,
                          profile_timeout, setup_file, "ekf", level, nSubsteps = nSubsteps,
                          lower = if (link == "identity") lower else rep(10, length(fit$par)))
    return(cbind(fit$par, ci))
  }
  
  terms_fun <- make_terms_fun(full_param_list, variable_indices, nSubsteps = nSubsteps)
  eng <- .ekf_ci_engine(fit$par, nll_NULL, terms_fun, method = ci_method,
                        level = level, hac_lag = hac_lag, link = link, lower = lower,
                        boundary_ci = boundary_ci)
  res <- cbind(fit$par, eng$ci)
  colnames(res)[1] <- ""
  res <- .attach_ci_info(res, eng, eng$method)
  attr(res, "nll") <- fit$value
  attr(res, "convergence") <- rf$convergence
  res
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
                              setup_file = "variableSetup.R",
                              link = c("log", "identity"), snap_tol = "lrt", boundary_ci = TRUE) {
  type <- match.arg(type)
  if (identical(snap_tol, "lrt")) snap_tol <- .lrt_cut(level)
  link <- match.arg(link)  
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
  
  if (snap_tol > 0) {
    sn <- .snap_to_bounds(function(th) { v <- nll_NULL(th); if (is.finite(v)) v else 1e10 },
                          fit$par, fit$value, rep_len(lower, length(fit$par)),
                          pmax(abs(initial_values), 1e-2), snap_tol)
    fit$par <- setNames(sn$par, names(fit$par)); fit$value <- sn$value
  }
  
  if (!calc.ci) return(cbind(fit$par))
  
  if (ci_method == "profile") {
    ci <- .legacy_profile(fit, solve(fit$hessian), full_param_list, variable_indices,
                          profile_timeout, setup_file, type, level, lower = lower)
    return(cbind(fit$par, ci))
  }
  
  terms_fun <- if (type == "wiener") make_terms_fun_wiener(full_param_list, variable_indices)
  else                  make_terms_fun_em(full_param_list, variable_indices)
  eng <- tryCatch(
    .ekf_ci_engine(fit$par, nll_NULL, terms_fun, method = ci_method,
                   level = level, hac_lag = hac_lag, link = link, lower = lower,
                   boundary_ci = boundary_ci),
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
                                     setup_file = "variableSetup.R",
                                     link = c("log", "identity")) {
  .simple_filter_ci("wiener", observation, unknown_parameters = unknown_parameters,
                    full_param_list = full_param_list, calc.ci = calc.ci,
                    ci_method = match.arg(ci_method), level = level, hac_lag = hac_lag,
                    profile_timeout = profile_timeout, setup_file = setup_file,
                    link = match.arg(link), ...)
}

em_inference_plus_ci <- function(observation, full_param_list = NULL, ...,
                                 unknown_parameters = c("C", "mu", "sigma"),
                                 calc.ci = TRUE,
                                 ci_method = c("sandwich", "adjprofile", "wald", "profile"),
                                 level = 0.95, hac_lag = NULL,
                                 profile_timeout = 3600,
                                 setup_file = "variableSetup.R",
                                 link = c("log", "identity")) {
  .simple_filter_ci("em", observation, unknown_parameters = unknown_parameters,
                    full_param_list = full_param_list, calc.ci = calc.ci,
                    ci_method = match.arg(ci_method), level = level, hac_lag = hac_lag,
                    profile_timeout = profile_timeout, setup_file = setup_file,
                    link = match.arg(link), ...)
}