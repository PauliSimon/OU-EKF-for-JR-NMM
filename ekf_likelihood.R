coef.fastmle <- function(object, ...) {
  object$coefficients
}

vcov.fastmle <- function(object, ...) {
  object$vcov
}
opt_control <- list(
  maxit = 10000,
  factr = 1e5,
  pgtol = 1e-10
)

.ekf_args <- function(full_param_list, nSubsteps = 1L, recompute_mats = TRUE) {
  if (nSubsteps > 1L || recompute_mats) {
    # transition mats must use INTERNAL step: Delta / nSubsteps
    delta      <- full_param_list[[6]] / nSubsteps
    dampingMat <- diag(c(full_param_list[[9]], full_param_list[[9]], full_param_list[[11]]))
    sigMat     <- diag(c(tau., full_param_list[[13]], tau.))
    full_param_list[[2]] <- expMt(dampingMat, delta)
    full_param_list[[3]] <- variableSigmaCovarMat(delta, sigMat, dampingMat)
  }
  unname(full_param_list[-13])
}

.ekf_burnin <- function(args) if (length(args) >= 16) as.integer(args[[16]]) else 2L

.ekf_nll <- function(args, nSubsteps = 1L) {
  if (nSubsteps == 1L) return(do.call(construct_ekf_likelihood, args))
  do.call(construct_ekf_combined_likelihood,
          c(args[1:15], list(nSubsteps = as.integer(nSubsteps), kStep = 1L,
                             kWeight = 1, burnin = .ekf_burnin(args))))
}

create_optimizer_function <- function(full_param_list, variable_indices, ...,
                                      link_fun = identity, nSubsteps = 1L){
  nSubsteps <- as.integer(nSubsteps)
  stopifnot(length(nSubsteps) == 1L, nSubsteps >= 1L)
  mats_vary <- any(c(9, 11, 13) %in% variable_indices)
  toOptim <- function(pars){
    pars <- link_fun(pars)
    for(idx in seq_along(variable_indices)){
      full_param_list[[variable_indices[idx]]] <- pars[idx]
    }
    if(nSubsteps>1){
      return(.ekf_nll(.ekf_args(full_param_list, nSubsteps, recompute_mats = mats_vary), nSubsteps))
    } else{
      return(do.call(construct_ekf_likelihood,.ekf_args(full_param_list, nSubsteps)[1:15]))
    }
  }
  return(toOptim)
}


optimize_parameters <- function(full_param_list,variable_indices,initials,...,link_fun=identity){
  toOptim <- create_optimizer_function(full_param_list,variable_indices)
  optim(link_fun(initials),toOptim,...)
}



library(bbmle)
# estimation_profile_ci <- function(unknown_parameters,observation,
#                             initial_values=NULL,
#                             nObs=length(observation),full_param_list=NULL,...,calc.ci=TRUE){
#     cPlace <- which("C"==unknown_parameters)
#     if(length(cPlace)>0)unknown_parameters[cPlace] <- "C1"
#     variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
#     if(is.null(initial_values))initial_values <- initialGuesses[unknown_parameters]
#     if(is.null(full_param_list)){
#       full_param_list <- usual_param_value
#       full_param_list[[1]] <- observation
#       full_param_list[[5]] <- nObs
#       stdString <- which(param_strings=="obsStd")
#       stdEstimated <- which(variable_indices==stdString)
#       
#       full_param_list[[stdString]] <- full_param_list[[stdString]]^2
#     }
#     
#     
#     nll_NULL <- create_optimizer_function(full_param_list,variable_indices)
#     loglik <- function(pars) - nll_NULL(pars)
#     
#     fit <- optim(
#       par = initial_values,
#       fn = nll_NULL,
#       method = "L-BFGS-B",lower=c(10,10,10),
#       control = opt_control,
#       hessian=TRUE
#     )
#     
#     V <- solve(fit$hessian)
#     dimnames(V) <- list(names(fit$par), names(fit$par))
#     
#     obj <- structure(
#       list(
#         coefficients = fit$par,
#         vcov = V
#       ),
#       class = "fastmle"
#     )
#     if(calc.ci){
#       ci <- profileCI::profileCI(
#         obj,
#         loglik = loglik,
#         faster = TRUE,
#         mult = 32,epsilon=0,
#         optim_args = list(
#           method = "L-BFGS-B",lower=c(10,10,10),
#           control = opt_control
#         )
#       )
#     } else {
#       ci <- NULL
#     }
#     
#     return(cbind(fit$par,ci))
#     
# # 
# #     f <- function(C1,mu,sigma){
# #       nll_NULL(c(C1,mu,sigma))
# #     }
# #     
# #     
# #     profile_fit <- function(){
# #       
# #       
# #       fit <- optim(
# #         par = c(C1 = 155, mu = 260, sigma = 2400),
# #         fn = function(p) f(p[1], p[2], p[3]),
# #         method = "L-BFGS-B",
# #         lower = c(0, 0, 0),
# #         upper = c(Inf, Inf, Inf),
# #         control = list(maxit = 5000)
# #       )
# #       fit <- mle2(f,as.list(fit$par),method="L-BFGS-B",lower=0,upper=Inf,
# #                   control = list(maxit = 50000,
# #                                  factr = 1e3,
# #                                  pgtol = 1e-10))
# #       
# #       ci <- try({
# #         prof <- profile(fit,maxsteps=1000)
# #         confint(prof,test="LRT")
# #       },silent=TRUE)
# #       
# #       if (inherits(ci, "try-error")) {
# #         H <- numDeriv::hessian(nll_NULL,fit@coef)
# #         se <- sqrt(diag(solve(H)))
# #         inherits(fit,"mle2")
# #         ci <- try({
# #           mapply(function(i)fit@coef+i*se,qnorm(c(.025,.975)))
# #           }, silent = TRUE)
# #         
# #       }
# #       if (inherits(ci, "try-error")) ci <- cbind(rep(0,2),Inf)
# #       
# #       res_mat <- cbind(fit@coef,ci)
# #       colnames(res_mat) <- c("Est","LRT_Lower","LRT_Upper")
# #       return(res_mat)
# #     }
# #     
# #     return(profile_fit())
# }


# estimation_profile_ci_obsd <- function(observation,
#                             nObs=length(observation),full_param_list=NULL,...){
#     unknown_parameters <- c("C","mu","sigma","obsStd")
#     cPlace <- which("C"==unknown_parameters)
#     if(length(cPlace)>0)unknown_parameters[cPlace] <- "C1"
#     variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
#     initial_values=initialGuesses[unknown_parameters]
#     
#     if(is.null(full_param_list)){
#       full_param_list <- usual_param_value
#       full_param_list[[1]] <- observation
#       full_param_list[[5]] <- nObs
#       stdString <- which(param_strings=="obsStd")
#       stdEstimated <- which(variable_indices==stdString)
#       
#       full_param_list[[stdString]] <- full_param_list[[stdString]]^2
#     }
#     
#     
#     nll_NULL <- create_optimizer_function(full_param_list,variable_indices,link_fun=identity)
# 
#     f <- function(C1,mu,sigma,obsStd){
#       nll_NULL(c(C1,mu,sigma,obsStd))
#     }
#     
#     
#     profile_fit <- function(){
#       fit <- optim(
#              par = c(C1 = 155, mu = 260, sigma = 2400, obsStd = 0.1),
#              fn = function(p) f(p[1], p[2], p[3], p[4]),
#              method = "L-BFGS-B",
#              lower = c(0, 0, 0, 0),
#              upper = c(Inf, Inf, Inf, Inf),
#              control = list(maxit = 5000)
#          )
#       fit <- mle2(f,as.list(fit$par),method="L-BFGS-B",lower=0,upper=Inf,
#                   control = list(maxit = 50000,
#                                  factr = 1e3,
#                                  pgtol = 1e-10))
#     
#       ci <- try({
#         prof <- profile(fit,maxsteps=1000)
#         confint(prof,test="LRT")
#       },silent=TRUE)
#       
#       if (inherits(ci, "try-error")) {
#         H <- numDeriv::hessian(nll_NULL,fit@coef)
#         se <- sqrt(diag(solve(H)))
#         inherits(fit,"mle2")
#         ci <- try({
#           mapply(function(i)fit@coef+i*se,qnorm(c(.025,.975)))
#           }, silent = TRUE)
#         
#       }
#       if (inherits(ci, "try-error")) ci <- cbind(rep(-Inf,2),Inf)
#       
#       res_mat <- cbind(fit@coef,ci)
#       colnames(res_mat) <- c("Est","LRT_Lower","LRT_Upper")
#       return(res_mat)
#     }
#     
#     return(profile_fit())
# }


# estimation_profile_ci_obsd <- function(observation,
#                                        nObs=length(observation),full_param_list=NULL,...,calc.ci=TRUE){
#   unknown_parameters <- c("C","mu","sigma","obsStd")
#   lb0      <- c(10, 10, 10, 0)
#   cPlace <- which("C"==unknown_parameters)
#   if(length(cPlace)>0)unknown_parameters[cPlace] <- "C1"
#   variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
#   initial_values=initialGuesses[unknown_parameters]
#   
#   if(is.null(full_param_list)){
#     full_param_list <- usual_param_value
#     full_param_list[[1]] <- observation
#     full_param_list[[5]] <- nObs
#   }
#   
#   
#   nll_NULL <- create_optimizer_function(full_param_list,variable_indices,link_fun=identity)
#   
#   loglik <- function(pars) - nll_NULL(pars)
#   
#   fit <- optim(
#     par = initial_values,
#     fn = nll_NULL,
#     method = "L-BFGS-B",lower=c(rep(10,3),0),
#     control = opt_control,
#     hessian=TRUE
#   )
#   
#   est <- fit$par; nm <- names(est) <- unknown_parameters
#   
#   ## scale-aware Hessian -> always-finite covariance
#   H  <- numDeriv::hessian(nll_NULL, est); H <- (H + t(H)) / 2
#   Dm <- diag(pmax(abs(est), 1e-3))                  # relative scaling
#   Hs <- Dm %*% H %*% Dm
#   ev <- eigen(Hs, symmetric = TRUE)
#   V  <- Dm %*% (ev$vectors %*% diag(1 / pmax(ev$values, 1e-8)) %*% t(ev$vectors)) %*% Dm
#   dimnames(V) <- list(nm, nm)
#   se <- sqrt(diag(V))
#   
#   ## Wald fallback (truncated at the bounds)
#   z    <- qnorm(0.975)
#   ci   <- cbind(pmax(est - z * se, lb0), est + z * se)
#   meth <- rep("wald", 4)
#   # dimnames(V) <- list(names(fit$par), names(fit$par))
#   
#   # obj <- structure(
#   #   list(
#   #     coefficients = fit$par,
#   #     vcov = V
#   #   ),
#   #   class = "fastmle"
#   # )
#   if(calc.ci){
#     obj <- structure(list(coefficients = est, vcov = V), class = "fastmle")
#     oa  <- list(method = "L-BFGS-B", lower = 1e-3,
#                 control = list(maxit = 100, factr = 1e7, pgtol = 1e-4))
#     ci <- profileCI::profileCI(obj,loglik,faster=TRUE,mult=32,epsilon=0,optim_args=oa)
#     for (j in 1:4) {
#       args <- list(obj, loglik = loglik, parm = nm[j], faster = TRUE,
#                    mult = 32, epsilon = 0, optim_args = oa)
#       if (est[j] > lb0[j]) args$lb <- lb0[j]
#       pj <- tryCatch(suppressWarnings(do.call(profileCI::profileCI, args)),
#                      error = function(e) NULL)
#       if (is.null(pj)) next
#       r <- unclass(pj)[1, 1:2]
#       if (nm[j] == "obsStd" && is.na(r[1]) && is.finite(r[2])) r[1] <- 0  # hit boundary
#       if (!anyNA(r)) { ci[j, ] <- r; meth[j] <- "profile" }
#     }
#     # ci <- profileCI::profileCI(
#     #   obj,
#     #   loglik = loglik,
#     #   faster = TRUE,
#     #   mult = 100,epsilon=0,
#     #   optim_args = list(
#     #     method = "L-BFGS-B",lower=c(0),
#     #     control = list(
#     #       maxit = 100,
#     #       factr = 1e7,
#     #       pgtol = 1e-4
#     #     )
#     #   )
#     
#     # )
#   } else {
#     ci <- NULL
#   }
#   
#   return(cbind(fit$par,ci))
# }

estimation_ci_split <- function(observation, nObs = length(observation),
                                initial_values = NULL, level = 0.95,profile_timeout=3600,
                                nSubsteps = 1L){
  par_names <- c("C1", "mu", "sigma", "obsStd")
  vi <- unname(sapply(par_names, function(i) which(param_strings %in% i)))
  if (is.null(initial_values)) initial_values <- initialGuesses[par_names]
  
  base <- usual_param_value
  base[[1]] <- observation
  base[[5]] <- nObs

  nll4 <- create_optimizer_function(base, vi, nSubsteps = nSubsteps)
  fit  <- optim(initial_values, nll4, method = "L-BFGS-B",
                lower = c(10, 10, 10, 0), control = opt_control)
  est  <- setNames(fit$par, c("C1", "mu", "sigma", "R"))
  f0   <- fit$value
  R_hat <- unname(est["R"])
  
  full3 <- base
  full3[[4]] <- R_hat          
  ci3 <- tryCatch(
    unclass(estimation_profile_ci(c("C", "mu", "sigma"), observation,
                                  initial_values = est[c("C1", "mu", "sigma")],
                                  nObs = nObs, full_param_list = full3,
                                  profile_timeout=profile_timeout,
                                  nSubsteps = nSubsteps))[, 1:3, drop = FALSE],
    error = function(e) { warning("3-param profile CI failed: ", conditionMessage(e))
      cbind(est[1:3], NA, NA) })
  
 
  rr <- ci_R_profile(nll4, est, f0, level)
  R_ci <- rr$ci; R_method <- rr$method
  
  out <- rbind(ci3, R = c(R_hat, R_ci))
  rownames(out) <- c("C", "mu", "sigma", "R")
  colnames(out) <- c("est", "2.5%", "97.5%")
  # attr(out, "R_method") <- R_method
  out
}


fast_wald_obsd <- function(observation, nObs = length(observation),
                           full_param_list = NULL, start = NULL, level = 0.95){
  par_names <- c("C1", "mu", "sigma", "obsStd")
  variable_indices <- unname(sapply(par_names, function(i) which(param_strings %in% i)))
  if (is.null(full_param_list)) {
    full_param_list <- usual_param_value
    full_param_list[[1]] <- observation
    full_param_list[[5]] <- nObs
  }
  nll_raw <- create_optimizer_function(full_param_list, variable_indices)
  nll <- function(lp) { v <- nll_raw(exp(lp)); if (is.finite(v)) v else 1e10 }
  
  if (is.null(start)) start <- log(unname(initialGuesses[par_names]))  
  fit <- optim(start, nll, method = "BFGS",
               control = list(maxit = 500, reltol = 1e-10))
  
  H  <- optimHess(fit$par, nll, control = list(ndeps = rep(1e-4, 4)))   # ca 40
  H  <- (H + t(H)) / 2
  ev <- eigen(H, symmetric = TRUE)
  pd <- all(ev$values > 1e-8)
  V  <- ev$vectors %*% diag(1 / pmax(ev$values, 1e-8)) %*% t(ev$vectors)
  se <- sqrt(diag(V))
  
  z  <- qnorm(1 - (1 - level) / 2)
  out <- cbind(est = exp(fit$par),
               lower = exp(fit$par - z * se),
               upper = exp(fit$par + z * se))
  rownames(out) <- c("C", "mu", "sigma", "R")  
  attr(out, "hessian_pd") <- pd
  attr(out, "par_log") <- fit$par                  
  out
}

ci_R_profile <- function(nll4, est, f0, level = 0.95){
  R_hat  <- unname(est[4]); th_hat <- unname(est[1:3])
  cut    <- qchisq(level, 1) / 2
  z      <- qnorm(1 - (1 - level) / 2)
  
  
  H  <- optimHess(unname(est), nll4,
                  control = list(parscale = pmax(abs(unname(est)), 1e-2), ndeps = rep(1e-4, 4)))
  H  <- (H + t(H)) / 2
  Dm <- diag(pmax(abs(unname(est)), 1e-2))
  ev <- eigen(Dm %*% H %*% Dm, symmetric = TRUE)
  V  <- Dm %*% (ev$vectors %*% diag(1 / pmax(ev$values, 1e-8)) %*% t(ev$vectors)) %*% Dm
  se_R <- sqrt(V[4, 4])
  wald <- c(max(R_hat - z * se_R, 0), R_hat + z * se_R)
  
  prof <- function(R) {
    o <- tryCatch(optim(th_hat, function(th) nll4(c(th, R)), method = "L-BFGS-B",
                        lower = c(10, 10, 10),
                        control = list(maxit = 100, factr = 1e7, pgtol = 1e-4)),
                  error = function(e) NULL)
    if (is.null(o) || !is.finite(o$value)) NA_real_ else o$value
  }
  g   <- function(R) { v <- prof(R); if (is.na(v)) NA_real_ else v - f0 - cut }
  tol <- 1e-3 * max(R_hat, 0.05)
  root <- function(lo, hi) tryCatch(uniroot(g, c(lo, hi), tol = tol)$root,
                                    error = function(e) NA_real_)
  
  hi <- max(wald[2], R_hat + 0.05); up <- NA_real_
  for (k in 1:8) {
    gv <- g(hi)
    if (is.na(gv)) break
    if (gv > 0) { up <- root(R_hat, hi); break }
    hi <- R_hat + 2 * (hi - R_hat)
  }
  if (R_hat <= 0) { lo <- 0 } else {
    g0 <- g(0)
    lo <- if (is.na(g0)) NA_real_ else if (g0 <= 0) 0 else root(0, R_hat)
  }
  
  ci <- c(lo, up); method <- "profile"
  if (anyNA(ci)) { ci[is.na(ci)] <- wald[is.na(ci)]; method <- "profile+wald" }
  list(ci = ci, method = method, se_joint = se_R)
}

estimation_profile_ci <- function(unknown_parameters, observation,
                                  initial_values = NULL,
                                  nObs = length(observation),
                                  full_param_list = NULL,
                                  ...,
                                  calc.ci = TRUE,
                                  profile_timeout = 3600,
                                  setup_file = "variableSetup.R",
                                  nSubsteps = 1L){#file.path(getwd(), "variableSetup.R")) {
  
  cPlace <- which("C" == unknown_parameters)
  if (length(cPlace) > 0) unknown_parameters[cPlace] <- "C1"
  
  variable_indices <- unname(
    sapply(unknown_parameters, function(i) which(param_strings %in% i))
  )
  
  if (is.null(initial_values)) {
    initial_values <- initialGuesses[unknown_parameters]
  }
  
  if (is.null(full_param_list)) {
    full_param_list <- usual_param_value
    full_param_list[[1]] <- observation
    full_param_list[[5]] <- nObs
    
    stdString <- which(param_strings == "obsStd")
    stdEstimated <- which(variable_indices == stdString)
    
    full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  }
  
  nll_NULL <- create_optimizer_function(
    full_param_list,
    variable_indices,
    nSubsteps = nSubsteps
  )
  
  fit <- optim(
    par = initial_values,
    fn = nll_NULL,
    method = "L-BFGS-B",
    lower = c(10, 10, 10),
    control = opt_control,
    hessian = TRUE
  )
  
  V <- solve(fit$hessian)
  dimnames(V) <- list(names(fit$par), names(fit$par))
  
  obj <- structure(
    list(
      coefficients = fit$par,
      vcov = V
    ),
    class = "fastmle"
  )
  
  if (calc.ci) {
    
    profile_result <- run_profile_ci_timeout(
      obj = obj,
      full_param_list = full_param_list,
      variable_indices = variable_indices,
      timeout = profile_timeout,
      setup_file = setup_file,
      nSubsteps = nSubsteps,
      optim_args = list(
        method = "L-BFGS-B",
        lower = c(10, 10, 10),
        control = opt_control
      )
    )
    
    if (profile_result$timed_out || profile_result$failed) {
      if (profile_result$timed_out) {
        warning(
          sprintf(
            "Profile CI exceeded %.1f minutes. Using Wald CI instead.",
            profile_timeout / 60
          ),
          call. = FALSE
        )
      } else {
        warning(
          paste(
            "Profile CI failed:",
            profile_result$error,
            "Using Wald CI instead."
          ),
          call. = FALSE
        )
      }
      
      # Wald backup
      z <- qnorm(0.975)
      se <- sqrt(pmax(diag(V), 0))
      
      ci <- cbind(
        `2.5%` = pmax(fit$par - z * se, 10),
        `97.5%` = fit$par + z * se
      )
      
      rownames(ci) <- names(fit$par)
      
    } else {
      
      ci <- profile_result$value
      
    }
    
  } else {
    ci <- NULL
  }
  
  return(cbind(fit$par, ci))
}

run_profile_ci_timeout <- function(obj, full_param_list,
                                   variable_indices,
                                   timeout = 3600,
                                   optim_args,
                                   setup_file = file.path(getwd(), "variableSetup.R")) {
  
  if (!requireNamespace("callr", quietly = TRUE)) {
    stop("Package 'callr' is required.")
  }
  
  setup_file <- normalizePath(setup_file, mustWork = TRUE)
  setup_dir <- dirname(setup_file)
  
  required_files <- c(
    "paramPrepper.R",
    "ekf_likelihood.R",
    "JansenAndRit.cpp",
    "WienerIntegrator.cpp",
    "EulerMaruyama.cpp"
  )
  
  missing_files <- required_files[
    !file.exists(file.path(setup_dir, required_files))
  ]
  
  if (length(missing_files)) {
    stop(
      "Missing project files in ", setup_dir, ": ",
      paste(missing_files, collapse = ", ")
    )
  }
  
  p <- callr::r_bg(
    
    func = function(obj, full_param_list, variable_indices,
                    optim_args, setup_file) {
      
      setup_dir <- dirname(setup_file)
      setwd(setup_dir)
      
      setup_code <- readLines(setup_file, warn = FALSE)
      
      setup_code <- setup_code[
        !grepl(
          "rstudioapi::getActiveDocumentContext",
          setup_code,
          fixed = TRUE
        )
      ]
      
      eval(parse(text = setup_code), envir = .GlobalEnv)
      
      nll_NULL <- create_optimizer_function(
        full_param_list,
        variable_indices
      )
      
      loglik <- function(pars) -nll_NULL(pars)
      
      profileCI::profileCI(
        obj,
        loglik = loglik,
        faster = TRUE,
        mult = 32,
        epsilon = 0,
        optim_args = optim_args
      )
    },
    
    args = list(
      obj = obj,
      full_param_list = full_param_list,
      variable_indices = variable_indices,
      optim_args = optim_args,
      setup_file = setup_file
    ),
    
    supervise = TRUE
  )
  
  on.exit({
    if (p$is_alive()) {
      try(p$kill(), silent = TRUE)
      try(p$wait(timeout = 5000), silent = TRUE)
    }
  }, add = TRUE)
  
  p$wait(timeout = timeout * 1000)
  
  if (p$is_alive()) {
    
    p$kill()
    p$wait(timeout = 5000)
    
    return(list(
      timed_out = TRUE,
      value = NULL
    ))
  }
  
  list(
    timed_out = FALSE,
    value = p$get_result()
  )
}


estimation_profile_ci_c <- function(unknown_parameters,observation,
                                  initial_values=initialGuesses[unknown_parameters],
                                  nObs=length(observation),full_param_list=NULL,...){
  cPlace <- which("C"==unknown_parameters)
  if(length(cPlace)>0)unknown_parameters[cPlace] <- "C1"
  variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
  
  if(is.null(full_param_list)){
    full_param_list <- usual_param_value
    full_param_list[[1]] <- observation
    full_param_list[[5]] <- nObs
  }
  
  stdString <- which(param_strings=="obsStd")
  stdEstimated <- which(variable_indices==stdString)
  
  full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  
  nll_NULL <- create_optimizer_function(full_param_list,variable_indices)
  
  f <- function(C1){
    nll_NULL(C1)
  }
  
  
  profile_fit <- function(){
    fit <- mle2(f,as.list(log(initial_values)),method="BFGS",
                control = list(maxit = 5000, reltol = 1e-12))
    
    ci <- try({
      prof <- profile(fit,maxsteps=1000)
      confint(prof,test="LRT")
    },silent=TRUE)
    
    if (inherits(ci, "try-error")) {
      H <- numDeriv::hessian(nll_NULL,fit@coef)
      se <- sqrt(diag(solve(H)))
      inherits(fit,"mle2")
      ci <- try({
        mapply(function(i)fit@coef+i*se,qnorm(c(.025,.975)))
      }, silent = TRUE)
      
    }
    if (inherits(ci, "try-error")) ci <- cbind(rep(-Inf,2),Inf)
    
    res_mat <- t(cbind(c(fit@coef,ci)))
    colnames(res_mat) <- c("Est","LRT_Lower","LRT_Upper")
    return(res_mat)
  }
  
  return(exp(profile_fit()))
}




get_hessian <- function(full_param_list,variable_indices,estimated_parameters,...){
  toOptim <- function(pars){
    for(idx in 1:length(variable_indices)){
      full_param_list[[variable_indices[idx]]] <- pars[idx]
    }
    if(any(c(9,11,13)%in%variable_indices)){
      dampingMat <- diag(c(full_param_list[[9]],full_param_list[[9]],full_param_list[[11]]))
      sigMat <- diag(c(10,full_param_list[[13]],10))
      full_param_list[[2]] <- expMt(dampingMat,full_param_list[[6]])
      full_param_list[[3]] <- variableSigmaCovarMat(full_param_list[[6]],sigMat,dampingMat)
    }
    full_param_list <- full_param_list[-13]
    # ekf <- do.call(ekf_ou,full_param_list)
    # -get_ekf_likelihood(ekf$mk,ekf$Pk,full_param_list[[1]],full_param_list[[5]],full_param_list[[4]])
    do.call(construct_ekf_likelihood,full_param_list)
  }
  
  numDeriv::hessian(toOptim,estimated_parameters)
}


# parameter_inference <- function(unknown_parameters,observation,stepMat,covMat,obsStd,nObs,Delta,mu,A,a,B,b,C1,...){
#   full_param_list <- list(observation,stepMat,covMat,obsStd^2,nObs,Delta,mu,A,a,B,b,C1)
#   variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
#   initials <- sapply(variable_indices,function(i)full_param_list[[i]])
#   res <- optimize_parameters(full_param_list,
#                              variable_indices,initials,...)
#   pars <- exp(res$par)
#   names(pars) <- unknown_parameters
#   return(c(pars,"Value"=res$value))
# }




# parameter_inference("C1",observation = observation,stepMat = x1StepMat.,covMat = covMatOU.,obsStd = obsErr,
#                     nObs = nObs,Delta = Delta.,mu = mu.,A = A,a = a,B = B,b = b,
#                     C1 = C1_guess,lower=120,upper=180)
# 
# 
# parameter_inference(c("A","B","mu"),observation = observation,stepMat = x1StepMat.,covMat = covMatOU.,obsStd = obsErr,
#                     nObs = nObs,Delta = Delta.,mu = mu_guess,A = A_guess,a = a,B = B_guess,b = b,
#                     C1 = C1)


jr_inference <- function(unknown_parameters,observation,initial_values=initialGuesses[unknown_parameters],
                         nObs=length(observation),full_param_list=NULL,...){
  cPlace <- which("C"==unknown_parameters)
  if(length(cPlace)>0)unknown_parameters[cPlace] <- "C1"
  variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
  
  if(is.null(full_param_list)){
    full_param_list <- usual_param_value
    full_param_list[[1]] <- observation
    full_param_list[[5]] <- nObs
  }
  
  stdString <- which(param_strings=="obsStd")
  stdEstimated <- which(variable_indices==stdString)
  
  full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  
  res <- optimize_parameters(full_param_list,
                             variable_indices,initial_values,method="BFGS",...)
  pars <- exp(res$par)
  if(length(stdEstimated)>0) pars[stdEstimated] <- sqrt(pars[stdEstimated])
  names(pars) <- unknown_parameters
  return(c(pars))
}



jr_inference_with_ci <- function(unknown_parameters,observation,initial_values=initialGuesses[unknown_parameters],
                         nObs=length(observation),full_param_list=NULL,method="BFGS",...){
  cPlace <- which("C"==unknown_parameters)
  if(length(cPlace)>0)unknown_parameters[cPlace] <- "C1"
  variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
  
  if(is.null(full_param_list)){
    full_param_list <- usual_param_value
    full_param_list[[1]] <- observation
    full_param_list[[5]] <- nObs
  }
  
  stdString <- which(param_strings=="obsStd")
  stdEstimated <- which(variable_indices==stdString)
  
  full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  
  res <- optimize_parameters(full_param_list,
                             variable_indices,initial_values,method=method,...)
  pars <- exp(res$par)
  if(length(stdEstimated)>0) pars[stdEstimated] <- sqrt(pars[stdEstimated])
  names(pars) <- unknown_parameters
  
  hess <- get_hessian(full_param_list,variable_indices,pars)
  #Symmetrise just to go sure
  hess <- 0.5 * (hess + t(hess))
  vcov <- tryCatch(
    solve(hess),
    error = function(e) {
      if (!use_pinv) stop("Hessian is not invertible.")
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' is required for generalized inverse.")
      }
      MASS::ginv(hess)
    }
  )
  
  se <- sqrt(pmax(diag(vcov), 0))
  
  #marginal wald
  z <- qnorm(1 - (1 - .95) / 2)
  ci <- cbind(
    estimate = pars,
    se = se,
    lower = pars - z * se,
    upper = pars + z * se
  )
  
  return(ci)
}




wiener_optimizer <- function(full_param_list,variable_indices,...){
  toOptim <- function(pars){
    for(idx in 1:length(variable_indices)){
      full_param_list[[variable_indices[idx]]] <- pars[idx]
    }
    # if(any(c(9,11,13)%in%variable_indices)){
    # dampingMat <- diag(c(full_param_list[[9]],full_param_list[[9]],full_param_list[[11]]))
    # full_param_list[[2]] <- expMt(dampingMat,full_param_list[[6]]/2)
    full_param_list[[3]] <- WIEcovMat(full_param_list[[13]])
    # }
    full_param_list <- full_param_list[-13]
    do.call(wi_ekf_nll,full_param_list)
  }
  return(toOptim)
}

wiener_inference_plus_ci <-   function(observation,full_param_list=NULL,...,calc.ci=TRUE){
  unknown_parameters <- c("C","mu","sigma")
  cPlace <- which("C"==unknown_parameters)
  if(length(cPlace)>0)unknown_parameters[cPlace] <- "C1"
  initial_values=initialGuesses[unknown_parameters]
  variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
  
  if(is.null(full_param_list)){
    full_param_list <- usual_param_value
    full_param_list[[2]] <- stepMatHalf.
    full_param_list[[1]] <- observation
    full_param_list[[3]] <- WIEcovMat(full_param_list[[13]])
    full_param_list[[5]] <- length(observation)
  }
  
  stdString <- which(param_strings=="obsStd")
  stdEstimated <- which(variable_indices==stdString)
  
  full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  
  nll_NULL <- wiener_optimizer(full_param_list,variable_indices)
  loglik <- function(pars) - nll_NULL(pars)
  
  fit <- optim(
    par = initial_values,
    fn = nll_NULL,
    method = "L-BFGS-B",lower=c(10,10,10),
    control = opt_control,
    hessian=TRUE
  )
  
  V <- solve(fit$hessian)
  dimnames(V) <- list(names(fit$par), names(fit$par))
  
  obj <- structure(
    list(
      coefficients = fit$par,
      vcov = V
    ),
    class = "fastmle"
  )
  if(calc.ci){
    ci <- profileCI::profileCI(
      obj,
      loglik = loglik,
      faster = TRUE,
      mult = 32,epsilon=0,
      optim_args = list(
        method = "L-BFGS-B",lower=c(10,10,10),
        control = opt_control
      )
    )
  } else {
    ci <- NULL
  }
  
  return(cbind(fit$par,ci))

}



wiener_inference_plus_ci_c <-   function(observation,full_param_list=NULL,...){
  unknown_parameters <- c("C")
  cPlace <- which("C"==unknown_parameters)
  if(length(cPlace)>0)unknown_parameters[cPlace] <- "C1"
  initial_values=initialGuesses[unknown_parameters]
  variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
  
  if(is.null(full_param_list)){
    full_param_list <- usual_param_value
    full_param_list[[2]] <- stepMatHalf.
    full_param_list[[1]] <- observation
    full_param_list[[5]] <- length(observation)
  }
  
  stdString <- which(param_strings=="obsStd")
  stdEstimated <- which(variable_indices==stdString)
  
  full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  
  nll_NULL <- wiener_optimizer(full_param_list,variable_indices)
  
  f <- function(C1){
    nll_NULL(c(C1))
  }
  
  # fit <- optim(initial_values,nll,method="BFGS",...)
  
  profile_fit <- function(){
    fit <- mle2(f,as.list(log(initial_values)),method="BFGS",
                control = list(maxit = 5000, reltol = 1e-12))
    
    ci <- try({
      prof <- profile(fit,maxsteps=1000)
      while(class(prof)[1]=="mle2"){
        prof <- profile(prof,maxsteps=1000)
      }
      confint(prof,test="LRT")
    },silent=TRUE)
    
    
    if (inherits(ci, "try-error")) {
      H <- numDeriv::hessian(nll_NULL,fit@coef)
      se <- sqrt(diag(solve(H)))
      ci <- try({
        mapply(function(i)fit@coef+i*se,qnorm(c(.025,.975)))
      }, silent = TRUE)
    }
    if (inherits(ci, "try-error")){
     warning("Could not compute confidence interval") 
     ci <- cbind(rep(-Inf,2),Inf)
    }
    res_mat <- t(cbind(c(fit@coef,ci)))
    colnames(res_mat) <- c("Est","LRT_Lower","LRT_Upper")
    return(res_mat)
  }
  
  return(exp(profile_fit()))
}


em_optimizer <- function(full_param_list,variable_indices,...){
  toOptim <- function(pars){
    for(idx in 1:length(variable_indices)){
      full_param_list[[variable_indices[idx]]] <- pars[idx]
    }
    # if(any(c(9,11,13)%in%variable_indices)){
    # dampingMat <- diag(c(full_param_list[[9]],full_param_list[[9]],full_param_list[[11]]))
    # full_param_list[[2]] <- expMt(dampingMat,full_param_list[[6]]/2)
    full_param_list[[3]] <- diag(c(rep(0,3),tau.,full_param_list[[13]],tau.)^2)*Delta.
    # }
    full_param_list <- full_param_list[-13]
    do.call(em_ekf_nll,full_param_list)
  }
  return(toOptim)
}
em_inference_plus_ci <-   function(observation,full_param_list=NULL,...,calc.ci=TRUE){
  unknown_parameters <- c("C","mu","sigma")
  cPlace <- which("C"==unknown_parameters)
  if(length(cPlace)>0)unknown_parameters[cPlace] <- "C1"
  initial_values=initialGuesses[unknown_parameters]
  variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
  
  if(is.null(full_param_list)){
    full_param_list <- usual_param_value
    full_param_list[[1]] <- observation
    full_param_list[[5]] <- length(observation)
  }
  
  stdString <- which(param_strings=="obsStd")
  stdEstimated <- which(variable_indices==stdString)
  
  full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  
  nll_NULL <- em_optimizer(full_param_list,variable_indices)
  loglik <- function(pars) - nll_NULL(pars)
  
  fit <- optim(
    par = initial_values,
    fn = nll_NULL,
    method = "L-BFGS-B",lower=c(10,10,10),
    control = opt_control,
    hessian=TRUE
  )
  
  V <- solve(fit$hessian)
  dimnames(V) <- list(names(fit$par), names(fit$par))
  
  obj <- structure(
    list(
      coefficients = fit$par,
      vcov = V
    ),
    class = "fastmle"
  )
  if(calc.ci){
    ci <- profileCI::profileCI(
      obj,
      loglik = loglik,
      faster = TRUE,
      mult = 32,
      optim_args = list(
        method = "L-BFGS-B",lower=c(10,10,10),
        control = opt_control
      )
    )
  } else {
    ci <- NULL
  }
  
  return(cbind(fit$par,ci))
}

em_inference_plus_ci_c <-   function(observation,full_param_list=NULL,...){
  unknown_parameters <- c("C")
  cPlace <- which("C"==unknown_parameters)
  if(length(cPlace)>0)unknown_parameters[cPlace] <- "C1"
  initial_values=initialGuesses[unknown_parameters]
  variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
  
  if(is.null(full_param_list)){
    full_param_list <- usual_param_value
    full_param_list[[1]] <- observation
    full_param_list[[5]] <- length(observation)
  }
  
  stdString <- which(param_strings=="obsStd")
  stdEstimated <- which(variable_indices==stdString)
  
  full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  
  nll_NULL <- em_optimizer(full_param_list,variable_indices)
  
  f <- function(C1){
    nll_NULL(c(C1))
  }
  
  # fit <- optim(initial_values,nll,method="BFGS",...)
  
  profile_fit <- function(){
    fit <- mle2(f,as.list(log(initial_values)),method="BFGS",
                control = list(maxit = 5000, reltol = 1e-12))
    
    ci <- try({
      prof <- profile(fit,maxsteps=1000)
      while(class(prof)[1]=="mle2"){
        prof <- profile(prof,maxsteps=1000)
      }
      confint(prof,test="LRT")
    },silent=TRUE)
    
    
    if (inherits(ci, "try-error")) {
      H <- numDeriv::hessian(nll_NULL,fit@coef)
      se <- sqrt(diag(solve(H)))
      ci <- try({
        mapply(function(i)fit@coef+i*se,qnorm(c(.025,.975)))
      }, silent = TRUE)
    }
    if (inherits(ci, "try-error")){
     warning("Could not compute confidence interval") 
     ci <- cbind(rep(-Inf,2),Inf)
    }
    res_mat <- t(cbind(c(fit@coef,ci)))
    colnames(res_mat) <- c("Est","LRT_Lower","LRT_Upper")
    return(res_mat)
  }
  
  return(exp(profile_fit()))
}



em_inference_plus_ci_nosig <-   function(observation,unknown_parameters=c("C","mu"),full_param_list=NULL,...){
  em_coptimizer <- function(full_param_list,variable_indices,...){
    toOptim <- function(pars){
      pars <- exp(pars)
      for(idx in 1:length(variable_indices)){
        full_param_list[[variable_indices[idx]]] <- pars[idx]
      }
      # if(any(c(9,11,13)%in%variable_indices)){
      # dampingMat <- diag(c(full_param_list[[9]],full_param_list[[9]],full_param_list[[11]]))
      # full_param_list[[2]] <- expMt(dampingMat,full_param_list[[6]]/2)
      # full_param_list[[3]] <- diag(c(rep(0,3),tau.,full_param_list[[13]],tau.))*Delta.
      # }
      full_param_list <- full_param_list[-13]
      do.call(em_ekf_nll,full_param_list)
    }
    return(toOptim)
  }
  cPlace <- which("C"==unknown_parameters)
  if(length(cPlace)>0)unknown_parameters[cPlace] <- "C1"
  initial_values=initialGuesses[unknown_parameters]
  variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
  
  if(is.null(full_param_list)){
    full_param_list <- usual_param_value
    full_param_list[[1]] <- observation
    full_param_list[[5]] <- length(observation)
  }
  
  stdString <- which(param_strings=="obsStd")
  stdEstimated <- which(variable_indices==stdString)
  
  full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  
  nll_NULL <- em_coptimizer(full_param_list,variable_indices)
  
  f <- function(C1,mu){
    nll_NULL(c(C1,mu))
  }
  
  # fit <- optim(initial_values,nll,method="BFGS",...)
  
  profile_fit <- function(){
    fit <- mle2(f,as.list(log(initial_values)),method="BFGS",
                control = list(maxit = 5000, reltol = 1e-12))
    
    ci <- try({
      prof <- profile(fit,maxsteps=1000)
      while(class(prof)[1]=="mle2"){
        prof <- profile(prof,maxsteps=1000)
      }
      confint(prof,test="LRT")
    },silent=TRUE)
    
    
    if (inherits(ci, "try-error")) {
      H <- numDeriv::hessian(nll_NULL,fit@coef)
      se <- sqrt(diag(solve(H)))
      ci <- try({
        mapply(function(i)fit@coef+i*se,qnorm(c(.025,.975)))
      }, silent = TRUE)
    }
    if (inherits(ci, "try-error")){
      warning("Could not compute confidence interval") 
      ci <- cbind(rep(-Inf,2),Inf)
    }
    res_mat <- if(length(unknown_parameters)==1) matrix(c(fit@coef,ci),1) else cbind(fit@coef,ci)
    # colnames(res_mat) <- c("Est","LRT_Lower","LRT_Upper")
    return(res_mat)
  }
  
  return(exp(profile_fit()))
}





# 
# 
# 
# jr_inference_with_ci <- function(unknown_parameters,observation,initial_values=initialGuesses[unknown_parameters],
#                                  nObs=length(observation),full_param_list=NULL,...){
#   cPlace <- which("C"==unknown_parameters)
#   if(length(cPlace)>0)unknown_parameters[cPlace] <- "C1"
#   variable_indices <- unname(sapply(unknown_parameters,function(i)which(param_strings%in%i)))
#   
#   if(is.null(full_param_list)){
#     full_param_list <- usual_param_value
#     full_param_list[[1]] <- observation
#     full_param_list[[5]] <- nObs
#   }
#   
#   stdString <- which(param_strings=="obsStd")
#   stdEstimated <- which(variable_indices==stdString)
#   
#   full_param_list[[stdString]] <- full_param_list[[stdString]]^2
#   
#   res <- optimize_parameters(full_param_list,
#                              variable_indices,initial_values,method="BFGS",hessian=TRUE,...)
#   pars <- exp(res$par)
#   if(length(stdEstimated)>0) pars[stdEstimated] <- sqrt(pars[stdEstimated])
#   names(pars) <- unknown_parameters
#   
#   hess <- get_hessian(full_param_list,variable_indices,pars)
#   #Symmetrise just to go sure
#   hess <- 0.5 * (hess + t(hess))
#   vcov <- tryCatch(
#     solve(hess),
#     error = function(e) {
#       if (!use_pinv) stop("Hessian is not invertible.")
#       if (!requireNamespace("MASS", quietly = TRUE)) {
#         stop("Package 'MASS' is required for generalized inverse.")
#       }
#       MASS::ginv(hess)
#     }
#   )
#   
#   #marginal wald
#   se <- sqrt(pmax(diag(vcov), 0))
#   
#   z <- qnorm(1 - (1 - .95) / 2)
#   
#   
#   
#   
#   #WALD ci log scale
#   log_se <- sqrt(pmax(diag(res$hessian),0))
#   
#   
#   
#   ci <- cbind(
#     estimate = pars,
#     se = se,
#     lower_wald = pars - z * se,
#     upper_wald = pars + z * se,
#     lower_logwald = exp(log(pars)-z*log_se),
#     upper_logwald = exp(log(pars)+z*log_se),
#     lower_profile
#   )
#   
#   return(ci)
# }


# nsim <- 500
# sigmas <- numeric(nsim)
# for(i in 1:nsim){
# observation <- simulate_jr()
# nObs <- length(observation)
# sigmas[i] <- round(jr_inference("sigma",800,observation,lower=50,upper=2000),2)
# }
# plot(density(sigmas))
# abline(v=sigma.,col="red",lty=2,lwd=3)
# 
# observation <- simulate_jr()
# nObs <- length(observation)
# observation <- observation+rnorm(nObs,sd=obsErr)
# 
# 
# round(jr_inference("sigma",800,observation,lower=50,upper=2000),2)
# sigma.
# 
# 
# 
# round(jr_inference(c("nu0","nuMax"),c(3,7),observation),2)
# c(nu0,nuMax)
# 
# 
# jr_inference("C1",120,observation,lower=50,upper=1000)
# C1
# 
# jr_inference(c("A","B","mu"),c(4,30,200),observation)
# c(A,B,mu.)
# 
# jr_inference(c("A","C1"),c(4,200),observation)
# c(A,C1)
# 
# jr_inference(c("A","mu","sigma","b"),c(4,200,800,40),observation)
# c(A,mu.,sigma.,b)
# 
# jr_inference(c("C1","mu","sigma"),c(100,200,800),observation)
# c(C1,mu.,sigma.)
# 
# jr_inference(c("C1","mu","sigma","A","b"),c(100,200,800,4,40),observation)
# c(C1,mu.,sigma.,A,b)
# 
# 
# 
# # jr_inference(c("A","B","C1"),c(4,30,200),observation) #This is hard af
# 
# 
# TODO MAKE EVEN EVEN EVEN SMOOTHER BY USING OOP. BASICALLY YOU CAN CREATE SOME
# var = VariableForInference("A",4,1,10)
# but also
# var2 = VariableForInference("C1",100)
# where then -Inf and Inf are defaults
# jr_inference(InfVar("A",4,1,10),InfVar("C1",200),InfVar("Sigma",2,1e-10))





run_profile_ci_timeout <- function(obj, full_param_list, variable_indices,
                                   timeout = 3600,
                                   optim_args,
                                   setup_file = file.path(getwd(), "variableSetup.R"),
                                   optimizer_type = c("ekf", "em", "wiener")) {
  
  optimizer_type <- match.arg(optimizer_type)
  
  if (!requireNamespace("callr", quietly = TRUE)) {
    stop("Package 'callr' is required.")
  }
  
  setup_file <- normalizePath(setup_file, mustWork = TRUE)
  setup_dir <- dirname(setup_file)
  
  required_files <- c(
    "paramPrepper.R",
    "ekf_likelihood.R",
    "JansenAndRit.cpp",
    "WienerIntegrator.cpp",
    "EulerMaruyama.cpp"
  )
  
  missing_files <- required_files[
    !file.exists(file.path(setup_dir, required_files))
  ]
  
  if (length(missing_files) > 0) {
    stop(
      "Missing project files: ",
      paste(missing_files, collapse = ", ")
    )
  }
  
  p <- callr::r_bg(
    
    func = function(obj, full_param_list, variable_indices,
                    optim_args, setup_file, optimizer_type) {
      
      setup_dir <- dirname(setup_file)
      setwd(setup_dir)
      
      setup_code <- readLines(setup_file, warn = FALSE)
      
      setup_code <- setup_code[
        !grepl(
          "rstudioapi::getActiveDocumentContext",
          setup_code,
          fixed = TRUE
        )
      ]
      
      eval(parse(text = setup_code), envir = .GlobalEnv)
      
      optimizer_builder <- switch(
        optimizer_type,
        ekf = create_optimizer_function,
        em = em_optimizer,
        wiener = wiener_optimizer
      )
      
      nll_NULL <- optimizer_builder(
        full_param_list,
        variable_indices
      )
      
      loglik <- function(pars) -nll_NULL(pars)
      
      profileCI::profileCI(
        obj,
        loglik = loglik,
        faster = TRUE,
        mult = 32,
        epsilon = 0,
        optim_args = optim_args
      )
    },
    
    args = list(
      obj = obj,
      full_param_list = full_param_list,
      variable_indices = variable_indices,
      optim_args = optim_args,
      setup_file = setup_file,
      optimizer_type = optimizer_type
    ),
    
    supervise = TRUE
  )
  
  on.exit({
    if (p$is_alive()) {
      try(p$kill(), silent = TRUE)
      try(p$wait(timeout = 5000), silent = TRUE)
    }
  }, add = TRUE)
  
  p$wait(timeout = timeout * 1000)
  
  if (p$is_alive()) {
    
    p$kill()
    p$wait(timeout = 5000)
    
    return(list(
      timed_out = TRUE,
      value = NULL
    ))
  }
  
  list(
    timed_out = FALSE,
    value = p$get_result()
  )
}


em_inference_plus_ci <- function(observation, full_param_list = NULL, ...,
                                 calc.ci = TRUE,
                                 profile_timeout = 3600,
                                 setup_file = file.path(getwd(), "variableSetup.R")) {
  
  unknown_parameters <- c("C1", "mu", "sigma")
  
  initial_values <- initialGuesses[unknown_parameters]
  
  variable_indices <- unname(
    sapply(unknown_parameters, function(i) which(param_strings %in% i))
  )
  
  if (is.null(full_param_list)) {
    full_param_list <- usual_param_value
    full_param_list[[1]] <- observation
    full_param_list[[5]] <- length(observation)
  }
  
  stdString <- which(param_strings == "obsStd")
  full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  
  
  nll_NULL <- em_optimizer(full_param_list, variable_indices)
  
  fit <- optim(
    par = initial_values,
    fn = nll_NULL,
    method = "L-BFGS-B",
    lower = c(10, 10, 10),
    control = opt_control,
    hessian = TRUE
  )
  
  V <- solve(fit$hessian)
  
  dimnames(V) <- list(names(fit$par), names(fit$par))
  
  obj <- structure(
    list(
      coefficients = fit$par,
      vcov = V
    ),
    class = "fastmle"
  )
  
  
  if (calc.ci) {
    
    profile_result <- run_profile_ci_timeout(
      obj = obj,
      full_param_list = full_param_list,
      variable_indices = variable_indices,
      timeout = profile_timeout,
      setup_file = setup_file,
      optimizer_type = "em",
      optim_args = list(
        method = "L-BFGS-B",
        lower = c(10, 10, 10),
        control = opt_control
      )
    )
    
    if (profile_result$timed_out || profile_result$failed) {
      if (profile_result$timed_out) {
        warning(
          sprintf(
            "Profile CI exceeded %.1f minutes. Using Wald CI instead.",
            profile_timeout / 60
          ),
          call. = FALSE
        )
      } else {
        warning(
          paste(
            "Profile CI failed:",
            profile_result$error,
            "Using Wald CI instead."
          ),
          call. = FALSE
        )
      }
      
      
      z <- qnorm(0.975)
      se <- sqrt(pmax(diag(V), 0))
      
      ci <- cbind(
        `2.5%` = pmax(fit$par - z * se, 10),
        `97.5%` = fit$par + z * se
      )
      
      rownames(ci) <- names(fit$par)
      
    } else {
      
      ci <- profile_result$value
      
    }
    
  } else {
    ci <- NULL
  }
  
  return(cbind(fit$par, ci))
}

wiener_inference_plus_ci <- function(observation, full_param_list = NULL, ...,
                                     calc.ci = TRUE,
                                     profile_timeout = 3600,
                                     setup_file = file.path(getwd(), "variableSetup.R")) {
  
  unknown_parameters <- c("C1", "mu", "sigma")
  
  initial_values <- initialGuesses[unknown_parameters]
  
  variable_indices <- unname(
    sapply(unknown_parameters, function(i) which(param_strings %in% i))
  )
  
  if (is.null(full_param_list)) {
    full_param_list <- usual_param_value
    full_param_list[[2]] <- stepMatHalf.
    full_param_list[[1]] <- observation
    full_param_list[[3]] <- WIEcovMat(full_param_list[[13]])
    full_param_list[[5]] <- length(observation)
  }
  
  stdString <- which(param_strings == "obsStd")
  full_param_list[[stdString]] <- full_param_list[[stdString]]^2
  
  
  nll_NULL <- wiener_optimizer(full_param_list, variable_indices)
  
  fit <- optim(
    par = initial_values,
    fn = nll_NULL,
    method = "L-BFGS-B",
    lower = c(10, 10, 10),
    control = opt_control,
    hessian = TRUE
  )
  
  V <- solve(fit$hessian)
  
  dimnames(V) <- list(names(fit$par), names(fit$par))
  
  obj <- structure(
    list(
      coefficients = fit$par,
      vcov = V
    ),
    class = "fastmle"
  )
  
  if (calc.ci) {
    
    profile_result <- run_profile_ci_timeout(
      obj = obj,
      full_param_list = full_param_list,
      variable_indices = variable_indices,
      timeout = profile_timeout,
      setup_file = setup_file,
      optimizer_type = "wiener",
      optim_args = list(
        method = "L-BFGS-B",
        lower = c(10, 10, 10),
        control = opt_control
      )
    )
    
    if (profile_result$timed_out || profile_result$failed) {
      if (profile_result$timed_out) {
        warning(
          sprintf(
            "Profile CI exceeded %.1f minutes. Using Wald CI instead.",
            profile_timeout / 60
          ),
          call. = FALSE
        )
      } else {
        warning(
          paste(
            "Profile CI failed:",
            profile_result$error,
            "Using Wald CI instead."
          ),
          call. = FALSE
        )
      }
      
      
      z <- qnorm(0.975)
      se <- sqrt(pmax(diag(V), 0))
      
      ci <- cbind(
        `2.5%` = pmax(fit$par - z * se, 10),
        `97.5%` = fit$par + z * se
      )
      
      rownames(ci) <- names(fit$par)
      
    } else {
      ci <- profile_result$value
    }
  } else {
    ci <- NULL
  }
  return(cbind(fit$par, ci))
}


run_profile_ci_timeout <- function(obj, full_param_list, variable_indices,
                                   timeout = 3600,
                                   optim_args,
                                   setup_file = NULL,
                                   optimizer_type = c("ekf", "em", "wiener"),
                                   nSubsteps = 1L) {
  
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    stop("Package 'R.utils' is required.")
  }
  
  optimizer_type <- match.arg(optimizer_type)
  
  timed_out <- FALSE
  failed <- FALSE
  error_message <- NULL
  
  value <- tryCatch({
    
    optimizer_builder <- switch(
      optimizer_type,
      ekf = function(fpl, vi) create_optimizer_function(fpl, vi, nSubsteps = nSubsteps),
      em = em_optimizer,
      wiener = wiener_optimizer
    )
    
    nll_NULL <- optimizer_builder(
      full_param_list,
      variable_indices
    )
    
    loglik <- function(pars) -nll_NULL(pars)
    
    R.utils::withTimeout(
      expr = {
        profileCI::profileCI(
          obj,
          loglik = loglik,
          faster = TRUE,
          mult = 32,
          epsilon = 0,
          optim_args = optim_args
        )
      },
      timeout = timeout,
      cpu = Inf,
      elapsed = timeout,
      onTimeout = "error"
    )
    
  },
  TimeoutException = function(e) {
    timed_out <<- TRUE
    NULL
  },
  error = function(e) {
    failed <<- TRUE
    error_message <<- conditionMessage(e)
    NULL
  })
  
  list(
    timed_out = timed_out,
    failed = failed,
    error = error_message,
    value = value
  )
}
