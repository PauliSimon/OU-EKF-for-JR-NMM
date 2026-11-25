# OU-EKF-for-JR-NMM

Code accompanying *"Splitting-based extended Kalman filtering for parameter inference in the stochastic Jansen–Rit neural mass model"* (Pauli, Tubikanec, Futschik).

The code estimates parameters of the stochastic JR-NMM from a one-dimensional, EEG-like signal `Y = X2 − X3`.

| Filter    | Discretisation                    | Main function                                    |
|-----------|-----------------------------------|--------------------------------------------------|
| **OU-EKF**  | Ornstein–Uhlenbeck Strang splitting | `estimation_profile_ci()`                        |
| **WIE-EKF** | Wiener Strang splitting             | `wiener_inference_plus_ci()`                     |
| **EM-EKF**  | Euler–Maruyama                      | `em_inference_plus_ci()`                         |

In the paper, OU-EKF performs best, so it is the recommended default.

---

## 0. Setup



**Start by sourcing `variableSetup.R`.** Make sure you have all dependencies listed in the file installed.

```r
install.packages(c("Rcpp","RcppArmadillo","expm","data.table","mvnfast","bbmle","profileCI","R.utils"))
```


`variableSetup.R` loads and initiates all required functions and defines the global model constants and simulators used below.

```r
setwd("path/to/this/folder")   # variableSetup.R uses relative paths
source("variableSetup.R")
```

`variableSetup.R` compiles all of the C++ files, so the first run takes a moment. You need a working C++ toolchain.

### Default settings after sourcing

These are the globals that `variableSetup.R` defines (they match the paper's simulation study):

| Object | Value | Meaning |
|--------|-------|---------|
| `C.` / `C1` | 135 | connectivity parameter C |
| `mu.` | 220 | mean of external input |
| `sigma.` | 2000 | noise intensity σ |
| `tau.` | 1 | noise intensity τ on the other two velocity components (fixed) |
| `Delta.` | 0.002 | observation step size Δ |
| `maxT.` | 20 | default time horizon, used by `simulate_jr()` |
| `obsStd` / `obsErr` | 0 | observation noise (see the note on `obsStd` below) |

---

## 1. Simulate some data

The paper generates data with a fine Euler–Maruyama scheme (step Δ/100) and keeps every 100th point, so the data-generating process does not favour any of the filters. `simulate_em_baseline()` does exactly this:

```r
set.seed(1)

# Δ = 0.002, T = 20 s  ->  10,000 observations
y <- simulate_em_baseline(subFac = 100, Delta = 0.002, maxT = 20,
                          sigma = 2000, obsSd = 0)
length(y)
plot(y, type = "l")
```

`obsSd` is a **standard deviation**. To simulate noisy observations with variance `R = 0.75`:

```r
y_noisy <- simulate_em_baseline(subFac = 100, Delta = 0.002, maxT = 20,
                                sigma = 2000, obsSd = sqrt(0.75))
```

Other generators: You can use `type = "wi"` / `"ou"` in `simulate_em_baseline()` for Wiener/OU-based data.

---

## 2. Estimate C, µ, σ with OU-EKF

`estimation_profile_ci()` is the main entry point.

```r
res <- estimation_profile_ci(unknown_parameters = c("C", "mu", "sigma"),
                             observation        = y)
res
```

The result is a matrix with one row per parameter: the point estimate and the lower/upper bound of the confidence interval (95% by default). Row names are `C1`, `mu`, `sigma`.

```
        (estimate)   2.5%   97.5%
C1         ...
mu         ...
sigma      ...
```

True values for this simulation: C = 135, µ = 220, σ = 2000.

### Useful arguments

```r
res <- estimation_profile_ci(
  unknown_parameters = c("C", "mu", "sigma"),
  observation        = y,
  initial_values     = c(C1 = 140, mu = 200, sigma = 2200),  # default: initialGuesses[...]
  ci_method          = "sandwich",   # "sandwich" (default), "adjprofile", "wald", "profile"
  level              = 0.95,
  hac_lag            = NULL,         # Newey–West lag; default ceiling(N^(1/3))
  nSubsteps          = 1L,           # see Section 5
  n_extra            = 2             # extra random restarts of the optimiser
)
```

- **Estimation scale.** Optimisation is done on the log scale (BFGS), with the best of `1 + n_extra` starts (original start plus jittered ones) kept. Estimates are therefore always positive.
- **`ci_method`**
  - `"sandwich"`: Hessian + Newey–West (HAC) score covariance on the log scale. This is the method described in the paper (Section 2.3).
  - `"adjprofile"`: adjusted profile-likelihood intervals; slower, falls back to sandwich limits where the profile never crosses the cut-off (with a warning).
  - `"wald"`: model-based Wald only (inverse Hessian, no HAC correction).
  - `"profile"`: legacy profile intervals via the `profileCI` package (needs `profileCI` and `R.utils`; times out after `profile_timeout` seconds and then falls back to Wald).
- **Point estimates only** (no confidence intervals, faster): `calc.ci = FALSE` returns a one-column matrix of estimates.

### Extra diagnostics attached to the result

For the non-`profile` methods the result carries attributes:

```r
attr(res, "ci_method")
attr(res, "ci_all")                  # list(wald, sandwich, adjprofile) intervals
attr(res, "se_log_model")            # model-based s.e. on the log scale
attr(res, "se_log_sandwich")         # sandwich s.e. on the log scale
attr(res, "se_ratio")                # sandwich / model; far from 1 hints at misspecification
attr(res, "innovation_diagnostics")  # mean, var, lag-1 ACF of standardised innovations (and of their squares)
```

If the Hessian is not positive definite at the optimum you will get a warning; check convergence and identifiability in that case.

### Estimating other parameters

Any subset of `param_strings` works, e.g. the synaptic gain A and µ:

```r
estimation_profile_ci(c("A","mu"), y,
                      initial_values = c(A = 4, mu = 200))
```

Parameters not listed stay at the defaults from `usual_param_value` (Table 1 of the paper).

---

## 3. EM-EKF and WIE-EKF

The three parameters considered in the paper (C,mu,sigma) can be estimated with the Euler–Maruyama or Wiener-splitting filters. The interface is similar, but note the differences below.

```r
res_wie <- wiener_inference_plus_ci(observation = y)
res_em  <- em_inference_plus_ci(observation = y)
```

The defaults: `unknown_parameters = c("C", "mu", "sigma")`, `ci_method = "sandwich"`, `level = 0.95`, are not variable here, as EM-EKF and WIE-EKF were used only as comparison to OU-EKF. The returned matrix has the same layout (estimate, lower, upper) as for OU-EKF.

Comparing all three on the same data:

```r
rbind(
  OU  = estimation_profile_ci(c("C", "mu", "sigma"), y)[, 1],
  WIE = wiener_inference_plus_ci(y)[, 1],
  EM  = em_inference_plus_ci(y)[, 1]
)
```

Differences to `estimation_profile_ci()`:

- There is **no `nSubsteps` argument**; substepping is only implemented for OU-EKF.
- If the confidence-interval computation fails, you get the estimates with `NA` intervals and a warning (the error message is stored in `attr(res, "ci_error")`).
- You cannot estimate different parameter sets

---

## 4. Also estimating the observation noise (`obsStd`)

Add `"obsStd"` to `unknown_parameters`.

> **Important:** despite its name (for what we will call legacy reasons), the `obsStd` slot of the EKF estimates the **observation-noise variance R**, not a standard deviation. The C++ likelihood takes `obs_var`, and `estimation_profile_ci()` squares the default `obsErr` when it builds the parameter list. So:
> - the estimate returned in the `obsStd` row is **R (variance)**;
> - `initialGuesses["obsStd"] = 0.5` is a starting value for a *variance*;
> - if you hand in your own `full_param_list`, slot 4 must contain the **variance**.

### OU-EKF

```r
y_noisy <- simulate_em_baseline(subFac = 100, Delta = 0.002, maxT = 20,
                                sigma = 2000, obsSd = sqrt(0.75))   # true R = 0.75

res_R <- estimation_profile_ci(c("C", "mu", "sigma", "obsStd"), y_noisy)
res_R
# the "obsStd" row estimates R = 0.75
```


### Fixing R to a positive value

Rather than estimating R, you can prescribe it (as done for the real-EEG sensitivity analysis in Section 6.1 of the paper). Build the parameter list yourself and set slot 4 to the **variance**:

```r
fpl <- usual_param_value        # list; see slot order below
fpl[[1]] <- y                   # observation
fpl[[5]] <- length(y)           # nObs
fpl[[4]] <- 2                   # fixed R = 2 (variance)

estimation_profile_ci(c("C", "mu", "sigma"), y, full_param_list = fpl)
```

`full_param_list` slot order (same as `param_strings`):

| Slot | Name |
|------|------|
| 1 | observation |
| 2 | stepMat |
| 3 | covMat |
| 4 | obsStd (= observation variance R) |
| 5 | nObs |
| 6 | Delta |
| 7 | mu |
| 8–11 | A, a, B, b |
| 12 | C1 |
| 13 | sigma |
| 14–16 | nuMax, nu0, r |
| 17 | burnin (default 500) |



---

## 5. Larger observation steps: internal substepping

If the observation interval is too large for a single numerical step (the paper's simulations show good coverage up to about Δ = 2·10⁻³), use `nSubsteps` with OU-EKF. The observations stay on the original grid; the filter propagates through `nSubsteps` finer internal steps of length `Δ / nSubsteps` between them.

```r
# e.g. Δ = 0.01 with internal step 0.001
y_coarse <- simulate_em_baseline(subFac = 500, Delta = 0.01, maxT = 50, sigma = 2000)

fpl <- usual_param_value
fpl[[1]] <- y_coarse
fpl[[5]] <- length(y_coarse)
fpl[[6]] <- 0.01                       # observation step Delta

estimation_profile_ci(c("C", "mu", "sigma"), y_coarse,
                      full_param_list = fpl,
                      nSubsteps = 10)
```

When `nSubsteps > 1`, the transition and covariance matrices are recomputed internally from `Delta / nSubsteps`, so you do not need to adjust `stepMat` and `covMat` yourself.

For the real EEG recording used in the paper (173.61 Hz sampling, i.e. Δ ≈ 0.00576) the authors used `nSubsteps = 10`.

---

## 6. Applying it to real data

```r
eeg <- scan("your_eeg_channel.txt")     # one-dimensional signal

fpl <- usual_param_value
fpl[[1]] <- eeg
fpl[[5]] <- length(eeg)
fpl[[6]] <- 1 / 173.61                  # sampling interval
fpl[[4]] <- 2                           # prescribed observation variance (paper's OU-EKF_unc uses R = 2)

fit <- estimation_profile_ci(c("C", "mu", "sigma"), eeg,
                             full_param_list = fpl,
                             nSubsteps = 10)
fit
```

The paper's real-data analysis considers three ways of coping with model discrepancy:

1. **Increased observation uncertainty**: fix `R > 0`. This is what the snippet above does.
2. **ACF-penalised likelihood** (penalty `λ·D_ACF(θ)` added to the approximate log-likelihood).
3. **Estimating the process-noise scale τ**. Add "tau" to the unknown_parameter vector.


---

## Troubleshooting

- **Always `source("variableSetup.R")` from the code folder.** Use `setwd()` to change the working directory to the folder in which `variableSetup.R` is placed. The `sourceCpp()` calls use relative paths.
- **`obsStd` is a variance** inside the filter (see Section 4).
- `full_param_list` overrides the defaults entirely. If you build it yourself, set the observation (slot 1), `nObs` (slot 5), and, for non-default step sizes, `Delta` (slot 6).
- If you change `Delta` but keep `sigma` fixed (not estimated) and use `nSubsteps = 1`, the transition matrices are *not* recomputed from the new `Delta`. Either estimate `sigma`, or use `nSubsteps > 1`, or rebuild `stepMat`/`covMat` with `expMt()` and `variableSigmaCovarMat()`.
- Some older functions in `ekf_likelihood.R` are overridden in `ciTests.R` as it is sourced afterwards, the versions in `ciTests.R` are the ones in effect for `estimation_profile_ci`, `wiener_inference_plus_ci` and `em_inference_plus_ci`.
- The sandwich intervals are asymptotic approximations; in the paper they give close to nominal coverage for OU-EKF when Δ is small enough, and deteriorate for large Δ without substepping.
