rm(list = ls())

# ---- packages ----
library(lme4)
library(lmerTest)
library(tidyverse)
library(MASS)   # for mvrnorm()


# --- Helper: simulate an AR(1) time series ---
r_ar1 <- function(n, s, r){
  innov_sd <- s * sqrt(1 - r^2)
  innov <- rnorm(n, mean = 0, sd = innov_sd)
  
  out <- numeric(n)
  out[1] <- rnorm(1, mean = 0, sd = s)
  for (t in 2:n) {
    out[t] <- r * out[t - 1] + innov[t]
  }
  out
}

# ---- core simulator ----
sim_diary <- function(
    n_child = 50,
    n_days  = 7,
    
    # Predictor X structure
    sd_between_x = 0.6,  # SD of each child's mean X across children (between-person variability)
    sd_within_x  = 0.8,  # SD of day-to-day deviations within a child (within-person variability)
    rho_x_ar1    = 0.2,  # AR(1) correlation in within-person deviations of X across consecutive days
    
    # Outcome Y structure
    beta0   = 0,         # intercept
    beta_w  = 0.20,      # within-person effect of x_cwc on Y  (in SD units of x_cwc)
    beta_b  = 0.20,      # between-person effect of x_bar on Y (in SD units of x_bar)
    tau0    = 0.60,      # SD of random intercepts
    tau1    = 0.15,      # SD of random slopes for x_cwc
    rho_u   = 0.0,       # correlation between intercept and slope random effects
    sigma_e = 0.80,      # residual SD (day-level noise), assumed independent over days
    
    weekday_fx = FALSE,  # include weekday fixed effects in the fitted model?
    
    # Missingness (day-level MCAR)
    mcar_rate = 0.05,    # probability any given day is missing (independent of everything)
    
    seed = NULL
){
  if (!is.null(seed)) set.seed(seed)
  
  # --- Random effects covariance matrix for (u0, u1) ---
  Sigma_u <- matrix(
    c(tau0^2,              rho_u * tau0 * tau1,
      rho_u * tau0 * tau1,        tau1^2),
    nrow = 2, byrow = TRUE
  )
  
  # --- Balanced panel: each child measured each day ---
  df <- expand.grid(id = 1:n_child, day = 1:n_days) %>%
    as_tibble() %>%
    arrange(id, day) %>%
    mutate(
      weekday = factor((day - 1) %% 7,
                       labels = c("Mon","Tue","Wed","Thu","Fri","Sat","Sun"))
    )
  
  # --- Person-level mean of X (between-person differences) ---
  x_bar <- rnorm(n_child, mean = 0, sd = sd_between_x)
  df <- df %>%
    left_join(tibble(id = 1:n_child, x_bar = x_bar), by = "id")
  
  # --- Within-person daily deviations of X with AR(1) structure ---
  df <- df %>%
    group_by(id) %>%
    mutate(x_within = r_ar1(n_days, sd_within_x, rho_x_ar1)) %>%
    ungroup() %>%
    mutate(
      x     = x_bar + x_within,
      x_cwc = x - x_bar
    )
  
  # --- Standardise predictors separately ---
  # IMPORTANT: beta_w and beta_b are interpreted in SD units of each component (x_cwc and x_bar).
  df <- df %>%
    mutate(
      x     = as.numeric(scale(x)),
      x_bar = as.numeric(scale(x_bar)),
      x_cwc = as.numeric(scale(x_cwc))
    )
  
  # --- Random effects per child (random intercept + random slope for x_cwc) ---
  RE <- MASS::mvrnorm(n_child, mu = c(0, 0), Sigma = Sigma_u)
  RE <- tibble(id = 1:n_child, u0 = RE[,1], u1 = RE[,2])
  df <- df %>% left_join(RE, by = "id")
  
  # --- Day-level residuals e (independent across days) ---
  # This keeps the fitted lmer() model correctly specified w.r.t. residual correlation.
  df <- df %>%
    mutate(e = rnorm(nrow(df), mean = 0, sd = sigma_e))
  
  # --- True weekday effect (set to zero here; weekday_fx just changes the fitted model) ---
  df <- df %>% mutate(weekday_eff = 0)
  
  # --- Outcome generation ---
  df <- df %>%
    mutate(
      y = beta0 + u0 +
        (beta_w + u1) * x_cwc +
        beta_b * x_bar +
        weekday_eff +
        e
    )
  
  # Optional: standardise outcome
  df <- df %>% mutate(y = as.numeric(scale(y)))
  
  # --- Missingness: day-level MCAR only ---
  if (mcar_rate > 0) {
    keep <- runif(nrow(df)) >= mcar_rate
    df <- df[keep, ]
  }
  
  df
}

# ---- one run: fit model and extract p-values/CIs ----
fit_once <- function(df, weekday_fx = TRUE) {
  
  if (weekday_fx) {
    f <- y ~ x_cwc + x_bar + weekday + (1 + x_cwc | id)
  } else {
    f <- y ~ x_cwc + x_bar + (1 + x_cwc | id)
  }
  
  m <- suppressMessages(lmer(f, data = df, REML = FALSE))
  
  broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE) %>%
    filter(term %in% c("x_cwc", "x_bar")) %>%
    transmute(
      term,
      est = estimate,
      lwr = conf.low,
      upr = conf.high,
      p   = p.value
    )
}

# ---- power wrapper over many simulations (FOR LOOP) ----
power_grid <- function( n_sims = 100, weekday_fx = FALSE, ...) {
  all_sims <- vector("list", n_sims)
  
  for (s in 1:n_sims) {
    df  <- sim_diary(seed = s, weekday_fx = weekday_fx, ...)
    est <- fit_once(df, weekday_fx = weekday_fx)
    est$sim <- s
    all_sims[[s]] <- est
  }
  
  res <- bind_rows(all_sims)
  
  res %>%
    group_by(term) %>%
    summarise(
      power_0.05 = mean(p < 0.05, na.rm = TRUE),
      median_est = median(est, na.rm = TRUE),
      mean_est   = mean(est, na.rm = TRUE),
      SE_est     = sd(est, na.rm = TRUE),
      q25        = quantile(est, 0.25, na.rm = TRUE),
      q75        = quantile(est, 0.75, na.rm = TRUE),
      .groups    = "drop"
    )
}

# ---- example run ----
quick_A <- power_grid(
  n_sims = 200,
  n_child = 50,
  n_days = 7,
  beta_w = 0.20,
  beta_b = 0.20,
  tau0 = 0.60,
  tau1 = 0.15,
  sigma_e = 0.80,
  rho_x_ar1 = 0.2,
  mcar_rate = 0.05
)

quick_A
