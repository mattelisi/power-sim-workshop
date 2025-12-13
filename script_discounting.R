
# delayed discounting

# ---- delayed discounting design (60 trials) ----
immediate <- c(
  rep(10, 15),
  13:27,
  rep(20, 15),
  seq(16, 44, length.out = 15)
)

delayed <- c(
  11:25,
  rep(30, 15),
  24:38,
  rep(45, 15)
)

design_dd <- data.frame(
  trial = 1:length(immediate),
  immediate = immediate,
  delayed = delayed
)

simulate_dd <- function(
    n_per_group = 30,
    groups = c("CL", "PD"),
    # ---- population-level parameters (baseline group) ----
    mu_pse   = 12.5,   # mean point of subjective equality (bias)
    mu_slope = 0.3,    # mean slope of psychometric function
    sd_pse   = 2.0,    # between-subject SD of PSE
    sd_slope = 0.15,   # between-subject SD of slope
    rho      = 0,      # correlation between PSE and slope
    # ---- group effect (small effect, approx d = 0.2) ----
    d_effect = 0.2,    # Cohen's d
    
    seed = NULL
){
  
  if (!is.null(seed)) set.seed(seed)
  
  # we scale tghe effect size according to the SD across people
  delta <- d_effect * sd_pse
  
  # ---- covariance matrix for (PSE, slope) ----
  Sigma <- matrix(
    c(sd_pse^2,
      rho * sd_pse * sd_slope,
      rho * sd_pse * sd_slope,
      sd_slope^2),
    nrow = 2, byrow = TRUE
  )
  
  dat <- list()
  row_i <- 1
  
  for (g in groups) {
    
    # group-specific mean shift
    mu_shift_pse   <- ifelse(g != groups[1], delta, 0)
    
    mu_vec <- c(mu_pse + mu_shift_pse, mu_slope)
    
    # ---- sample subject-level parameters ----
    theta <- MASS::mvrnorm(n_per_group, mu = mu_vec, Sigma = Sigma)
    colnames(theta) <- c("pse", "slope")
    
    for (s in 1:n_per_group) {
      
      d_i <- design_dd[sample(nrow(design_dd)), ]
      
      diff <- d_i$delayed - d_i$immediate
      
      # psychometric (logistic) choice rule
      p_delayed <- 1 / (1 + exp(-theta[s, "slope"] *
                                  (diff - theta[s, "pse"])))
      
      choice <- rbinom(nrow(d_i), size = 1, prob = p_delayed)
      
      dat[[row_i]] <- d_i %>%
        mutate(
          choose_delayed = choice,
          p_delayed = p_delayed,
          subjectID = paste0(g, "_", s),
          group = g,
          pse = theta[s, "pse"],
          slope = theta[s, "slope"]
        )
      
      row_i <- row_i + 1
    }
  }
  
  bind_rows(dat)
}


### model fitting function

# Fit a multilevel logistic model for delayed discounting data
# (always includes the group-by-diff interaction AND random slopes).
#
# Expected columns in dat:
#   choose_delayed : 0/1 outcome
#   diff           : delayed - immediate (numeric)
#   group          : factor/character
#   subjectID      : participant ID
#
# Returns a tidy-ish data.frame with estimates and p-values for:
#   - group main effect
#   - diff main effect
#   - group:diff interaction
# plus convergence/singularity flags (useful in power sims).

fit_dd_glmer <- function(dat) {
  
  # --- model: interaction + random slope ---
  dat$diff <- dat$delayed - dat$immediate
  dat$group <- ifelse(dat$group=="PD",1,0)
  fit <-  lme4::glmer(choose_delayed ~ group * diff + (1 + diff | subjectID),
                      data = dat,
                      family = binomial(link = "logit"),
                      control = lme4::glmerControl(optimizer = "bobyqa"))
  
  
  # --- extract Wald tests from model summary (simple & standard) ---
  coefs <- summary(fit)$coefficients
  out <- coefs[rownames(coefs)=="group",]
  names(out) <- c("estimate", "se", "z", "p")
  
  # add other info
  out <- as.data.frame(t(out), stringsAsFactors = FALSE)
  out$N <- length(unique(dat$subjectID))
  
  # add effect on the PSE scale
  slope_CL <- coefs[rownames(coefs)=="diff",1]
  slope_PD <- coefs[rownames(coefs)=="diff",1]
  PSE_CL <- - coefs[rownames(coefs)=="(Intercept)",1] / slope_CL
  PSE_PD <- - (coefs[rownames(coefs)=="(Intercept)",1]+coefs[rownames(coefs)=="group",1]) / slope_PD
  out$delta_PSE <- PSE_PD - PSE_CL
  
  return(out)
}

# dat <- simulate_dd(n_per_group = 400, d_effect=0.2)
# fit_dd_glmer(dat)


# run actual simulations
n_sim <- 100
n_participants <- seq(25, 350, by=25)

total_runs <- length(n_participants) * n_sim
pb <- txtProgressBar(min = 0, max = total_runs, style = 3)
run_counter <- 0

sim_res <- {}
for(n in n_participants){
  for(i in 1:n_sim){
    
    run_counter <- run_counter + 1
    setTxtProgressBar(pb, run_counter)
    
    sim_i <- fit_dd_glmer(simulate_dd(n_per_group = n)) 
    sim_res <- rbind(sim_res, sim_i)
  }
}

write_csv(sim_res, "dd_sim.csv")

