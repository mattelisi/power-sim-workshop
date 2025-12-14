
# delayed discounting

library(tidyverse)
library(lme4)

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
    sd_pse   = 3.0,    # between-subject SD of PSE
    sd_slope = 0.1,    # between-subject SD of slope
    rho      = -0.3,   # correlation between PSE and slope
    # ---- group effect (small effect, approx d = 0.2) ----
    d_effect = 0.5,    # conceptually similar to Cohen's d
    
    seed = NULL
){
  
  if (!is.null(seed)) set.seed(seed)
  
  # we scale the effect size according to total 
  # latent variability (between subject SD + logistic noise)
  delta <- d_effect * sqrt(sd_pse^2 + ((pi/sqrt(3))/mu_slope)^2)
  
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
      
      # here set the experimental design
      d_i <- design_dd
      
      diff <- d_i$delayed - d_i$immediate
      
      # psychometric (logistic) choice function
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

fit_dd_glmer <- function(dat) {
  
  # --- model: interaction + random slope ---
  dat$diff <- dat$delayed - dat$immediate
  dat$group <- ifelse(dat$group=="PD",1,0)
  fit <-  lme4::glmer(choose_delayed ~ group * diff + (1 + diff | subjectID),
                      data = dat,
                      family = binomial(link = "logit"))
  
  
  # --- extract Wald tests from model summary (simple & standard) ---
  coefs <- summary(fit)$coefficients
  out <- coefs[rownames(coefs)=="group",]
  names(out) <- c("estimate", "se", "z", "p")
  
  # add other info
  out <- as.data.frame(t(out), stringsAsFactors = FALSE)
  out$N <- length(unique(dat$subjectID))
  
  # add effect on the PSE scale
  slope_CL <- coefs[rownames(coefs)=="diff",1]
  slope_PD <- slope_CL + coefs[rownames(coefs)=="group:diff",1] 
  PSE_CL <- - coefs[rownames(coefs)=="(Intercept)",1] / slope_CL
  PSE_PD <- - (coefs[rownames(coefs)=="(Intercept)",1]+coefs[rownames(coefs)=="group",1]) / slope_PD
  out$delta_PSE <- PSE_PD - PSE_CL
  
  return(out)
}

# dat <- simulate_dd(n_per_group = 150, d_effect=0.5)
# fit_dd_glmer(dat)


# run actual simulations
n_sim <- 100
n_participants <- seq(25, 150, by=25)

# # this bit if you want a progress bar
# total_runs <- length(n_participants) * n_sim
# pb <- txtProgressBar(min = 0, max = total_runs, style = 3)
# run_counter <- 0
# 
# # iterate over sample sizes
# sim_res <- {}
# for(n in n_participants){
#   for(i in 1:n_sim){
#     sim_i <- fit_dd_glmer(simulate_dd(n_per_group = n,
#                                           d_effect = 0.5)) 
#     
#     # store output
#     sim_i$n_per_group <- n
#     sim_i$sim <- i
#     sim_res <- rbind(sim_res, sim_i)
#     sim_i$ES <- 0.5
#   }
# }
# 
# write_csv(sim_res, "dd_sim_ES2.csv")


# parallel version

library(parallel)

# Set up parallel processing
n_cores <- 10
cl <- makeCluster(n_cores)

# Load required libraries on each worker
clusterEvalQ(cl, {
  library(tidyverse)
  library(lme4)
})

# Export necessary functions and objects to the cluster
clusterExport(cl, c("fit_dd_glmer", "simulate_dd", "n_sim", "design_dd"))

# Create a function that runs simulations for a single sample size
run_sim_for_n <- function(n) {
  results <- list()
  for(i in 1:n_sim) {
    sim_i <- fit_dd_glmer(simulate_dd(n_per_group = n,
                                          d_effect = 0.5))
    sim_i$n_per_group <- n
    sim_i$sim <- i
    sim_i$ES <- 0.5
    results[[i]] <- sim_i
  }
  do.call(rbind, results)
}

# Run simulations in parallel across sample sizes
sim_res_list <- parLapply(cl, n_participants, run_sim_for_n)

# Stop the cluster
stopCluster(cl)

# Combine all results
sim_res <- do.call(rbind, sim_res_list)
write_csv(sim_res, "dd_sim.csv")


