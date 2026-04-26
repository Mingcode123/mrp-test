library(tidyverse)
library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = detectCores(logical = FALSE))

set.seed(123)

# Settings
B = 10
rep_start = 1
rep_end = B

chains = 2
iter = 2000
warmup = 250

truth_type = "weighted"

if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("results/ht")) dir.create("results/ht", recursive = TRUE)
if (!dir.exists("results/nsb")) dir.create("results/nsb", recursive = TRUE)
if (!dir.exists("results/sb")) dir.create("results/sb", recursive = TRUE)

# Data
hps_population = read_rds("data/data_clean/hps_week1_population.rds")
poststrat_cells = read_rds("data/data_clean/hps_week1_poststrat_cells.rds")
hps_sample_ids = read_rds("results/hps_sample_ids.rds")

states = levels(hps_population$state)
age_levels = levels(hps_population$age_group)
sex_levels = levels(hps_population$sex)
race_levels = levels(hps_population$race_eth)

hps_population = hps_population %>%
  mutate(
    state = factor(state, levels = states),
    age_group = factor(age_group, levels = age_levels),
    sex = factor(sex, levels = sex_levels),
    race_eth = factor(race_eth, levels = race_levels)
  )

poststrat_cells = poststrat_cells %>%
  mutate(
    state = factor(state, levels = states),
    age_group = factor(age_group, levels = age_levels),
    sex = factor(sex, levels = sex_levels),
    race_eth = factor(race_eth, levels = race_levels),
    state_id = as.integer(state)
  )

if (truth_type=="weighted") {
  poststrat_cells = poststrat_cells %>%
    mutate(N_target = N_cell_w)
} else {
  poststrat_cells = poststrat_cells %>%
    mutate(N_target = N_cell)
}

poststrat_cells = poststrat_cells %>%
  filter(!is.na(N_target), N_target>0)

X_poststrat = model.matrix(~ age_group+sex+race_eth, data = poststrat_cells)

# Stan models
nsb_model = stan_model("stan/Binomial.stan")
sb_model = stan_model("stan/Binomial_ICAR.stan")

# State adjacency for spatial model
state_edges = read_csv(
  "data/state_adjacency.csv",
  col_types = cols(
    state1 = col_character(),
    state2 = col_character()
  )
)

icar_edges = state_edges %>%
  filter(state1 %in% states, state2 %in% states) %>%
  mutate(
    node1 = as.integer(factor(state1, levels = states)),
    node2 = as.integer(factor(state2, levels = states)),
    edge_start = pmin(node1,node2),
    edge_end = pmax(node1,node2)
  ) %>%
  distinct(edge_start, edge_end) %>%
  transmute(
    node1 = edge_start,
    node2 = edge_end
  )

target_reps = 1:B
reps = target_reps[target_reps>=rep_start & target_reps<=rep_end]

# HT
for (b in reps) {
  
  if (!(b %in% hps_sample_ids$replication)) {
    cat("Skipping HT replication", b, ": sample ids not found\n")
    next
  }
  
  if (file.exists("results/ht/estimates.rds")) {
    old_ht_estimates = read_rds("results/ht/estimates.rds")
    
    if (b %in% old_ht_estimates$replication) {
      cat("Skipping HT replication", b, ": already in estimates\n")
      next
    }
  }
  
  cat("HT replication", b, "of", B, "\n")
  
  ids_b = hps_sample_ids %>%
    filter(replication==b) %>%
    select(SCRAM)
  
  hps_sample = hps_population %>%
    semi_join(ids_b, by = "SCRAM") %>%
    mutate(
      state = factor(state, levels = states),
      age_group = factor(age_group, levels = age_levels),
      sex = factor(sex, levels = sex_levels),
      race_eth = factor(race_eth, levels = race_levels)
    )
  
  if (truth_type=="weighted") {
    ht_state_estimates = hps_sample %>%
      group_by(state) %>%
      summarise(
        theta_hat = sum(PWEIGHT*y, na.rm = TRUE)/sum(PWEIGHT, na.rm = TRUE),
        lower = NA_real_,
        upper = NA_real_,
        .groups = "drop"
      )
  } else {
    ht_state_estimates = hps_sample %>%
      group_by(state) %>%
      summarise(
        theta_hat = mean(y),
        lower = NA_real_,
        upper = NA_real_,
        .groups = "drop"
      )
  }
  
  ht_estimates_b = ht_state_estimates %>%
    mutate(
      method = "HT",
      replication = b
    ) %>%
    select(replication, method, state, theta_hat, lower, upper)
  
  if (file.exists("results/ht/estimates.rds")) {
    old_ht_estimates = read_rds("results/ht/estimates.rds")
    
    ht_estimates = old_ht_estimates %>%
      filter(replication!=b) %>%
      bind_rows(ht_estimates_b) %>%
      arrange(replication, state)
  } else {
    ht_estimates = ht_estimates_b
  }
  
  write_rds(ht_estimates, "results/ht/estimates.rds")
  write_csv(ht_estimates, "results/ht/estimates.csv")
}

# NSB
for (b in reps) {
  
  if (!(b %in% hps_sample_ids$replication)) {
    cat("Skipping NSB replication", b, ": sample ids not found\n")
    next
  }
  
  if (file.exists("results/nsb/estimates.rds")) {
    old_nsb_estimates = read_rds("results/nsb/estimates.rds")
    
    if (b %in% old_nsb_estimates$replication) {
      cat("Skipping NSB replication", b, ": already in estimates\n")
      next
    }
  }
  
  cat("NSB replication", b, "of", B, "\n")
  
  ids_b = hps_sample_ids %>%
    filter(replication==b) %>%
    select(SCRAM)
  
  hps_sample = hps_population %>%
    semi_join(ids_b, by = "SCRAM") %>%
    mutate(
      state = factor(state, levels = states),
      age_group = factor(age_group, levels = age_levels),
      sex = factor(sex, levels = sex_levels),
      race_eth = factor(race_eth, levels = race_levels),
      state_id = as.integer(state)
    )
  
  if (truth_type=="weighted") {
    hps_sample = hps_sample %>%
      mutate(model_w = PWEIGHT/mean(PWEIGHT, na.rm = TRUE))
  } else {
    hps_sample = hps_sample %>%
      mutate(model_w = 1)
  }
  
  X_sample = model.matrix(~ age_group+sex+race_eth, data = hps_sample)
  
  nsb_data = list(
    D = length(states),
    N = nrow(hps_sample),
    P = ncol(X_sample),
    Y = hps_sample$y,
    wgt = hps_sample$model_w,
    cty = hps_sample$state_id,
    X = X_sample
  )
  
  nsb_fit = sampling(
    nsb_model,
    data = nsb_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    refresh = 50,
    seed = 123+b
  )
  
  nsb_post = extract(nsb_fit)
  nsb_beta = nsb_post$beta
  nsb_mu = nsb_post$mu
  
  nsb_state_estimates = tibble()
  
  for (s in states) {
    idx = poststrat_cells$state==s
    state_id = as.integer(factor(s, levels = states))
    X_state = X_poststrat[idx,,drop = FALSE]
    N_state = poststrat_cells$N_target[idx]
    
    eta = X_state %*% t(nsb_beta)
    eta = sweep(eta,2,nsb_mu[,state_id],"+")
    p = 1/(1+exp(-eta))
    theta_draw = as.numeric(t(N_state/sum(N_state)) %*% p)
    
    nsb_state_estimates = bind_rows(
      nsb_state_estimates,
      tibble(
        state = factor(s, levels = states),
        theta_hat = mean(theta_draw),
        lower = quantile(theta_draw,0.025),
        upper = quantile(theta_draw,0.975)
      )
    )
  }
  
  nsb_estimates_b = nsb_state_estimates %>%
    mutate(
      method = "NSB",
      replication = b
    ) %>%
    select(replication, method, state, theta_hat, lower, upper)
  
  saveRDS(nsb_fit, paste0("results/nsb/fit_rep", b, ".rds"))
  
  if (file.exists("results/nsb/estimates.rds")) {
    old_nsb_estimates = read_rds("results/nsb/estimates.rds")
    
    nsb_estimates = old_nsb_estimates %>%
      filter(replication!=b) %>%
      bind_rows(nsb_estimates_b) %>%
      arrange(replication, state)
  } else {
    nsb_estimates = nsb_estimates_b
  }
  
  write_rds(nsb_estimates, "results/nsb/estimates.rds")
  write_csv(nsb_estimates, "results/nsb/estimates.csv")
}

# SB
for (b in reps) {
  
  if (!(b %in% hps_sample_ids$replication)) {
    cat("Skipping SB replication", b, ": sample ids not found\n")
    next
  }
  
  if (file.exists("results/sb/estimates.rds")) {
    old_sb_estimates = read_rds("results/sb/estimates.rds")
    
    if (b %in% old_sb_estimates$replication) {
      cat("Skipping SB replication", b, ": already in estimates\n")
      next
    }
  }
  
  cat("SB replication", b, "of", B, "\n")
  
  ids_b = hps_sample_ids %>%
    filter(replication==b) %>%
    select(SCRAM)
  
  hps_sample = hps_population %>%
    semi_join(ids_b, by = "SCRAM") %>%
    mutate(
      state = factor(state, levels = states),
      age_group = factor(age_group, levels = age_levels),
      sex = factor(sex, levels = sex_levels),
      race_eth = factor(race_eth, levels = race_levels),
      state_id = as.integer(state)
    )
  
  if (truth_type=="weighted") {
    hps_sample = hps_sample %>%
      mutate(model_w = PWEIGHT/mean(PWEIGHT, na.rm = TRUE))
  } else {
    hps_sample = hps_sample %>%
      mutate(model_w = 1)
  }
  
  X_sample = model.matrix(~ age_group+sex+race_eth, data = hps_sample)
  
  sb_data = list(
    D = length(states),
    D_edges = nrow(icar_edges),
    node1 = icar_edges$node1,
    node2 = icar_edges$node2,
    N = nrow(hps_sample),
    P = ncol(X_sample),
    Y = hps_sample$y,
    wgt = hps_sample$model_w,
    cty = hps_sample$state_id,
    X = X_sample
  )
  
  sb_fit = sampling(
    sb_model,
    data = sb_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    refresh = 50,
    seed = 456+b
  )
  
  sb_post = extract(sb_fit)
  sb_beta = sb_post$beta
  sb_mu = sb_post$mu
  sb_phi = sb_post$phi
  sb_sigma_phi = sb_post$sigmaPhi
  
  sb_state_estimates = tibble()
  
  for (s in states) {
    idx = poststrat_cells$state==s
    state_id = as.integer(factor(s, levels = states))
    X_state = X_poststrat[idx,,drop = FALSE]
    N_state = poststrat_cells$N_target[idx]
    
    eta = X_state %*% t(sb_beta)
    eta = sweep(eta,2,sb_mu[,state_id]+sb_sigma_phi*sb_phi[,state_id],"+")
    p = 1/(1+exp(-eta))
    theta_draw = as.numeric(t(N_state/sum(N_state)) %*% p)
    
    sb_state_estimates = bind_rows(
      sb_state_estimates,
      tibble(
        state = factor(s, levels = states),
        theta_hat = mean(theta_draw),
        lower = quantile(theta_draw,0.025),
        upper = quantile(theta_draw,0.975)
      )
    )
  }
  
  sb_estimates_b = sb_state_estimates %>%
    mutate(
      method = "SB",
      replication = b
    ) %>%
    select(replication, method, state, theta_hat, lower, upper)
  
  saveRDS(sb_fit, paste0("results/sb/fit_rep", b, ".rds"))
  
  if (file.exists("results/sb/estimates.rds")) {
    old_sb_estimates = read_rds("results/sb/estimates.rds")
    
    sb_estimates = old_sb_estimates %>%
      filter(replication!=b) %>%
      bind_rows(sb_estimates_b) %>%
      arrange(replication, state)
  } else {
    sb_estimates = sb_estimates_b
  }
  
  write_rds(sb_estimates, "results/sb/estimates.rds")
  write_csv(sb_estimates, "results/sb/estimates.csv")
}