library(tidyverse)
library(rstan)
library(rstanarm)
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
n_draws = 1000

truth_type = "weighted"

if (!dir.exists("results/mrp")) dir.create("results/mrp", recursive = TRUE)
if (!dir.exists("results/mrp2")) dir.create("results/mrp2", recursive = TRUE)
if (!dir.exists("results/mrp3")) dir.create("results/mrp3", recursive = TRUE)

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
    race_eth = factor(race_eth, levels = race_levels)
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

hps_population = hps_population %>%
  mutate(
    age_race = interaction(age_group,race_eth, sep = ":", drop = FALSE),
    state_race = interaction(state,race_eth, sep = ":", drop = FALSE),
    sex_race = interaction(sex,race_eth, sep = ":", drop = FALSE)
  )

poststrat_cells = poststrat_cells %>%
  mutate(
    age_race = interaction(age_group,race_eth, sep = ":", drop = FALSE),
    state_race = interaction(state,race_eth, sep = ":", drop = FALSE),
    sex_race = interaction(sex,race_eth, sep = ":", drop = FALSE)
  )

age_race_levels = levels(poststrat_cells$age_race)
state_race_levels = levels(poststrat_cells$state_race)
sex_race_levels = levels(poststrat_cells$sex_race)

hps_population = hps_population %>%
  mutate(
    age_race = factor(age_race, levels = age_race_levels),
    state_race = factor(state_race, levels = state_race_levels),
    sex_race = factor(sex_race, levels = sex_race_levels)
  )

poststrat_cells = poststrat_cells %>%
  mutate(
    age_race = factor(age_race, levels = age_race_levels),
    state_race = factor(state_race, levels = state_race_levels),
    sex_race = factor(sex_race, levels = sex_race_levels)
  )

target_reps = 1:B
reps = target_reps[target_reps>=rep_start & target_reps<=rep_end]

# MRP
for (b in reps) {
  
  if (!(b %in% hps_sample_ids$replication)) {
    cat("Skipping MRP replication", b, ": sample ids not found\n")
    next
  }
  
  if (file.exists("results/mrp/estimates.rds")) {
    old_mrp_estimates = read_rds("results/mrp/estimates.rds")
    
    if (b %in% old_mrp_estimates$replication) {
      cat("Skipping MRP replication", b, ": already in estimates\n")
      next
    }
  }
  
  cat("MRP replication", b, "\n")
  
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
      age_race = factor(age_race, levels = age_race_levels),
      state_race = factor(state_race, levels = state_race_levels),
      sex_race = factor(sex_race, levels = sex_race_levels)
    )
  
  if (truth_type=="weighted") {
    hps_sample = hps_sample %>%
      mutate(model_w = PWEIGHT/mean(PWEIGHT, na.rm = TRUE))
  } else {
    hps_sample = hps_sample %>%
      mutate(model_w = 1)
  }
  
  mrp_fit = stan_glmer(
    y ~ (1|age_group)+(1|sex)+(1|race_eth)+(1|state),
    family = binomial(link = "logit"),
    data = hps_sample,
    weights = model_w,
    chains = chains,
    iter = iter,
    warmup = warmup,
    refresh = 50,
    seed = 789+b
  )
  
  mrp_draws = posterior_epred(
    mrp_fit,
    newdata = poststrat_cells,
    draws = n_draws,
    allow_new_levels = TRUE
  )
  
  mrp_state_estimates = tibble()
  
  for (s in states) {
    idx = poststrat_cells$state==s
    N_state = poststrat_cells$N_target[idx]
    theta_draw = as.numeric(mrp_draws[,idx,drop = FALSE] %*% (N_state/sum(N_state)))
    
    mrp_state_estimates = bind_rows(
      mrp_state_estimates,
      tibble(
        state = factor(s, levels = states),
        theta_hat = mean(theta_draw),
        lower = quantile(theta_draw,0.025),
        upper = quantile(theta_draw,0.975)
      )
    )
  }
  
  mrp_estimates_b = mrp_state_estimates %>%
    mutate(
      method = "MRP",
      replication = b
    ) %>%
    select(replication, method, state, theta_hat, lower, upper)
  
  saveRDS(mrp_fit, paste0("results/mrp/fit_rep", b, ".rds"))
  
  if (file.exists("results/mrp/estimates.rds")) {
    old_mrp_estimates = read_rds("results/mrp/estimates.rds")
    
    mrp_estimates = old_mrp_estimates %>%
      filter(replication!=b) %>%
      bind_rows(mrp_estimates_b) %>%
      arrange(replication, state)
  } else {
    mrp_estimates = mrp_estimates_b
  }
  
  write_rds(mrp_estimates, "results/mrp/estimates.rds")
  write_csv(mrp_estimates, "results/mrp/estimates.csv")
}

# MRP2
for (b in reps) {
  
  if (!(b %in% hps_sample_ids$replication)) {
    cat("Skipping MRP2 replication", b, ": sample ids not found\n")
    next
  }
  
  if (file.exists("results/mrp2/estimates.rds")) {
    old_mrp2_estimates = read_rds("results/mrp2/estimates.rds")
    
    if (b %in% old_mrp2_estimates$replication) {
      cat("Skipping MRP2 replication", b, ": already in estimates\n")
      next
    }
  }
  
  cat("MRP2 replication", b, "\n")
  
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
      age_race = factor(age_race, levels = age_race_levels),
      state_race = factor(state_race, levels = state_race_levels),
      sex_race = factor(sex_race, levels = sex_race_levels)
    )
  
  if (truth_type=="weighted") {
    hps_sample = hps_sample %>%
      mutate(model_w = PWEIGHT/mean(PWEIGHT, na.rm = TRUE))
  } else {
    hps_sample = hps_sample %>%
      mutate(model_w = 1)
  }
  
  mrp2_fit = stan_glmer(
    y ~ (1|age_group)+(1|sex)+(1|race_eth)+(1|state)+(1|age_race)+(1|state_race),
    family = binomial(link = "logit"),
    data = hps_sample,
    weights = model_w,
    chains = chains,
    iter = iter,
    warmup = warmup,
    refresh = 50,
    seed = 101112+b
  )
  
  mrp2_draws = posterior_epred(
    mrp2_fit,
    newdata = poststrat_cells,
    draws = n_draws,
    allow_new_levels = TRUE
  )
  
  mrp2_state_estimates = tibble()
  
  for (s in states) {
    idx = poststrat_cells$state==s
    N_state = poststrat_cells$N_target[idx]
    theta_draw = as.numeric(mrp2_draws[,idx,drop = FALSE] %*% (N_state/sum(N_state)))
    
    mrp2_state_estimates = bind_rows(
      mrp2_state_estimates,
      tibble(
        state = factor(s, levels = states),
        theta_hat = mean(theta_draw),
        lower = quantile(theta_draw,0.025),
        upper = quantile(theta_draw,0.975)
      )
    )
  }
  
  mrp2_estimates_b = mrp2_state_estimates %>%
    mutate(
      method = "MRP2",
      replication = b
    ) %>%
    select(replication, method, state, theta_hat, lower, upper)
  
  saveRDS(mrp2_fit, paste0("results/mrp2/fit_rep", b, ".rds"))
  
  if (file.exists("results/mrp2/estimates.rds")) {
    old_mrp2_estimates = read_rds("results/mrp2/estimates.rds")
    
    mrp2_estimates = old_mrp2_estimates %>%
      filter(replication!=b) %>%
      bind_rows(mrp2_estimates_b) %>%
      arrange(replication, state)
  } else {
    mrp2_estimates = mrp2_estimates_b
  }
  
  write_rds(mrp2_estimates, "results/mrp2/estimates.rds")
  write_csv(mrp2_estimates, "results/mrp2/estimates.csv")
}

# MRP3
for (b in reps) {
  
  if (!(b %in% hps_sample_ids$replication)) {
    cat("Skipping MRP3 replication", b, ": sample ids not found\n")
    next
  }
  
  if (file.exists("results/mrp3/estimates.rds")) {
    old_mrp3_estimates = read_rds("results/mrp3/estimates.rds")
    
    if (b %in% old_mrp3_estimates$replication) {
      cat("Skipping MRP3 replication", b, ": already in estimates\n")
      next
    }
  }
  
  cat("MRP3 replication", b, "\n")
  
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
      age_race = factor(age_race, levels = age_race_levels),
      state_race = factor(state_race, levels = state_race_levels),
      sex_race = factor(sex_race, levels = sex_race_levels)
    )
  
  if (truth_type=="weighted") {
    hps_sample = hps_sample %>%
      mutate(model_w = PWEIGHT/mean(PWEIGHT, na.rm = TRUE))
  } else {
    hps_sample = hps_sample %>%
      mutate(model_w = 1)
  }
  
  mrp3_fit = stan_glmer(
    y ~ (1|age_group)+(1|sex)+(1|race_eth)+(1|state)+(1|age_race)+(1|sex_race),
    family = binomial(link = "logit"),
    data = hps_sample,
    weights = model_w,
    chains = chains,
    iter = iter,
    warmup = warmup,
    refresh = 50,
    seed = 131415+b
  )
  
  mrp3_draws = posterior_epred(
    mrp3_fit,
    newdata = poststrat_cells,
    draws = n_draws,
    allow_new_levels = TRUE
  )
  
  mrp3_state_estimates = tibble()
  
  for (s in states) {
    idx = poststrat_cells$state==s
    N_state = poststrat_cells$N_target[idx]
    theta_draw = as.numeric(mrp3_draws[,idx,drop = FALSE] %*% (N_state/sum(N_state)))
    
    mrp3_state_estimates = bind_rows(
      mrp3_state_estimates,
      tibble(
        state = factor(s, levels = states),
        theta_hat = mean(theta_draw),
        lower = quantile(theta_draw,0.025),
        upper = quantile(theta_draw,0.975)
      )
    )
  }
  
  mrp3_estimates_b = mrp3_state_estimates %>%
    mutate(
      method = "MRP3",
      replication = b
    ) %>%
    select(replication, method, state, theta_hat, lower, upper)
  
  saveRDS(mrp3_fit, paste0("results/mrp3/fit_rep", b, ".rds"))
  
  if (file.exists("results/mrp3/estimates.rds")) {
    old_mrp3_estimates = read_rds("results/mrp3/estimates.rds")
    
    mrp3_estimates = old_mrp3_estimates %>%
      filter(replication!=b) %>%
      bind_rows(mrp3_estimates_b) %>%
      arrange(replication, state)
  } else {
    mrp3_estimates = mrp3_estimates_b
  }
  
  write_rds(mrp3_estimates, "results/mrp3/estimates.rds")
  write_csv(mrp3_estimates, "results/mrp3/estimates.csv")
}