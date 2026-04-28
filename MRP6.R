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

if (!dir.exists("results/mrp6")) dir.create("results/mrp6", recursive = TRUE)

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
    N_target = N_cell
  ) %>%
  filter(!is.na(N_target), N_target>0)

hps_population = hps_population %>%
  mutate(
    age_race = interaction(age_group,race_eth, sep = ":", drop = FALSE),
    state_race = interaction(state,race_eth, sep = ":", drop = FALSE),
    sex_race = interaction(sex,race_eth, sep = ":", drop = FALSE),
    age_sex = interaction(age_group,sex, sep = ":", drop = FALSE)
  )

poststrat_cells = poststrat_cells %>%
  mutate(
    age_race = interaction(age_group,race_eth, sep = ":", drop = FALSE),
    state_race = interaction(state,race_eth, sep = ":", drop = FALSE),
    sex_race = interaction(sex,race_eth, sep = ":", drop = FALSE),
    age_sex = interaction(age_group,sex, sep = ":", drop = FALSE)
  )

age_race_levels = levels(poststrat_cells$age_race)
state_race_levels = levels(poststrat_cells$state_race)
sex_race_levels = levels(poststrat_cells$sex_race)
age_sex_levels = levels(poststrat_cells$age_sex)

hps_population = hps_population %>%
  mutate(
    age_race = factor(age_race, levels = age_race_levels),
    state_race = factor(state_race, levels = state_race_levels),
    sex_race = factor(sex_race, levels = sex_race_levels),
    age_sex = factor(age_sex, levels = age_sex_levels)
  )

poststrat_cells = poststrat_cells %>%
  mutate(
    age_race = factor(age_race, levels = age_race_levels),
    state_race = factor(state_race, levels = state_race_levels),
    sex_race = factor(sex_race, levels = sex_race_levels),
    age_sex = factor(age_sex, levels = age_sex_levels)
  )

target_reps = 1:B
reps = target_reps[target_reps>=rep_start & target_reps<=rep_end]

# MRP6: old MRP2 structure under new Sun-style weighting
for (b in reps) {
  
  if (!(b %in% hps_sample_ids$replication)) {
    cat("Skipping MRP6 replication", b, ": sample ids not found\n")
    next
  }
  
  if (file.exists("results/mrp6/estimates.rds")) {
    old_mrp6_estimates = read_rds("results/mrp6/estimates.rds")
    
    if (b %in% old_mrp6_estimates$replication) {
      cat("Skipping MRP6 replication", b, ": already in estimates\n")
      next
    }
  }
  
  cat("MRP6 replication", b, "\n")
  
  ids_b = hps_sample_ids %>%
    filter(replication==b) %>%
    select(SCRAM, pik, sun_w)
  
  hps_sample = hps_population %>%
    inner_join(ids_b, by = "SCRAM") %>%
    mutate(
      state = factor(state, levels = states),
      age_group = factor(age_group, levels = age_levels),
      sex = factor(sex, levels = sex_levels),
      race_eth = factor(race_eth, levels = race_levels),
      age_race = factor(age_race, levels = age_race_levels),
      state_race = factor(state_race, levels = state_race_levels),
      sex_race = factor(sex_race, levels = sex_race_levels),
      age_sex = factor(age_sex, levels = age_sex_levels),
      model_w = length(sun_w)*sun_w/sum(sun_w)
    )
  
  mrp6_fit = stan_glmer(
    y ~ (1|age_group)+(1|sex)+(1|race_eth)+(1|state)+(1|age_race)+(1|state_race),
    family = binomial(link = "logit"),
    data = hps_sample,
    weights = model_w,
    chains = chains,
    iter = iter,
    warmup = warmup,
    refresh = 50,
    seed = 222324+b
  )
  
  mrp6_draws = posterior_epred(
    mrp6_fit,
    newdata = poststrat_cells,
    draws = n_draws,
    allow_new_levels = TRUE
  )
  
  mrp6_state_estimates = tibble()
  
  for (s in states) {
    idx = poststrat_cells$state==s
    N_state = poststrat_cells$N_target[idx]
    theta_draw = as.numeric(mrp6_draws[,idx,drop = FALSE] %*% (N_state/sum(N_state)))
    
    mrp6_state_estimates = bind_rows(
      mrp6_state_estimates,
      tibble(
        state = factor(s, levels = states),
        theta_hat = mean(theta_draw),
        lower = quantile(theta_draw,0.025),
        upper = quantile(theta_draw,0.975)
      )
    )
  }
  
  mrp6_estimates_b = mrp6_state_estimates %>%
    mutate(
      method = "MRP6",
      replication = b
    ) %>%
    select(replication, method, state, theta_hat, lower, upper)
  
  saveRDS(mrp6_fit, paste0("results/mrp6/fit_rep", b, ".rds"))
  
  if (file.exists("results/mrp6/estimates.rds")) {
    old_mrp6_estimates = read_rds("results/mrp6/estimates.rds")
    
    mrp6_estimates = old_mrp6_estimates %>%
      filter(replication!=b) %>%
      bind_rows(mrp6_estimates_b) %>%
      arrange(replication, state)
  } else {
    mrp6_estimates = mrp6_estimates_b
  }
  
  write_rds(mrp6_estimates, "results/mrp6/estimates.rds")
  write_csv(mrp6_estimates, "results/mrp6/estimates.csv")
}