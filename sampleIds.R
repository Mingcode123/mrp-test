library(tidyverse)

set.seed(123)

# Settings
B = 100
sample_prop = 0.10

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# Data
hps_population = read_rds("data/data_clean/hps_week1_population.rds")

states = levels(hps_population$state)

hps_population = hps_population %>%
  mutate(state = factor(state, levels = states))

# Sample ids
hps_sample_ids = list()

for (b in 1:B) {
  
  cat("Sampling replication", b, "\n")
  
  set.seed(123+b)
  
  hps_sample = hps_population %>%
    group_by(state) %>%
    slice_sample(prop = sample_prop) %>%
    ungroup()
  
  hps_sample_ids[[b]] = hps_sample %>%
    transmute(
      replication = b,
      SCRAM
    )
}

hps_sample_ids = bind_rows(hps_sample_ids)

write_rds(hps_sample_ids, "results/hps_sample_ids.rds")
write_csv(hps_sample_ids, "results/hps_sample_ids.csv")