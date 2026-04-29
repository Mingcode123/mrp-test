library(tidyverse)
library(sampling)

set.seed(123)

B = 100

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

hps_population = read_rds("data/data_clean/hps_week1_population.rds")

states = levels(hps_population$state)
if (is.null(states)) states = sort(unique(as.character(hps_population$state)))

hps_population = hps_population %>%
  mutate(
    state = factor(state, levels = states),
    pps_size = log(PWEIGHT)+2*(y==0)
  )

if (any(is.na(hps_population$pps_size))) stop("Missing pps_size")
if (any(hps_population$pps_size<=0)) stop("Non-positive pps_size")

n_sample = round(nrow(hps_population)/15)

pik = inclusionprobabilities(hps_population$pps_size, n_sample)

hps_population = hps_population %>%
  mutate(
    pik = pik,
    sun_w = 1/pik
  )

hps_sample_ids = list()

for (b in 1:B) {
  
  cat("Sampling replication", b, "\n")
  
  set.seed(123+b)
  
  sample_ind = UPmidzuno(hps_population$pik)
  
  hps_sample_ids[[b]] = hps_population %>%
    filter(sample_ind==1) %>%
    transmute(
      replication = b,
      SCRAM,
      pik,
      sun_w,
      pps_size
    )
}

hps_sample_ids = bind_rows(hps_sample_ids)

write_rds(hps_sample_ids, "results/hps_sample_ids.rds")
write_csv(hps_sample_ids, "results/hps_sample_ids.csv")