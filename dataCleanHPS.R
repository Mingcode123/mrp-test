library(tidyverse)

out_dir = "data/data_clean"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

hps = read_csv(
  "data/HPS_Week01_PUF_CSV/pulse2020_puf_01.csv",
  show_col_types = FALSE
)

repweights = read_csv(
  "data/HPS_Week01_PUF_CSV/pulse2020_repwgt_puf_01.csv",
  show_col_types = FALSE
)

hps_population = hps %>%
  filter(EXPCTLOSS %in% c(1,2)) %>%
  mutate(
    y = if_else(EXPCTLOSS==1,1,0),
    state = str_pad(as.character(EST_ST),2,pad = "0"),
    age = 2020-TBIRTH_YEAR,
    age_group = case_when(
      age>=18 & age<=24 ~ "18-24",
      age>=25 & age<=39 ~ "25-39",
      age>=40 & age<=54 ~ "40-54",
      age>=55 & age<=64 ~ "55-64",
      age>=65 ~ "65+",
      TRUE ~ NA_character_
    ),
    sex = case_when(
      EGENDER==1 ~ "Male",
      EGENDER==2 ~ "Female",
      TRUE ~ NA_character_
    ),
    race_eth = case_when(
      RHISPANIC==2 ~ "Hispanic",
      RHISPANIC==1 & RRACE==1 ~ "Non-Hispanic White",
      RHISPANIC==1 & RRACE==2 ~ "Non-Hispanic Black",
      RHISPANIC==1 & RRACE==3 ~ "Non-Hispanic Asian",
      RHISPANIC==1 & RRACE==4 ~ "Non-Hispanic Other",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !state %in% c("02","15"),
    !is.na(age_group),
    !is.na(sex),
    !is.na(race_eth)
  ) %>%
  mutate(
    state = factor(state),
    age_group = factor(age_group, levels = c("18-24","25-39","40-54","55-64","65+")),
    sex = factor(sex, levels = c("Male","Female")),
    race_eth = factor(
      race_eth,
      levels = c(
        "Hispanic",
        "Non-Hispanic White",
        "Non-Hispanic Black",
        "Non-Hispanic Asian",
        "Non-Hispanic Other"
      )
    ),
    cell_id = interaction(state, age_group, sex, race_eth, drop = TRUE)
  ) %>%
  select(
    SCRAM,
    WEEK,
    state,
    PWEIGHT,
    y,
    EXPCTLOSS,
    age,
    age_group,
    sex,
    race_eth,
    cell_id,
    everything()
  )

hps_repweights = repweights %>%
  semi_join(hps_population, by = "SCRAM")

poststrat_cells = hps_population %>%
  count(state, age_group, sex, race_eth, name = "N_cell") %>%
  group_by(state) %>%
  mutate(
    N_state = sum(N_cell),
    cell_prop_state = N_cell/N_state
  ) %>%
  ungroup()

poststrat_cells_w = hps_population %>%
  group_by(state, age_group, sex, race_eth) %>%
  summarise(
    N_cell_w = sum(PWEIGHT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(state) %>%
  mutate(
    N_state_w = sum(N_cell_w),
    cell_prop_state_w = N_cell_w/N_state_w
  ) %>%
  ungroup()

poststrat_cells = poststrat_cells %>%
  left_join(
    poststrat_cells_w,
    by = c("state","age_group","sex","race_eth")
  )

state_truth = hps_population %>%
  group_by(state) %>%
  summarise(
    n_state = n(),
    y_sum = sum(y),
    theta_pop = mean(y),
    theta_pop_w = sum(PWEIGHT*y, na.rm = TRUE)/sum(PWEIGHT, na.rm = TRUE),
    .groups = "drop"
  )

write_rds(hps_population, file.path(out_dir, "hps_week1_population.rds"))
write_rds(hps_repweights, file.path(out_dir, "hps_week1_repweights.rds"))
write_rds(poststrat_cells, file.path(out_dir, "hps_week1_poststrat_cells.rds"))
write_rds(state_truth, file.path(out_dir, "hps_week1_state_truth.rds"))
write_csv(state_truth, file.path(out_dir, "hps_week1_state_truth.csv"))