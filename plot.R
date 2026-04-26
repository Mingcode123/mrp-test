library(tidyverse)
library(maps)

truth_type = "weighted"

if (!dir.exists("results/plot")) dir.create("results/plot", recursive = TRUE)

hps_population = read_rds("data/data_clean/hps_week1_population.rds")
state_truth = read_rds("data/data_clean/hps_week1_state_truth.rds")
hps_sample_ids = read_rds("results/hps_sample_ids.rds")

ht_estimates = read_rds("results/ht/estimates.rds")
nsb_estimates = read_rds("results/nsb/estimates.rds")
sb_estimates = read_rds("results/sb/estimates.rds")
mrp_estimates = read_rds("results/mrp/estimates.rds")

if (file.exists("results/mrp2/estimates.rds")) {
  mrp2_estimates = read_rds("results/mrp2/estimates.rds")
} else {
  mrp2_estimates = tibble()
}

if (file.exists("results/mrp3/estimates.rds")) {
  mrp3_estimates = read_rds("results/mrp3/estimates.rds")
} else {
  mrp3_estimates = tibble()
}

states = levels(hps_population$state)
method_levels = c("HT","NSB","SB","MRP","MRP2","MRP3")

state_labels = c(
  "01"="AL","02"="AK","04"="AZ","05"="AR","06"="CA","08"="CO","09"="CT",
  "10"="DE","11"="DC","12"="FL","13"="GA","15"="HI","16"="ID","17"="IL",
  "18"="IN","19"="IA","20"="KS","21"="KY","22"="LA","23"="ME","24"="MD",
  "25"="MA","26"="MI","27"="MN","28"="MS","29"="MO","30"="MT","31"="NE",
  "32"="NV","33"="NH","34"="NJ","35"="NM","36"="NY","37"="NC","38"="ND",
  "39"="OH","40"="OK","41"="OR","42"="PA","44"="RI","45"="SC","46"="SD",
  "47"="TN","48"="TX","49"="UT","50"="VT","51"="VA","53"="WA","54"="WV",
  "55"="WI","56"="WY"
)

state_region = tibble(
  state_abbr = state.abb,
  region = tolower(state.name)
) %>%
  bind_rows(
    tibble(state_abbr = "DC", region = "district of columbia")
  )

hps_population = hps_population %>%
  mutate(state = factor(state, levels = states))

state_truth = state_truth %>%
  mutate(state = factor(state, levels = states))

if (truth_type=="weighted") {
  state_truth = state_truth %>%
    mutate(theta_true = theta_pop_w)
} else {
  state_truth = state_truth %>%
    mutate(theta_true = theta_pop)
}

state_truth = state_truth %>%
  select(state, theta_true)

sample_state_size = hps_sample_ids %>%
  left_join(
    hps_population %>% select(SCRAM, state),
    by = "SCRAM"
  ) %>%
  count(replication, state, name = "n_state_sample")

all_methods_estimates = bind_rows(
  ht_estimates,
  nsb_estimates,
  sb_estimates,
  mrp_estimates,
  mrp2_estimates,
  mrp3_estimates
) %>%
  mutate(
    method = factor(method, levels = method_levels),
    state = factor(state, levels = states)
  ) %>%
  left_join(state_truth, by = "state") %>%
  left_join(sample_state_size, by = c("replication","state")) %>%
  mutate(
    error = theta_hat-theta_true,
    abs_error = abs(error),
    sq_error = error^2,
    covered = if_else(is.na(lower) | is.na(upper),NA,theta_true>=lower & theta_true<=upper),
    state_abbr = recode(as.character(state), !!!state_labels)
  ) %>%
  filter(!is.na(method))

state_metrics = all_methods_estimates %>%
  group_by(method, state) %>%
  summarise(
    mse = mean(sq_error, na.rm = TRUE),
    rmse = sqrt(mse),
    bias = mean(error, na.rm = TRUE),
    squared_bias = bias^2,
    coverage = mean(covered, na.rm = TRUE),
    theta_hat_mean = mean(theta_hat, na.rm = TRUE),
    theta_true = first(theta_true),
    n_state_sample_mean = mean(n_state_sample, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    coverage = if_else(is.nan(coverage),NA_real_,coverage),
    state_abbr = recode(as.character(state), !!!state_labels),
    error_mean = theta_hat_mean-theta_true
  )

overall_metrics = state_metrics %>%
  group_by(method) %>%
  summarise(
    mse = mean(mse, na.rm = TRUE),
    rmse = sqrt(mse),
    squared_bias = mean(squared_bias, na.rm = TRUE),
    coverage = mean(coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    coverage = if_else(is.nan(coverage),NA_real_,coverage)
  )

ht_rmse = overall_metrics %>%
  filter(method=="HT") %>%
  pull(rmse)

summary_metrics = overall_metrics %>%
  mutate(
    rmse_reduction = (ht_rmse-rmse)/ht_rmse
  )

write_rds(all_methods_estimates, "results/plot/all_methods_estimates.rds")
write_csv(all_methods_estimates, "results/plot/all_methods_estimates.csv")

write_rds(state_metrics, "results/plot/state_metrics.rds")
write_csv(state_metrics, "results/plot/state_metrics.csv")

write_rds(overall_metrics, "results/plot/overall_metrics.rds")
write_csv(overall_metrics, "results/plot/overall_metrics.csv")

write_rds(summary_metrics, "results/plot/summary_metrics.rds")
write_csv(summary_metrics, "results/plot/summary_metrics.csv")

# Figure 1: estimates against truth
axis_range = range(
  c(all_methods_estimates$theta_true,all_methods_estimates$theta_hat),
  na.rm = TRUE
)

axis_min = max(0, axis_range[1]-0.02)
axis_max = min(1, axis_range[2]+0.02)

p_truth = ggplot(all_methods_estimates, aes(x = theta_true, y = theta_hat)) +
  geom_point(size = 0.5, alpha = 0.65) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal(xlim = c(axis_min,axis_max), ylim = c(axis_min,axis_max)) +
  facet_wrap(~method) +
  labs(
    x = "Full Week 1 truth",
    y = "Estimated state rate",
    title = "State estimates against full Week 1 truth"
  ) +
  theme_minimal()

ggsave(
  "results/plot/figure1_estimate_truth.png",
  p_truth,
  width = 8,
  height = 6
)

# Figure 2: state-level estimation error, ordered by sample size
state_order_tbl = state_metrics %>%
  distinct(state, state_abbr, n_state_sample_mean) %>%
  arrange(n_state_sample_mean) %>%
  mutate(order_id = row_number())

error_plot_data = state_metrics %>%
  left_join(
    state_order_tbl %>% select(state, state_abbr, n_state_sample_mean, order_id),
    by = "state"
  )

max_abs_error_line = max(abs(error_plot_data$error_mean), na.rm = TRUE)

p_error = ggplot(
  error_plot_data,
  aes(x = order_id, y = error_mean, color = method, group = method)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_point(size = 1.2, alpha = 0.9) +
  scale_x_continuous(
    breaks = state_order_tbl$order_id,
    labels = state_order_tbl$state_abbr
  ) +
  coord_cartesian(ylim = c(-max_abs_error_line,max_abs_error_line)) +
  labs(
    x = "State, ordered by sample size",
    y = "Estimation error",
    color = "Method",
    title = "State-level estimation error, ordered by sample size"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ggsave(
  "results/plot/figure2_error_by_state_size.png",
  p_error,
  width = 11,
  height = 5
)

# Map data
state_lookup = tibble(
  state = names(state_labels),
  state_abbr = unname(state_labels)
) %>%
  filter(state %in% states) %>%
  left_join(state_region, by = "state_abbr")

state_map = map_data("state") %>%
  as_tibble() %>%
  filter(region %in% state_lookup$region)

map_state_metrics = state_metrics %>%
  left_join(
    state_lookup %>% select(state, region),
    by = "state"
  )

mrp_map_data = map_state_metrics %>%
  filter(method=="MRP") %>%
  select(region, error_mean)

sb_map_data = map_state_metrics %>%
  filter(method=="SB") %>%
  select(region, error_mean)

mrp_sb_map_data = map_state_metrics %>%
  filter(method %in% c("MRP","SB")) %>%
  select(method, state, region, theta_hat_mean) %>%
  pivot_wider(
    names_from = method,
    values_from = theta_hat_mean
  ) %>%
  mutate(diff_mrp_sb = MRP-SB) %>%
  select(region, diff_mrp_sb)

mrp_map = state_map %>%
  left_join(mrp_map_data, by = "region")

sb_map = state_map %>%
  left_join(sb_map_data, by = "region")

mrp_sb_map = state_map %>%
  left_join(mrp_sb_map_data, by = "region")

max_abs_map_error = max(
  abs(c(mrp_map$error_mean,sb_map$error_mean)),
  na.rm = TRUE
)

max_abs_mrp_sb = max(
  abs(mrp_sb_map$diff_mrp_sb),
  na.rm = TRUE
)

# Figure 3: MRP error map
p_map_mrp = ggplot(mrp_map, aes(x = long, y = lat, group = group, fill = error_mean)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.3) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-max_abs_map_error,max_abs_map_error),
    name = "MRP - truth"
  ) +
  labs(
    title = "MRP error against full Week 1 truth"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  "results/plot/figure3_map_mrp_minus_truth.png",
  p_map_mrp,
  width = 8,
  height = 5
)

# Figure 4: SB error map
p_map_sb = ggplot(sb_map, aes(x = long, y = lat, group = group, fill = error_mean)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.3) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-max_abs_map_error,max_abs_map_error),
    name = "SB - truth"
  ) +
  labs(
    title = "SB error against full Week 1 truth"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  "results/plot/figure4_map_sb_minus_truth.png",
  p_map_sb,
  width = 8,
  height = 5
)

# Figure 5: MRP minus SB map
p_map_mrp_sb = ggplot(mrp_sb_map, aes(x = long, y = lat, group = group, fill = diff_mrp_sb)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.3) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-max_abs_mrp_sb,max_abs_mrp_sb),
    name = "MRP - SB"
  ) +
  labs(
    title = "Difference between MRP and SB"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  "results/plot/figure5_map_mrp_minus_sb.png",
  p_map_mrp_sb,
  width = 8,
  height = 5
)

print(summary_metrics)