library(tidyverse)
library(maps)

if (!dir.exists("results/plot")) dir.create("results/plot", recursive = TRUE)

B_eval = 10
eval_reps = 1:B_eval

method_levels_plot = c("HT","NSB","SB","MRP","MRP2","MRP2.1","MRP3","MRP4","MRP5")

methods_all = c("HT","NSB","SB","MRP","MRP2","MRP2.1","MRP3","MRP4","MRP5")
methods_focus = c("HT","NSB","SB","MRP2.1","MRP4")

method_colors = c(
  "HT" = "#3C3C3C",
  "NSB" = "#1F77B4",
  "SB" = "#17BECF",
  "MRP" = "#FF7F0E",
  "MRP2" = "#FDB863",
  "MRP2.1" = "#8C564B",
  "MRP3" = "#F46D43",
  "MRP4" = "#D62728",
  "MRP5" = "#9467BD"
)

theme_set(
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
)

read_estimates = function(path) {
  if (file.exists(path)) {
    read_rds(path)
  } else {
    tibble(
      replication = integer(),
      method = character(),
      state = character(),
      theta_hat = numeric(),
      lower = numeric(),
      upper = numeric()
    )
  }
}

hps_population = read_rds("data/data_clean/hps_week1_population.rds")
state_truth = read_rds("data/data_clean/hps_week1_state_truth.rds")
hps_sample_ids = read_rds("results/hps_sample_ids.rds")

states = levels(hps_population$state)
if (is.null(states)) states = sort(unique(as.character(hps_population$state)))

hps_population = hps_population %>%
  mutate(state = factor(state, levels = states))

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

state_truth = state_truth %>%
  mutate(
    state = factor(state, levels = states),
    theta_true = theta_pop
  ) %>%
  select(state, theta_true)

ht_estimates = read_estimates("results/ht/estimates.rds")
nsb_estimates = read_estimates("results/nsb/estimates.rds")
sb_estimates = read_estimates("results/sb/estimates.rds")
mrp_estimates = read_estimates("results/mrp/estimates.rds")
mrp2_estimates = read_estimates("results/mrp2/estimates.rds")
mrp2_1_estimates = read_estimates("results/mrp2.1/estimates.rds")
mrp3_estimates = read_estimates("results/mrp3/estimates.rds")
mrp4_estimates = read_estimates("results/mrp4/estimates.rds")
mrp5_estimates = read_estimates("results/mrp5/estimates.rds")

all_estimates = bind_rows(
  ht_estimates,
  nsb_estimates,
  sb_estimates,
  mrp_estimates,
  mrp2_estimates,
  mrp2_1_estimates,
  mrp3_estimates,
  mrp4_estimates,
  mrp5_estimates
) %>%
  filter(replication %in% eval_reps) %>%
  mutate(
    method = as.character(method),
    method = factor(method, levels = method_levels_plot),
    state = factor(state, levels = states)
  ) %>%
  filter(!is.na(method)) %>%
  left_join(state_truth, by = "state") %>%
  mutate(
    error = theta_hat-theta_true,
    sq_error = error^2,
    covered = if_else(is.na(lower) | is.na(upper), NA, theta_true>=lower & theta_true<=upper),
    state_abbr = recode(as.character(state), !!!state_labels, .default = as.character(state))
  )

sample_state_size = hps_sample_ids %>%
  filter(replication %in% eval_reps) %>%
  left_join(
    hps_population %>% select(SCRAM, state),
    by = "SCRAM"
  ) %>%
  mutate(state = factor(state, levels = states)) %>%
  group_by(replication, state) %>%
  summarise(
    n_state_sample = n(),
    sun_w_state_total = sum(sun_w, na.rm = TRUE),
    .groups = "drop"
  )

all_estimates = all_estimates %>%
  left_join(sample_state_size, by = c("replication","state"))

runtime_chain = read_csv("results/method_chain_runtime.csv", show_col_types = FALSE)

if (!"status" %in% names(runtime_chain)) runtime_chain$status = "ok"
if (!"total_seconds" %in% names(runtime_chain)) runtime_chain$total_seconds = NA_real_
if (!"warmup_seconds" %in% names(runtime_chain)) runtime_chain$warmup_seconds = NA_real_
if (!"sampling_seconds" %in% names(runtime_chain)) runtime_chain$sampling_seconds = NA_real_
if (!"chain" %in% names(runtime_chain)) runtime_chain$chain = NA_integer_

runtime_chain = runtime_chain %>%
  filter(replication %in% eval_reps) %>%
  filter(status=="ok") %>%
  mutate(
    method = as.character(method),
    method = factor(method, levels = method_levels_plot),
    total_seconds = as.numeric(total_seconds),
    warmup_seconds = as.numeric(warmup_seconds),
    sampling_seconds = as.numeric(sampling_seconds),
    total_seconds_clean = if_else(
      is.na(total_seconds),
      warmup_seconds+sampling_seconds,
      total_seconds
    )
  ) %>%
  filter(!is.na(method))

runtime_rep = runtime_chain %>%
  group_by(method, replication) %>%
  summarise(
    runtime_minutes = max(total_seconds_clean, na.rm = TRUE)/60,
    n_chains = n_distinct(chain),
    .groups = "drop"
  )

runtime_summary = runtime_rep %>%
  group_by(method) %>%
  summarise(
    runtime_mean_minutes = mean(runtime_minutes, na.rm = TRUE),
    runtime_median_minutes = median(runtime_minutes, na.rm = TRUE),
    runtime_min_minutes = min(runtime_minutes, na.rm = TRUE),
    runtime_max_minutes = max(runtime_minutes, na.rm = TRUE),
    .groups = "drop"
  )

write_rds(runtime_summary, "results/plot/runtime_summary.rds")
write_csv(runtime_summary, "results/plot/runtime_summary.csv")

method_rep_check = all_estimates %>%
  count(method, replication, name = "n_state_rows") %>%
  arrange(method, replication)

method_summary_check = all_estimates %>%
  group_by(method) %>%
  summarise(
    n_rep = n_distinct(replication),
    n_rows = n(),
    .groups = "drop"
  )

write_rds(method_summary_check, "results/plot/method_summary_check.rds")
write_csv(method_summary_check, "results/plot/method_summary_check.csv")

overall_metrics = all_estimates %>%
  group_by(method) %>%
  summarise(
    mse = mean(sq_error, na.rm = TRUE),
    bias = mean(error, na.rm = TRUE),
    coverage = mean(covered, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    coverage = if_else(is.nan(coverage), NA_real_, coverage)
  ) %>%
  right_join(
    tibble(method = factor(method_levels_plot, levels = method_levels_plot)),
    by = "method"
  ) %>%
  left_join(runtime_summary, by = "method") %>%
  arrange(method)

write_rds(overall_metrics, "results/plot/overall_metrics.rds")
write_csv(overall_metrics, "results/plot/overall_metrics.csv")

state_metrics = all_estimates %>%
  group_by(method, state, state_abbr) %>%
  summarise(
    mse = mean(sq_error, na.rm = TRUE),
    bias = mean(error, na.rm = TRUE),
    coverage = mean(covered, na.rm = TRUE),
    theta_hat_mean = mean(theta_hat, na.rm = TRUE),
    theta_true = first(theta_true),
    n_state_sample_mean = mean(n_state_sample, na.rm = TRUE),
    sun_w_state_total_mean = mean(sun_w_state_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    coverage = if_else(is.nan(coverage), NA_real_, coverage),
    error_mean = theta_hat_mean-theta_true
  )

write_rds(all_estimates, "results/plot/all_methods_estimates.rds")
write_csv(all_estimates, "results/plot/all_methods_estimates.csv")

write_rds(state_metrics, "results/plot/state_metrics.rds")
write_csv(state_metrics, "results/plot/state_metrics.csv")

state_order_tbl = state_metrics %>%
  group_by(state, state_abbr) %>%
  summarise(
    n_state_sample_order = mean(n_state_sample_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(n_state_sample_order) %>%
  mutate(order_id = row_number())

state_metrics = state_metrics %>%
  left_join(
    state_order_tbl %>% select(state, state_abbr, order_id),
    by = c("state","state_abbr")
  )

# Figure 1: Overall metrics for all 9 methods
methods_fig1 = c("HT","NSB","SB","MRP","MRP2","MRP2.1","MRP3","MRP4","MRP5")

format_label_fig1 = function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    abs(x)>=10 ~ sprintf("%.1f", x),
    abs(x)>=1 ~ sprintf("%.2f", x),
    abs(x)>=0.1 ~ sprintf("%.3f", x),
    abs(x)>=0.001 ~ sprintf("%.5f", x),
    TRUE ~ sprintf("%.6f", x)
  )
}

overall_long = overall_metrics %>%
  filter(method %in% methods_fig1) %>%
  select(method, mse, bias, coverage, runtime_mean_minutes) %>%
  pivot_longer(
    cols = c(mse, bias, coverage, runtime_mean_minutes),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    method = factor(method, levels = methods_fig1),
    metric = factor(
      metric,
      levels = c("mse","bias","coverage","runtime_mean_minutes"),
      labels = c("MSE","Bias","Coverage","Runtime (minutes)")
    )
  )

label_data_fig1 = overall_long %>%
  filter(!is.na(value)) %>%
  mutate(
    label = format_label_fig1(value),
    label_vjust = if_else(value>=0, -0.35, 1.25)
  )

p1 = ggplot(overall_long, aes(x = method, y = value, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_col(width = 0.75, show.legend = FALSE, na.rm = TRUE) +
  geom_text(
    data = label_data_fig1,
    aes(label = label, vjust = label_vjust),
    size = 2.8,
    show.legend = FALSE
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = method_colors[methods_fig1], drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.15,0.18))) +
  labs(
    x = NULL,
    y = NULL,
    title = "Figure 1.1. Overall metrics for all methods"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = "results/plot/figure1_1_overall_metrics.png",
  plot = p1,
  width = 13,
  height = 8,
  dpi = 300
)

# Figure 1.2: Replication-level MSE and Bias by method
methods_fig1_2 = c("HT","NSB","SB","MRP","MRP2","MRP2.1","MRP3","MRP4","MRP5")

replication_metrics_fig1_2 = all_estimates %>%
  filter(method %in% methods_fig1_2) %>%
  group_by(method, replication) %>%
  summarise(
    mse = mean((theta_hat-theta_true)^2, na.rm = TRUE),
    bias = mean(theta_hat-theta_true, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(mse, bias),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    method = factor(method, levels = methods_fig1_2),
    metric = factor(metric, levels = c("mse","bias"), labels = c("MSE","Bias"))
  )

p1_2 = ggplot(replication_metrics_fig1_2, aes(x = method, y = value, color = method)) +
  geom_hline(
    data = tibble(metric = factor(c("MSE","Bias"), levels = c("MSE","Bias")), yint = c(0,0)),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.3,
    color = "grey40"
  ) +
  geom_point(
    size = 2,
    alpha = 0.8,
    position = position_jitter(width = 0.12, height = 0)
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.45,
    color = "black"
  ) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(values = method_colors[methods_fig1_2], guide = "none", drop = FALSE) +
  labs(
    x = NULL,
    y = NULL,
    title = "Figure 1.2. Replication-level MSE and Bias by method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = "results/plot/figure1_2_replication_mse_bias.png",
  plot = p1_2,
  width = 11,
  height = 8,
  dpi = 300
)

# Figure 2.1: Estimates against truth for all 9 methods
methods_fig2_1 = c("HT","NSB","SB","MRP","MRP2","MRP2.1","MRP3","MRP4","MRP5")

plot_data_fig2_1 = all_estimates %>%
  filter(method %in% methods_fig2_1) %>%
  mutate(method = factor(method, levels = methods_fig2_1))

axis_range = range(
  c(plot_data_fig2_1$theta_true, plot_data_fig2_1$theta_hat),
  na.rm = TRUE
)

axis_min = max(0, axis_range[1]-0.02)
axis_max = min(1, axis_range[2]+0.02)

p2_1 = ggplot(plot_data_fig2_1, aes(x = theta_true, y = theta_hat, color = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "black") +
  geom_point(size = 0.5, alpha = 0.65) +
  facet_wrap(~method, ncol = 3) +
  coord_equal(xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max)) +
  scale_color_manual(values = method_colors[methods_fig2_1], guide = "none", drop = FALSE) +
  labs(
    x = "Full Week 1 truth",
    y = "Estimated state rate",
    title = "Figure 2.1. State estimates against full Week 1 truth"
  )

ggsave(
  filename = "results/plot/figure2_1_estimates_against_truth.png",
  plot = p2_1,
  width = 12,
  height = 10,
  dpi = 300
)

# Figure 2.2: State-level signed error ordered by average sample size, all 9 methods
methods_fig2_2 = c("HT","NSB","SB","MRP","MRP2","MRP2.1","MRP3","MRP4","MRP5")

plot_data_fig2_2 = state_metrics %>%
  filter(method %in% methods_fig2_2) %>%
  mutate(
    method = factor(method, levels = methods_fig2_2)
  )

max_abs_error_fig2_2 = max(abs(plot_data_fig2_2$error_mean), na.rm = TRUE)

p2_2 = ggplot(
  plot_data_fig2_2,
  aes(x = order_id, y = error_mean, color = method, group = method)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_line(linewidth = 0.45, alpha = 0.75) +
  geom_point(size = 1.2, alpha = 0.85) +
  scale_x_continuous(
    breaks = state_order_tbl$order_id,
    labels = state_order_tbl$state_abbr
  ) +
  coord_cartesian(
    ylim = c(-max_abs_error_fig2_2, max_abs_error_fig2_2)
  ) +
  scale_color_manual(values = method_colors[methods_fig2_2], drop = FALSE) +
  labs(
    x = "State, ordered by average state sample size",
    y = "Mean estimation error",
    color = "Method",
    title = "Figure 2.2. State-level signed error ordered by sample size"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )

ggsave(
  filename = "results/plot/figure2_2_state_signed_error_ordered_by_sample_size.png",
  plot = p2_2,
  width = 13,
  height = 6,
  dpi = 300
)

# Figure 3: Runtime by replication
methods_fig3 = c("NSB","SB","MRP","MRP2","MRP2.1","MRP3","MRP4","MRP5")

plot_data_fig3 = runtime_rep %>%
  filter(method %in% methods_fig3) %>%
  mutate(method = factor(method, levels = methods_fig3))

p3 = ggplot(plot_data_fig3, aes(x = method, y = runtime_minutes, color = method)) +
  geom_point(size = 2, alpha = 0.8, position = position_jitter(width = 0.12, height = 0)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, color = "black") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(values = method_colors[methods_fig3], guide = "none", drop = FALSE) +
  labs(
    x = NULL,
    y = "Runtime (minutes)",
    title = "Figure 3. Runtime by replication"
  )

ggsave(
  filename = "results/plot/figure3_runtime_by_replication.png",
  plot = p3,
  width = 10,
  height = 5,
  dpi = 300
)

# Figure 4.1: 9 methods against truth
methods_fig4_1 = c("HT","NSB","SB","MRP","MRP2","MRP2.1","MRP3","MRP4","MRP5")

state_region_map = tibble(
  state_abbr = state.abb,
  region = tolower(state.name)
) %>%
  bind_rows(
    tibble(state_abbr = "DC", region = "district of columbia")
  )

state_lookup_map = tibble(
  state = names(state_labels),
  state_abbr = unname(state_labels)
) %>%
  filter(state %in% states) %>%
  left_join(state_region_map, by = "state_abbr")

us_map_base = map_data("state") %>%
  as_tibble() %>%
  filter(region %in% state_lookup_map$region) %>%
  left_join(
    state_lookup_map %>% select(state, state_abbr, region),
    by = "region"
  )

truth_state_map = state_metrics %>%
  distinct(state, state_abbr, theta_true) %>%
  mutate(state = as.character(state))

map_diff_truth = state_metrics %>%
  filter(method %in% methods_fig4_1) %>%
  mutate(
    state = as.character(state),
    method = factor(method, levels = methods_fig4_1)
  ) %>%
  select(method, state, state_abbr, theta_hat_mean) %>%
  left_join(
    truth_state_map %>% select(state, theta_true),
    by = "state"
  ) %>%
  mutate(
    diff_from_true = theta_hat_mean-theta_true
  )

map_plot_truth = us_map_base %>%
  mutate(state = as.character(state)) %>%
  left_join(
    map_diff_truth,
    by = "state",
    relationship = "many-to-many"
  ) %>%
  mutate(
    map_group = interaction(method, group, drop = TRUE)
  )

max_abs_diff_truth = max(abs(map_plot_truth$diff_from_true), na.rm = TRUE)

if (!is.finite(max_abs_diff_truth) || max_abs_diff_truth==0) {
  max_abs_diff_truth = 0.001
}

p4_1 = ggplot(
  map_plot_truth,
  aes(x = long, y = lat, group = map_group, fill = diff_from_true)
) +
  geom_polygon(color = "white", linewidth = 0.12) +
  coord_fixed(1.3, clip = "off") +
  facet_wrap(~method, ncol = 3) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs_diff_truth, max_abs_diff_truth),
    na.value = "grey90",
    name = "Estimate - True"
  ) +
  labs(
    title = "Figure 4.1. State maps relative to true values"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = "results/plot/figure4_1_maps_against_truth.png",
  plot = p4_1,
  width = 13,
  height = 10,
  dpi = 300
)

# Figure 4.2: Pairwise method comparison maps
methods_fig4_2_row = c("HT","NSB","SB","MRP","MRP2","MRP2.1","MRP3","MRP4","MRP5")
methods_fig4_2_col = c("HT","NSB","SB","MRP","MRP2","MRP2.1","MRP3","MRP4","MRP5")

estimates_fig4_2 = state_metrics %>%
  filter(method %in% unique(c(methods_fig4_2_row, methods_fig4_2_col))) %>%
  mutate(
    state = as.character(state),
    method = as.character(method)
  ) %>%
  select(method, state, theta_hat_mean)

map_diff_pair = expand_grid(
  method_row = methods_fig4_2_row,
  method_col = methods_fig4_2_col,
  state = as.character(states)
) %>%
  left_join(
    estimates_fig4_2 %>%
      rename(method_row = method, theta_row = theta_hat_mean),
    by = c("method_row","state")
  ) %>%
  left_join(
    estimates_fig4_2 %>%
      rename(method_col = method, theta_col = theta_hat_mean),
    by = c("method_col","state")
  ) %>%
  mutate(
    diff_value = theta_row-theta_col,
    panel_row = factor(method_row, levels = methods_fig4_2_row),
    panel_col = factor(method_col, levels = methods_fig4_2_col)
  )

map_plot_pair = us_map_base %>%
  mutate(state = as.character(state)) %>%
  left_join(
    map_diff_pair,
    by = "state",
    relationship = "many-to-many"
  ) %>%
  mutate(
    map_group = interaction(method_row, method_col, group, drop = TRUE)
  )

max_abs_diff_pair = max(abs(map_plot_pair$diff_value), na.rm = TRUE)

if (!is.finite(max_abs_diff_pair) || max_abs_diff_pair==0) {
  max_abs_diff_pair = 0.001
}

p4_2 = ggplot(
  map_plot_pair,
  aes(x = long, y = lat, group = map_group, fill = diff_value)
) +
  geom_polygon(color = "white", linewidth = 0.08) +
  coord_fixed(1.3, clip = "off") +
  facet_grid(panel_row~panel_col, switch = "y") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs_diff_pair, max_abs_diff_pair),
    na.value = "grey90",
    name = "Row - Column"
  ) +
  labs(
    title = "Figure 4.2. Pairwise state map comparisons"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8, face = "bold"),
    strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),
    strip.placement = "outside",
    legend.position = "right"
  )

ggsave(
  filename = "results/plot/figure4_2_pairwise_maps.png",
  plot = p4_2,
  width = 16,
  height = 16,
  dpi = 300
)

# Figure 5.1: Focus methods, state-level MSE and signed bias against sample size
methods_fig5_1 = c("HT","NSB","SB","MRP2.1","MRP4")
state_label_every_fig5_1 = 1

plot_data_fig5_1 = state_metrics %>%
  filter(method %in% methods_fig5_1) %>%
  mutate(method = factor(method, levels = methods_fig5_1)) %>%
  select(method, state, state_abbr, n_state_sample_mean, mse, bias) %>%
  pivot_longer(
    cols = c(mse, bias),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("mse","bias"),
      labels = c("State-level MSE","State-level bias")
    )
  )

fig5_1_range = plot_data_fig5_1 %>%
  group_by(metric) %>%
  summarise(
    y_min = min(value, na.rm = TRUE),
    y_max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_rng = pmax(y_max-y_min, abs(y_max)*0.10, 1e-6),
    label_y = y_min-0.25*y_rng,
    tick_y0 = y_min-0.08*y_rng,
    tick_y1 = y_min-0.02*y_rng
  )

fig5_1_metric_tbl = tibble(
  metric = factor(levels(plot_data_fig5_1$metric), levels = levels(plot_data_fig5_1$metric))
)

state_label_positions_fig5_1 = state_metrics %>%
  distinct(state, state_abbr, n_state_sample_mean) %>%
  arrange(n_state_sample_mean) %>%
  mutate(label_id = row_number()) %>%
  filter((label_id-1)%%state_label_every_fig5_1==0) %>%
  crossing(fig5_1_metric_tbl) %>%
  left_join(fig5_1_range, by = "metric")

p5_1 = ggplot(
  plot_data_fig5_1,
  aes(x = n_state_sample_mean, y = value, color = method)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_point(size = 1.4, alpha = 0.7) +
  geom_smooth(se = FALSE, linewidth = 0.7, method = "loess", formula = y ~ x) +
  geom_segment(
    data = state_label_positions_fig5_1,
    aes(
      x = n_state_sample_mean,
      xend = n_state_sample_mean,
      y = tick_y0,
      yend = tick_y1
    ),
    inherit.aes = FALSE,
    color = "grey40",
    linewidth = 0.25
  ) +
  geom_text(
    data = state_label_positions_fig5_1,
    aes(
      x = n_state_sample_mean,
      y = label_y,
      label = state_abbr
    ),
    inherit.aes = FALSE,
    angle = 90,
    size = 2.1,
    vjust = 0.5,
    color = "grey20"
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = method_colors[methods_fig5_1], drop = FALSE) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Average state sample size",
    y = NULL,
    color = "Method",
    title = "Figure 5.1. State-level MSE and bias against sample size",
    subtitle = "Bias is signed"
  ) +
  theme(
    plot.margin = margin(10,10,45,10)
  )

ggsave(
  filename = "results/plot/figure5_1_state_mse_bias_against_sample_size.png",
  plot = p5_1,
  width = 11,
  height = 8,
  dpi = 300
)

# Figure 5.2: MSE decomposition by method
methods_fig5_2 = c("HT","NSB","SB","MRP2.1","MRP4")

mse_decomp_state = all_estimates %>%
  filter(method %in% methods_fig5_2) %>%
  mutate(method = factor(method, levels = methods_fig5_2)) %>%
  group_by(method, state, state_abbr) %>%
  summarise(
    mse = mean(sq_error, na.rm = TRUE),
    bias = mean(error, na.rm = TRUE),
    squared_bias = bias^2,
    variance_component = mse-squared_bias,
    n_state_sample_mean = mean(n_state_sample, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    variance_component = pmax(variance_component, 0)
  )

mse_decomp_method = mse_decomp_state %>%
  group_by(method) %>%
  summarise(
    variance_component = mean(variance_component, na.rm = TRUE),
    squared_bias = mean(squared_bias, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(variance_component, squared_bias),
    names_to = "component",
    values_to = "value"
  ) %>%
  mutate(
    method = factor(method, levels = methods_fig5_2),
    component = factor(
      component,
      levels = c("variance_component","squared_bias"),
      labels = c("Variance component","Squared bias")
    )
  )

write_rds(mse_decomp_method, "results/plot/mse_decomp_method.rds")
write_csv(mse_decomp_method, "results/plot/mse_decomp_method.csv")

p5_2 = ggplot(
  mse_decomp_method,
  aes(x = method, y = value, fill = component)
) +
  geom_col(width = 0.75) +
  labs(
    x = NULL,
    y = "Average state-level MSE component",
    fill = NULL,
    title = "Figure 5.2. MSE decomposition by method"
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = "results/plot/figure5_2_mse_decomposition_by_method.png",
  plot = p5_2,
  width = 9,
  height = 6,
  dpi = 300
)

# Figure 5.3: MSE decomposition against sample size
methods_fig5_3 = c("HT","NSB","SB","MRP2.1","MRP4")

plot_data_fig5_3 = mse_decomp_state %>%
  filter(method %in% methods_fig5_3) %>%
  mutate(method = factor(method, levels = methods_fig5_3)) %>%
  select(method, state, state_abbr, n_state_sample_mean, variance_component, squared_bias) %>%
  pivot_longer(
    cols = c(variance_component, squared_bias),
    names_to = "component",
    values_to = "value"
  ) %>%
  mutate(
    component = factor(
      component,
      levels = c("variance_component","squared_bias"),
      labels = c("Variance component","Squared bias")
    )
  )

p5_3 = ggplot(
  plot_data_fig5_3,
  aes(x = n_state_sample_mean, y = value, color = method)
) +
  geom_point(size = 1.3, alpha = 0.65) +
  geom_smooth(se = FALSE, linewidth = 0.7, method = "loess", formula = y ~ x) +
  facet_wrap(~component, scales = "free_y", ncol = 1) +
  scale_color_manual(values = method_colors[methods_fig5_3], drop = FALSE) +
  labs(
    x = "Average state sample size",
    y = NULL,
    color = "Method",
    title = "Figure 5.3. MSE decomposition against sample size"
  )

ggsave(
  filename = "results/plot/figure5_3_mse_decomposition_against_sample_size.png",
  plot = p5_3,
  width = 10,
  height = 8,
  dpi = 300
)

# Figure 5.4: State maps for MRP2.1 and MRP4 against NSB and SB
comparison_pairs_fig5_4 = tribble(
  ~panel,            ~method_a, ~method_b,
  "MRP2.1 vs NSB",   "MRP2.1",  "NSB",
  "MRP2.1 vs SB",    "MRP2.1",  "SB",
  "MRP4 vs NSB",     "MRP4",    "NSB",
  "MRP4 vs SB",      "MRP4",    "SB"
)

state_region_fig5_4 = tibble(
  state_abbr = state.abb,
  region = tolower(state.name)
) %>%
  bind_rows(
    tibble(state_abbr = "DC", region = "district of columbia")
  )

state_lookup_fig5_4 = tibble(
  state = names(state_labels),
  state_abbr = unname(state_labels)
) %>%
  filter(state %in% states) %>%
  left_join(state_region_fig5_4, by = "state_abbr")

state_map_fig5_4 = map_data("state") %>%
  as_tibble() %>%
  filter(region %in% state_lookup_fig5_4$region) %>%
  left_join(
    state_lookup_fig5_4 %>% select(state, state_abbr, region),
    by = "region"
  )

estimates_fig5_4 = state_metrics %>%
  filter(method %in% c("NSB","SB","MRP2.1","MRP4")) %>%
  mutate(
    state = as.character(state),
    method = as.character(method)
  ) %>%
  select(method, state, theta_hat_mean)

map_diff_fig5_4 = comparison_pairs_fig5_4 %>%
  crossing(state = as.character(states)) %>%
  left_join(
    estimates_fig5_4 %>%
      rename(method_a = method, theta_a = theta_hat_mean),
    by = c("method_a","state")
  ) %>%
  left_join(
    estimates_fig5_4 %>%
      rename(method_b = method, theta_b = theta_hat_mean),
    by = c("method_b","state")
  ) %>%
  mutate(
    diff_value = theta_a-theta_b,
    panel = factor(
      panel,
      levels = c(
        "MRP2.1 vs NSB",
        "MRP2.1 vs SB",
        "MRP4 vs NSB",
        "MRP4 vs SB"
      )
    )
  )

map_plot_fig5_4 = state_map_fig5_4 %>%
  mutate(state = as.character(state)) %>%
  left_join(
    map_diff_fig5_4,
    by = "state",
    relationship = "many-to-many"
  ) %>%
  mutate(
    map_group = interaction(panel, group, drop = TRUE)
  )

max_abs_diff_fig5_4 = max(abs(map_plot_fig5_4$diff_value), na.rm = TRUE)

if (!is.finite(max_abs_diff_fig5_4) || max_abs_diff_fig5_4==0) {
  max_abs_diff_fig5_4 = 0.001
}

p5_4 = ggplot(
  map_plot_fig5_4,
  aes(x = long, y = lat, group = map_group, fill = diff_value)
) +
  geom_polygon(color = "white", linewidth = 0.15) +
  coord_fixed(1.3, clip = "off") +
  facet_wrap(~panel, ncol = 2) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs_diff_fig5_4, max_abs_diff_fig5_4),
    na.value = "grey90",
    name = "Method A minus Method B"
  ) +
  labs(
    title = "Figure 5.4. State map differences for MRP2.1 and MRP4 against NSB and SB"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = "results/plot/figure5_4_map_mrp21_mrp4_vs_nsb_sb.png",
  plot = p5_4,
  width = 10,
  height = 8,
  dpi = 300
)

print(method_summary_check)
print(overall_metrics)