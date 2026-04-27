library(tidyverse)

states = c(
  "01","04","05","06","08","09","10","11","12","13",
  "16","17","18","19","20","21","22","23","24","25",
  "26","27","28","29","30","31","32","33","34","35",
  "36","37","38","39","40","41","42","44","45","46",
  "47","48","49","50","51","53","54","55","56"
)

county_adjacency = read_delim(
  "data/county_adjacency.txt",
  delim = "|",
  col_types = cols(.default = col_character())
)

state_adjacency = county_adjacency %>%
  transmute(
    state1 = substr(`County GEOID`,1,2),
    state2 = substr(`Neighbor GEOID`,1,2)
  ) %>%
  filter(
    state1 %in% states,
    state2 %in% states,
    state1!=state2
  ) %>%
  mutate(
    edge_start = pmin(state1,state2),
    edge_end = pmax(state1,state2)
  ) %>%
  distinct(edge_start, edge_end) %>%
  transmute(
    state1 = edge_start,
    state2 = edge_end
  ) %>%
  arrange(state1, state2)

write_csv(state_adjacency, "data/state_adjacency.csv")