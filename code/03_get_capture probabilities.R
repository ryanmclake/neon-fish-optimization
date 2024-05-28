library(tidyverse)
library(ubms)
library(brms)
library(janitor)
library(tidybayes)

three_pass_data_wide = read_csv("input/three_pass_data_wide_total_fish.csv")
three_pass_model = readRDS(file = "input/three_pass_model.rds")

reach_info = three_pass_data_wide |> ungroup() |>
  mutate(site_int = parse_number(as.character(site_int))) |>
  distinct(site_int, reach_id, year) |>
  separate(reach_id, into = c("site_id", "date", "reach"), remove = F)

per_pass_det = getP(three_pass_model) 

colnames(per_pass_det) = c("first", "second", "third")

per_pass_posts = as_tibble(per_pass_det) |>
  mutate(site_int = row_number()) |>
  pivot_longer(cols = -site_int) |> 
  separate(name, into = c("pass", ".draw"))

total_prob_summary = per_pass_posts |> 
  left_join(reach_info) |>
  group_by(site_int, site_id, reach_id, date, reach, year, .draw) |>
  summarize(value = sum(value)) |>
  group_by(site_int, site_id, year, reach_id, reach, date) |>
  median_qi(value) |>
  mutate(date = ymd(date),
         month = month(date)) |>
  rename(reachID = reach_id)
