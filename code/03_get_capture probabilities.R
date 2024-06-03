# clear your environment and memory
rm(list=ls())
gc()

# read the exported files from 01_download_stack_save_NEON_data.R script
joined_fish <- readr::read_csv("./input/joined_fish_tables.csv")


# put passes in a matrix (for ubms)
three_pass_matrix = joined_fish %>% 
  ungroup() %>%
  dplyr::select(`1`,`2`,`3`) %>% 
  as.matrix()

# assign covariates (for ubms)
three_pass_frame <- unmarkedFrameMPois(three_pass_matrix,
                                       siteCovs=as.data.frame(joined_fish %>% ungroup() %>% dplyr::select(reachID,year)),
                                       type = "removal")
saveRDS(three_pass_frame, file = "./input/three_pass_frame.rds")

# fit model
# this estimates population size and capture efficiency for each reach_id (called site_int here)
three_pass_frame = readRDS(file =  "./input/three_pass_frame.rds")

three_pass_model = stan_multinomPois(formula = ~ reachID ~ reachID + (1|year),
                                     data = three_pass_frame,
                                     chains = 3, iter = 1500)
# save the model
saveRDS(three_pass_model, file = "./input/three_pass_model.rds")

# read the model
three_pass_model = readRDS(file = "input/three_pass_model.rds")

# obtain the reach information
reach_info = joined_fish |> ungroup() |>
  dplyr::mutate(site_int = row_number()) |>
  distinct(site_int, reachID, year) |>
  separate(reachID, into = c("siteID", "date", "reach"), remove = F)

# get the three pass model detection probabilities
per_pass_det = getP(three_pass_model) 

# change column names to match the pass number
colnames(per_pass_det) = c("first", "second", "third")

# pivot the DF
per_pass_posts = as_tibble(per_pass_det) |>
  mutate(site_int = row_number()) |>
  pivot_longer(cols = -site_int) |> 
  separate(name, into = c("pass", ".draw")) 

# summarize to get the median cpture probabilities and the 95% CIs
total_prob_summary = per_pass_posts |> 
  left_join(reach_info) |>
  group_by(site_int, siteID, reachID, date, reach, year, .draw) |>
  summarize(value = sum(value)) |>
  group_by(site_int, siteID, year, reachID, reach, date) |>
  median_qi(value) |>
  mutate(date = ymd(date),
         month = month(date))

readr::write_csv(total_prob_summary, "./input/total_capt_probability.csv")
