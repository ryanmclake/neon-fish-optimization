
# put passes in a matrix (for ubms)
three_pass_matrix = EF_joined %>% 
  ungroup() %>%
  dplyr::select(`1`,`2`,`3`) %>% 
  as.matrix()

# assign covariates (for ubms)
three_pass_frame <- unmarkedFrameMPois(three_pass_matrix,
                                       siteCovs=as.data.frame(EF_joined %>% ungroup() %>% dplyr::select(reachID,year)),
                                       type = "removal")
saveRDS(three_pass_frame, file = "./input/three_pass_frame.rds")

# fit model
# this estimates population size and capture efficiency for each reach_id (called site_int here)
three_pass_frame = readRDS(file =  "./input/three_pass_frame.rds")

three_pass_model = stan_multinomPois(formula = ~ reachID ~ reachID + (1|year),
                                     data = three_pass_frame,
                                     chains = 3, iter = 1500)

saveRDS(three_pass_model, file = "./input/three_pass_model.rds")

three_pass_model = readRDS(file = "input/three_pass_model.rds")

reach_info = EF_joined |> ungroup() |>
  dplyr::mutate(site_int = row_number()) |>
  distinct(site_int, reachID, year) |>
  separate(reachID, into = c("siteID", "date", "reach"), remove = F)

per_pass_det = getP(three_pass_model) 

colnames(per_pass_det) = c("first", "second", "third")

per_pass_posts = as_tibble(per_pass_det) |>
  mutate(site_int = row_number()) |>
  pivot_longer(cols = -site_int) |> 
  separate(name, into = c("pass", ".draw")) 

total_prob_summary = per_pass_posts |> 
  left_join(reach_info) |>
  group_by(site_int, siteID, reachID, date, reach, year, .draw) |>
  summarize(value = sum(value)) |>
  group_by(site_int, siteID, year, reachID, reach, date) |>
  median_qi(value) |>
  mutate(date = ymd(date),
         month = month(date))
