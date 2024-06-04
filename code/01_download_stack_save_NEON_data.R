
# clear your environment and memory
rm(list=ls())
gc()

# install packages needed for the analysis
# Install pacman if you haven't already
if (!"pacman" %in% installed.packages()) install.packages("pacman")

# use pacman to obtain packages
pacman::p_load(sf,brms,ubms,fuzzyjoin,arm,reshape2,doBy,lme4,janitor,tidyr,
               readr,dplyr,neonUtilities,minpack.lm,FSA, ggplot2, lubridate,
               tidybayes)

# setup the exclusion "in" function
'%!in%' <- function(x,y)!('%in%'(x,y))

# identify lake sites to exclud in this analyses because we're focusing on streams
lakesites <- c("PRPO", "PRLA", "CRAM", "LIRO", "TOOK")

if(!file.exists("./input/joined_fish_tables.csv")){

  # download the fish data
fish <- neonUtilities::loadByProduct(
  dpID='DP1.20107.001',
  check.size=F,
  site = "all",
  package='expanded',
  startdate = "2015-01-01",
  enddate = as.character(Sys.Date()),
  token = Sys.getenv('NEON_PAT'))

# unpack the fish data
list2env(fish,envir=.GlobalEnv)

rea <- neonUtilities::loadByProduct(
  dpID='DP1.20190.001',
  check.size=F,
  site = "all",
  package='expanded',
  startdate = "2015-01-01",
  enddate = as.character(Sys.Date()),
  token = Sys.getenv('NEON_PAT'))

list2env(rea,envir=.GlobalEnv)

mean_wetted_width = rea$rea_widthFieldData %>%
  clean_names() %>% 
  dplyr::select(site_id, collect_date, wetted_width) %>% 
  mutate(year = year(collect_date),
         month = month(collect_date),
         year_month = paste(year,month, sep = "_")) %>% 
  group_by(site_id) %>% 
  summarize(mean_wetted_width_m = mean(wetted_width, na.rm = T),
            sd_wetted_width_m = sd(wetted_width, na.rm = T))

readr::write_csv(mean_wetted_width, "./input/mean_wetted_width.csv")

# identify the reach lengths for each reach eithin each stream site
reach_lengths = fish$fsh_fieldData |> 
  distinct(reachID, measuredReachLength) |>
  group_by(reachID) |>
  add_tally() |>
  filter(n == 1)

# identify the individual fishing events to join later
EF_events <- fsh_fieldData |>
  dplyr::filter(siteID %!in% lakesites) |>
  dplyr::filter(grepl("fixed",fixedRandomReach)) |>
  dplyr::select(reachID, fixedRandomReach)

# group the bulk count fish data to join with the perFish table
EF_count <- fsh_bulkCount |>
  dplyr::filter(siteID %!in% lakesites) |>
  dplyr::filter(grepl("e-fisher",eventID)) |>
  dplyr::select(eventID, passNumber, bulkFishCount) |>
  dplyr::group_by(eventID, passNumber) |>
  dplyr::summarize(totalFish = sum(bulkFishCount))

# tally up the perFish data to join with the bulkCount data
EF_fish <- fsh_perFish |>
  dplyr::filter(siteID %!in% lakesites) |>
  dplyr::filter(grepl("e-fisher",eventID)) |>
  dplyr::select(eventID, passNumber, taxonID) |>
  dplyr::group_by(eventID, passNumber) |>
  tally()

# join up the events, bulk count, and per fish data tables
EF_joined <- EF_fish |>
  dplyr::left_join(EF_count, by = c("eventID", "passNumber")) |>
  dplyr::mutate(totalFish = ifelse(is.na(totalFish), 0, totalFish),
                summed_fish_per_pass = n + totalFish) |>
  separate(eventID, into = c("siteID", "date", "reach"), remove = F) |>
  dplyr::mutate(reachID = paste0(siteID,".",date,".",reach)) |>
  dplyr::mutate(date = lubridate::ymd(date)) %>%
  dplyr::left_join(., EF_events, by = "reachID", relationship = "many-to-many") %>%
  dplyr::group_by(reachID) |>
  dplyr::filter(n() == 3) %>%
  dplyr::select(-eventID, -n,  -totalFish) |>
  pivot_wider(names_from = passNumber, values_from = summed_fish_per_pass) |>
  na.omit() |>
  dplyr::mutate(increased = case_when(`3` >= `1` ~ "no depletion", 
                                      TRUE ~ "depletion")) %>%
  filter(increased == "depletion") %>%
  left_join(., reach_lengths, by = "reachID") %>%
  dplyr::mutate(year = lubridate::year(date),
                month = lubridate::month(date))

# Write the joined fish table so a user doesn't need to re-download the NEON data
readr::write_csv(EF_joined, "./input/joined_fish_tables.csv")

} else {
  cat("You have already downloaded the NEON fish data and saved the joined file, go to script 02_estimate_abundance.R")
}




