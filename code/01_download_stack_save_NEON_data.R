library(FSA)
library(tidyverse)
library(minpack.lm)
library(neonUtilities)
library(dplyr)
library(readr)
library(tidyr)
library(janitor)
library(lme4)
library(doBy, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
library(arm, warn.conflicts = FALSE)
library(fuzzyjoin, warn.conflicts = FALSE)
library(ubms)
library(brms)
library(viridis)

predict_function <- function(x, beta, alpha){
  est <- (alpha*x)/(beta+x)
  return(est)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

lakesites <- c("PRPO", "PRLA", "CRAM", "LIRO", "TOOK")

fish <- neonUtilities::loadByProduct(
  dpID='DP1.20107.001',
  check.size=F,
  site = "all",
  package='expanded',
  startdate = "2015-01-01",
  enddate = as.character(Sys.Date()),
  token = Sys.getenv('NEON_PAT'))

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

mean_wetted_widths = rea$rea_widthFieldData |>
  dplyr::select(siteID, collectDate, wettedWidth) |>
  mutate(year = year(collectDate),
         month = month(collectDate),
         year_month = paste(year,month, sep = "_")) |>
  group_by(siteID, year, month) |>
  summarize(mean_wetted_width_m = mean(wettedWidth, na.rm = T),
            sd_wetted_width_m = sd(wettedWidth, na.rm = T))

reach_lengths = fish$fsh_fieldData |> 
  distinct(reachID, measuredReachLength) |>
  group_by(reachID) |>
  add_tally() |>
  filter(n == 1)

EF_events <- fsh_fieldData |>
  dplyr::filter(siteID %!in% lakesites) |>
  dplyr::filter(grepl("fixed",fixedRandomReach)) |>
  dplyr::select(reachID, fixedRandomReach)

EF_count <- fsh_bulkCount |>
  dplyr::filter(siteID %!in% lakesites) |>
  dplyr::filter(grepl("e-fisher",eventID)) |>
  dplyr::select(eventID, passNumber, bulkFishCount) |>
  dplyr::group_by(eventID, passNumber) |>
  dplyr::summarize(totalFish = sum(bulkFishCount))

EF_fish <- fsh_perFish |>
  dplyr::filter(siteID %!in% lakesites) |>
  dplyr::filter(grepl("e-fisher",eventID)) |>
  dplyr::select(eventID, passNumber, taxonID) |>
  dplyr::group_by(eventID, passNumber) |>
  tally()

EF_joined <- EF_fish |>
  dplyr::left_join(EF_count, by = c("eventID", "passNumber")) |>
  dplyr::mutate(totalFish = ifelse(is.na(totalFish), 0, totalFish),
         summed_fish_per_pass = n + totalFish) |>
  separate(eventID, into = c("siteID", "date", "reach"), remove = F) |>
  dplyr::mutate(reachID = paste0(siteID,".",date,".",reach)) |>
  dplyr::mutate(date = ymd(date)) %>%
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
  

