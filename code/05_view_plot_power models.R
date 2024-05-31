
files <- list.files("./output/1_site", recursive = T)

(f <- file.path("./output/1_site", c(files)))

site <- c("ARIK", "BIGC", "BLDE", "CUPE", "HOPB", "KING", "LECO", "LEWI", "MART", "MAYF", "MCDI", "MCRA",
           "POSE", "PRIN", "SYCA", "WALK", "WLOU")

d1 <- lapply(f, read_csv)

for( i in seq_along(d1)){
  d1[[i]]$sites <- rep(site[i],nrow(d1[[i]]))
}

d2 <- purrr::map_df(d1, ~as.data.frame(.)) %>%
  dplyr::select(-`...1`) %>%
  filter(rhat == -0.05)

d2 %>%
  ggplot(.,aes(n.years, removal.05, color = as.character(sites)))+
  geom_line(lwd = 1)+
  geom_hline(yintercept = 0.7)+
  ylim(c(0,1))+
  ylab("Power")+
  xlab("Years into sampling")+
  labs(title="5% trend with 1 fixed reach")+
  theme_classic()


d5 <- data.frame(n.years = seq(from = 5, to = 30, by  =1)) %>%
  tibble()

yearly <- list()

for(i in 1:length(site)){
  h <- d2 %>% filter(sites == site[i])
  
  j <- left_join(d5, h, by = "n.years") %>%
    dplyr::select(n.years, removal.05) %>%
    mutate(removal.05 = imputeTS::na_interpolation(removal.05, method = "linear")) %>%
    mutate(site = sites[i]) %>%
    dplyr::filter(abs(removal.05 - 0.7) == min(abs(removal.05 - 0.7)))
  
  yearly[[i]] <- j
  
}

one_reach <- do.call("rbind", yearly) %>%
  mutate(type = "one_reach") %>%
  summarize(years = mean(n.years),
            sd_years = sd(n.years))


by_site <- do.call("rbind", yearly) %>%
  mutate(type = "one_reach") %>%
  dplyr::select(-removal.05) %>%
  rename(n.years.1 = n.years)


files <- list.files("./output/3_site", recursive = T)

(f <- file.path("./output/3_site", c(files)))

d1 <- lapply(f, read_csv)

for( i in seq_along(d1)){
  d1[[i]]$sites <- rep(site[i],nrow(d1[[i]]))
}

d2 <- purrr::map_df(d1, ~as.data.frame(.)) %>%
  dplyr::select(-`...1`) %>%
  filter(rhat == -0.05)

d2 %>%
  ggplot(.,aes(n.years, removal.05, color = as.character(sites)))+
  geom_line(lwd = 1)+
  geom_hline(yintercept = 0.7)+
  ylim(c(0,1))+
  ylab("Power")+
  xlab("Years into sampling")+
  labs(title="5% trend with 3 fixed reaches")+
  theme_classic()


d5 <- data.frame(n.years = seq(from = 5, to = 30, by  =1)) %>%
  tibble()

yearly <- list()

for(i in 1:length(site)){
  h <- d2 %>% ungroup() %>% filter(sites == site[i])
  
  j <- left_join(d5, h, by = "n.years") %>%
    dplyr::select(n.years, removal.05) %>%
    mutate(removal.05 = imputeTS::na_interpolation(removal.05, method = "linear")) %>%
    mutate(site = sites[i]) %>%
    dplyr::filter(abs(removal.05 - 0.7) == min(abs(removal.05 - 0.7)))
  
  yearly[[i]] <- j
  
}

three_reach <- do.call("rbind", yearly) %>%
  mutate(type = "three_reach") %>%
  summarize(years = mean(n.years),
            sd_years = sd(n.years))

by_site2 <- do.call("rbind", yearly) %>%
  mutate(type = "three_reach") %>%
  dplyr::select(-removal.05)%>%
  rename(n.years.3 = n.years)

compare <- left_join(by_site, by_site2, by = "site") %>%
  mutate(`Time Difference` = n.years.1 - n.years.3) %>%
  filter(`Time Difference` > 0) %>%
  rename(siteID = site)

sites_in_analysis <- c(compare$siteID)

sites_w_coords <- fsh_fieldData %>%
  dplyr::select(siteID, decimalLatitude, decimalLongitude) %>%
  group_by(siteID) %>%
  slice(1) %>%
  filter(siteID %in% sites_in_analysis) %>%
  left_join(., compare, by = "siteID") %>%
  left_join(., species_richness, by = "siteID") %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")


usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

ggplot(usa) +
  geom_sf(color = "#2b2b2b", fill = "white", size=0.125) +
  geom_sf(data = sites_w_coords, aes(size = `Time Difference`))+
  coord_sf(crs = st_crs("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"), datum = NA) +
  ggthemes::theme_map()

