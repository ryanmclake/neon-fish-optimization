
files <- list.files("./output/1_site", recursive = T)

(f <- file.path("./output/1_site", c(files)))

sites <- c(unique(EF_joined$siteID))

d1 <- lapply(f, read_csv)

for( i in seq_along(d1)){
  d1[[i]]$sites <- rep(sites[i],nrow(d1[[i]]))
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

for(i in 1:length(sites)){
  h <- d2 %>% filter(sites == sites[i])
  
  j <- left_join(d5, h, by = "n.years") %>%
    dplyr::select(n.years, removal.05) %>%
    mutate(removal.05 = imputeTS::na_interpolation(removal.05, method = "linear")) %>%
    mutate(site = sites[i]) %>%
    dplyr::filter(abs(removal.05 - 0.7) == min(abs(removal.05 - 0.7)))
  
  yearly[[i]] <- j
  
}

one_reach <- do.call("rbind", yearly) %>%
  mutate(type = "one_reach") %>%
  summarize(years = median(n.years),
            sd_years = sd(n.years))





files <- list.files("./output/3_site", recursive = T)

(f <- file.path("./output/3_site", c(files)))

sites <- c(unique(EF_joined$siteID))

d1 <- lapply(f, read_csv)

for( i in seq_along(d1)){
  d1[[i]]$sites <- rep(sites[i],nrow(d1[[i]]))
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

for(i in 1:length(sites)){
  h <- d2 %>% ungroup() %>% filter(sites == sites[i])
  
  j <- left_join(d5, h, by = "n.years") %>%
    dplyr::select(n.years, removal.05) %>%
    mutate(removal.05 = imputeTS::na_interpolation(removal.05, method = "linear")) %>%
    mutate(site = sites[i]) %>%
    dplyr::filter(abs(removal.05 - 0.7) == min(abs(removal.05 - 0.7)))
  
  yearly[[i]] <- j
  
}

three_reach <- do.call("rbind", yearly) %>%
  mutate(type = "three_reach") %>%
  group_by(site) %>%
  summarize(years = mean(n.years),
            sd_years = sd(n.years))

