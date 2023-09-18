# make some maps showing d/mtd/mntd at different time points
# pre-alpha
# during alpha
# during delta
library(sf)
load("data/model_dfs.RData")

# reformat data
all_dates <- as.character(seq(as.Date("2020-01-01"), as.Date("2021-11-01"), by = "month"))
fra_regions <- unique(french_disp$Area)
ger_regions <- unique(ger_disp$Area)
spa_regions <- unique(spain_disp$Area)
ita_regions <- unique(italy_disp$Area)

fra_df <- data.frame(month = rep(all_dates, each = length(fra_regions)), Area = rep(fra_regions, length(all_dates)))
ger_df <- data.frame(month = rep(all_dates, each = length(ger_regions)), Area = rep(ger_regions, length(all_dates)))
spa_df <- data.frame(month = rep(all_dates, each = length(spa_regions)), Area = rep(spa_regions, length(all_dates)))
ita_df <- data.frame(month = rep(all_dates, each = length(ita_regions)), Area = rep(ita_regions, length(all_dates)))

french_disp <- left_join(fra_df, french_disp, by = c("month", "Area"))
ger_disp <- left_join(ger_df, ger_disp, by = c("month", "Area"))
spain_disp <- left_join(spa_df, spain_disp, by = c("month", "Area"))
italy_disp <- left_join(ita_df, italy_disp, by = c("month", "Area"))

euro_disp$region <- "europe"; french_disp$region <- "france"; ger_disp$region <- "germany"; spain_disp$region <- "spain"
italy_disp$region <- "italy"
#french_disp$month_group <- spain_disp$month_group <- NULL
data <- rbind(french_disp, ger_disp, spain_disp, italy_disp)
data$doy <- as.numeric(as.Date(data$month))
names(data)[24:26] <- c("d","d1","d0")
# data <- data[abs(data$d) < 5,] # nope - keep all the data!
data$date <- as.Date(data$month)
# data <- na.omit(data) # keep the NAs!!

# Add in the variant stages as binary variables
data$variant <- "baseline"
data$variant[data$date>="2020-12-01"] <- "alpha"
data$variant[data$date>="2021-04-01"] <- "delta"
data$variant <- factor(data$variant)
data$variant <- relevel(data$variant, "baseline")

# Make region factor (for GAMs)
data$f.r <- factor(data$region)


# get some map data to match up with these regions (GID1)

states <- st_read("~/Documents/areadata/data/gis/gadm-states.shp")

# can we just merge the states frame into the data?
# nope, there are some name mismatches (damn)
# fix it
# metricsdata_fra <- unique(data_noeuro[data_noeuro$region == "france",]$Area)
# statesnames_fra <- states[states$GID_0 == "FRA",]$NAME_1
# 
# metricsdata_ger <- unique(data_noeuro[data_noeuro$region == "germany",]$Area)
# statesnames_ger <- states[states$GID_0 == "DEU",]$NAME_1
# 
# metricsdata_spa <- unique(data_noeuro[data_noeuro$region == "spain",]$Area)
# statesnames_spa <- states[states$GID_0 == "ESP",]$NAME_1
# 
# metricsdata_ita <- unique(data[data$region == "italy",]$Area)
# statesnames_ita <- states[states$GID_0 == "ITA",]$NAME_1

# mismatch betwen anglicised names and official names in the shapefile
# need to create a file that deals with this

# metricsnames <- c(metricsdata_fra, metricsdata_ger, metricsdata_spa, NA) # the NA because one less name here for some reason
# statesnames <- c(statesnames_fra, statesnames_ger, statesnames_spa)

# tempdat <- data.frame(metricsnames, statesnames)
# write.csv(tempdat, "data/name-matching.csv", row.names = FALSE)
# mess around to align them

name_matches <- read.csv("data/name-matching.csv")
names(name_matches)[1] <- "Area"

states_small <- states[states$NAME_1 %in% name_matches$statesnames,]
states_small <- merge(states_small, name_matches, by.x = "NAME_1", by.y = "statesnames")
states_small <- states_small[states_small$GID_0 %in% c("FRA", "DEU", "ESP", "ITA"),]

# need to fix this such that we get every region for every date,
# even if that region/date combo didn't have a metric reading

# then merge with the data
combined_data <- left_join(states_small, data, by = "Area")

# exclude the canary islands because that's silly
combined_data <- combined_data[combined_data$Area != "Canary Islands",]

# pull out three interesting dates
plotting_data <- combined_data[combined_data$month %in% c("2020-03-01", "2021-01-01", "2021-06-01"),]

# make some arbitrary places to cutoff the outlier points so that plotting looks sensible
hist(plotting_data$d)
hist(plotting_data$ses.mntd.mntd.obs.z)
hist(plotting_data$ses.mpd.mpd.obs.z)

plotting_data[!is.na(plotting_data$d) &
                plotting_data$d > 10,]$d <- NA
plotting_data[!is.na(plotting_data$d) &
                plotting_data$d < -4,]$d <- NA

plotting_data[!is.na(plotting_data$ses.mpd.mpd.obs.z) &
                plotting_data$ses.mpd.mpd.obs.z < -15,]$ses.mpd.mpd.obs.z <- NA

plotting_theme <- theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())

# plot the sweeps of the variants
alpha_plot <- ggplot(plotting_data) +
  geom_sf(aes(fill = VOC_prop), size = 0.2) +
  scale_fill_viridis_c(name = "Alpha") +
  facet_wrap(~variant) +
  plotting_theme
alpha_plot

delta_plot <- ggplot(plotting_data) +
  geom_sf(aes(fill = VOC_D_prop), size = 0.2) +
  scale_fill_viridis_c(name = "Delta") +
  facet_wrap(~variant) +
  plotting_theme
delta_plot

# plot the community phylogenetic metrics
d_plot <- ggplot(plotting_data) +
  geom_sf(aes(fill = d), size = 0.2) +
  scale_fill_gradientn(name = "d", colours = c("red", "yellow", "white", "blue"), values = c(0, 0.5, 0.6, 1)) +
  facet_wrap(~variant) +
  plotting_theme
d_plot

mntd_plot <- ggplot(plotting_data) +
  geom_sf(aes(fill = ses.mntd.mntd.obs.z), size = 0.2) +
  scale_fill_gradient2(name = "MNTD", low = "red", high = "blue", mid = "white", midpoint = 0) +
  facet_wrap(~variant) +
  plotting_theme
mntd_plot

mpd_plot <- ggplot(plotting_data) +
  geom_sf(aes(fill = ses.mpd.mpd.obs.z), size = 0.2) +
  scale_fill_gradient2(name = "MPD", low = "red", high = "blue", mid = "white", midpoint = 0) +
  facet_wrap(~variant) +
  plotting_theme
mpd_plot

ggsave("results/alpha_plot.svg", alpha_plot, width = 9, height = 3)
ggsave("results/delta_plot.svg", delta_plot, width = 9, height = 3)
ggsave("results/d_plot.svg", d_plot, width = 9, height = 3)
ggsave("results/mntd_plot.svg", mntd_plot, width = 9, height = 3)
ggsave("results/mpd_plot.svg", mpd_plot, width = 9, height = 3)
