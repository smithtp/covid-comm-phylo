##########################################
# Script to plot diversity metrics
# from dated covid phylogenies through time
#
# Tom Smith 2021
##########################################

library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)

setwd("~/Documents/covid-comm-phylo/")

# ggplotting theme:
main_theme <- theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.text.y = element_text(size = 16),
                    axis.text.x = element_text(size = 16, angle = 90),
                    axis.title.y = element_text(size = 18),
                    axis.title.x = element_blank())

euro_disp <- read.csv("results/europe_dispersion.csv")
ger_disp <- read.csv("results/germany_dispersion.csv")

# .obs = observed metric in community
# .rand.mean = mean of null communities (obtained by randomising - can set different ways of producing null model though)
# .obs.z = Standardised effect size (mpd.z = -NRI; mntd.z = -NTI)
# NRI = "Net related index"
# NTI = "Nearest taxon index"
#
# Cooper 2008 Proc B:
# If competition affects community membership, then species in a community will be more distantly related
# than expected by chance (phylogenetic overdispersion; Webb et al. 2002). If, conversely, community membership 
# is determined by habitat filtering, the species within a community will be more closely related than expected 
# by chance (phylogenetic clustering; Webb et al. 2002). Finally, if community assembly is not strongly 
# influenced by phylogeny, or if multiple factors oppose and nullify each other, community lists will be randomly
# assembled with respect to phylogeny (Helmus et al. 2007).
#
# clustering = more species per genus than expected = positive NRI & NTI (so negative ses.mntd.obs.z)
# overdispersed = less species per genus than expected = negative NRI & NTI (so positive ses.mntd.obs.z)

###################
# Make the plots!
###################

europe_mntd <- ggplot(euro_disp[euro_disp$area != "United Kingdom",], aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa)) +
  geom_point(fill = "black", shape = 21, alpha = 0.5) +
  geom_point(data = euro_disp[euro_disp$area == "United Kingdom",], fill = "red", shape = 21, alpha= 0.6) +
  scale_x_date(date_breaks = "3 months") +
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
europe_mntd


europe_mpd <- ggplot(euro_disp[euro_disp$area != "United Kingdom",], aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mntd.ntaxa)) +
  geom_point(fill = "black", shape = 21, alpha = 0.5) +
  geom_point(data = euro_disp[euro_disp$area == "United Kingdom",], fill = "red", shape = 21, alpha= 0.6) +
  scale_x_date(date_breaks = "3 months") +
  labs(size = "number of sequences", y = "Net Related Index (NRI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
europe_mpd


ger_mntd <- ggplot(ger_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa)) +
  geom_point(fill = "black", shape = 21, alpha = 0.5) +
  scale_x_date(date_breaks = "3 months") +
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
ger_mntd

ger_mpd <- ggplot(ger_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa)) +
  geom_point(fill = "black", shape = 21, alpha = 0.5) +
  scale_x_date(date_breaks = "3 months") +
  labs(size = "number of sequences", y = "Net Related Index (NRI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
ger_mpd

# save the plots

ggsave("results/europe_mntd_z.svg", europe_mntd)
ggsave("results/europe_mpd_z.svg", europe_mpd)
ggsave("results/ger_mntd_z.svg", ger_mntd)
ggsave("results/ger_mpd_z.svg", ger_mpd)