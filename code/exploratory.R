## EXPLORATORY PLOTS ##

## Script to make exploratory plots from the dispersion metrics calculated in 
## the phylo_shape.R and data-wrangling.R scripts. Split by country/dataset
## Plots can be seen in the EcoCovid_Notes doc
## Jess Hodge 2021
## Updated Tom Smith 2023

#### SET UP ####
rm(list = ls())
library(ggplot2)
library(tidyr)
## load in altered meta data
load("data/model_dfs.RData")

## set up plots
main_theme <- theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.text.y = element_text(size = 16),
                    axis.text.x = element_text(size = 16, angle = 90),
                    axis.title.y = element_text(size = 18),
                    axis.title.x = element_blank())

main_theme1 <- theme(axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(),
                     axis.text.y = element_text(size = 16),
                     axis.text.x = element_text(size = 16),
                     axis.title.y = element_text(size = 18),
                     axis.title.x = element_text(size = 18))


###############
#### EURO: DATA ####
## EU data is from next strain, the quality of the phylo is lower than those 
  ## provided by Manon

## The non-finite value warnings are just from areas/months that have too few 
  ## strains to calc dispersion
## Plot NRI
ggplot(euro_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)", title = "Europe NRI") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
## NTI
ggplot(euro_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)", title = "Europe NTI") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
## D, only samples larger than 2
ggplot(euro_disp[euro_disp$ses.mntd.ntaxa > 2,], aes(x = as.Date(month), y = d.D, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Phylogenetic signal (D)", title = "Europe D") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme


#### EURO: REGIONAL ####
## Hard to see regions with the EU data - lots of countries, consider filtering for only countries with more than 5 samples 
## NRI
ggplot(euro_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1+
  facet_wrap(~Area)
## NTI
ggplot(euro_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme+
  facet_wrap(~Area)
## D, only samples bigger than 2
ggplot(euro_disp[euro_disp$ses.mntd.ntaxa > 2,], aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1+
  facet_wrap(~Area)


#### EURO: VOC VS DISP ####

## Proportion of alpha variant against the dispersion metric
## NRI coloured by date without lines
ggplot(euro_disp, aes(y = -ses.mpd.mpd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(y = "Net Related Index (NRI)", x = "Prop of Variant of Concern", title = "EU NRI vs VOC") +
  main_theme1
## NRI Separate by date with lines
ggplot(euro_disp, aes(y = -ses.mpd.mpd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(y = "Net Related Index (NRI)", x = "Prop of Variant of Concern", title = "Euro NRI vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1

## NTI - x axis and colour both voc
ggplot(euro_disp, aes(y = -ses.mntd.mntd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(x = "Prop of Variant of Concern", y = "Nearest Taxon Index", title = "Europe NTI vs VOC by Date") +
  main_theme1
## NTI separated by date
ggplot(euro_disp, aes(y = -ses.mntd.mntd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(x = "Prop of Variant of Concern", y = "Nearest Taxon Index", title = "Europe NTI vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1

## D  
ggplot(euro_disp[euro_disp$ses.mpd.ntaxa > 2,], aes(y = d.D, x = VOC_prop, size = ses.mpd.ntaxa, fill = month)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(xlab = "Prop of Variant of Concern", y = "Phylogenetic signal (D)") +
  main_theme1
## D with lines
ggplot(euro_disp[euro_disp$ses.mpd.ntaxa > 2,], aes(y = d.D, x = VOC_prop, size = ses.mpd.ntaxa, fill = month)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(xlab = "Prop of Variant of Concern", y = "Phylogenetic signal (D)", title = "Euro D vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  ylim(-2.5,3)+
  main_theme1


#### EURO: NUMBER OF LINEAGES ####
## Want to see how the number of lineages affected the dispersion metrics 
## NRI
ggplot(euro_disp, aes(y = -ses.mpd.mpd.obs.z, x = nlin, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c() +
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## NTI
ggplot(euro_disp, aes(y = -ses.mntd.mntd.obs.z, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## D
ggplot(euro_disp[euro_disp$ses.mntd.ntaxa > 2,], aes(y = d.D, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1

#### EURO: SPLIT BY VOC OR DATE ####
## Just useful for trying to see what the trend may have been before this voc emerged
  ## and how it changed after
## NRI PRES/ABS
ggplot(euro_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme1
## NRI PRE/POST
ggplot(euro_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

## NTI PRES/ABS
ggplot(euro_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme1
## NTI PRE/POST
ggplot(euro_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

## D PRES/ABS
ggplot(euro_disp[euro_disp$ses.mpd.ntaxa > 2,], aes(x = as.Date(month), y = d.D, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Phylo Signal (D)", title = "Europe D VOC pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme
## VOC PRE/POST
ggplot(euro_disp[euro_disp$ses.mpd.ntaxa > 2,], aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Phylo Signal (D)", title = "Europe D pre/post emergence") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1



#### EURO: DATE SHIFTING ####
## Try shifting the dates so that the areas match up with timing of VOC, so 0 is 
  ## when the VOC first emerges
## Try plotting that out
ggplot(euro_disp, aes(x = shift_month, y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Net Related Index (NRI)", x = "Shifted month", title = 'Euro NRI shifted') +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## vs old
ggplot(euro_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Net Related Index (NRI)", x = "month", title = 'Euro NRI') +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## Okay so there is a difference


###############
#### GERMAN (UPDATED): DATA ####
## German data produced by Manon
## Plot NRI
ggplot(ger_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c() +
  labs(size = "number of seqs", y = "Net Related Index (NRI)", 
       title = "German NRI", fill = "Alpha prop") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
## NTI
ggplot(ger_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c() +
  labs(size = "number of seqs", y = "Nearest Taxon Index (NTI)", 
       title = "German NTI", fill = "Alpha prop") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
## D, only samples larger than 5
ggplot(ger_disp[ger_disp$ses.mntd.ntaxa > 10,], aes(x = as.Date(month), y = d.D, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Phylogenetic signal (D)", 
       title = "German D", fill = "Alpha prop") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme

#### GERMAN: REGIONAL ####
## Hard to see regions  
## NRI
ggplot(ger_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1+
  facet_wrap(~Area)
## NTI
ggplot(ger_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme+
  facet_wrap(~Area)
## D, only samples bigger than 5
ggplot(ger_disp[ger_disp$ses.mntd.ntaxa > 5,], aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1+
  facet_wrap(~Area)

#### GERMAN: VOC VS DISP ####

## Proportion of alpha variant against the dispersion metric coloured by date
## Need to group months for lines
## NRI
ggplot(ger_disp, aes(y = -ses.mpd.mpd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(y = "Net Related Index (NRI)", x = "Prop of Variant of Concern", title = "German NRI vs VOC") +
  main_theme1
## NRI with lines 
ggplot(ger_disp, aes(y = -ses.mpd.mpd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(y = "Net Related Index (NRI)", x = "Prop of Variant of Concern", title = "German NRI vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## hard to see any patterns here, think I would need to include most recent VOC 

## NTI 
ggplot(ger_disp, aes(y = -ses.mntd.mntd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(x = "Prop of Variant of Concern", y = "Nearest Taxon Index", title = "German NTI vs VOC") +
  main_theme1
## NTI with lines
ggplot(ger_disp, aes(y = -ses.mntd.mntd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(x = "Prop of Variant of Concern", y = "Nearest Taxon Index", title = "German NTI vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## D
ggplot(ger_disp[ger_disp$ses.mpd.ntaxa > 5,], aes(y = d.D, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(xlab = "Prop of Variant of Concern", y = "Phylogenetic signal (D)", title = "German D vs VOC") +
  main_theme1
## D with lines
ggplot(ger_disp[ger_disp$ses.mpd.ntaxa > 5,], aes(y = d.D, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(xlab = "Prop of Variant of Concern", y = "Phylogenetic signal (D)", title = "German D vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  ylim(-2.5,3)+
  main_theme1

#### GERMAN: NUMBER OF LINEAGES ####
## Want to see how the number of lineages affected the dispersion metrics 
## NRI
ggplot(ger_disp, aes(y = -ses.mpd.mpd.obs.z, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## NTI
ggplot(ger_disp, aes(y = -ses.mntd.mntd.obs.z, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## NRI
ggplot(ger_disp[ger_disp$ses.mntd.ntaxa > 5,], aes(y = d.D, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1

#### GERMAN: SPLIT BY VOC OR DATE ####
## Just useful for trying to see what the trend may have been before this voc emerged
## and how it changed after
## NRI PRES/ABS
ggplot(ger_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)", title = "German NRI pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme1
## NRI PRE/POST
ggplot(ger_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)", title = "German D pre/post") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

## NTI PRES/ABS
ggplot(ger_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)", title = "German NTI pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme1
## NTI PRE/POST
ggplot(ger_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)", title = "German D pre/post") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

## D PRES/ABS
ggplot(ger_disp[ger_disp$ses.mpd.ntaxa > 5,], aes(x = as.Date(month), y = d.D, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Phylo Signal (D)", title = "German D pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme
## VOC PRE/POST
ggplot(ger_disp[ger_disp$ses.mpd.ntaxa > 5,], aes(x = as.Date(month), y = d.D, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Phylo Signal (D)", title = "German D pre/post") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

#### GERMAN: DATE SHIFTING ####
ggplot(ger_disp, aes(x = shift_month, y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Net Related Index (NRI)", x = "Shifted month", title = 'German NRI shifted') +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## vs old
ggplot(ger_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Net Related Index (NRI)", x = "month", title = 'German NRI') +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1


###############
#### FRENCH: DATA ####
## French data produced by Manon

## Checked some of the more extreme dispersion values by plotting on the tree,
## all seem genuine values so kept them in the data

#### FRENCH: DISP ~ DATE COLOURED BY VOC ####
## Plot NRI
ggplot(french_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Net Related Index (NRI)", title = "French NRI") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
## NTI
ggplot(french_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Nearest Taxon Index (NTI)", title = "French NTI") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
## D, only samples larger than 6
ggplot(french_disp[french_disp$ses.mntd.ntaxa > 6,], aes(x = as.Date(month), y = d.D, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Phylogenetic signal (D)", title = "French D") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme

#### FRENCH: REGIONAL ####
## NRI
ggplot(french_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1+
  facet_wrap(~Area)
## NTI
ggplot(french_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme+
  facet_wrap(~Area)
## D, only samples bigger than 5
ggplot(french_disp[french_disp$ses.mntd.ntaxa > 6,], aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1+
  facet_wrap(~Area)

#### FRENCH: VOC VS DISP ####
## For the french data this bit works better in three month groups 
## Proportion of alpha variant against the dispersion metric coloured by date
## NRI
ggplot(french_disp, aes(y = -ses.mpd.mpd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(y = "Net Related Index (NRI)", x = "Prop of Variant of Concern", title = "French NRI vs VOC") +
  main_theme1
## NRI with lines 
ggplot(french_disp, aes(y = -ses.mpd.mpd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(y = "Net Related Index (NRI)", x = "Prop of Variant of Concern", title = "French NRI vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1

## NTI 
ggplot(french_disp, aes(y = -ses.mntd.mntd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(x = "Prop of Variant of Concern", y = "Nearest Taxon Index", title = "French NTI vs VOC") +
  main_theme1
## NTI with lines
ggplot(french_disp, aes(y = -ses.mntd.mntd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(x = "Prop of Variant of Concern", y = "Nearest Taxon Index", title = "French NTI vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1

## D
ggplot(french_disp[french_disp$ses.mpd.ntaxa > 6,], aes(y = d.D, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(xlab = "Prop of Variant of Concern", y = "Phylogenetic signal (D)", title = "French D vs VOC") +
  main_theme1
## D with lines
ggplot(french_disp[french_disp$ses.mpd.ntaxa > 6,], aes(y = d.D, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(xlab = "Prop of Variant of Concern", y = "Phylogenetic signal (D)", title = "French D vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  ylim(-2.5,3)+
  main_theme1

#### FRENCH: NUMBER OF LINEAGES ####
## Want to see how the number of lineages affected the dispersion metrics 
## NRI
ggplot(french_disp, aes(y = -ses.mpd.mpd.obs.z, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## NTI
ggplot(french_disp, aes(y = -ses.mntd.mntd.obs.z, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## NRI
ggplot(french_disp[french_disp$ses.mntd.ntaxa > 6,], aes(y = d.D, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1

#### FRENCH: SPLIT BY VOC OR DATE ####
## Just useful for trying to see what the trend may have been before this voc emerged
## and how it changed after
## NRI PRES/ABS
ggplot(french_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)", title = "French NRI pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme1
## NRI PRE/POST
ggplot(french_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)", title = "French D pre/post") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

## NTI PRES/ABS
ggplot(french_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)", title = "French NTI pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme1
## NTI PRE/POST
ggplot(french_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)", title = "French D pre/post") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

## D PRES/ABS
ggplot(french_disp[french_disp$ses.mpd.ntaxa > 6,], aes(x = as.Date(month), y = d.D, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Phylo Signal (D)", title = "French D pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme
## VOC PRE/POST
ggplot(french_disp[french_disp$ses.mpd.ntaxa > 6,], aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Phylo Signal (D)", title = "French D pre/post") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

#### FRENCH: DATE SHIFTING ####
ggplot(french_disp, aes(x = shift_month, y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Net Related Index (NRI)", x = "Shifted month", title = 'German NRI shifted') +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## vs old
ggplot(french_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Net Related Index (NRI)", x = "month", title = 'German NRI') +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1


###############
#### SPAIN: DATA ####
## Spain data and phylo provided by Manon
## Checked some of the more extreme dispersion values by plotting on the tree,
## all seem genuine values so kept them in the data

#### SPAIN: DISP ~ DATE COLOURED BY VOC ####
ggplot(spain_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)", title = "Spain NRI") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
## NTI
ggplot(spain_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)", title = "Spain NTI") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
## D, only samples larger than 5
ggplot(spain_disp[spain_disp$ses.mntd.ntaxa > 5,], aes(x = as.Date(month), y = d.D, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Phylogenetic signal (D)", title = "Spain D") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme

#### SPAIN: REGIONAL ####
## NRI
ggplot(spain_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1+
  facet_wrap(~Area)
## NTI
ggplot(spain_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme+
  facet_wrap(~Area)
## D, only samples bigger than 5
ggplot(spain_disp[spain_disp$ses.mntd.ntaxa > 5,], aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1+
  facet_wrap(~Area)

#### SPAIN: VOC VS DISP ####
## month groups not working for spain ... 
## Proportion of alpha variant against the dispersion metric coloured by date
## NRI
ggplot(spain_disp, aes(y = -ses.mpd.mpd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(y = "Net Related Index (NRI)", x = "Prop of Variant of Concern", title = "Spain NRI vs VOC") +
  main_theme1
## NRI with lines 
ggplot(spain_disp, aes(y = -ses.mpd.mpd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(y = "Net Related Index (NRI)", x = "Prop of Variant of Concern", title = "Spain NRI vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1

## NTI 
ggplot(spain_disp, aes(y = -ses.mntd.mntd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(x = "Prop of Variant of Concern", y = "Nearest Taxon Index", title = "Spain NTI vs VOC") +
  main_theme1
## NTI with lines
ggplot(spain_disp, aes(y = -ses.mntd.mntd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(x = "Prop of Variant of Concern", y = "Nearest Taxon Index", title = "Spain NTI vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## D
ggplot(spain_disp[spain_disp$ses.mpd.ntaxa > 5,], aes(y = d.D, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(xlab = "Prop of Variant of Concern", y = "Phylogenetic signal (D)", title = "Spain D vs VOC") +
  main_theme1
## D with lines
ggplot(spain_disp[spain_disp$ses.mpd.ntaxa > 5,], aes(y = d.D, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(xlab = "Prop of Variant of Concern", y = "Phylogenetic signal (D)", title = "Spain D vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  ylim(-2.5,3)+
  main_theme1

#### SPAIN: NUMBER OF LINEAGES ####
## Want to see how the number of lineages affected the dispersion metrics 
## NRI
ggplot(spain_disp, aes(y = -ses.mpd.mpd.obs.z, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## NTI
ggplot(spain_disp, aes(y = -ses.mntd.mntd.obs.z, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## NRI
ggplot(spain_disp[spain_disp$ses.mntd.ntaxa > 5,], aes(y = d.D, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1

#### SPAIN: SPLIT BY VOC OR DATE ####
## NRI PRES/ABS
ggplot(spain_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)", title = "Spain NRI pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme1
## NRI PRE/POST
ggplot(spain_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)", title = "Spain D pre/post") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

## NTI PRES/ABS
ggplot(spain_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)", title = "Spain NTI pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme1
## NTI PRE/POST
ggplot(spain_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)", title = "Spain D pre/post") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

## D PRES/ABS
ggplot(spain_disp[spain_disp$ses.mpd.ntaxa > 5,], aes(x = as.Date(month), y = d.D, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Phylo Signal (D)", title = "Spain D pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme
## VOC PRE/POST
ggplot(spain_disp[spain_disp$ses.mpd.ntaxa > 5,], aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Phylo Signal (D)", title = "Spain D pre/post") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

#### SPAIN: DATE SHIFTING ####
ggplot(spain_disp, aes(x = shift_month, y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Net Related Index (NRI)", x = "Shifted month", title = 'Spain NRI shifted') +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## vs old
ggplot(spain_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Net Related Index (NRI)", x = "month", title = 'Spain NRI') +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1


###############
#### ITALY: DATA ####
## italy data and phylo provided by Manon
## Checked some of the more extreme dispersion values by plotting on the tree,
## all seem genuine values so kept them in the data

#### ITALY: DISP ~ DATE COLOURED BY VOC ####
## Plot NRI
ggplot(italy_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Net Related Index (NRI)", title = "Italy NRI") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
## NTI
ggplot(italy_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Nearest Taxon Index (NTI)", title = "Italy NTI") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme
## D, only samples larger than 5
ggplot(italy_disp[italy_disp$ses.mntd.ntaxa > 5,], aes(x = as.Date(month), y = d.D, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Phylogenetic signal (D)", title = "Italy D") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme

#### ITALY: REGIONAL ####
## NRI
ggplot(italy_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1+
  facet_wrap(~Area)
## NTI
ggplot(italy_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme+
  facet_wrap(~Area)
## D, only samples bigger than 5
ggplot(italy_disp[italy_disp$ses.mntd.ntaxa > 5,], aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mntd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1+
  facet_wrap(~Area)

#### ITALY: VOC VS DISP ####
## month groups not working for italy ... 
## Proportion of alpha variant against the dispersion metric coloured by date
## NRI
ggplot(italy_disp, aes(y = -ses.mpd.mpd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(y = "Net Related Index (NRI)", x = "Prop of Variant of Concern", title = "Italy NRI vs VOC") +
  main_theme1
## NRI with lines - lines dont work cos not enough points for stat_smooth
ggplot(italy_disp, aes(y = -ses.mpd.mpd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(y = "Net Related Index (NRI)", x = "Prop of Variant of Concern", title = "Italy NRI vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## NTI 
ggplot(italy_disp, aes(y = -ses.mntd.mntd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(x = "Prop of Variant of Concern", y = "Nearest Taxon Index", title = "italy NTI vs VOC") +
  main_theme1
## NTI with lines
ggplot(italy_disp, aes(y = -ses.mntd.mntd.obs.z, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(x = "Prop of Variant of Concern", y = "Nearest Taxon Index", title = "Italy NTI vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1

## D
ggplot(italy_disp[italy_disp$ses.mpd.ntaxa > 5,], aes(y = d.D, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(xlab = "Prop of Variant of Concern", y = "Phylogenetic signal (D)", title = "Italy D vs VOC") +
  main_theme1
## D with lines
ggplot(italy_disp[italy_disp$ses.mpd.ntaxa > 5,], aes(y = d.D, x = VOC_prop, size = ses.mpd.ntaxa, fill = month_group)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(xlab = "Prop of Variant of Concern", y = "Phylogenetic signal (D)", title = "Italy D vs VOC by Date") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  ylim(-2.5,3)+
  main_theme1

#### ITALY: NUMBER OF LINEAGES ####
## NRI
ggplot(italy_disp, aes(y = -ses.mpd.mpd.obs.z, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## NTI
ggplot(italy_disp, aes(y = -ses.mntd.mntd.obs.z, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## NRI
ggplot(italy_disp[italy_disp$ses.mntd.ntaxa > 5,], aes(y = d.D, x = nlin, size = ses.mpd.ntaxa,fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  ylab( "Net Related Index (NRI)")+ 
  xlab("Number of lineages present") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1

#### ITALY: SPLIT BY VOC OR DATE ####
## Just useful for trying to see what the trend may have been before this voc emerged
## and how it changed after
## NRI PRES/ABS
ggplot(italy_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)", title = "italy NRI pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme1
## NRI PRE/POST
ggplot(italy_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Net Related Index (NRI)", title = "italy D pre/post") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

## NTI PRES/ABS
ggplot(italy_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)", title = "italy NTI pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme1
## NTI PRE/POST
ggplot(italy_disp, aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Nearest Taxon Index (NTI)", title = "italy D pre/post") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

## D PRES/ABS
ggplot(italy_disp[italy_disp$ses.mpd.ntaxa > 5,], aes(x = as.Date(month), y = d.D, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Phylo Signal (D)", title = "italy D pres/abs") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~V_pres, nrow = 2)+
  main_theme
## VOC PRE/POST
ggplot(italy_disp[italy_disp$ses.mpd.ntaxa > 5,], aes(x = as.Date(month), y = -ses.mntd.mntd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis_c()+
  labs(size = "number of sequences", y = "Phylo Signal (D)", title = "italy D pre/post") +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  facet_wrap(~prepost, nrow = 2)+
  main_theme1

#### ITALY: DATE SHIFTING ####
ggplot(italy_disp, aes(x = shift_month, y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Net Related Index (NRI)", x = "Shifted month", title = 'Italy NRI shifted') +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1
## vs old
ggplot(italy_disp, aes(x = as.Date(month), y = -ses.mpd.mpd.obs.z, size = ses.mpd.ntaxa, fill = VOC_prop)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_fill_viridis_c()+
  labs(size = "number of seqs", y = "Net Related Index (NRI)", x = "month", title = 'Italy NRI') +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  main_theme1

