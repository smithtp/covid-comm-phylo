##########################################
# Script to reconfigure the diversity metrics
# calculated from phylo_shape.R, ready for
# analysis
#
# Tom Smith 2021 (updated 2023)
# Updated based on JH's exploratory analysis
##########################################

#### SET UP ####
rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(phytools)
library(ape)
library(viridis)
library(zoo) # messing with dates
library(gdata)
## load in altered meta data
load("data/metas.RData")


###############
#### EURO: DATA ####
## EU data is from next strain, the quality of the phylo is lower than those 
## provided by Manon

euro_disp <- read.csv("data/europe_dispersion.csv")
## Checked some of the more extreme dispersion values by plotting on the tree,
## all seem genuine values so kepe them in the data

#### EURO: DISP ~ DATE COLOURED BY VOC ####
## format for merging
euro_meta$month <- format(as.Date(euro_meta$Collection.Data), "%Y-%m-01")
## Calculate the proportion of alpha and delta, AY.43 not a VOC but common strain included for interest
euro_area_VOC_prop <- euro_meta%>%group_by(Area, month)%>%summarise(VOC_prop = sum(PANGO.Lineage=='B.1.1.7')/length(PANGO.Lineage), 
                                                                    VOC_D_prop = sum(PANGO.Lineage=='B.1.617.2')/length(PANGO.Lineage),
                                                                    VOC_AY43_prop = sum(PANGO.Lineage=='AY.43')/length(PANGO.Lineage))
## add to disp df
euro_disp <- merge(euro_area_VOC_prop, euro_disp, by.x = c("Area", "month"), by.y = c("area", "month"))
rm(euro_area_VOC_prop)

#### EURO: NUMBER OF LINEAGES ####
## Want to see how the number of lineages affected the dispersion metrics 
euro_meta <- euro_meta%>%group_by(Area,month)%>%mutate(nlineages = length(unique(PANGO.Lineage)))
euro_nlin <- euro_meta%>%group_by(Area, month)%>%summarise(nlin = length(unique(PANGO.Lineage)))
euro_disp <- merge(euro_nlin, euro_disp, by = c("Area", "month"))
rm(euro_nlin)

#### EURO: SPLIT BY VOC OR DATE ####
## Just useful for trying to see what the trend may have been before this voc emerged
## and how it changed after
euro_disp$V_pres <- ifelse(euro_disp$VOC_prop > 0, 1,0)
euro_disp$prepost <- ifelse(as.Date(euro_disp$month) >= "2020-12-01", 1,0)

#### EURO: DATE SHIFTING ####
## Try shifting the dates so that the areas match up with timing of VOC, so 0 is 
## when the VOC first emerges
euro_disp[with(euro_disp, order("VOC_prop"))]
## find the earliest occurence of VOC
first <- euro_disp%>%group_by(Area)%>%filter(VOC_prop>0)
## Now I need to number the rest
first <- first%>%select(Area, month)%>%group_by(Area)%>%mutate(shift_month = 1:n())
## Want the first one to be 0 
first$shift_month <- first$shift_month -1 
euro_disp <- merge(euro_disp, first, by = c("Area", "month"), all.x = T, na.rm = T)
## Then I want to add negative numbers 
last <- euro_disp%>%group_by(Area)%>%filter(VOC_prop<=0)
last <- last%>%select(Area, month)%>%group_by(Area)%>%mutate(shift_month = n():1)
last$shift_month <- last$shift_month*-1
## merge
euro_disp <- merge(euro_disp, last, by = c("Area", "month"), all.x = T, na.rm = T)
# euro_disp1 <- euro_disp1%>%mutate(shift_month = shift_month.x + shift_month.y)
euro_disp$shift_month <- rowSums(euro_disp[,c("shift_month.x", "shift_month.y")], na.rm = T)
euro_disp <- euro_disp%>%select(-c(shift_month.x, shift_month.y))
rm(first, last)

###############
#### GERMAN (UPDATED): DATA ####
## German data produced by Manon

ger_disp <- read.csv("data/germany_dispersion2.csv")

## Checked some of the more extreme dispersion values by plotting on the tree,
## all seem genuine values so kept them in the data

#### GERMAN: DISP ~ DATE COLOURED BY VOC ####
## format for merging
ger_meta$month <- format(as.Date(ger_meta$Collection_date), "%Y-%m-01")
## Calculate the proportion of alpha and delta, AY.43 not a VOC but common strain included for interest
ger_area_VOC_prop <- ger_meta%>%group_by(Area, month)%>%summarise(VOC_prop = sum(Pango_lineage=='B.1.1.7')/length(Pango_lineage), 
                                                                  VOC_D_prop = sum(Pango_lineage=='B.1.617.2')/length(Pango_lineage),
                                                                  VOC_AY43_prop = sum(Pango_lineage=='AY.43')/length(Pango_lineage))
## add to disp df
ger_disp <- merge(ger_area_VOC_prop, ger_disp, by.x = c("Area", "month"), by.y = c("area", "month"))
rm(ger_area_VOC_prop)

#### GERMAN: VOC VS DISP ####

## Proportion of alpha variant against the dispersion metric coloured by date
## Need to group months for lines
ger_disp$month_group <- as.Date(ger_disp$month)
ger_disp$month_group <- cut.Date(ger_disp$month_group, "3 months", include.lowest = T)

#### GERMAN: NUMBER OF LINEAGES ####
## Want to see how the number of lineages affected the dispersion metrics 
ger_meta <- ger_meta%>%group_by(Area,month)%>%mutate(nlineages = length(unique(Pango_lineage)))
ger_nlin <- ger_meta%>%group_by(Area, month)%>%summarise(nlin = length(unique(Pango_lineage)))
ger_disp <- merge(ger_nlin, ger_disp, by = c("Area", "month"))
rm(ger_nlin)

#### GERMAN: SPLIT BY VOC OR DATE ####
## Just useful for trying to see what the trend may have been before this voc emerged
## and how it changed after
ger_disp$V_pres <- ifelse(ger_disp$VOC_prop > 0, 1,0)
ger_disp$prepost <- ifelse(as.Date(ger_disp$month) >= "2020-12-01", 1,0)

#### GERMAN: DATE SHIFTING ####
ger_disp[with(ger_disp, order("VOC_prop"))]
## find the earliest occurence of VOC
first <- ger_disp%>%group_by(Area)%>%filter(VOC_prop>0)
## Now I need to number the rest
first <- first%>%select(Area, month)%>%group_by(Area)%>%mutate(shift_month = 1:n())
## Want the first one to be 0 
first$shift_month <- first$shift_month -1 
ger_disp <- merge(ger_disp, first, by = c("Area", "month"), all.x = T, na.rm = T)
## Then I want to add negative numbers 
last <- ger_disp%>%group_by(Area)%>%filter(VOC_prop<=0)
last <- last%>%select(Area, month)%>%group_by(Area)%>%mutate(shift_month = n():1)
last$shift_month <- last$shift_month*-1
## merge
ger_disp <- merge(ger_disp, last, by = c("Area", "month"), all.x = T, na.rm = T)
# ger_disp1 <- ger_disp1%>%mutate(shift_month = shift_month.x + shift_month.y)
ger_disp$shift_month <- rowSums(ger_disp[,c("shift_month.x", "shift_month.y")], na.rm = T)
ger_disp <- ger_disp%>%select(-c(shift_month.x, shift_month.y))
rm(first, last)

###############
#### FRENCH: DATA ####
## French data produced by Manon

french_disp <- read.csv("data/french_dispersion.csv")

## Checked some of the more extreme dispersion values by plotting on the tree,
## all seem genuine values so kept them in the data

#### FRENCH: DISP ~ DATE COLOURED BY VOC ####
## format for merging
french_meta$month <- format(as.Date(french_meta$Collection_date), "%Y-%m-01")
## Calculate the proportion of alpha and delta, AY.43 not a VOC but common strain included for interest
french_area_VOC_prop <- french_meta%>%group_by(Area, month)%>%summarise(VOC_prop = sum(Pango_lineage=='B.1.1.7')/length(Pango_lineage), 
                                                                        VOC_D_prop = sum(Pango_lineage=='B.1.617.2')/length(Pango_lineage),
                                                                        VOC_AY43_prop = sum(Pango_lineage=='AY.43')/length(Pango_lineage))
## add to disp df
french_disp <- merge(french_area_VOC_prop, french_disp, by.x = c("Area", "month"), by.y = c("area", "month"))
rm(french_area_VOC_prop)

#### FRENCH: VOC VS DISP ####
## For the french data this bit works better in three month groups 
french_disp$month_group <- as.Date(french_disp$month)
french_disp$month_group <- cut.Date(french_disp$month_group, "3 months", include.lowest = T)

#### FRENCH: NUMBER OF LINEAGES ####
## Want to see how the number of lineages affected the dispersion metrics 
french_meta <- french_meta%>%group_by(Area,month)%>%mutate(nlineages = length(unique(Pango_lineage)))
french_nlin <- french_meta%>%group_by(Area, month)%>%summarise(nlin = length(unique(Pango_lineage)))
french_disp <- merge(french_nlin, french_disp, by = c("Area", "month"))
rm(french_nlin)

#### FRENCH: SPLIT BY VOC OR DATE ####
## Just useful for trying to see what the trend may have been before this voc emerged
## and how it changed after
french_disp$V_pres <- ifelse(french_disp$VOC_prop > 0, 1,0)
french_disp$prepost <- ifelse(as.Date(french_disp$month) >= "2020-12-01", 1,0)

#### FRENCH: DATE SHIFTING ####
french_disp[with(french_disp, order("VOC_prop"))]
## find the earliest occurence of VOC
first <- french_disp%>%group_by(Area)%>%filter(VOC_prop>0)
## Now I need to number the rest
first <- first%>%select(Area, month)%>%group_by(Area)%>%mutate(shift_month = 1:n())
## Want the first one to be 0 
first$shift_month <- first$shift_month -1 
french_disp <- merge(french_disp, first, by = c("Area", "month"), all.x = T, na.rm = T)
## Then I want to add negative numbers 
last <- french_disp%>%group_by(Area)%>%filter(VOC_prop<=0)
last <- last%>%select(Area, month)%>%group_by(Area)%>%mutate(shift_month = n():1)
last$shift_month <- last$shift_month*-1
## merge
french_disp <- merge(french_disp, last, by = c("Area", "month"), all.x = T, na.rm = T)
# french_disp1 <- french_disp1%>%mutate(shift_month = shift_month.x + shift_month.y)
french_disp$shift_month <- rowSums(french_disp[,c("shift_month.x", "shift_month.y")], na.rm = T)
french_disp <- french_disp%>%select(-c(shift_month.x, shift_month.y))
rm(first, last)


###############
#### SPAIN: DATA ####
## Spain data and phylo provided by Manon

spain_disp <- read.csv("data/spain_dispersion.csv")

## Checked some of the more extreme dispersion values by plotting on the tree,
## all seem genuine values so kept them in the data

#### SPAIN: DISP ~ DATE COLOURED BY VOC ####
## format for merging
spain_meta$month <- as.Date(spain_meta$Collection_date, format = "%d/%m/%Y")
spain_meta$month <- format(as.Date(spain_meta$month), "%Y-%m-01")
## Calculate the proportion of alpha and delta, AY.43 not a VOC but common strain included for interest
spain_area_VOC_prop <- spain_meta%>%group_by(Area, month)%>%summarise(VOC_prop = sum(Pango_lineage=='B.1.1.7')/length(Pango_lineage), 
                                                                      VOC_D_prop = sum(Pango_lineage=='B.1.617.2')/length(Pango_lineage),
                                                                      VOC_AY43_prop = sum(Pango_lineage=='AY.43')/length(Pango_lineage))
## add to disp df
spain_disp$month <- as.character(spain_disp$month)
spain_disp <- merge(spain_area_VOC_prop, spain_disp, by.x = c("Area", "month"), by.y = c("area", "month"))
rm(spain_area_VOC_prop)

#### SPAIN: VOC VS DISP ####
## Similar to french, works better in groups 
spain_disp$month_group <- as.Date(spain_disp$month)
spain_disp$month_group <- cut.Date(spain_disp$month_group, "3 months", include.lowest = T)

#### SPAIN: NUMBER OF LINEAGES ####
## Want to see how the number of lineages affected the dispersion metrics 
spain_meta <- spain_meta%>%group_by(Area,month)%>%mutate(nlineages = length(unique(Pango_lineage)))
spain_nlin <- spain_meta%>%group_by(Area, month)%>%summarise(nlin = length(unique(Pango_lineage)))
spain_disp <- merge(spain_nlin, spain_disp, by = c("Area", "month"))
rm(spain_nlin)

#### SPAIN: SPLIT BY VOC OR DATE ####
## Just useful for trying to see what the trend may have been before this voc emerged
## and how it changed after
spain_disp$V_pres <- ifelse(spain_disp$VOC_prop > 0, 1,0)
spain_disp$prepost <- ifelse(as.Date(spain_disp$month) >= "2020-12-01", 1,0)

#### SPAIN: DATE SHIFTING ####
spain_disp[with(spain_disp, order("VOC_prop"))]
## find the earliest occurence of VOC
first <- spain_disp%>%group_by(Area)%>%filter(VOC_prop>0)
## Now I need to number the rest
first <- first%>%select(Area, month)%>%group_by(Area)%>%mutate(shift_month = 1:n())
## Want the first one to be 0 
first$shift_month <- first$shift_month -1 
spain_disp <- merge(spain_disp, first, by = c("Area", "month"), all.x = T, na.rm = T)
## Then I want to add negative numbers 
last <- spain_disp%>%group_by(Area)%>%filter(VOC_prop<=0)
last <- last%>%select(Area, month)%>%group_by(Area)%>%mutate(shift_month = n():1)
last$shift_month <- last$shift_month*-1
## merge
spain_disp <- merge(spain_disp, last, by = c("Area", "month"), all.x = T, na.rm = T)
# spain_disp1 <- spain_disp1%>%mutate(shift_month = shift_month.x + shift_month.y)
spain_disp$shift_month <- rowSums(spain_disp[,c("shift_month.x", "shift_month.y")], na.rm = T)
spain_disp <- spain_disp%>%select(-c(shift_month.x, shift_month.y))
rm(first, last)

###############
#### ITALY: DATA ####
## italy data and phylo provided by Manon

italy_disp <- read.csv("data/italy_dispersion.csv")

## Checked some of the more extreme dispersion values by plotting on the tree,
## all seem genuine values so kept them in the data

#### ITALY: DISP ~ DATE COLOURED BY VOC ####
## format for merging
italy_meta$month <- as.Date(italy_meta$Collection_date, format = "%Y-%m-%d")
italy_meta$month <- format(as.Date(italy_meta$month), "%Y-%m-01")
## Calculate the proportion of alpha and delta, AY.43 not a VOC but common strain included for interest
italy_area_VOC_prop <- italy_meta%>%group_by(Area, month)%>%summarise(VOC_prop = sum(Pango_lineage=='B.1.1.7')/length(Pango_lineage), 
                                                                      VOC_D_prop = sum(Pango_lineage=='B.1.617.2')/length(Pango_lineage),
                                                                      VOC_AY43_prop = sum(Pango_lineage=='AY.43')/length(Pango_lineage))
## add to disp df
italy_disp <- merge(italy_area_VOC_prop, italy_disp, by.x = c("Area", "month"), by.y = c("area", "month"))
rm(italy_area_VOC_prop)

#### ITALY: VOC VS DISP ####
## Similar to french, works better in groups 
italy_disp$month_group <- as.Date(italy_disp$month)
italy_disp$month_group <- cut.Date(italy_disp$month_group, "2 months", include.lowest = T)

#### ITALY: NUMBER OF LINEAGES ####
## Want to see how the number of lineages affected the dispersion metrics 
italy_meta <- italy_meta%>%group_by(Area,month)%>%mutate(nlineages = length(unique(Pango_lineage)))
italy_nlin <- italy_meta%>%group_by(Area, month)%>%summarise(nlin = length(unique(Pango_lineage)))
italy_disp <- merge(italy_nlin, italy_disp, by = c("Area", "month"))
rm(italy_nlin)

#### ITALY: SPLIT BY VOC OR DATE ####
## Just useful for trying to see what the trend may have been before this voc emerged
## and how it changed after
italy_disp$V_pres <- ifelse(italy_disp$VOC_prop > 0, 1,0)
italy_disp$prepost <- ifelse(as.Date(italy_disp$month) >= "2020-12-01", 1,0)

#### ITALY: DATE SHIFTING ####
italy_disp[with(italy_disp, order("VOC_prop"))]
## find the earliest occurence of VOC
first <- italy_disp%>%group_by(Area)%>%filter(VOC_prop>0)
## Now I need to number the rest
first <- first%>%select(Area, month)%>%group_by(Area)%>%mutate(shift_month = 1:n())
## Want the first one to be 0 
first$shift_month <- first$shift_month -1 
italy_disp <- merge(italy_disp, first, by = c("Area", "month"), all.x = T, na.rm = T)
## Then I want to add negative numbers 
last <- italy_disp%>%group_by(Area)%>%filter(VOC_prop<=0)
last <- last%>%select(Area, month)%>%group_by(Area)%>%mutate(shift_month = n():1)
last$shift_month <- last$shift_month*-1
## merge
italy_disp <- merge(italy_disp, last, by = c("Area", "month"), all.x = T, na.rm = T)
# italy_disp1 <- italy_disp1%>%mutate(shift_month = shift_month.x + shift_month.y)
italy_disp$shift_month <- rowSums(italy_disp[,c("shift_month.x", "shift_month.y")], na.rm = T)
italy_disp <- italy_disp%>%select(-c(shift_month.x, shift_month.y))
rm(first, last)

###############
#### SAVE FOR MODELS ####

keep(euro_disp, ger_disp, french_disp, spain_disp, italy_disp,
     sure = T)
save.image("model_dfs.RData")
