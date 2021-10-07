##########################################
# Script to calculate diversity metrics
# from dated covid phylogenies
#
# Tom Smith 2021
##########################################

library(pez)
library(phytools)
library(zoo) # messing with dates
library(ggplot2)
library(directlabels)
library(tidyr)
library(dplyr)
library(lubridate)

### Function to calculate and return shape/dispersion metrics
### for a given month
### (now just returning dispersion metric)
monthly_metrics <- function(tree, meta, date){
  # cut to the month needed
  meta_subs <- meta[meta$month == date,]
  # make a matrix of strains by country
  subs_matrix <- as.matrix(table(meta_subs[,c("Area", "Strain")]))
  # cut the tree to just the tips for that month
  subs_tree <- keep.tip(tree, meta_subs$Strain)
  # resolve polytomies
  multi2tree <- multi2di(subs_tree)
  # create community comparative object
  compcomm <- comparative.comm(multi2tree, subs_matrix)
  
  # generate the phylo metrics
  metrics_gamma <- .gamma(compcomm)
  metrics_mntd <- .ses.mntd(compcomm)
  metrics_mpd <- .ses.mpd(compcomm)
  
  # add prefixes to the column names to align with pez outputs
  metrics_mntd <- metrics_mntd %>% 
    rename_with( ~ paste("ses.mntd.", .x, sep = ""))
  
  metrics_mpd <- metrics_mpd %>% 
    rename_with( ~ paste("ses.mpd.", .x, sep = ""))
  
  # bind together into one dataframe
  phylo_metrics <- cbind(metrics_mntd, metrics_mpd)
  phylo_metrics$gamma <- metrics_gamma
  
  # add further data fields and remove row names
  phylo_metrics$area <- row.names(phylo_metrics)
  rownames(phylo_metrics) <- c()
  
  phylo_metrics$month <- as.Date(as.yearmon(date, "%m-%Y"))
  
  # add nunber of tips in tree
  phylo_metrics$tree_size <- length(subs_tree$tip.label)

  return(phylo_metrics)
}

#############################
# European NextStrain tree
#############################

# Load data (note problems with NextStrain branch lengths)
euro_tree <- read.tree("data/nextstrain_ncov_europe_timetree_20210426.tre")
euro_meta <- read.delim("data/nextstrain_ncov_europe_metadata_20210426.tsv", as.is=TRUE, quote="")

# Bin European phylogeny by date, then calculate ses,mntd by country
euro_meta$month <- format(as.Date(euro_meta$Collection.Data), "%m-%Y")
euro_meta$Area <- euro_meta$Country

# try to run a list of months together
euro_month_list <- unique(euro_meta$month)
euro_disp <- do.call(rbind, lapply(euro_month_list, function(x) monthly_metrics(euro_tree, euro_meta, x)))


#############################
# Germany tree
#############################

ger_tree <- read.nexus("data/Germany11kb.nwk.result.date_fixed0.0008.nexus")
ger_meta <- read.csv("data/Germany11kmeta.csv")

# create a column in the meta which matches the tree tip labels
# a concatenation of the strain name and the decimal date?
ger_meta$Strain <- paste(ger_meta$Accession_ID, ger_meta$decimalDate, sep = "_")

# why do we have more accession numbers in the metadata than tips in the tree?
# better cut the metadata down to the tips in the tree and check they match up after
ger_meta <- ger_meta[ger_meta$Strain %in% ger_tree$tip.label,]

ger_meta$month <- format(as.Date(ger_meta$Collection_date), "%m-%Y")

# now need to get places so they're all at the same province level (not sub-divisions)
# can do it by nesting a couple of gsubs

#sub(".*Europe / Germany / ", "", unique(ger_meta$Location))
#unique(sub(" /.*", "", sub(".*Europe / Germany / ", "", unique(ger_meta$Location))))
ger_meta$Area <- sub(" /.*", "", sub(".*Europe / Germany / ", "", ger_meta$Location))

ger_month_list <- unique(ger_meta$month)
ger_disp <- do.call(rbind, lapply(ger_month_list[2:18], function(x) monthly_metrics(ger_tree, ger_meta, x)))


####################################
# write these results out for later
####################################


write.csv(euro_disp, "results/europe_dispersion.csv", row.names = FALSE)
write.csv(ger_disp, "results/germany_dispersion.csv", row.names = FALSE)
