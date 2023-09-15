##########################################
# Script to calculate diversity metrics
# from dated covid phylogenies
#
# Tom Smith 2021
# Altered by JMH to include more countries & D
##########################################

library(pez)
library(phytools)
library(zoo) # messing with dates
library(ggplot2)
library(directlabels)
library(tidyr)
library(dplyr)
library(lubridate)
library(gdata)

### Function to calculate and return shape/dispersion metrics
### for a given month
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
  ## Fix the nodel labels and negative length
  compcomm$phy$node.label <- apply(expand.grid(letters, letters, letters), 1, paste, collapse = "")[seq_len(compcomm$phy$Nnode)]
  compcomm$phy$edge.length[compcomm$phy$edge.length <=0] <- .0001
  
  # generate the phylo metrics
  metrics_gamma <- .gamma(compcomm)
  metrics_mntd <- .ses.mntd(compcomm)
  metrics_mpd <- .ses.mpd(compcomm)
  metrics_d <- .d(compcomm)
  
  # add prefixes to the column names to align with pez outputs
  metrics_mntd <- metrics_mntd %>% 
    rename_with( ~ paste("ses.mntd.", .x, sep = ""))
  
  metrics_mpd <- metrics_mpd %>% 
    rename_with( ~ paste("ses.mpd.", .x, sep = ""))
  
  # bind together into one dataframe
  phylo_metrics <- cbind(metrics_mntd, metrics_mpd)
  phylo_metrics$gamma <- metrics_gamma
  phylo_metrics$d <- metrics_d
  
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
## This line takes a while to run - if need to rerun then best on roy
euro_disp <- do.call(rbind, lapply(euro_month_list, function(x) monthly_metrics(euro_tree, euro_meta, x)))
## save
write.csv(euro_disp, "europe_dispersion.csv", row.names = FALSE)

#############################
# Germany tree
#############################

ger_tree <- read.nexus("data/Germany2_10k.nwk.treefile.result.date.nexus")
ger_meta <- read.csv("data/Germany2_meta.csv")

# create a column in the meta which matches the tree tip labels
# a concatenation of the strain name and the decimal date?
ger_meta$Strain <- ger_meta$Accession_ID

# why do we have more accession numbers in the metadata than tips in the tree?
# better cut the metadata down to the tips in the tree and check they match up after
ger_meta <- ger_meta[ger_meta$Strain %in% ger_tree$tip.label,]

ger_meta$month <- format(as.Date(ger_meta$Collection_date), "%m-%Y")

# now need to get places so they're all at the same province level (not sub-divisions)
# can do it by nesting a couple of gsubs

#sub(".*Europe / Germany / ", "", unique(ger_meta$Location))
#unique(sub(" /.*", "", sub(".*Europe / Germany / ", "", unique(ger_meta$Location))))
ger_meta$Area <- sub(" /.*", "", sub(".*Europe / Germany / ", "", ger_meta$Location))
## check regions 
unique(ger_meta$Area)
ger_month_list <- unique(ger_meta$month)
## check months - 01-2020 only has 6 - remove it
table(ger_meta$month)
ger_meta <- filter(ger_meta, month != "01-2020")
ger_month_list <- unique(ger_meta$month)
## save image
save.image("data/ger_base.RData")
## Run the below on roy if needed
ger_disp <- do.call(rbind, lapply(ger_month_list, function(x) monthly_metrics(ger_tree, ger_meta, x)))
## save
write.csv(ger_disp, "data/germany_dispersion2.csv", row.names = FALSE)


#############################
# France tree
#############################

french_tree <- read.nexus("data/France_10k_clean.nwk.treefile.result.date.nexus")
french_meta <- read.csv("data/France_10k_clean_meta.csv")

# create a column in the meta which matches the tree tip labels
# a concatenation of the strain name and the decimal date?
french_meta$Strain <- french_meta$Accession_ID

# why do we have more accession numbers in the metadata than tips in the tree?
# better cut the metadata down to the tips in the tree and check they match up after
french_meta <- french_meta[french_meta$Strain %in% french_tree$tip.label,]

french_meta$month <- format(as.Date(french_meta$Collection_date), "%m-%Y")

# now need to get places so they're all at the same province level (not sub-divisions)
french_meta$Area <- sub(" /.*", "", sub(".*Europe / France / ", "", french_meta$Location))
## unfortunately some misnames
french_meta$Area <- gsub('Marseille',"Provence-Alpes-Cote d'Azur", french_meta$Area)
french_meta$Area <- gsub('Pays de La Loire',"Pays de la Loire", french_meta$Area)
french_meta$Area <- gsub('Haute-Garonne',"Occitanie", french_meta$Area)


french_month_list <- unique(french_meta$month)
## check the months 
table(french_meta$month)
## so dont want to include "01-2020" because only one sample
french_month_list <- french_month_list[french_month_list != "01-2020"] 
## Run on roy
french_disp <- do.call(rbind, lapply(french_month_list, function(x) monthly_metrics(french_tree, french_meta, x)))
## save
write.csv(french_disp, "data/french_dispersion.csv", row.names = FALSE)


#############################
# Spain tree
#############################

spain_tree <- read.nexus("data/Spain_9k_reduced.nwk.treefile.result.date.nexus")
spain_meta <- read.csv("data/Spain_meta.csv")

# create a column in the meta which matches the tree tip labels
spain_meta$Strain <- spain_meta$Accession_ID

## cut the metadata down to the tips in the tree and check they match up after
spain_meta <- spain_meta[spain_meta$Strain %in% spain_tree$tip.label,]
## add month column
spain_meta$month <- as.Date(spain_meta$Collection_date, format = "%d/%m/%Y")
spain_meta$month <- format(as.Date(spain_meta$month), "%m-%Y")

# now need to get places so they're all at the same province level (not sub-divisions)
spain_meta$Area <- sub(" /.*", "", sub(".*Europe / Spain / ", "", spain_meta$Location))
## check names 
unique(spain_meta$Area)
## So Logrono is a city in La rioja 
spain_meta$Area <- gsub("Logrono","La Rioja", spain_meta$Area)
## Melilla is a city in africa under spanish - keep for now. dont have extremadura
## Keep the canaries for now 
spain_meta$Area <- gsub('Africa',"Canary Islands", spain_meta$Area)

spain_month_list <- unique(spain_meta$month)
## check months 
table(spain_meta$month) ## should be able to use all months
## save image then run in roy
spain_disp <- do.call(rbind, lapply(spain_month_list, function(x) monthly_metrics(spain_tree, spain_meta, x)))
## save 
write.csv(spain_disp, "data/spain_dispersion.csv", row.names = FALSE)


#############################
# Italy tree
#############################

italy_tree <- read.nexus("data/Italy2_10k.nwk.treefile.result.date2.nexus")
italy_meta <- read.csv("data/Italy_meta.csv")

# create a column in the meta which matches the tree tip labels
italy_meta$Strain <- italy_meta$Accession_ID

## cut the metadata down to the tips in the tree and check they match up after
italy_meta <- italy_meta[italy_meta$Strain %in% italy_tree$tip.label,]
## add month column
italy_meta$month <- as.Date(italy_meta$Collection_date, format = "%Y-%m-%d")
italy_meta$month <- format(as.Date(italy_meta$month), "%m-%Y")

# now need to get places so they're all at the same province level (not sub-divisions)
italy_meta$Area <- sub(" /.*", "", sub(".*Europe / Italy / ", "", italy_meta$Location))
## check names 
unique(italy_meta$Area)
## Miss-spelt
italy_meta$Area <- gsub("Trentino Alto Adige","Trentino-Alto Adige", italy_meta$Area)
italy_meta$Area <- gsub('Lombardy',"Lombardia", italy_meta$Area)
italy_meta$Area <- gsub('Sicily',"Sicialia", italy_meta$Area)
italy_meta$Area <- gsub('Sardinia',"Sardegna", italy_meta$Area)
italy_meta$Area <- gsub('Piedmont',"Piemonte", italy_meta$Area)
italy_meta$Area <- gsub('Friuli Venezia Giulia',"Friuli-Venezia Giulia", italy_meta$Area)
italy_meta$Area <- gsub('Emilia - Romagna',"Emilia-Romagna", italy_meta$Area)
## theres no data for valle d'aosta
italy_month_list <- unique(italy_meta$month)
## check months 
table(italy_meta$month) ## first month only has 3, remove it
italy_meta <- filter(italy_meta, month != "01-2020")
italy_month_list <- unique(italy_meta$month)
## save image then run in roy
save.image("data/italy_base.RData")
## run in roy
italy_disp <- do.call(rbind, lapply(italy_month_list, function(x) monthly_metrics(italy_tree, italy_meta, x)))
## save 
write.csv(italy_disp, "data/italy_dispersion.csv", row.names = FALSE)


#### SAVE ####
## Save the altered meta data 
keep(euro_meta, ger_meta, french_meta, spain_meta, italy_meta, sure = T)
save.image("data/metas.RData")
