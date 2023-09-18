##########################################
# Eco-phylogenetic modelling
# Will Pearse - 2022-02-09
# Edited by Jess Hodge and Tom Smith 2022
# Revised by Tom Smith 2023
##########################################

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

########################################
# Headers and data #####################
########################################
# Packages
library(mgcv)

# Load data
load("data/model_dfs.RData")

# Reformat a bit
euro_disp$region <- "europe"; french_disp$region <- "france"; ger_disp$region <- "germany"; spain_disp$region <- "spain"; italy_disp$region <- "italy"
french_disp$month_group <- spain_disp$month_group <- italy_disp$month_group <- ger_disp$month_group<- NULL
data <- rbind(euro_disp, french_disp, ger_disp, spain_disp, italy_disp)
data$doy <- as.numeric(as.Date(data$month))
names(data)[24:26] <- c("d","d1","d0")
data <- data[abs(data$d) < 5,]
data$date <- as.Date(data$month)
data <- na.omit(data)

# Add in the variant stages as binary variables
data$variant <- "baseline"
data$variant[data$date>="2020-12-01"] <- "alpha"
data$variant[data$date>="2021-04-01"] <- "delta"
data$variant <- factor(data$variant)
data$variant <- relevel(data$variant, "baseline")

# Make region factor (for GAMs)
data$f.r <- factor(data$region)

# filter out the "Europe" tree, best to use just the
# specific country ones IMO

data <- data %>% filter(f.r != "europe")

# Load in Colless' Index data 
colless <- readRDS("data/colless.RDS")

########################################
# Through time #########################
########################################
# Make a plot of D through time
pdf("results/figure-D-overview1.pdf")
par(mar=c(2.1,4.1,0.1,2.1))
.wrap.loess <- function(x){
  model <- loess(d ~ doy, data=data[data$region==x,])
  pred <- expand.grid(doy=seq(min(data$doy),max(data$doy),.1))
  pred$pred <- predict(model, pred)
  return(pred)
}
region.cols <- c("france"="red", "germany"="black", "spain"="orange", "italy"="blue")
with(data, plot(d ~ date, type="n", axes=FALSE, xlab="", ylab=""))
axis(1, data$date, format(data$date, "%b %y"), gap.axis=2, ylab="D")
mtext("D", 2, 2, font=2, cex=1.5)
axis(2)
abline(h=0, lwd=3, col="grey50", lty=2)
mtext("Brownian", 2, -2, at=-.15, las=1, padj=1, col="grey50", font=2)
mtext("Brownian", 4, -2, at=-.15, las=1, padj=1, col="grey50", font=2)
abline(h=1, lwd=3, col="grey50", lty=2)
mtext("Random", 2, -2, at=1.35, las=1, padj=1, col="grey50", font=2)
mtext("Random", 4, -2, at=1.35, las=1, padj=1, col="grey50", font=2)
for(reg in names(region.cols))
  lines(.wrap.loess(reg), col=region.cols[reg], lwd=5)
with(data, points(d ~ jitter(as.numeric(as.Date(month))), pch=20, cex=.5, col=region.cols[region]))
arrows(as.numeric(as.Date("2020-12-01")), -3.5, as.numeric(as.Date("2021-12-01")), code=3, angle=90)
arrows(as.numeric(as.Date("2021-04-01")), -3, as.numeric(as.Date("2021-12-01")), code=3, angle=90)
text(as.numeric(as.Date("2021-01-01")), -3.35, "Alpha")
text(as.numeric(as.Date("2021-05-01")), -2.85, "Delta")
legend("top", c("France","Germany","Spain","Italy"), col=region.cols, lwd=3, bty="n", horiz=TRUE, pch=20)
dev.off()

# Make a plot of SESmpd through time
pdf("results/figure-SESmpd-overview1.pdf")
par(mar=c(2.1,4.1,0.1,2.1))
.wrap.loess <- function(x){
  model <- loess(ses.mpd.mpd.obs.z ~ doy, data=data[data$region==x,])
  pred <- expand.grid(doy=seq(min(data$doy),max(data$doy),.1))
  pred$pred <- predict(model, pred)
  return(pred)
}
region.cols <- c("france"="red", "germany"="black", "spain"="orange", "italy"="blue")
with(data, plot(ses.mpd.mpd.obs.z ~ date, type="n", axes=FALSE, xlab="", ylab="", ylim=c(-10,5)))
axis(1, data$date, format(data$date, "%b %y"), gap.axis=2, ylab="D")
mtext(expression(SES[MPD]), 2, 2, font=2, cex=1.5)
axis(2)
abline(h=0, lwd=3, col="grey50", lty=2)
mtext("closely\nrelated", 2, -1, at=-1, las=1, padj=1, col="grey50", font=2)
mtext("distantly\nrelated  ", 2, -1, at=1, las=1, padj=1, col="grey50", font=2)
for(reg in names(region.cols))
  lines(.wrap.loess(reg), col=region.cols[reg], lwd=5)
with(data, points(ses.mpd.mpd.obs.z ~ jitter(as.numeric(as.Date(month))), pch=20, cex=.5, col=region.cols[region]))
arrows(as.numeric(as.Date("2020-12-01")), -9, as.numeric(as.Date("2021-12-01")), code=3, angle=90)
arrows(as.numeric(as.Date("2021-04-01")), -8, as.numeric(as.Date("2021-12-01")), code=3, angle=90)
text(as.numeric(as.Date("2021-01-01")), -8.85, "Alpha")
text(as.numeric(as.Date("2021-05-01")), -7.85, "Delta")
legend("top", c("France","Germany","Spain","Italy"), col=region.cols, lwd=3, bty="n", horiz=TRUE, pch=20)
dev.off()

# Make a plot of SESmntd through time
pdf("results/figure-SESmntd-overview1.pdf")
par(mar=c(2.1,4.1,0.1,2.1))
.wrap.loess <- function(x){
  model <- loess(ses.mntd.mntd.obs.z ~ doy, data=data[data$region==x,])
  pred <- expand.grid(doy=seq(min(data$doy),max(data$doy),.1))
  pred$pred <- predict(model, pred)
  return(pred)
}
region.cols <- c("france"="red", "germany"="black", "spain"="orange", "italy"="blue")
with(data, plot(ses.mntd.mntd.obs.z ~ date, type="n", axes=FALSE, xlab="", ylab="", ylim=c(-10,5)))
axis(1, data$date, format(data$date, "%b %y"), gap.axis=2, ylab="D")
mtext(expression(SES[MNTD]), 2, 2, font=2, cex=1.5)
axis(2)
abline(h=0, lwd=3, col="grey50", lty=2)
mtext("closely\nrelated", 2, -1, at=-1, las=1, padj=1, col="grey50", font=2)
mtext("distantly\nrelated  ", 2, -1, at=1, las=1, padj=1, col="grey50", font=2)
for(reg in names(region.cols))
  lines(.wrap.loess(reg), col=region.cols[reg], lwd=5)
with(data, points(ses.mntd.mntd.obs.z ~ jitter(as.numeric(as.Date(month))), pch=20, cex=.5, col=region.cols[region]))
arrows(as.numeric(as.Date("2020-12-01")), -9, as.numeric(as.Date("2021-12-01")), code=3, angle=90)
arrows(as.numeric(as.Date("2021-04-01")), -8, as.numeric(as.Date("2021-12-01")), code=3, angle=90)
text(as.numeric(as.Date("2021-01-01")), -8.85, "Alpha")
text(as.numeric(as.Date("2021-05-01")), -7.85, "Delta")
legend("top", c("France","Germany","Spain","Italy"), col=region.cols, lwd=3, bty="n", horiz=TRUE, pch=20)
dev.off()

# Make a plot of Colless' Index through time
pdf("results/figure-colless1.pdf")
par(mar=c(2.1,4.1,0.8,2.1))
region.cols <- c("france"="red", "germany"="black", "spain"="orange", "italy"="blue")
with(colless, plot(coll ~ date, type="n", axes=FALSE, xlab="", ylab=""))
axis(1, colless$date, format(colless$date, "%b %y"), gap.axis=2, ylab="Colless' Index")
mtext("Colless' index", 2, 2, font=2, cex=1.5)
axis(2)
for(reg in names(region.cols))
  lines(colless$doy[colless$region==reg], colless$coll[colless$region == reg], col = region.cols[reg], type = "l", lwd = 5)
arrows(as.numeric(as.Date("2020-12-01")), -600, as.numeric(as.Date("2021-12-01")), code=3, angle=90)
arrows(as.numeric(as.Date("2021-04-01")), 2000, as.numeric(as.Date("2021-12-01")), code=3, angle=90)
text(as.numeric(as.Date("2021-01-01")), 200, "Alpha")
text(as.numeric(as.Date("2021-05-01")), 2800, "Delta")
legend("top", c("France","Germany","Spain","Italy"), col=region.cols, lwd=3, bty="n", horiz=TRUE, pch=20)
dev.off()

# Some GAM visulations
par(mfrow=c(2,2))
plot(gam(d ~ s(as.numeric(date), by=f.r), data=data))
plot(gam(ses.mpd.mpd.obs.z ~ s(as.numeric(date), by=f.r)+f.r, data=data))
plot(gam(ses.mntd.mntd.obs.z ~ s(as.numeric(date), by=f.r)+f.r, data=data))

# ... We need some kind of model here?...

########################################
# Impact of variants ###################
########################################
d.gam <- gam(d ~ s(as.numeric(date), by=f.r)+f.r+variant, data=data)
mpd.gam <- gam(ses.mpd.mpd.obs.z ~ s(as.numeric(date), by=f.r)+f.r+variant, data=data)
mntd.gam <- gam(ses.mntd.mntd.obs.z ~ s(as.numeric(date), by=f.r)+f.r+variant, data=data)

summary(d.gam)
summary(mpd.gam)
summary(mntd.gam)
#... smaller change in D, but lots of predictability in SESmntd (32%!)
#... with no significant effect of the variants for SESmntd or D. So there's something about close structure that is different here...

# TS:
# What if instead of adding alpha/delta into the dataframe as present/absent
# we use the proportion of those sequences in the population directly?

d.gam.ts <- gam(d ~ s(as.numeric(date), by=f.r)+f.r+VOC_prop+VOC_D_prop, data=data)
plot(mpd.gam.ts <- gam(ses.mpd.mpd.obs.z ~ s(as.numeric(date), by=f.r)+f.r+VOC_prop+VOC_D_prop, data=data))
mntd.gam.ts <- gam(ses.mntd.mntd.obs.z ~ s(as.numeric(date), by=f.r)+f.r+VOC_prop+VOC_D_prop, data=data)

summary(d.gam.ts)
summary(mpd.gam.ts)
summary(mntd.gam.ts)

# -- Now we do see significant effects of the variants in SESmntd and SESmpd, the the predictability of mpd has gone up