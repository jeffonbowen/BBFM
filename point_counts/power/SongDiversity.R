## Diversity Analysis of the Songbird Data (Top) ====

# These are analysis categories:

  # Total diversity (richness, Shannon and Simpson diversity) over time. 
  # Mean station diversity (richness, Shannon and Simpson diversity) over time
  # Mean station density over time.
  # Mean species density per station before/after.

## Set up analysis environment ----

library(tidyverse)
library(lubridate)
library(vegan)
library(dismo)
library(BiodiversityR)
library(ggplot2)
library(lme4)
library(MuMIn)                # model selection
library(effects)              # for extracting effects terms. Use effects or Alleffects
library(emmeans)              # for estiamting marginal means
library(glmmTMB)
library(MASS)
library(DHARMa)               # Analysis of residuals for mixed models. By simulation.

library(sjplot)
library(iNEXT)
library(pscl)                 # zero-inflated models
library(emmeans)              # for estiamting marginal means
library(PerformanceAnalytics)
library(lattice)
library(gridExtra)

rm(list=ls()) # Clear the workspace

setwd("C:/Users/jeff.matheson/Documents/R Working/SiteC/Song/Diversity")
getwd()


## Import and make species and env matrices ----

dat <- read_csv("C:/Users/jeff.matheson/Documents/R Working/SiteC/Song/songdat.csv")
dat <- filter(dat, ValleyUpstPine == "Y")
dat <- dat %>% rename(BHC20 = BHC20_100)
dat <- dat %>% rename(Station = "Sample Station Label")
dat <- mutate(dat, year = year(dat$Date))
dat$year <- factor(dat$year, levels = c("2006", "2008", "2011", "2012", "2016", "2017", "2018", "2019"))
table(dat$year)
dat <- mutate(dat, MonYear = year(dat$Date)-2016)
dat <- mutate(dat, xMonYear = scale(MonYear))
dat$JDAY <- yday(dat$Date)/365
dat <- mutate(dat, xEasting = scale(Easting))
#
# Set up habitat data
dat$BHC16 <- dat$BHC20
dat$BHC16[dat$BHC16 %in% "FBT"] <- "FBS/FBT"
dat$BHC16[dat$BHC16 %in% "FBS"] <- "FBS/FBT"
dat$BHC16[dat$BHC16 %in% "DSG"] <- "DSG/DSS"
dat$BHC16[dat$BHC16 %in% "DSS"] <- "DSG/DSS"
# Remove ANT
dat <- filter(dat, BHC16 != "ANT")
dat$BHC16 <- factor(dat$BHC16, levels = c("CSH","CYF","CMF","DSH","DYF","DMF","RSH","RYF","RMF",
                                          "FBS/FBT","WGR","WSH","WRI","DSG/DSS","CUL","NVE"))
table(dat$BHC16)
# Habitats BHC7
dat$BHC7 <- dat$BHC20
dat$BHC7[dat$BHC7 %in% "CSH"]	<- "Conif"
dat$BHC7[dat$BHC7 %in% "CYF"]	<- "Conif"
dat$BHC7[dat$BHC7 %in% "CMF"]	<- "Conif"
dat$BHC7[dat$BHC7 %in% "DSH"]	<- "Decid"
dat$BHC7[dat$BHC7 %in% "DYF"]	<- "Decid"
dat$BHC7[dat$BHC7 %in% "DMF"]	<- "Decid"
dat$BHC7[dat$BHC7 %in% "RSH"]	<- "RipForest"
dat$BHC7[dat$BHC7 %in% "RYF"]	<- "RipForest"
dat$BHC7[dat$BHC7 %in% "RMF"]	<- "RipForest"
dat$BHC7[dat$BHC7 %in% "FBS"]	<- "Wetland"
dat$BHC7[dat$BHC7 %in% "FBT"]	<- "Conif"
dat$BHC7[dat$BHC7 %in% "WGR"]	<- "Wetland"
dat$BHC7[dat$BHC7 %in% "WSH"]	<- "Wetland"
dat$BHC7[dat$BHC7 %in% "WRI"]	<- "Wetland"
dat$BHC7[dat$BHC7 %in% "DSG"]	<- "Dry slopes"
dat$BHC7[dat$BHC7 %in% "DSS"]	<- "Dry slopes"
dat$BHC7[dat$BHC7 %in% "CUL"]	<- "Cultiv"
dat$BHC7[dat$BHC7 %in% "NVE"]	<- "Non-veg"
#
## Remove ANT
dat <- filter(dat, BHC7 != "ANT")
dat$BHC7 <- factor(dat$BHC7, levels = c("Conif", "Decid", 'RipForest', "Wetland", "Dry slopes", "Cultiv", "Non-veg"))
table(dat$BHC7)
#
dat$Sample <- paste(dat$Station, dat$Visit, sep = "-")

# Use only data from last 3 years as that is what will be used in the power analysis.
dat <- filter(dat, year %in% c(2017, 2018, 2019))

# Read in list of songbird species
song <- read_csv("C:/Users/jeff.matheson/Documents/R Working/SiteC/Song/songspecies.csv")

#library(tidyselect) # is this needed?
sp <- dplyr::select(dat, Station, Visit, Sample, one_of(song$Code4))
env <- dplyr::select(dat, Station:dur, Sample, year, JDAY, BHC7, BHC20, MonYear, xEasting, xMonYear)

sp.v <- as.data.frame(sp)            # .v stands for vegan format.  
row.names(sp.v) <- sp.v$Sample
sp.v$Station <- NULL
sp.v$Visit <- NULL
sp.v$Sample <- NULL

#Check for dupliates
dups <- dat %>% 
  group_by(Sample) %>% 
  filter(n()>1)


## Basic summary ----

# Place for basic data calcs


## Total diversity (richness, Shannon and Simpson diversity) over time. ----

# Will use iNext I think. Possibly BiodiversityR. Take this on later.


## Mean station diversity (richness, Shannon and Simpson diversity) over time ----

# Won't be able to use densities from QPAD approach since not all species done. 
# Will have to account for survey-level variability as random effects. 

library(BiodiversityR)

rich <- diversityresult(sp.v, index = "richness", method = "each site", sortit = FALSE, digits = 6)
shan <- diversityresult(sp.v, index = "Shannon", method = "each site", sortit = FALSE, digits = 6)
simp <- diversityresult(sp.v, index = "Simpson", method = "each site", sortit = FALSE, digits = 6)
divdat <- bind_cols(env, rich, shan, simp)

hist(divdat$richness)
hist(divdat$Shannon)
hist(divdat$Simpson)

#or
ggplot(divdat, aes(x=richness)) + geom_histogram(binwidth = 1)

library(vcd)
fit <- goodfit(divdat$richness)
fit
summary(fit)
plot(fit)
rootogram(divdat$richness)

#gaus1 <- lmer(richness ~ BHC7 + MonYear + TSSR + JDAY + (1|Station), data = divdat)
#summary(gaus1)

pois0 <- glmmTMB(richness ~ 1 + (1|Station), data = divdat, family = (poisson("log")))
summary(pois0)

pois1 <- glmmTMB(richness ~ MonYear + (1|Station), data = divdat, family = (poisson("log")))
summary(pois1)

pois2 <- glmmTMB(richness ~ BHC7 + (1|Station), data = divdat, family = (poisson("log")))
summary(pois2)

pois3 <- glmmTMB(richness ~ xEasting + (1|Station), data = divdat, family = (poisson("log")))
summary(pois3)

pois4 <- glmer(richness ~ MonYear + BHC7 + (1|Station), data = divdat, family = (poisson("log")))
summary(pois4)

pois5 <- glmmTMB(richness ~ MonYear + xEasting + (1|Station), data = divdat, family = (poisson("log")))
summary(pois5)

pois6 <- glmmTMB(richness ~ MonYear + xEasting + BHC7 + (1|Station), data = divdat, family = (poisson("log")))
summary(pois6)

pois7 <- glmmTMB(richness ~ MonYear + (1|TSSR) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois7)

pois8 <- glmmTMB(richness ~ BHC7 + (1|TSSR) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois8)

pois9 <- glmmTMB(richness ~ xEasting + (1|TSSR) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois9)

pois10 <- glmmTMB(richness ~ MonYear + BHC7 + (1|TSSR) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois10)

pois11 <- glmmTMB(richness ~ MonYear + xEasting + (1|TSSR) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois11)

pois12 <- glmmTMB(richness ~ MonYear + xEasting + BHC7 + (1|TSSR) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois12)

pois13 <- glmmTMB(richness ~ MonYear + (1|JDAY) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois13)

pois14 <- glmmTMB(richness ~ BHC7 + (1|JDAY) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois14)

pois15 <- glmmTMB(richness ~ xEasting + (1|JDAY) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois15)

pois16 <- glmmTMB(richness ~ MonYear + BHC7 + (1|JDAY) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois16)

pois17 <- glmmTMB(richness ~ MonYear + xEasting + (1|JDAY) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois17)

pois18 <- glmmTMB(richness ~ MonYear + xEasting + BHC7 + (1|JDAY) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois18)

pois19 <- glmmTMB(richness ~ MonYear + xEasting + BHC7 + (1|JDAY) + (1|TSSR) + (1|Station), data = divdat, family = (poisson("log")))
summary(pois19)

model.sel(pois0, pois1, pois2, pois3, pois4, pois5, pois6, pois7, pois8, pois9, pois10, pois11, pois12, 
          pois13, pois14, pois15, pois16, pois17, pois18, pois19)

bestmod <- pois19

library(visreg)
visreg(bestmod)

# Plot all effects
ae <- allEffects(bestmod)
plot(ae, residuals="TRUE")
ae

# Or plot individual effects
e <- predictorEffect("BHC7", bestmod)
plot(e)

# Examine Residuals for Gaussian
qqnorm(residuals(bestmod))
qqline(residuals(bestmod))
plot(bestmod,MonYear~resid(.))

#For other model distributions, use DHARMAa (simulation)
res <- simulateResiduals(fittedModel = bestmod)
plot(res)

testDispersion(res)

# Test for overdispersion
deviance(bestmod)/df.residual(bestmod)
# or
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(bestmod)

## For Power Analysis ----

library(simr)

pois6 <- glmer(richness ~ MonYear + xEasting + BHC7 + (1|Station), data = divdat, family = (poisson("log")), 
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(pois6)

ps <- powerSim(pois6, nsim = 20)
ps  


# Extend to 20 years of monitoring
# alpha set at 0.1 to make test more liberal. Would rather have false positive over false negative. 
# Effect set at 20% negative change per year. This arbitrary but would suggest realtively low, 

msim <- pois6
summary(msim)

ps <- powerSim(msim, nsim = 20)
ps  

# Set effect size at -20% over 10 years. Annual change of -2.21% (corrected for compounding)
# Slope for the annual effect size on log scale is -0.0224
fixef(msim)["MonYear"]
fixef(msim)["MonYear"]<- -0.0224
summary(msim)

# Determine how many years of monitoring required to detect 20% change
msim.extend <- extend(msim, along = "MonYear", n=10)
summary(msim.extend)
years.pc <- powerCurve(msim.extend, nsim = 100, alpha = 0.05, breaks = c(3,4,5,6,8,10))
years.out <- summary(years.pc)
years.out

write_csv(years.out, path = paste("C:/Users/jeff.matheson/Documents/R Working/SiteC/Song/PowerOut/rich-power-yearsout-20-05.csv", 
                               sep=""))

years.plot <- ggplot(years.out, aes(x=nlevels, y=mean))+
  geom_line()+ geom_point()+ xlab(label = "Monitoring Year")+ ylab(label = "Power")+
  labs(subtitle="Mean Station Richness")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05)+
  theme_bw()+ 
  scale_y_continuous(breaks = seq(0, 1.0, by=0.2))+
  scale_x_continuous(breaks = seq(0, 10, by=1))+
  theme(axis.title=element_text(size=10))+
  theme(legend.position = "none")+
  annotate(geom="text", x=7, y=0.35, label="Power to detect 2.2% change per year",
         color="black")+
  annotate(geom="text", x=7, y=0.30, label="",
         color="black")
years.plot

ggsave(paste("C:/Users/jeff.matheson/Documents/R Working/SiteC/Song/PowerOut/rich-power-year-20-05.png", sep=""), 
       years.plot)

# Assuming 10 years of monitoring, what is min number of stations needed to detect 20% change over 10 year?
msim10 <- extend(msim, along = "MonYear", n=10)
summary(msim10)
stations.pc <- powerCurve(msim.extend, along = "Station", nsim = 75, alpha = 0.05)
stations.out <- summary(stations.pc)
stations.out

write_csv(stations.out, path = paste("C:/Users/jeff.matheson/Documents/R Working/SiteC/Song/PowerOut/rich-power-stationsout-20-05.csv", sep=""))

stations.plot <- ggplot(stations.out, aes(x=nlevels, y=mean))+
  geom_line()+ geom_point()+ xlab(label = "Number of Stations")+ ylab(label = "Power")+
  labs(subtitle="Mean Station Richness")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05)+
  scale_y_continuous(breaks = seq(0, 1.0, by=0.2))+
  scale_x_continuous(breaks = seq(0, 300, by=50))+
  theme_bw()+ 
  #theme(axis.title=element_text(size=10))+
  theme(legend.position = "none")+
  annotate(geom="text", x=175, y=0.35, label="Power to detect 20% change in 10 years",
         color="black")+
  annotate(geom="text", x=175, y=0.30, label="at varying number of survey stations",
           color="black")
stations.plot

ggsave(paste("C:/Users/jeff.matheson/Documents/R Working/SiteC/Song/PowerOut/rich-power-stations-50-05.png", sep=""), stations.plot)


