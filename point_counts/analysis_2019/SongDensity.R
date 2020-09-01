###################################
## Site C Songbird Data Analysis ##
###################################

rm(list = ls())      

library(tidyverse)
library(ggplot2)
library(lubridate)
library(lme4)                 # GLMMs
library(emmeans)              # for estiamting marginal means
library(DHARMa)               # Analysis of residuals for mixed models. By simulation.
library(effects)
library(glmmTMB)              # for nb and zi
library(MuMIn)                # model selection
library(vcd)                  # Fits a discrete (count data) distribution for goodness-of-fit tests.
library(broom.mixed)          # Tidy regression tables

setwd("C:/Users/jeff.matheson/Documents/R Working/SiteC/Song")
dat <- read_csv("songdat.csv")
names(dat)
dat <- filter(dat, ValleyUpstPine == "Y")
dat <- dat %>% rename(BHC20 = BHC20_100)
dat <- dat %>% rename(Station = "Sample Station Label")
table(dat$BHC20)

dat <- mutate(dat, year = year(dat$Date))
dat$year <- factor(dat$year, levels = c("2006", "2008", "2011", "2012", "2016", "2017", "2018", "2019"))
table(dat$year)

dat$JDAY <- yday(dat$Date)/365
jdayspring <- yday(as_date("2019-05-17"))
dat <- mutate(dat, surveyjday = yday(dat$Date))
dat <- mutate(dat, DSLS = (yday(dat$Date)-jdayspring)/365)
summary(dat$DSLS)

## Duration and Dist
dat$dist <- 100

## Habitats BHC16
dat$BHC16 <- dat$BHC20
dat$BHC16[dat$BHC16 %in% "FBT"] <- "FBS/FBT"
dat$BHC16[dat$BHC16 %in% "FBS"] <- "FBS/FBT"
dat$BHC16[dat$BHC16 %in% "DSG"] <- "DSG/DSS"
dat$BHC16[dat$BHC16 %in% "DSS"] <- "DSG/DSS"
## Remove ANT
dat <- filter(dat, BHC16 != "ANT")

dat$BHC16 <- factor(dat$BHC16, levels = c("CSH","CYF","CMF","DSH","DYF","DMF","RSH","RYF","RMF",
                                      "FBS/FBT","WGR","WSH","WRI","DSG/DSS","CUL","NVE"))
table(dat$BHC16)
levels(dat$BHC16)

## Habitats BHC7
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
          ## Remove ANT
dat <- filter(dat, BHC7 != "ANT")

dat$BHC7 <- factor(dat$BHC7, levels = c("Conif", "Decid", 'RipForest', "Wetland", "Dry slopes", "Cultiv", "Non-veg"))
table(dat$BHC7)
levels(dat$BHC7)

names(dat)

##########
## QPAD ##
##########

library(mefa4)
library(QPAD)
QPAD::load_BAM_QPAD(3)
"BAM Version="
getBAMversion()

offdat <- dat

################################################################################
#########
spp <- "YEWA"
##########
##############

#offdat$count <- dat[, as.character(spp)]
names(offdat)
offdat$MAXDUR <- offdat$dur
offdat$MAXDIS <- offdat$dist
offdat$JDAY2 <- offdat$JDAY^2
offdat$TSSR2 <- offdat$TSSR^2
offdat$DSLS2 <- offdat$DSLS^2
offdat$LCC4 <- as.character(offdat$Label)
offdat$LCC4[offdat$LCC4 %in% c("Decid", "Mixed")] <- "DecidMixed"
offdat$LCC4[offdat$LCC4 %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- "Open"
offdat$LCC4 <- factor(offdat$LCC4,
                      c("DecidMixed", "Conif", "Open", "Wet"))
offdat$LCC2 <- as.character(offdat$LCC4)
offdat$LCC2[offdat$LCC2 %in% c("DecidMixed", "Conif")] <- "Forest"
offdat$LCC2[offdat$LCC2 %in% c("Open", "Wet")] <- "OpenWet"
offdat$LCC2 <- factor(offdat$LCC2, c("Forest", "OpenWet"))
table(offdat$LCC4, offdat$LCC2)
Xp <- cbind("(Intercept)"=1, as.matrix(offdat[,c("TSSR","JDAY","DSLS","TSSR2","JDAY2","DSLS2")]))
Xq <- cbind("(Intercept)"=1, TREE=offdat$TREE,
            LCC2OpenWet=ifelse(offdat$LCC2=="OpenWet", 1, 0),
            LCC4Conif=ifelse(offdat$LCC4=="Conif", 1, 0),
            LCC4Open=ifelse(offdat$LCC4=="Open", 1, 0),
            LCC4Wet=ifelse(offdat$LCC4=="Wet", 1, 0))
p <- rep(NA, nrow(offdat))
A <- q <- p
cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
mi <- bestmodelBAMspecies(spp, type="BIC")
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
OKp <- rowSums(is.na(Xp2)) == 0
Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
OKq <- rowSums(is.na(Xq2)) == 0
p[!OKp] <- sra_fun(offdat$MAXDUR[!OKp], cf0[1])
unlim <- ifelse(offdat$MAXDIS[!OKq] == Inf, TRUE, FALSE)
A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * offdat$MAXDIS[!OKq]^2)
q[!OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[!OKq], cf0[2]))
phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
p[OKp] <- sra_fun(offdat$MAXDUR[OKp], phi1)
unlim <- ifelse(offdat$MAXDIS[OKq] == Inf, TRUE, FALSE)
A[OKq] <- ifelse(unlim, pi * tau1^2, pi * offdat$MAXDIS[OKq]^2)
q[OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[OKq], tau1))
ii <- which(p == 0)
p[ii] <- sra_fun(offdat$MAXDUR[ii], cf0[1])
offdat$p <- p
offdat$A <- A
offdat$q <- q
offdat$C <- (p*A*q)
#offdat$D <- offdat$count/offdat$C
offdat$offset <- log(p) + log(A) + log(q)
"Coefficients"
summary(offdat$C)


########################
## Density Estimation ##
########################

########STOP###################
offdat$count <- 0
offdat$count <- dat$"YEWA"
###############################

#offdat$count <- dat[, as.character(spp)]   
#### BHC
offdat$BHC <- offdat$BHC7

#names(offdat)
str(offdat$count)
summary(offdat$count)
fh <- (ggplot(offdat, aes(x=count))
       +geom_histogram(binwidth = NULL, bins = 50)
       +stat_bin(binwidth=0.5, bins=5)
       +labs(title=spp) 
       +xlab(label="Number of individuals detected in sample")
       +ylab(label="Number of samples")
)
fh
fh1 <- (ggplot(offdat, aes(count))
        +geom_histogram(binwidth = NULL, bins = 50)
        +stat_bin(binwidth=0.5, bins=5)
        +facet_wrap(BHC~., ncol=3)
        +labs(title=spp) 
        +xlab(label="Number of individuals detected in sample")
        +ylab(label="Number of samples")
)
fh1

## Model abundance
## Fit either poisson or neg binom distribution and test fit
fit <- goodfit(offdat$count)
fit
summary(fit)

m0 <- glmmTMB(count ~ (1|Station), data = offdat, offset = offset, family = (poisson("log")))
summary(m0)
m1 <- glmmTMB(count ~ BHC + (1|Station), data = offdat, offset = offset, family = (poisson("log")))
summary(m1)
m2 <- glmmTMB(count ~ BHC + year + (1|Station), data = offdat, offset = offset, family = (poisson("log")))
summary(m2)
mod.sel <- model.sel(m0, m1, m2)
mod.sel

# zip1 <- glmmTMB(count ~ BHC + (1|Station), data = offdat, offset = offset, family = poisson(link="log"), zi=~1)
# summary(zip1)
# zip2 <- glmmTMB(count ~ BHC + year + (1|Station), data = offdat, offset = offset, family = poisson(link="log"), zi=~1)
# summary(zip2)
# nb1 <- glmmTMB(count ~ BHC + (1|Station), data = offdat, offset = offset, family = nbinom2)
# summary(nb1)
# nb2 <- glmmTMB(count ~ BHC + year + (1|Station), data = offdat, offset = offset, family = nbinom2)
# summary(nb2)
# mod.sel <- model.sel(m0, m1, m2, zip1, zip2, nb1, nb2)
# mod.sel

#########
bestmod <- m1
summary(bestmod)

## EMMs for BHC
ref_grid(bestmod)
emmH <- as.data.frame(emmeans(bestmod, ~BHC))
emmH
emmH.out <- data.frame(emmH)
emmH.out$mean <- exp(emmH$emmean)
emmH.out$LCL <- exp(emmH$emmean-(1.96*emmH$SE))
emmH.out$UCL <- exp(emmH$emmean+(1.96*emmH$SE))
emmH.out$species <- paste(spp)
emmH.out
emmH.out$UCL[emmH.out$UCL == Inf] <- 0
g0 <- ggplot(emmH.out, aes(x=BHC, y=mean, fill=BHC))+
  geom_bar(stat = "identity", width = 0.5)+
  xlab(label = "Bird Habitat Class") +
  ylab(label = "Density (males/ha) with 95% CI")+
#  coord_cartesian(ylim=c(0,0.02))+
  labs(subtitle=spp)+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.05)+
  theme_bw()+
  theme(axis.title=element_text(size=10))+
  theme(legend.position = "none")
g0

## EMMs for Year
ref_grid(bestmod)
emmY <- as.data.frame(emmeans(bestmod, ~year))
emmY
emmY.out <- data.frame(emmY)
emmY.out$mean <- exp(emmY$emmean)
emmY.out$LCL <- exp(emmY$emmean-(1.96*emmY$SE))
emmY.out$UCL <- exp(emmY$emmean+(1.96*emmY$SE))
emmY.out$species <- paste(spp)
emmY.out
emmY.out$UCL[emmY.out$UCL == Inf] <- 0
g1 <- ggplot(emmY.out, aes(x=year, y=mean, fill=year))+
  geom_bar(stat = "identity", width = 0.5)+
  xlab(label = "Year") +
  ylab(label = "Density (males/ha) with 95% CI")+
#  coord_cartesian(ylim=c(0,0.2))+
  labs(subtitle=spp)+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.05)+
  theme_bw()+
  theme(axis.title=element_text(size=10))+
  theme(legend.position = "none")
g1

#####################
## Output when finalized for each species
#####################
# Make model selection table
mod.selout <- mod.sel %>% as_tibble(rownames="Model") %>% select(Model, AICc)
mod.selout$Species <- spp
write_csv(mod.selout, "densmodelsel.csv", col_names = TRUE, append = TRUE)
# Make best model parameter table
mout <- tidy(m1)
mout$Species <- spp
mout
write_csv(mout, "densmodelparams.csv", col_names = TRUE, append = TRUE)
emmH.out <- select(emmH.out, species, BHC, mean, LCL, UCL)
emmH.out
# not needed
# emmH.longer <- emmH.out %>% pivot_longer(-c(species, BHC), names_to = "estimate", values_to = "value")
write_csv(emmH.out, "densestimates.csv", col_names = TRUE, append = TRUE)
write_csv(emmY.out, "densestimatesyear.csv", col_names = TRUE, append = TRUE)
# Write the figures
ggsave(paste(spp,"-BHC.png", sep=""), g0, width = 8, height = 4)
ggsave(paste(spp,"-Year.png", sep=""), g1, width = 4, height = 4)

#########
## End ##
#########

