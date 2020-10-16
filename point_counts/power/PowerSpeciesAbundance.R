# Power Analysis Using Existing Data

# This is a replacement for the Rmd version. 

# Use the data from the last 3 years as pilot data and run power simulation for extending along years.
# The base simulation would 122 stations per year (366 over 3 years). Year will be continuous to simulate trend.
# Would also like to do simulation with 75 and 100 statons per year. Need to figure that out.

## Generate Offset ----

library(tidyverse)
library(ggplot2)
library(lubridate)
library(lme4)                 # GLMMs
library(emmeans)              # for estiamting marginal means
setwd("C:/Users/jeff.matheson/Documents/R Working/SiteC/Song")
dat <- read_csv("songdat.csv")
dat <- filter(dat, ValleyUpstPine == "Y")
dat <- dat %>% rename(BHC20 = BHC20_100)
dat <- dat %>% rename(Station = "Sample Station Label")
table(dat$BHC20)
dat <- mutate(dat, year = year(dat$Date))
dat$year <- factor(dat$year, levels = c("2006", "2008", "2011", "2012", "2016", "2017", "2018", "2019"))
table(dat$year)
dat <- mutate(dat, MonYear = year(dat$Date)-2016)
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
spp <- "OSFL"
##########
##############
#offdat$count <- dat[, as.character(spp)]
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

## Model ----

offdat <- filter(offdat, MonYear > 0)

########STOP###################
offdat$count <- 0
offdat$count <- offdat$"OSFL"
###############################

#offdat$count <- dat[, as.character(spp)]   
#### BHC
offdat$BHC <- offdat$BHC7

# offdat$count[offdat$count>1] <- 1  # for logistic

m <- glmer(count ~ MonYear + (1|BHC) + (1|Station), data = offdat, offset = offset, family = (poisson("log")))
summary(m)

## Simulate ----

library(simr)

# Extend to 10 yearts of monitoring
# alpha set at p=0.05.  
# Effect set at 50% negative change per year which is equal to about -6.7% per year.

msim <- m
summary(msim)

# Set effect size at 50%
fixef(msim)["MonYear"]
fixef(msim)["MonYear"]<- -0.0693
summary(msim)

# Assuming 10 years of monitoring, what is min number of stations to detect 50% over 10 years?
msim10 <- extend(msim, along = "MonYear", n=10)
summary(msim10)

c2.pc <- powerCurve(msim10, along = "Station", nsim = 50, alpha = 0.05)
c2.out <- summary(c2.pc)
c2.out

write_csv(c2.out, path = paste("./PowerOut/",spp,"-power-c1out.csv", sep=""))

c2.plot <- ggplot(c2.out, aes(x=nlevels, y=mean))+
  geom_line()+ geom_point()+ xlab(label = "Number of Stations")+ ylab(label = "Power")+
  labs(subtitle=spp)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05)+
  scale_y_continuous(breaks = seq(0, 1.0, by=0.2))+
  scale_x_continuous(breaks = seq(0, 300, by=50))+
  theme_bw()+ 
  coord_cartesian(ylim=c(0,1))+
  #theme(axis.title=element_text(size=10))+
  theme(legend.position = "none")+
  annotate(geom="text", x=175, y=0.3, label="Power to detect 50% change over 10 years",
           color="black")+
  annotate(geom="text", x=175, y=0.25, label="at vaying number of survey stations",
           color="black")
c2.plot

ggsave(paste("./PowerOut/",spp,"-power.png", sep=""), c2.plot, width = 7, height = 5)

