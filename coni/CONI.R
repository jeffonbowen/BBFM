###################
## CONI Analysis ##
###################

rm(list = ls())      
      
library(tidyverse)

CONIobs <- read_csv("CONIobs.csv", col_types = cols(
                    DistClass = col_character()))
table(CONIobs$DistClass)

CONIobs$DistClass <- factor(CONIobs$DistClass, labels = c("0-50m", "51-100m", "100+m"))

fct_count(CONIobs$DistClass)

CONIheard <- filter(CONIobs, Detect_Type2 == "H")

g1 <- ggplot(data=CONIheard, aes(x=DistClass))+
          geom_bar(stat="count", width = 0.3)+
          theme_bw()+
          labs(x="Distance From Observer", y="Count of CONI Heard")
g1

ggsave("CONIdistance.png", g1, width = 5, height = 3.5)



###################
## CONI Hab Dat  ##
###################

CONIhab <- read_csv("CONI-BHC.csv")
CONIhab$Rowid_ <- NULL
CONIhabl <- CONIhab %>% pivot_longer(-STATIONID, names_to = "BirdHabCat", values_to = "Area")

result <- CONIhabl %>% 
  group_by(STATIONID) %>%
  filter(Area == max(Area))

table(result$BirdHabCat)

# alternate
result <- CONIhabl %>% 
  group_by(STATIONID) %>%
  top_n(n=1)


############
## Paired ##
############

rm(list = ls())      

library(tidyverse)
library(lme4)                 # GLMMs
library(emmeans)              # for estiamting marginal means
library(DHARMa)               # Analysis of residuals for mixed models. By simulation.
library(glmmTMB)              # for nb and zi

paired <- read_csv("PairedData.csv")
paired$Year <- factor(paired$Year)

paired

# lubridate::yday(x)   Not needed right now

g2 <- ggplot(paired, aes(Count,colour=Method))+
       geom_histogram(stat="count", position="dodge", fill="white")

g2

paired %>% 
  group_by(Method) %>%
  tally()

SimpleMean <- paired %>% 
  group_by(Method) %>%
  summarize(mean=mean(Count))
SimpleMean

m4 <- glmer(Count ~ Method + (1|StationID), data = paired, family = (poisson("log")))
summary(m4)

qqnorm(residuals(m4))
qqline(residuals(m4))
plot(m4,Habitat~resid(.))

#For other model distributions, use DHARMAa (simulation)
res <- simulateResiduals(fittedModel = m4)
plot(res)

# However, for fixed effects, need to recalculate predictions by hand - see help ?predict.glmmTMB.
# To compute population-level predictions for a given grouping variable (i.e., setting all random effects for that grouping variable to zero), set the group value to NA.
newdata=paired
newdata$group = NA
pred = predict(m4, newdata = newdata)
# Now perform the plot
plotResiduals(pred, res$scaledResiduals)

ref_grid(m4)

emmeans(m4, ~Method)
pairs(emmeans(m4, ~Method))

emm <- as.data.frame(emmeans(m4, ~Method))
emm


emm$mean <- exp(emm$emmean)
emm$LCL <- exp(emm$emmean-(1.96*emm$SE))
emm$UCL <- exp(emm$emmean+(1.96*emm$SE))
emm

write_csv(as.data.frame(emm), "pairedout.csv")


"Correction Factor ARU:PC"
emm[1,7]/emm[2,7]

g3<- ggplot(emm, aes(x=Method, y=mean, fill=Method))+
  geom_bar(stat = "identity", width = 0.5)+
  xlab(label = "Survey Method") +
  ylab(label = "Mean Relative Abundance with 95% CI")+
  theme_bw()+
  theme(legend.position = "none")+
    geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.05)
g3

ggsave("CONIpaired.png", g3, width = 5, height = 3.5)


########################
## Density by Habitat ##
########################

rm(list = ls())      

library(tidyverse)
library(lubridate)
library(lme4)                 # GLMMs
library(emmeans)              # for estiamting marginal means
library(DHARMa)               # Analysis of residuals for mixed models. By simulation.
library(glmmTMB)              # for nb and zi
library(MuMIn)                # model selection
library(vcd)                 # Fits a discrete (count data) distribution for goodness-of-fit tests.

## Read in in the ARU data

arudat <- read_csv("ARUHL.csv")

table(arudat$BHC)
arudat$BHC[arudat$BHC %in% "WSH"] <- "WGR"
arudat$BHC[arudat$BHC %in% "DSG"] <- "DSS"
arudat$BHC <- factor(arudat$BHC, levels = c("DSH", "RSH", "WGR", "WRI", "DSS", "CUL", "NVE", "ANT"))
table(arudat$BHC)
levels(arudat$BHC)

arudat <- arudat %>% mutate(year = as.character(year(arudat$RecDate))) 
arudat$year <- factor(arudat$year, levels = c("2018", "2019"))
table(arudat$year)

arudat <- mutate(arudat, day = yday(arudat$RecDate)/365)
arudat <- mutate(arudat, hour = hour(arudat$RecTime)/24)

fh <- ggplot(arudat, aes(Count))+
       geom_histogram(binwidth = NULL, bins = 50)+
       stat_bin(binwidth=0.5, bins=5)+
       xlab(label="Number of individuals detected in sample")+
       ylab(label="Number of samples")

fh

## Model abundance

## Fit either poisson or neg binom distribution and test fit
fit <- goodfit(arudat$Count)
fit
summary(fit)

m0 <- glmmTMB(Count ~ (1|Station), data = arudat, family = (poisson("log")))
summary(m0)
m1 <- glmmTMB(Count ~ BHC + (1|Station), data = arudat, family = (poisson("log")))
summary(m1)
m2 <- glmmTMB(Count ~ year + (1|Station), data = arudat, family = (poisson("log")))
summary(m2)
m3 <- glmmTMB(Count ~ day + (1|Station), data = arudat, family = (poisson("log")))
summary(m3)
m4 <- glmmTMB(Count ~ hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m4)
m5 <- glmmTMB(Count ~ BHC + year + (1|Station), data = arudat, family = (poisson("log")))
summary(m5)
m6 <- glmmTMB(Count ~ BHC + day + (1|Station), data = arudat, family = (poisson("log")))
summary(m6)
m7 <- glmmTMB(Count ~ BHC + hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m7)
m8 <- glmmTMB(Count ~ BHC + year + day + (1|Station), data = arudat, family = (poisson("log")))
summary(m8)
m9 <- glmmTMB(Count ~ BHC + day + hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m9)
m10 <- glmmTMB(Count ~ BHC + year + hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m10)
m11 <- glmmTMB(Count ~ BHC + hour + day + year + (1|Station), data = arudat, family = (poisson("log")))
summary(m11)
m12 <- glmmTMB(Count ~ year + day + (1|Station), data = arudat, family = (poisson("log")))
summary(m12)
m13 <- glmmTMB(Count ~ year + hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m13)
m14 <- glmmTMB(Count ~ year + day + hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m14)

mod.sel <- model.sel(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14)
mod.sel

write_csv(mod.sel, "modesel1.csv")

bestmod <- m10

## EMMs for BHC
ref_grid(bestmod)
emmH <- as.data.frame(emmeans(bestmod, ~BHC))
emmH
emmH.out <- data.frame(emmH)

emmH.out$mean <- exp(emmH$emmean)
emmH.out$LCL <- exp(emmH$emmean-(1.96*emmH$SE))
emmH.out$UCL <- exp(emmH$emmean+(1.96*emmH$SE))
emmH.out

write_csv(as.data.frame(emmH.out), "aruout.csv")

g0 <- ggplot(emmH.out, aes(x=BHC, y=mean, fill=BHC))+
  geom_bar(stat = "identity", width = 0.5)+
#  scale_fill_manual(values=c("#1D4D18", "#5BC355", "#ABEB22","#EBB57C"))+
  xlab(label = "Bird Habitat Class") +
  ylab(label = "Relative Abundance (CONI/station) with 95% CI")+
  theme_bw()+
  theme(legend.position = "none")+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.05)+
  theme(axis.title=element_text(size=10))
g0

ggsave("aruabundance.png", g0, width = 8, height = 5)

########################
## CONI Booming ##
########################

rm(list = ls())      

library(tidyverse)
library(lubridate)
library(lme4)                 # GLMMs
library(emmeans)              # for estiamting marginal means
library(DHARMa)               # Analysis of residuals for mixed models. By simulation.
library(glmmTMB)              # for nb and zi
library(MuMIn)                # model selection
library(vcd)                 # Fits a discrete (count data) distribution for goodness-of-fit tests.

## Read in in the ARU data

arudat <- read_csv("ARUHL.csv")

table(arudat$BHC)
arudat$BHC[arudat$BHC %in% "WSH"] <- "WGR"
arudat$BHC[arudat$BHC %in% "DSG"] <- "DSS"
arudat$BHC <- factor(arudat$BHC, levels = c("DSH", "RSH", "WGR", "WRI", "DSS", "CUL", "NVE", "ANT"))
table(arudat$BHC)
levels(arudat$BHC)

arudat <- arudat %>% mutate(year = as.character(year(arudat$RecDate))) 
arudat$year <- factor(arudat$year, levels = c("2018", "2019"))
table(arudat$year)

arudat <- mutate(arudat, day = yday(arudat$RecDate)/365)
arudat <- mutate(arudat, hour = hour(arudat$RecTime)/24)

arudat$BoomCount <- ifelse(arudat$BoomCount > 0, 1, 0)

fh <- ggplot(arudat, aes(Count))+
  geom_histogram(binwidth = NULL, bins = 50)+
  stat_bin(binwidth=0.5, bins=5)+
  xlab(label="Number of individuals detected in sample")+
  ylab(label="Number of samples")

fh

## Model presence

## Fit either poisson or neg binom distribution and test fit
fit <- goodfit(arudat$BoomCount)
fit
summary(fit)

m0 <- glmmTMB(BoomCount ~ (1|Station), data = arudat, family = binomial)
summary(m0)
m1 <- glmmTMB(BoomCount ~ BHC + (1|Station), data = arudat, family = (poisson("log")))
summary(m1)
m2 <- glmmTMB(BoomCount ~ year + (1|Station), data = arudat, family = (poisson("log")))
summary(m2)
m3 <- glmmTMB(BoomCount ~ day + (1|Station), data = arudat, family = (poisson("log")))
summary(m3)
m4 <- glmmTMB(BoomCount ~ hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m4)
m5 <- glmmTMB(BoomCount ~ BHC + year + (1|Station), data = arudat, family = (poisson("log")))
summary(m5)
m6 <- glmmTMB(BoomCount ~ BHC + day + (1|Station), data = arudat, family = (poisson("log")))
summary(m6)
m7 <- glmmTMB(BoomCount ~ BHC + hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m7)
m8 <- glmmTMB(BoomCount ~ BHC + year + day + (1|Station), data = arudat, family = (poisson("log")))
summary(m8)
m9 <- glmmTMB(BoomCount ~ BHC + day + hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m9)
m10 <- glmmTMB(BoomCount ~ BHC + year + hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m10)
m11 <- glmmTMB(BoomCount ~ BHC + hour + day + year + (1|Station), data = arudat, family = (poisson("log")))
summary(m11)
m12 <- glmmTMB(BoomCount ~ year + day + (1|Station), data = arudat, family = (poisson("log")))
summary(m12)
m13 <- glmmTMB(BoomCount ~ year + hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m13)
m14 <- glmmTMB(BoomCount ~ year + day + hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m14)

mod.sel <- model.sel(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14)
mod.sel

write_csv(mod.sel, "modesel1.csv")

bestmod <- m10

## EMMs for BHC
ref_grid(bestmod)
emmH <- as.data.frame(emmeans(bestmod, ~BHC))
emmH
emmH.out <- data.frame(emmH)

emmH.out$mean <- exp(emmH$emmean)
emmH.out$LCL <- exp(emmH$emmean-(1.96*emmH$SE))
emmH.out$UCL <- exp(emmH$emmean+(1.96*emmH$SE))
emmH.out

write_csv(as.data.frame(emmH.out), "aruout.csv")

g0 <- ggplot(emmH.out, aes(x=BHC, y=mean, fill=BHC))+
  geom_bar(stat = "identity", width = 0.5)+
  #  scale_fill_manual(values=c("#1D4D18", "#5BC355", "#ABEB22","#EBB57C"))+
  xlab(label = "Bird Habitat Class") +
  ylab(label = "Relative Abundance (CONI/station) with 95% CI")+
  theme_bw()+
  theme(legend.position = "none")+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.05)+
  theme(axis.title=element_text(size=10))
g0

ggsave("aruabundance.png", g0, width = 8, height = 5)

###############
## occupancy ##
###############

rm(list = ls())      

library(unmarked)
library(tidyverse)
library(lubridate)
library(glmmTMB)              # for nb and zi
library(MuMIn)                # model selection
library(emmeans)

arudat <- read_csv("ARUHLOCCU.csv")

table(arudat$BHC)

arudat$BHC[arudat$BHC %in% "WSH"] <- "WGR"
arudat$BHC[arudat$BHC %in% "DSG"] <- "DSS"
arudat$BHC <- factor(arudat$BHC, levels = c("DSH", "RSH", "WGR", "WRI", "DSS", "CUL", "NVE", "ANT"))
table(arudat$BHC)
levels(arudat$BHC)

arudat$Year <- factor(arudat$Year, levels = c("2018", "2019"))
table(arudat$Year)

arudat <- mutate(arudat, day1 = yday(arudat$Date1)/365)
arudat <- mutate(arudat, day2 = yday(arudat$Date2)/365)

arudat <- mutate(arudat, hour1 = hour(arudat$Time1)/24)
arudat <- mutate(arudat, hour2 = hour(arudat$Time2)/24)


obs <- select(arudat, Boom1, Boom2)

obsCovs <- list(day=arudat[,c("day1", "day2")], hour=arudat[,c("hour1", "hour2")])

siteCovs <- arudat[,c("Year", "BHC")]

unmarkeddat <- unmarkedFrameOccu(y = obs, siteCovs = siteCovs, obsCovs = obsCovs)
summary(unmarkeddat)

fm0 <- occu(~1 ~1, unmarkeddat)
fm0
fm1 <- occu(~ BHC ~ 1, unmarkeddat)
fm1
fm2 <- occu(~ Year ~ 1, unmarkeddat)
fm2
fm3 <- occu(~ day ~ 1, unmarkeddat)
fm3
fm4 <- occu(~ hour ~ 1, unmarkeddat)
fm4
fm5 <- occu(~ BHC + Year ~ 1, unmarkeddat)
fm5
fm6 <- occu(~ BHC + day ~ 1, unmarkeddat)
fm6
fm7 <- occu(~ BHC + hour ~ 1, unmarkeddat)
fm7
fm8 <- occu(~ BHC + Year + day ~ 1, unmarkeddat)
fm8
fm9 <- occu(~ BHC + day + hour ~ 1, unmarkeddat)
fm9
fm10 <- occu(~ BHC + Year + hour ~ 1, unmarkeddat)
fm10
fm11 <- occu(~ BHC + Year + day + hour ~ 1, unmarkeddat)
fm11
fm12 <- occu(~ Year + day ~ 1, unmarkeddat)
fm12
fm13 <- occu(~ Year + hour ~ 1, unmarkeddat)
fm13
fm14 <- occu(~ Year + day + hour ~ 1, unmarkeddat)
fm14

mod.sel <- model.sel(fm0, fm1, fm2, fm3, fm4, fm5, fm6, fm7, fm8, fm9, fm10, fm11, fm12, fm13, fm14)
mod.sel

bestmod <- fm10

write_csv(mod.sel, "modeseloccu.csv")

backTransform(fm10, 'state')


backTransform(linearComb(fm10, coefficients = c(1,0,0,0,0), type = 'det'))



grid <- expand.grid(
  BHC = c("DSH", "RSH", "WGR", "WRI", "DSS", "CUL", "NVE", "ANT"), 
  Year = c("2018","2019"),
  hour = 0.934
  )
grid

p <- predict(fm10, type = 'det', newdata = grid, appendData=TRUE)
p

write_csv(p, "p.csv")

g1 <- ggplot(p, aes(x=BHC, y=Predicted, fill=Year))+
  geom_bar(stat = "identity", width = 0.5, position=position_dodge())+
  scale_fill_manual(values=c("dodgerblue3", "orangered3"))+
  xlab(label = "Bird Habitat Class") +
  ylab(label = "Probability of Detection with 95% CI")+
  theme_bw()+
  theme(legend.position = c(0.7,0.9))+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, position=position_dodge(width = 0.5))+
  theme(axis.title=element_text(size=10))
g1

ggsave("aruoccupancy.png", g1, width = 8, height = 5)

citation("unmarked")

###########################
## Mitigation Properties ##
###########################

rm(list = ls())      

library(tidyverse)
library(lubridate)
library(lme4)                 # GLMMs
library(emmeans)              # for estiamting marginal means
library(DHARMa)               # Analysis of residuals for mixed models. By simulation.
library(glmmTMB)              # for nb and zi
library(MuMIn)                # model selection
library(vcd)                 # Fits a discrete (count data) distribution for goodness-of-fit tests.

## Read in in the ARU data

arudat <- read_csv("ARUHL.csv")

arudat <- filter(arudat, MitProp %in% c("Marl", "Rutledge", "Wilder"))

table(arudat$MitProp)

table(arudat$BHC)
arudat$BHC[arudat$BHC %in% "WSH"] <- "WGR"
arudat$BHC[arudat$BHC %in% "DSG"] <- "DSS"
arudat$BHC <- factor(arudat$BHC, levels = c("DSH", "RSH", "WGR", "WRI", "DSS", "CUL", "NVE", "ANT"))
table(arudat$BHC)
levels(arudat$BHC)

arudat <- arudat %>% mutate(year = as.character(year(arudat$RecDate))) 
arudat$year <- factor(arudat$year, levels = c("2018", "2019"))
table(arudat$year)

arudat <- mutate(arudat, day = yday(arudat$RecDate)/365)
arudat <- mutate(arudat, hour = hour(arudat$RecTime)/24)

arudat <- filter(arudat, MitProp == c("Marl", "Rutledge", "Wilder"))

fh <- ggplot(arudat, aes(Count))+
  geom_histogram(binwidth = NULL, bins = 50)+
  stat_bin(binwidth=0.5, bins=5)+
  xlab(label="Number of individuals detected in sample")+
  ylab(label="Number of samples")

fh

## Model abundance

## Fit either poisson or neg binom distribution and test fit
fit <- goodfit(arudat$Count)
fit
summary(fit)

m0 <- glmmTMB(Count ~ (1|Station), data = arudat, family = (poisson("log")))
summary(m0)
m1 <- glmmTMB(Count ~ MitProp + (1|Station), data = arudat, family = (poisson("log")))
summary(m1)
m3 <- glmmTMB(Count ~ day + (1|Station), data = arudat, family = (poisson("log")))
summary(m3)
m4 <- glmmTMB(Count ~ hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m4)
m6 <- glmmTMB(Count ~ MitProp + day + (1|Station), data = arudat, family = (poisson("log")))
summary(m6)
m7 <- glmmTMB(Count ~ MitProp + hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m7)
m9 <- glmmTMB(Count ~ MitProp + day + hour + (1|Station), data = arudat, family = (poisson("log")))
summary(m9)

mod.sel <- model.sel(m0, m1, m3, m4, m6, m7, m9)
mod.sel

write_csv(mod.sel, "modesel1.csv")

bestmod <- m1

## EMMs for BHC
ref_grid(bestmod)
emmH <- as.data.frame(emmeans(bestmod, ~MitProp))
emmH
emmH.out <- data.frame(emmH)

emmH.out$mean <- exp(emmH$emmean)
emmH.out$LCL <- exp(emmH$emmean-(1.96*emmH$SE))
emmH.out$UCL <- exp(emmH$emmean+(1.96*emmH$SE))
emmH.out

write_csv(as.data.frame(emmH.out), "aruout.csv")

g0 <- ggplot(emmH.out, aes(x=MitProp, y=mean, fill=MitProp))+
  geom_bar(stat = "identity", width = 0.5)+
  #  scale_fill_manual(values=c("#1D4D18", "#5BC355", "#ABEB22","#EBB57C"))+
  xlab(label = "Bird Habitat Class") +
  ylab(label = "Relative Abundance (CONI/station) with 95% CI")+
  theme_bw()+
  theme(legend.position = "none")+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.05)+
  theme(axis.title=element_text(size=10))
g0

ggsave("aruabundance.png", g0, width = 8, height = 5)
