#########################################
## Site C Songbird Occupancy Analysis ##
#########################################


library(unmarked)
library(tidyverse)
library(lubridate)
library(MuMIn)                # model selection
library(emmeans)

rm(list = ls())     

setwd("C:/Users/jeff.matheson/Documents/R Working/SiteC/Song/out")

dat <- read_csv("CAWA.csv")
         spp <- "CAWA"
##############
names(dat)

dat <- dat %>% rename(BHC20 = BHC20_100)
dat <- dat %>% rename(Station = "Sample Station Label")
table(dat$BHC20)
summary(dat$DSLS)
## Habitats BHC7
dat$BHC7 <- dat$BHC20
dat$BHC7[dat$BHC7 %in% "CSH"]	<- "Conif"
dat$BHC7[dat$BHC7 %in% "CYF"]	<- "Conif"
dat$BHC7[dat$BHC7 %in% "CMF"]	<- "Conif"
dat$BHC7[dat$BHC7 %in% "DSH"]	<- "Decid"
dat$BHC7[dat$BHC7 %in% "DYF"]	<- "Decid"
dat$BHC7[dat$BHC7 %in% "DMF"]	<- "Decid"
dat$BHC7[dat$BHC7 %in% "RSH"]	<- "Ripforest"
dat$BHC7[dat$BHC7 %in% "RYF"]	<- "Ripforest"
dat$BHC7[dat$BHC7 %in% "RMF"]	<- "Ripforest"
dat$BHC7[dat$BHC7 %in% "BSF"]	<- "Wetland"
dat$BHC7[dat$BHC7 %in% "BTF"]	<- "Conif"
dat$BHC7[dat$BHC7 %in% "WGR"]	<- "Wetland"
dat$BHC7[dat$BHC7 %in% "WSH"]	<- "Wetland"
dat$BHC7[dat$BHC7 %in% "WRI"]	<- "Wetland"
dat$BHC7[dat$BHC7 %in% "DSG"]	<- "Dry slopes"
dat$BHC7[dat$BHC7 %in% "DSS"]	<- "Dry slopes"
dat$BHC7[dat$BHC7 %in% "CUL"]	<- "Cultiv"
dat$BHC7[dat$BHC7 %in% "NVE"]	<- "Non-veg"
          ## Remove ANT
dat <- filter(dat, BHC7 != "ANT")
dat$BHC <- dat$BHC7
dat$BHC <- factor(dat$BHC, levels = c("Conif", "Decid", 'Ripforest', "Wetland", "Dry slopes", "Cultiv", "Non-veg"))
table(dat$BHC)
levels(dat$BHC)
names(dat)
dat$xYear <- as.numeric(scale(dat$Year))
dat$Year <- factor(dat$Year, levels = c("2006", "2008", '2011', "2012", "2016", "2017", "2018", "2019"))

summary(dat$xYear)
table(dat$Year)

obs <- select(dat, visit1, visit2, visit3, visit4)
obsCovs <- list(day=dat[,c("day1", "day2", "day3", "day4")], tssr=dat[,c("tssr1", "tssr2", "tssr3", "tssr4")], 
                duration=dat[,c("dur1", "dur2", "dur3", "dur4")])
siteCovs <- dat[,c("Year","BHC", "xYear")]
unmarkeddat <- unmarkedFrameOccu(y = obs, siteCovs = siteCovs, obsCovs = obsCovs)

summary(unmarkeddat)

m0 <- occu(~1 ~1, unmarkeddat)
m0
m1 <- occu(~ BHC ~ 1, unmarkeddat)
m1
m2 <- occu(~ tssr ~ 1, unmarkeddat)
m2
m3 <- occu(~ BHC + day ~ 1, unmarkeddat)
m3
m4 <- occu(~ BHC + tssr ~ 1, unmarkeddat)
m4
m5 <- occu(~ BHC + day + tssr ~ 1, unmarkeddat)
m5
m7 <- occu(~ xYear + tssr ~ 1, unmarkeddat)
m7
m8 <- occu(~ xYear + day + tssr ~ 1, unmarkeddat)
m8
m9 <- occu(~ BHC + xYear + day + tssr ~ 1, unmarkeddat)
m9

mod.sel <- model.sel(m0, m1, m2, m3, m4, m5, m7, m8, m9)
mod.sel

bestmod <- m9

# m9
grid <- expand.grid(
  BHC = c("Conif", "Decid", 'Ripforest', "Wetland", "Dry slopes", "Cultiv", "Non-veg"), 
  tssr = 0.061,
  day = 0.4436,
  xYear = 0
)
grid

#
grid <- expand.grid(
  BHC = c("Conif", "Decid", 'Ripforest', "Wetland", "Dry slopes", "Cultiv", "Non-veg") 
,tssr = 0.061
#  ,day = 0.4436
)
grid

p <- predict(bestmod, type = 'det', newdata = grid, appendData=TRUE)
p

### Manual removal o CIs
p$upper[p$upper == 1.0000000] <- 0
p

#p$upper[1] <- 0
#p$upper[7] <- 0
p

g1 <- ggplot(p, aes(x=BHC, y=Predicted, fill = BHC))+
  geom_bar(stat = "identity", width = 0.5, position=position_dodge())+
  labs(subtitle=spp)+
  xlab(label = "Bird Habitat Class") +
  ylab(label = "Probability of Detection with 95% CI")+
  theme_bw()+
  theme(legend.position = "none")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, position=position_dodge(width = 0.5))+
  theme(axis.title=element_text(size=10))
#  ylim(0, 0.5)
g1

#####################
## Output when finalized
#####################
# Make model selection table
mod.selout <- mod.sel %>% as_tibble(rownames="Model") %>% select(Model, AICc)
mod.selout$Species <- spp
write_csv(mod.selout, "occumodelsel.csv", col_names = TRUE, append = TRUE)
# Make best model parameter table
mout <- as.data.frame(coef(bestmod))
mout$species <- spp
mout <- rownames_to_column(mout, var = "coef")
mout
write_csv(mout, "occumodelparams.csv", col_names = TRUE, append = TRUE)
#Export estimates
p$Species <- spp
p
write_csv(p, "occuestimates.csv", col_names = TRUE, append = TRUE)
# Write the figure
ggsave(paste(spp,"-occu.png", sep=""), g1, width = 8, height = 4)


citation("unmarked")


#########
## End ##
#########


