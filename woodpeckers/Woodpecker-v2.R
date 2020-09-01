#####################################
## Site C Woodpecker Data Analysis ##
#####################################

rm(list = ls())      

library(tidyverse)
library(ggplot2)
library(lubridate)
library(lme4)                 # GLMMs
library(emmeans)              # for estiamting marginal means
library(DHARMa)               # Analysis of residuals for mixed models. By simulation.
library(glmmTMB)              # for nb and zi
library(MuMIn)                # model selection
library(vcd)                  # Fits a discrete (count data) distribution for goodness-of-fit tests.
library(broom.mixed)

dat <- read_csv("woodat.csv")
dat <- dat %>% rename("BHC" = "BHC20_100")
table(dat$BHC)

dat$Year <- factor(dat$Year, levels = c("2018", "2019"))
table(dat$Year)
dat <- mutate(dat, day = yday(dat$Date)/365)
dat <- mutate(dat, hour = hour(dat$Time)/24)

#dat$BHC <- factor(dat$BHC, levels = c("ANT", "CMF","CUL", "CYF", "DMF", "DSH", "DYF", "FBS", 
#                                      "FBT", "RMF", "RSH", "RYF", "WRI", "WSH"))
table(dat$BHC)
levels(dat$BHC)

## Lump and remove ANT and CUL

dat <- filter(dat, BHC != "ANT")
dat <- filter(dat, BHC != "CUL")
dat$BHC[dat$BHC %in% "FBT"] <- "FBT/FBS"
dat$BHC[dat$BHC %in% "FBS"] <- "FBT/FBS"
dat$BHC[dat$BHC %in% "WSH"] <- "WRI/WSH"
dat$BHC[dat$BHC %in% "WRI"] <- "WRI/WSH"

dat$BHC <- factor(dat$BHC, levels = c("CYF", "CMF", "DSH", "DYF", "DMF",  
                                      "RSH", "RYF", "RMF", "FBT/FBS", "WRI/WSH"))


table(dat$BHC)
levels(dat$BHC)

spp <- "DOWO"
spp
dat$count <- 0
dat$count <- dat$YBSA

#dat$count <- dat[, as.character(spp)]     Not working
#colnames(dat)[17] <- "count"

names(dat)
str(dat)
summary(dat$count)

fh <- (ggplot(dat, aes(x=count))
       +geom_histogram(binwidth = NULL, bins = 50)
       +stat_bin(binwidth=0.5, bins=5)
       +labs(title=spp) 
       +xlab(label="Number of individuals detected in sample")
       +ylab(label="Number of samples")
)
fh
fh1 <- (ggplot(dat, aes(count))
        +geom_histogram(binwidth = NULL, bins = 50)
        +stat_bin(binwidth=0.5, bins=5)
        +facet_wrap(BHC~., ncol=2)
        +labs(title=spp) 
        +xlab(label="Number of individuals detected in sample")
        +ylab(label="Number of samples")
)
fh1

## Model abundance
## Fit either poisson or neg binom distribution and test fit
fit <- goodfit(dat$count)
fit
summary(fit)

m0 <- glmmTMB(count ~ (1|Station), data = dat, family = (poisson("log")))
summary(m0)
m1 <- glmmTMB(count ~ BHC + (1|Station), data = dat, family = (poisson("log")))
summary(m1)
m2 <- glmmTMB(count ~ BHC + Year + (1|Station), data = dat, family = (poisson("log")))
summary(m2)
m3 <- glmmTMB(count ~ BHC + day + (1|Station), data = dat, family = (poisson("log")))
summary(m3)
m4 <- glmmTMB(count ~ BHC + hour + (1|Station), data = dat, family = (poisson("log")))
summary(m4)
m5 <- glmmTMB(count ~ BHC + Year + day + (1|Station), data = dat, family = (poisson("log")))
summary(m5)
m6 <- glmmTMB(count ~ BHC + day + hour + (1|Station), data = dat, family = (poisson("log")))
summary(m6)
m7 <- glmmTMB(count ~ BHC + Year + hour + (1|Station), data = dat, family = (poisson("log")))
summary(m7)
m8 <- glmmTMB(count ~ BHC + hour + day + Year + (1|Station), data = dat, family = (poisson("log")))
summary(m8)
mod.sel <- model.sel(m0, m1, m2, m3, m4, m5, m6, m7, m8)
mod.sel

bestmod <- m4
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
  #  scale_fill_manual(values=c("#1D4D18", "#5BC355", "#ABEB22","#EBB57C"))+
  xlab(label = "Bird Habitat Class") +
  ylab(label = "Relative Abundance (#/station) with 95% CI")+
  coord_cartesian(ylim=c(0,0.8))+
  theme_bw()+
  theme(legend.position = "none")+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.05)+
  theme(axis.title=element_text(size=10))
g0

#####################
## Output when finalized for each species
#####################

# Make model selection table
mod.selout <- mod.sel %>% as_tibble(rownames="Model") %>% select(Model, AICc)
mod.selout$Species <- spp
write_csv(mod.selout, "modelsel.csv", col_names = TRUE, append = TRUE)
# Make best model parameter table
mout <- tidy(m1)
mout$Species <- spp
mout
write_csv(mout, "modelparams.csv", col_names = TRUE, append = TRUE)
emmH.out <- select(emmH.out, species, BHC, mean, LCL, UCL)
emmH.out
# not needed
# emmH.longer <- emmH.out %>% pivot_longer(-c(species, BHC), names_to = "estimate", values_to = "value")
write_csv(emmH.out, "estimates.csv", col_names = TRUE, append = TRUE)
# Write the figure
ggsave(paste(spp,"g0.png", sep=""), g0, width = 8, height = 4)


