#### Snag Data ----

library(readxl)
library(tidyverse)
library(vcd)
library(MASS)
library(pscl)
library(visreg)               #regression visualization
library(MuMIn)                #model selection
library(effects)              #for extracting effects terms. Use effects or Alleffects
library(emmeans)              #Estimate marginal means
library(gridExtra)
library(here)


### Some Analysis ---
snag_dat <- read_xlsx(here("cavity_trees", "SnagData_2020-08-03.xlsx"), 
                      sheet = "Snag_Dat_Unified")

snag_dat$Structural_Stage <- as_factor(snag_dat$Structural_Stage)

table(snag_dat$Bird_Habitat_Class)
table(snag_dat$Structural_Stage)

hist(snag_dat$SN50_Total)

## SN11_Total
g_box <- ggplot(snag_dat) + 
  geom_boxplot(aes(Structural_Stage, SN11_Total))
g_box

snag_dat <- snag_dat %>% 
  filter(SN11_Total < 20)

fit.nb <- goodfit(snag_dat$SN11_Total, type = "nbinomial")
summary(fit.nb)
rootogram(fit.nb)

# Negative binomial is best fit.
m.nb <- glm.nb(SN11_Total ~ Structural_Stage, data = snag_dat)
summary(m.nb)

# best to simulate residuals to check model fit. 
simulationOutput <- simulateResiduals(fittedModel = m.nb, plot = T)

# Quick check
ae <- allEffects(m.nb)
plot(ae)

# Predict means
# using emmeans
means_tot <- emmeans(m.nb, specs = "Structural_Stage", type = "response")
summary(means_tot)
means_tot <- tidy(summary(means_tot))
means_tot

g_means_tot <- means_tot %>% 
  ggplot(aes(x=Structural_Stage, y=response))+
  geom_bar(stat = "identity", width = 0.5)+
  xlab(label = "Structural Stage") +
  ylab(label = "Snag Density (# / 400m2) with 95% CI")+
  labs(subtitle = "Total Snags")+
  coord_cartesian(ylim=c(0,5))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05)+
  theme_bw()+
  theme(axis.title=element_text(size=10))+
  theme(legend.position = "none")
g_means_tot
  
## SN11_Conif
snag_dat <- read_xlsx("SnagData_2020-08-03.xlsx", sheet = "Snag_Dat_Unified")
snag_dat$Structural_Stage <- as_factor(snag_dat$Structural_Stage)

g_box <- ggplot(snag_dat) + 
  geom_boxplot(aes(Structural_Stage, SN11_Conif_Tot))
g_box
snag_dat <- snag_dat %>% 
  filter(SN11_Conif_Tot < 15)

fit.nb <- goodfit(snag_dat$SN11_Conif_Tot, type = "nbinomial")
summary(fit.nb)
rootogram(fit.nb)

# Negative binomial is best fit.
m.nb <- glm.nb(SN11_Conif_Tot ~ Structural_Stage, data = snag_dat)
summary(m.nb)

# best to simulate residuals to check model fit. 
simulationOutput <- simulateResiduals(fittedModel = m.nb, plot = T)

# Quick check
ae <- allEffects(m.nb)
plot(ae)

# Predict means
# using emmeans
means_conif <- emmeans(m.nb, specs = "Structural_Stage", type = "response")
summary(means_conif)
means_conif <- tidy(summary(means_conif))
means_conif

g_means_conif <- means_conif %>% 
  ggplot(aes(x=Structural_Stage, y=response))+
  geom_bar(stat = "identity", width = 0.5, fill = "blue")+
  xlab(label = "Structural Stage") +
  ylab(label = "Snag Density (# / 400m2) with 95% CI")+
  labs(subtitle = "Coniferous Snags")+
  coord_cartesian(ylim=c(0,5))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05)+
  theme_bw()+
  theme(axis.title=element_text(size=10))+
  theme(legend.position = "none")
g_means_conif

## SN11_decid
snag_dat <- read_xlsx("SnagData_2020-08-03.xlsx", sheet = "Snag_Dat_Unified")
snag_dat$Structural_Stage <- as_factor(snag_dat$Structural_Stage)

g_box <- ggplot(snag_dat) + 
  geom_boxplot(aes(Structural_Stage, SN11_Decid_Tot))
g_box

fit.nb <- goodfit(snag_dat$SN_11_Decid_Tot, type = "nbinomial")
summary(fit.nb)
rootogram(fit.nb)

# Negative binomial is best fit.
m.nb <- glm.nb(SN11_Decid_Tot ~ Structural_Stage, data = snag_dat)
summary(m.nb)

# best to simulate residuals to check model fit. 
simulationOutput <- simulateResiduals(fittedModel = m.nb, plot = T)

# Quick check
ae <- allEffects(m.nb)
plot(ae)

# Predict means
# using emmeans
means_decid <- emmeans(m.nb, specs = "Structural_Stage", type = "response")
summary(means_decid)
means_decid <- tidy(summary(means_decid))
means_decid

g_means_decid <- means_decid %>% 
  ggplot(aes(x=Structural_Stage, y=response))+
  geom_bar(stat = "identity", width = 0.5, fill = "green")+
  labs(subtitle = "Deciduous Snags")+
  xlab(label = "Structural Stage") +
  ylab(label = "Snag Density (# / 400m2) with 95% CI")+
  coord_cartesian(ylim=c(0,5))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05)+
  theme_bw()+
  theme(axis.title=element_text(size=10))+
  theme(legend.position = "none")
g_means_decid

g_all <- grid.arrange(g_means_conif, g_means_decid, g_means_tot, nrow = 1)

write_csv(means_conif, "snag_estiamtes.csv", append = TRUE)
write_csv(means_decid, "snag_estiamtes.csv", append = TRUE)
write_csv(means_tot, "snag_estiamtes.csv", append = TRUE)

###---###

## Generic Function

# We have already explored the data so know that negbinom is fit. No
# need included any EDA. 

snag_dat_raw <- read_xlsx("SnagData_2020-08-03.xlsx", sheet = "Snag_Dat_Unified")
snag_dat <- snag_dat_raw
snag_dat$Structural_Stage <- as_factor(snag_dat$Structural_Stage)

# Fit negbinom

model <- SN11_Decid_Tot ~ Structural_Stage

m.nb <- glm.nb(model, data = snag_dat)
summary(m.nb)
# best to simulate residuals to check model fit. 
simulationOutput <- simulateResiduals(fittedModel = m.nb, plot = T)

# Quick check
ae <- allEffects(m.nb)
plot(ae)

# Predict means
# using emmeans
means_decid <- emmeans(m.nb, specs = "Structural_Stage", type = "response")
summary(means_decid)
means_decid <- tidy(summary(means_decid))
means_decid

g_means_decid <- means_decid %>% 
  ggplot(aes(x=Structural_Stage, y=response))+
  geom_bar(stat = "identity", width = 0.5, fill = "green")+
  labs(subtitle = "Deciduous Snags")+
  xlab(label = "Structural Stage") +
  ylab(label = "Snag Density (# / 400m2) with 95% CI")+
  coord_cartesian(ylim=c(0,5))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05)+
  theme_bw()+
  theme(axis.title=element_text(size=10))+
  theme(legend.position = "none")
g_means_decid
