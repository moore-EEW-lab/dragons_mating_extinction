
library(car)
library(MASS)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)
library(phytools)
library(MCMCglmm)
library(visreg)

all.dat <- read.csv('nalley_moore_revised_ext.csv')
head(all.dat)


############ analyses in main text


### test if the local persistence/extinction of perchers vs fliers depends on local changes in temperature and local wildfire activity. First, fit overall model
new.a.mod00_RS <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.all.burned + Repro_strategy:z.temp.change + Repro_strategy:z.all.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                        family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)
summary(new.a.mod00_RS)

# estimate slopes of relationship between local persistence/extinction and regional warming for perchers vs fliers
emtrends(new.a.mod00_RS, 'Repro_strategy', var = 'z.temp.change')
# Repro_strategy z.temp.change.trend   SE  df asymp.LCL asymp.UCL
# f                           -0.283 0.153 Inf    -0.583    0.0165
# p                           -0.323 0.133 Inf    -0.584   -0.0616

# fliers have higher extinction risk than perchers although perchers may be more negatively affected by regional warming, but no real difference in magnitude of slopes


# estimate slopes of relationship between local persistence/extinction and wildfire activity for perchers vs fliers
emtrends(new.a.mod00_RS, 'Repro_strategy', var = 'z.all.burned')
# Repro_strategy z.all.burned.trend    SE  df asymp.LCL asymp.UCL
# f                          -0.392 0.160 Inf    -0.706   -0.0791
# p                          -0.311 0.143 Inf    -0.592   -0.0306


# both perchers and fliers are negatively affected by wildfire activity


### test if the local persistence/extinction of ornamented vs non-ornamented species depends on local changes in temperature and local wildfire activity. First, fit overall model
new.a.mod00 <- glmer(Persistence ~ M_wing_color + z.temp.change + z.all.burned + M_wing_color:z.temp.change + M_wing_color:z.all.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                     family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)
summary(new.a.mod00)

# estimate slopes of relationship between local persistence/extinction and regional warming for ornamented vs non-ornamented species
emtrends(new.a.mod00, 'M_wing_color', var = 'z.temp.change')
# M_wing_color z.temp.change.trend   SE  df asymp.LCL asymp.UCL
# n                         -0.136 0.13 Inf    -0.392     0.119
# y                         -0.449 0.13 Inf    -0.704    -0.194

## Increases in extinction risk with increases in temperature for ornamented species but not non-ornamented species (based on 95% CIs)


# estimate slopes of relationship between local persistence/extinction and wildfire activity for ornamented vs non-ornamented species
emtrends(new.a.mod00, 'M_wing_color', var = 'z.all.burned')
# M_wing_color z.all.burned.trend    SE  df asymp.LCL asymp.UCL
# n                         -0.20 0.148 Inf    -0.490    0.0904
# y                         -0.45 0.132 Inf    -0.709   -0.1911

## Increases in extinction risk with increases in wildfire activity for ornamented species but not non-ornamented species (based on 95% CIs)



#### Hypothesis tests to compare models with and without the interaction between mating traits all-severity wildfires 


# global model with interaction between reproductive strategy & all burned 
new.a.mod00_RS <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.all.burned + 
                          Repro_strategy:z.temp.change + Repro_strategy:z.all.burned + 
                          (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# model w/o interaction between reproductive strategy & temperature change 
new.a.mod00_RSb <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.all.burned + 
                           Repro_strategy:z.all.burned + 
                           (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# model w/o interaction between reproductive strategy & all burned
new.a.mod00_RSc <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.all.burned + 
                           Repro_strategy:z.temp.change + 
                           (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)


# Anova likelihood ratio test for models with and without reproductive strategy and temperature change interaction
anova(new.a.mod00_RS, new.a.mod00_RSb) 
#                 npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# new.a.mod00_RSb    9 1299.9 1348.4 -640.97   1281.9                     
# new.a.mod00_RS    10 1301.9 1355.7 -640.95   1281.9 0.0489  1      0.825

# both models are about the same and the effect of the interaction is not significant (p > 0.05)


# Anova likelihood ratio test for models with and without reproductive strategy and all severity wildfire activity interaction
anova(new.a.mod00_RS, new.a.mod00_RSc) 
#                 npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# new.a.mod00_RSc    9 1300.1 1348.5 -641.05   1282.1                     
# new.a.mod00_RS    10 1301.9 1355.7 -640.95   1281.9 0.2121  1     0.6451

# both models are about the same and the effect of the interaction is not significant (p > 0.05)



# global model with interaction between wing color & all burned
new.a.mod00 <- glmer(Persistence ~ M_wing_color + z.temp.change + z.all.burned + 
                       M_wing_color:z.temp.change + M_wing_color:z.all.burned + 
                       (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# model w/o interaction between wing color & temperature change
new.a.mod00b <- glmer(Persistence ~ M_wing_color + z.temp.change + z.all.burned + 
                        M_wing_color:z.all.burned + 
                        (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# model w/o interaction between wing color & all burned 
new.a.mod00c <- glmer(Persistence ~ M_wing_color + z.temp.change + z.all.burned + M_wing_color:z.temp.change + 
                        (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)


# Anova likelihood ratio test for models with and without wing color and temperature change interaction
anova(new.a.mod00, new.a.mod00b) 
#              npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# new.a.mod00b    9 1301.3 1349.7 -641.64   1283.3                       
# new.a.mod00    10 1299.9 1353.7 -639.94   1279.9 3.3849  1    0.06579 .

# model with with the interaction is marginally better, and the effect of the interaction is marginally significant


# Anova likelihood ratio test for models with and without wing color and all severity wildfire activity interaction
anova(new.a.mod00, new.a.mod00c) 
#              npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# new.a.mod00c    9 1300.0 1348.5 -641.03   1282.0                     
# new.a.mod00    10 1299.9 1353.7 -639.94   1279.9 2.1644  1     0.1412



#### Variance inflation factor to account for multicollearnarity between temperature and fire activity

# Wing ornamentation + temperature change + all wildfire severity
vif(new.a.mod00)
# M_wing_color              z.temp.change               z.all.burned M_wing_color:z.temp.change  M_wing_color:z.all.burned 
# 1.003829                   1.651430                   1.717610                   1.662162                   1.730428 

# All VIF < 1.8


# Wing ornamentation + temperature change + high wildfire severity
vif(new.h.mod00)
# M_wing_color              z.temp.change              z.high.burned M_wing_color:z.temp.change M_wing_color:z.high.burned 
# 1.002764                   1.645929                   1.662879                   1.636682                   1.653744 

# All VIF < 1.8



################## analyses in supporting information

##### Table S1
#### Examining sampling intensity 

## we would expect it to most strongly affect the variance explained by the random intercept for FID, since sampling intensity is a function of the grid cell and, by including it at all, we're already accounting for the possibility that the probability of re-observation will be different from one grid cell to the next anyways

# model that includes sampling intensity
new.a.mod00b <- glmer(Persistence ~ M_wing_color + z.temp.change + z.all.burned + M_wing_color:z.temp.change + M_wing_color:z.all.burned + z.per.obs + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)


summary(new.a.mod00)
# Random effects:
# Groups        Name        Variance Std.Dev.
# FID           (Intercept) 1.7056   1.3060  
# Genus_species (Intercept) 0.6524   0.8077  
# Genus         (Intercept) 0.7383   0.8592  
# Family        (Intercept) 0.3913   0.6255  

summary(new.a.mod00b)
# Random effects:
# Groups        Name        Variance Std.Dev.
# FID           (Intercept) 1.4122   1.1884  
# Genus_species (Intercept) 0.6988   0.8360  
# Genus         (Intercept) 0.7577   0.8705  
# Family        (Intercept) 0.4241   0.6512  

confint(new.a.mod00, parm = "M_wing_colory:z.temp.change", method = 'Wald')
# M_wing_colory:z.temp.change -0.6309089 0.006079071

confint(new.a.mod00, parm = "M_wing_colory:z.all.burned", method = 'Wald')
# M_wing_colory:z.all.burned -0.5685327 0.06825472

confint(new.a.mod00b, parm = "M_wing_colory:z.temp.change", method = 'Wald')
# M_wing_colory:z.temp.change -0.6598017 0.01456937

confint(new.a.mod00b, parm = "M_wing_colory:z.all.burned", method = 'Wald')
# M_wing_colory:z.all.burned -0.5885825 0.105667


emtrends(new.a.mod00, 'M_wing_color', var = 'z.temp.change')
# M_wing_color z.temp.change.trend   SE  df asymp.LCL asymp.UCL
# n                         -0.136 0.13 Inf    -0.392     0.119
# y                         -0.449 0.13 Inf    -0.704    -0.194

emtrends(new.a.mod00, 'M_wing_color', var = 'z.all.burned')
# M_wing_color z.all.burned.trend    SE  df asymp.LCL asymp.UCL
# n                         -0.20 0.148 Inf    -0.490    0.0904
# y                         -0.45 0.132 Inf    -0.709   -0.1911

emtrends(new.a.mod00b, 'M_wing_color', var = 'z.temp.change')
# M_wing_color z.temp.change.trend    SE  df asymp.LCL asymp.UCL
# n                        -0.0984 0.133 Inf    -0.359     0.163
# y                        -0.4210 0.134 Inf    -0.684    -0.158

emtrends(new.a.mod00b, 'M_wing_color', var = 'z.all.burned')
# M_wing_color z.all.burned.trend    SE  df asymp.LCL asymp.UCL
# n                        -0.110 0.155 Inf    -0.413    0.1930
# y                        -0.351 0.137 Inf    -0.620   -0.0825



##### Table S2
#### Species’ range limits and local extinction

FID_geometry <- read.csv('FIDcoordinates_states.csv')
head(FID_geometry)


count_data <- merge(all.dat, FID_geometry, by = "FID")
head(count_data)
unique(count_data$Genus_species)

### Table with observed latitudinal range and northernmost and southermost extinction points

summary_table <- count_data %>%
  group_by(Genus_species) %>%
  summarise(
    Overall_Northern_Latitude = max(Latitude, na.rm = TRUE),  # Overall maximum latitude for each species
    Northern_Latitude_Extinct = if (any(Persistence == 0)) max(Latitude[Persistence == 0], na.rm = TRUE) else NA_real_,  # Northern latitude of extinction, NA if none extinct
    Overall_Southern_Latitude = min(Latitude, na.rm = TRUE),  # Overall minimum latitude for each species
    Southern_Latitude_Extinct = if (any(Persistence == 0)) min(Latitude[Persistence == 0], na.rm = TRUE) else NA_real_  # Southern latitude of extinction, NA if none extinct
  )



### Table S3
##### Assessing if full topology of phylogeny should be directly modeled

library(MCMCglmm)
drags <- read.tree('mating.extinction.phylo.tre')
summary(drags)


all.dat$animal <- all.dat$Genus_species
drags$node.label <- seq(1:59)

## fit models for ornamentation, compare DIC. Easier to do it in MCMCglmm, but it's a Bayesian framework, so fit flat non-informative priors. 
# w/ phylogeny
prior.1 <- list(R = list(V=1, fix=1), G = list(G1 = list(V = 1, nu=0.0001), G2 = list(V = 1, nu=0.0001), G3 = list(V = 1, nu=0.0001)))

p.mod01 <- MCMCglmm(as.factor(Persistence) ~ M_wing_color + z.temp.change + z.high.burned + M_wing_color:z.temp.change + M_wing_color:z.high.burned, family = 'categorical', random =~Genus_species + FID + animal, prior = prior.1, pedigree = drags, node = 'tips', burnin = 50000, nitt = 200000, thin = 200, data = all.dat, pr = TRUE, verbose = FALSE)
summary(p.mod01) # DIC = 1124.558. note that DIC values might jump around a bit because of the Bayesian MCMC of it all


# w/ nested random effects
prior.2 <- list(R = list(V=1, fix =1), G = list(G1 = list(V = 1, nu = 0.0001), G2 = list(V = 1, nu=0.0001), G3 = list(V = 1, nu=0.0001), G4 = list(V = 1, nu=0.0001)))

p.mod02 <- MCMCglmm(as.factor(Persistence) ~ M_wing_color + z.temp.change + z.high.burned + M_wing_color:z.temp.change + M_wing_color:z.high.burned, family = 'categorical', random =~Genus_species + FID + Genus + Family, prior = prior.2, burnin = 50000, nitt = 200000, thin = 200, data = all.dat, pr = TRUE, verbose = FALSE)
summary(p.mod02) # DIC 1124.275



##### Table S4
#### Distribution of dragonfly species reproductive traits

head(all.dat)

# Data frame with only unique values of Genus_species so each species is represented once
all.dat.unique <- all.dat[!duplicated(all.dat$Genus_species), ] 

# Contigency table 
contingency_traits <- table(all.dat.unique$Repro_strategy, all.dat.unique$M_wing_color)
head(contingency_traits)

colnames(contingency_traits) <- c("Non-Ornamented", "Ornamented")
rownames(contingency_traits) <- c("Flier", "Percher")

# Chi-Square test of independence
chisq.test(contingency_traits)
# data:  contingency_traits
# X-squared = 3.1715, df = 1, p-value = 0.07493

# p > 0.05 suggests each of the two trait groups independently influence extinction risk rather than their interaction



##### Table S5
#### Parameter estimates for the effect of high-severity wildfires on extinction risk

### test if the local persistence/extinction of perchers vs fliers depends on local high-severity wildfire activity
new.h.mod00_RS <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.high.burned + Repro_strategy:z.temp.change + Repro_strategy:z.high.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# estimate slopes of relationship between local persistence/extinction and wildfire activity for perchers vs fliers
emtrends(new.h.mod00_RS, 'Repro_strategy', var = 'z.high.burned')
# Repro_strategy z.high.burned.trend    SE  df asymp.LCL asymp.UCL
# f                          0.00867 0.169 Inf    -0.323    0.3403
# p                         -0.26178 0.161 Inf    -0.577    0.0532

# neither perchers nor fliers are negatively affected by regional warming


### test if the local persistence/extinction of non-ornamented vs ornamented species depends on local high-severity wildfire activity
new.h.mod00 <- glmer(Persistence ~ M_wing_color + z.temp.change + z.high.burned + M_wing_color:z.temp.change + M_wing_color:z.high.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# estimate slopes of relationship between local persistence/extinction and wildfire activity for non-ornamented vs ornamented species
emtrends(new.h.mod00, 'M_wing_color', var = 'z.high.burned')
# M_wing_color z.high.burned.trend    SE  df asymp.LCL asymp.UCL
# n                         0.0327 0.157 Inf    -0.276     0.341
# y                        -0.3006 0.147 Inf    -0.589    -0.012

# Increases in extinction risk with increases in high severity wildfire activity for ornamented species but not non-ornamented species (based on 95% CIs)


#### Hypothesis tests to compare models with and without the interaction between mating traits high-severity wildfires 

# global model with interaction between wing color & high burned 
new.h.mod00 <- glmer(Persistence ~ M_wing_color + z.temp.change + z.high.burned + 
                       M_wing_color:z.temp.change + M_wing_color:z.high.burned + 
                       (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)


# model w/o interaction between wing color & high burned 
new.h.mod00c <- glmer(Persistence ~ M_wing_color + z.temp.change + z.high.burned + 
                        M_wing_color:z.temp.change + 
                        (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)


# Anova likelihood ratio test for models with and without wing color and high severity wildfire activity interaction
anova(new.h.mod00, new.h.mod00c) 
#             npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# new.h.mod00c    9 1307.1 1355.5 -644.55   1289.1                       
# new.h.mod00    10 1306.0 1359.8 -643.02   1286.0 3.0661  1    0.07994 . 

# model with the interaction between wing color and high severity wildfire is marginally better, although not quite significant (p > 0.05)


# global model with interaction between reproductive strategy & high burned
new.h.mod00_RS <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.high.burned + 
                          Repro_strategy:z.temp.change + Repro_strategy:z.high.burned + 
                          (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# model w/o interaction between reproductive strategy  & high burned 
new.h.mod00_RSc <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.high.burned + 
                           Repro_strategy:z.temp.change + 
                           (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)


# Anova likelihood ratio test for models with and without reproductive strategy and high severity wildfire activity interaction
anova(new.h.mod00_RS, new.h.mod00_RSc) 
#                 npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# new.h.mod00_RSc    9 1306.8 1355.3 -644.43   1288.8                     
# new.h.mod00_RS    10 1306.9 1360.7 -643.47   1286.9 1.9156  1     0.1663

# both models are about the same and the effect of the interaction is not significant (p > 0.05)



##### Table S6
#### Comparison of persistence models with and without a grid cell's historical temperature included as a covariate


# global ornamentation model with all burn severity 
new.a.mod00 <- glmer(Persistence ~ M_wing_color + z.temp.change + z.all.burned + M_wing_color:z.temp.change + M_wing_color:z.all.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                     family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# ornamentation model with all burn severity + historical temp
new.a.mod01 <- glmer(Persistence ~ M_wing_color + z.temp.change + z.all.burned + z.hist.temp + M_wing_color:z.temp.change + M_wing_color:z.all.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                     family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)
AICc(new.a.mod00, new.a.mod01) 
# df     AICc
# new.a.mod00 10 1300.085
# new.a.mod01 11 1297.814

# marginally improves ornamentation + all burn severity model, but really just barely


# global ornamentation model with high burn severity
new.h.mod00 <- glmer(Persistence ~ M_wing_color + z.temp.change + z.high.burned + M_wing_color:z.temp.change + M_wing_color:z.high.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                     family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# ornamentation model with high burn severity + historical temp
new.h.mod01 <- glmer(Persistence ~ M_wing_color + z.temp.change + z.high.burned + z.hist.temp + M_wing_color:z.temp.change + M_wing_color:z.high.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                     family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)
AICc(new.h.mod00, new.h.mod01) 
# df     AICc
# new.h.mod00 10 1306.176
# new.h.mod01 11 1306.521

# no evidence that adding historical temp improves ornamentation + high burn severity model


# global reproductive strategy model with all burn severity 
new.a.mod00_RS <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.all.burned + Repro_strategy:z.temp.change + Repro_strategy:z.all.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                        family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# reproductive strategy model with all burn severity + historical temp
new.a.mod01_RS <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.all.burned + z.hist.temp + Repro_strategy:z.temp.change + Repro_strategy:z.all.burned +
                          (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

AICc(new.a.mod00_RS, new.a.mod01_RS) 
# df     AICc
# new.a.mod00_RS 10 1302.030
# new.a.mod01_RS 11 1301.294

# no evidence that adding historical temp improves the reproductive strategy + all burn severity model


# global reproductive strategy model with high burn severity
new.h.mod00_RS <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.high.burned + Repro_strategy:z.temp.change + Repro_strategy:z.high.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                        family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# reproductive strategy model with high burned severity + historical temp
new.h.mod01_RS <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.high.burned + z.hist.temp + Repro_strategy:z.temp.change + Repro_strategy:z.high.burned + 
                          (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

AICc(new.h.mod00_RS, new.h.mod01_RS) 
# df     AICc
# new.h.mod00_RS 10 1307.073
# new.h.mod01_RS 11 1307.706

# no evidence that adding historical temp improves reproductive strategy + high burn severity model



##### Table S7
#### Parameter estimates ± SE and 95% confidence intervals for all severity wildfire persistence models with a grid cell's historical temperature included as a covariate

# ornamentation model with all burn severity + historical temp
new.a.mod01 <- glmer(Persistence ~ M_wing_color + z.temp.change + z.all.burned + z.hist.temp + M_wing_color:z.temp.change + M_wing_color:z.all.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                     family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)
emtrends(new.a.mod01, 'M_wing_color', var = 'z.temp.change')
# M_wing_color z.temp.change.trend    SE  df asymp.LCL asymp.UCL
# n                         -0.180 0.143 Inf     -0.46     0.100
# y                         -0.508 0.144 Inf     -0.79    -0.226

emtrends(new.a.mod01, 'M_wing_color', var = 'z.all.burned')
# M_wing_color z.all.burned.trend    SE  df asymp.LCL asymp.UCL
# n                        -0.273 0.162 Inf    -0.590    0.0442
# y                        -0.539 0.147 Inf    -0.827   -0.2510


# reproductive strategy model with all burn severity + historical temp
new.a.mod01_RS <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.all.burned + z.hist.temp + Repro_strategy:z.temp.change + Repro_strategy:z.all.burned +
                          (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)
emtrends(new.a.mod01_RS, 'Repro_strategy', var = 'z.temp.change')
# Repro_strategy z.temp.change.trend        SE  df asymp.LCL asymp.UCL
# f                          -0.3162 0.0006107 Inf   -0.3174   -0.3150
# p                          -0.3733 0.0007483 Inf   -0.3748   -0.3719

emtrends(new.a.mod01_RS, 'Repro_strategy', var = 'z.all.burned')
# Repro_strategy z.all.burned.trend        SE  df asymp.LCL asymp.UCL
# f                         -0.4413 0.0006105 Inf   -0.4425   -0.4401
# p                         -0.3887 0.0007480 Inf   -0.3902   -0.3873



##### Table S8
#### Comparison of persistence models with ornamentation quantified as variation in darkness of color (2=black, 1=brown/yellow, 0=no ornamented) vs. binary (ornamented vs non-ornamented)


# Make dark score a 3-level factor instead of a pseudo-continuous variable
all.dat$m.dark.score <- as.factor(all.dat$m.dark.score) 


# Binary color + all severity burn 
new.a.mod00 <- glmer(Persistence ~ M_wing_color + z.temp.change + z.all.burned + M_wing_color:z.temp.change + M_wing_color:z.all.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                     family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# Color variation + all severity burn
new.a.mod02 <- glmer(Persistence ~ m.dark.score + z.temp.change + z.all.burned + m.dark.score:z.temp.change + m.dark.score:z.all.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                     family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)
AICc(new.a.mod00, new.a.mod02) 
# df     AICc
# new.a.mod00 10 1300.085
# new.a.mod02 13 1304.537

# Color variation model decidedly worse than binary ornamented vs non-ornamented for all wildfire severity models


# Binary color + high severity burn
new.h.mod00 <- glmer(Persistence ~ M_wing_color + z.temp.change + z.high.burned + M_wing_color:z.temp.change + M_wing_color:z.high.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                     family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)

# Color variation + high severity burn
new.h.mod02 <- glmer(Persistence ~ m.dark.score + z.temp.change + z.high.burned + m.dark.score:z.temp.change + m.dark.score:z.high.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)
AICc(new.h.mod00, new.h.mod02) # decidedly worse than ornamented vs non-ornamented
# df     AICc
# new.h.mod00 10 1306.176
# new.h.mod02 13 1310.263

# Color variation model also decidedly worse than binary ornamented vs non-ornamented for high wildfire severity models



##### Table S9
#### Evaluating persistence results when different classifications of dragonfly wing ornamentation was used


## Color variation + all severity wildfire burned area
new.a.mod02 <- glmer(Persistence ~ m.dark.score + z.temp.change + z.all.burned + m.dark.score:z.temp.change + m.dark.score:z.all.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), 
                     family = 'binomial', data = all.dat, control = glmerControl('bobyqa'), na.action = na.fail)
emtrends(new.a.mod02, 'm.dark.score', var = 'z.temp.change')
# m.dark.score z.temp.change.trend    SE  df asymp.LCL asymp.UCL
# 0                         -0.136 0.130 Inf    -0.391    0.1186
# 1                         -0.615 0.199 Inf    -1.005   -0.2246
# 2                         -0.356 0.155 Inf    -0.660   -0.0517

emtrends(new.a.mod02, 'm.dark.score', var = 'z.all.burned')
# m.dark.score z.temp.change.trend    SE  df asymp.LCL asymp.UCL
# 0                        -0.199 0.148 Inf    -0.488    0.0907
# 1                        -0.381 0.202 Inf    -0.777    0.0142
# 2                        -0.489 0.152 Inf    -0.787   -0.1919



############### figures


####### fig 2a

## make dummy dataframe to predict values onto
Repro_strategy <- c(rep('f', 1000), rep('p', 1000))
z.temp.change <- rep(seq(min(all.dat$z.temp.change), max(all.dat$z.temp.change), length.out = 1000), 2)
z.all.burned <- rep(mean(all.dat$z.all.burned), 2000)
col.frame <- data.frame(Repro_strategy, z.temp.change, z.all.burned)
col.frame$pred.persist <- predict(new.a.mod00_RS, newdata = col.frame, re.form = NA, type = 'response')

## function for bootstrapping the 95% confidence intervals
myFunc2 <- function(mm) {
  predict(mm, newdata = col.frame, re.form = NA, type = 'response')
}

col.boot <- bootMer(new.a.mod00_RS, FUN = myFunc2, use.u = TRUE, type = 'parametric', nsim = 500) # need to use use.u = TRUE, otherwise it spits out prediction intervals instead of confidence intervals (which we don't want)

pred.col.boot <- t(apply(col.boot$t, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)))
col.frame$lci <- pred.col.boot [,1] # lower bound of the confidence band (allegedly) - this looks sensible to me based on the model output
col.frame$uci <- pred.col.boot [,2] # upper bound of the confidence band (allegedly) - this looks sensible to me based on the model output
col.frame$se <- apply(col.boot$t, MARGIN = 2, FUN = sd) # allegedly this is the standard error


### make columns that allow us to plot the raw data onto the figures
all.dat$persist.temp.percher <- ifelse(all.dat$Repro_strategy == 'p' & all.dat$Persistence == 1, all.dat$z.temp.change, NA)
all.dat$not.persist.temp.percher <- ifelse(all.dat$Repro_strategy == 'p' & all.dat$Persistence == 0, all.dat$z.temp.change, NA)
all.dat$persist.temp.flier <- ifelse(all.dat$Repro_strategy == 'f' & all.dat$Persistence == 1, all.dat$z.temp.change, NA)
all.dat$not.persist.temp.flier <- ifelse(all.dat$Repro_strategy == 'f' & all.dat$Persistence == 0, all.dat$z.temp.change, NA)


### make figure that depicts local extinction/persistence for perchers vs fliers in reponse to regional warming
rs.a.persist.plot <- 
  ggplot(data = col.frame, aes(x = z.temp.change, group = Repro_strategy)) +
  geom_rug(data = all.dat, aes(x = persist.temp.flier), color = '#5D9D9F', sides = 't', outside = FALSE, length = unit(0.4, 'lines')) +
  geom_rug(data = all.dat, aes(x = not.persist.temp.flier), color = '#5D9D9F', sides = 'b', outside = FALSE, length = unit(0.4, 'lines')) +
  geom_rug(data = all.dat, aes(x = persist.temp.percher), color = '#881213', sides = 't', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
  geom_rug(data = all.dat, aes(x = not.persist.temp.percher), color = '#881213', sides = 'b', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
  geom_ribbon(data = col.frame, aes(x = z.temp.change, ymin = lci, ymax = uci, fill = Repro_strategy), alpha = 0.32, size = 1) +
  geom_line(data = col.frame, aes(x = z.temp.change, y = pred.persist, color = Repro_strategy, linetype = Repro_strategy), linewidth = 1) +
  scale_color_manual(values = c('#5D9D9F', '#881213'), guide = 'none') + 
  scale_fill_manual(values = c('#609EA0', '#DCBFC1'), guide = 'none') +
  scale_linetype_manual(values = c(2, 1), guide = 'none') +
  coord_cartesian(ylim = c(0.29, 1)) +
  scale_y_continuous(breaks = c(0.3, 0.5, 0.7, 0.9)) +
  scale_x_continuous(breaks = c(-2.5, -1.25, 0, 1.25, 2.5)) +
  labs(y = 'Probability of Persistence', x = 'Increase in Air Temperature\n(z-transformed)')

bgrd01 =
  theme(axis.text = element_text(color="Black"),
        axis.title.x = element_text(face = "plain", size = 12, margin = margin(t = 10)), 
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(face = "plain", size = 12, margin = margin(r = 10)), 
        axis.text.y = element_text(size = 10),
        panel.background = element_rect(fill = "White"),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        axis.line.x = element_line(linetype = 'solid', color = 'black', linewidth = 0.8),
        axis.line.y = element_line(linetype = 'solid', color = 'black', linewidth = 0.8),
        legend.key = element_rect(fill = 'white')) 

rs.a.temp.persist.fig <- rs.a.persist.plot + bgrd01
png('revised.repro_strat.temp.persist.png', width = 3.7, height = 5.28, units = 'in', res = 600)
print(rs.a.temp.persist.fig)
dev.off()


#### fig 2b

### dummy dataset to predict values onto
Repro_strategy <- c(rep('f', 1000), rep('p', 1000))
z.temp.change <- rep(mean(all.dat$z.temp.change), 2000)
z.all.burned <- rep(seq(min(all.dat$z.all.burned), max(all.dat$z.all.burned), length.out = 1000), 2)
burn.frame <- data.frame(Repro_strategy, z.temp.change, z.all.burned)
burn.frame$pred.persist <- predict(new.a.mod00_RS, newdata = burn.frame, re.form = NA, type = 'response')


## function for bootstrapping the 95% confidence intervals
myFunc3 <- function(mm) {
  predict(mm, newdata = burn.frame, re.form = NA, type = 'response')
}

wf.boot <- bootMer(new.a.mod00_RS, FUN = myFunc3, use.u = TRUE, type = 'parametric', nsim = 500) # need to use use.u = TRUE, otherwise it spits out prediction intervals instead of confidence intervals (which we don't want)

pred.burn.boot <- t(apply(wf.boot$t, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)))
burn.frame$lci <- pred.burn.boot [,1] # lower bound of the confidence band (allegedly) - this looks sensible to me based on the model output
burn.frame$uci <- pred.burn.boot [,2] # upper bound of the confidence band (allegedly) - this looks sensible to me based on the model output
burn.frame$se <- apply(wf.boot$t, MARGIN = 2, FUN = sd) # allegedly this is the standard error


### make columns that allow us to plot the raw data onto the figures
all.dat$persist.burn.percher <- ifelse(all.dat$Repro_strategy == 'p' & all.dat$Persistence == 1, all.dat$z.all.burned, NA)
all.dat$not.persist.burn.percher <- ifelse(all.dat$Repro_strategy == 'p' & all.dat$Persistence == 0, all.dat$z.all.burned, NA)
all.dat$persist.burn.flier <- ifelse(all.dat$Repro_strategy == 'f' & all.dat$Persistence == 1, all.dat$z.all.burned, NA)
all.dat$not.persist.burn.flier <- ifelse(all.dat$Repro_strategy == 'f' & all.dat$Persistence == 0, all.dat$z.all.burned, NA)


## make figure that depicts relationship between local extinction/persistence for percher vs flier species in response to wildfire activity
rs.b.persist.plot <- 
  ggplot(data = burn.frame, aes(x = z.all.burned, group = Repro_strategy)) +
  geom_rug(data = all.dat, aes(x = persist.burn.flier), color = '#5D9D9F', sides = 't', outside = FALSE, length = unit(0.4, 'lines')) +
  geom_rug(data = all.dat, aes(x = not.persist.burn.flier), color = '#5D9D9F', sides = 'b', outside = FALSE, length = unit(0.4, 'lines')) +
  geom_rug(data = all.dat, aes(x = persist.burn.percher), color = '#881213', sides = 't', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
  geom_rug(data = all.dat, aes(x = not.persist.burn.percher), color = '#881213', sides = 'b', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
  geom_ribbon(data = burn.frame, aes(x = z.all.burned, ymin = lci, ymax = uci, fill = Repro_strategy), alpha = 0.32, size = 1) +
  geom_line(data = burn.frame, aes(x = z.all.burned, y = pred.persist, color = Repro_strategy, linetype = Repro_strategy), size = 1) +
  scale_color_manual(values = c('#5D9D9F', '#881213'), guide = 'none') +
  scale_fill_manual(values = c('#609EA0', '#DCBFC1'), guide = 'none') +
  scale_linetype_manual(values = c(1, 1), guide = 'none') +
  coord_cartesian(ylim = c(0.29, 1)) +
  scale_y_continuous(breaks = c(0.3, 0.5, 0.7, 0.9)) +
  scale_x_continuous(breaks = c(-1.50, -0.75, 0, 0.75, 1.50)) +
  labs(y = 'Probability of Persistence', x = 'Total Hectares Burned\n(z-transformed)')

bgrd01 =
  theme(axis.text = element_text(color="Black"),
        axis.title.x = element_text(face = "plain", size = 12, margin = margin(t = 10)), 
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(face = "plain", size = 12, margin = margin(r = 10)), 
        axis.text.y = element_text(size = 10),
        panel.background = element_rect(fill = "White"),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        axis.line.x = element_line(linetype = 'solid', color = 'black', size = 0.8),
        axis.line.y = element_line(linetype = 'solid', color = 'black', size = 0.8),
        legend.key = element_rect(fill = 'white')) 

rs.b.burn.persist.fig <- rs.b.persist.plot + bgrd01
png('revised.repro_strat.burn.persist.png', width = 3.7, height = 5.28, units = 'in', res = 600)
print(rs.burn.persist.fig)
dev.off()


######### ornament figures

#### fig 3a

## make dummy dataframe to predict values onto
M_wing_color <- c(rep('n', 1000), rep('y', 1000))
z.temp.change <- rep(seq(min(all.dat$z.temp.change), max(all.dat$z.temp.change), length.out = 1000), 2)
z.all.burned <- rep(mean(all.dat$z.all.burned), 2000)
col.frame <- data.frame(M_wing_color, z.temp.change, z.all.burned)
col.frame$pred.persist <- predict(new.a.mod00, newdata = col.frame, re.form = NA, type = 'response')

## function for bootstrapping the 95% confidence intervals
myFunc2 <- function(mm) {
  predict(mm, newdata = col.frame, re.form = NA, type = 'response')
}

col.boot <- bootMer(new.a.mod00, FUN = myFunc2, use.u = TRUE, type = 'parametric', nsim = 500) # need to use use.u = TRUE, otherwise it spits out prediction intervals instead of confidence intervals (which we don't want)

pred.col.boot <- t(apply(col.boot$t, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)))
col.frame$lci <- pred.col.boot [,1] # lower bound of the confidence band (allegedly) - this looks sensible to me based on the model output
col.frame$uci <- pred.col.boot [,2] # upper bound of the confidence band (allegedly) - this looks sensible to me based on the model output
col.frame$se <- apply(col.boot$t, MARGIN = 2, FUN = sd) # allegedly this is the standard error


### make columns that allow us to plot the raw data onto the figures
all.dat$persist.temp.orn <- ifelse(all.dat$M_wing_color == 'y' & all.dat$Persistence == 1, all.dat$z.temp.change, NA)
all.dat$not.persist.temp.orn <- ifelse(all.dat$M_wing_color == 'y' & all.dat$Persistence == 0, all.dat$z.temp.change, NA)
all.dat$persist.temp.non.orn <- ifelse(all.dat$M_wing_color == 'n' & all.dat$Persistence == 1, all.dat$z.temp.change, NA)
all.dat$not.persist.temp.non.orn <- ifelse(all.dat$M_wing_color == 'n' & all.dat$Persistence == 0, all.dat$z.temp.change, NA)

### make figure for local persistence/extinction for ornamented vs non-ornamented species in response to regional warming 
o.persist.plot <- 
  ggplot(data = col.frame, aes(x = z.temp.change, group = M_wing_color)) +
  geom_rug(data = all.dat, aes(x = persist.temp.non.orn), color = 'slategray', sides = 't', outside = FALSE, length = unit(0.4, 'lines')) +
  geom_rug(data = all.dat, aes(x = not.persist.temp.non.orn), color = 'slategray', sides = 'b', outside = FALSE, length = unit(0.4, 'lines')) +
  geom_rug(data = all.dat, aes(x = persist.temp.orn), color = 'sienna', sides = 't', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
  geom_rug(data = all.dat, aes(x = not.persist.temp.orn), color = 'sienna', sides = 'b', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
  geom_ribbon(data = col.frame, aes(x = z.temp.change, ymin = lci, ymax = uci, fill = M_wing_color), alpha = 0.32, size = 1) +
  geom_line(data = col.frame, aes(x = z.temp.change, y = pred.persist, color = M_wing_color, linetype = M_wing_color), size = 1) +
  scale_color_manual(values = c('slategray4', 'sienna'), guide = 'none') +
  scale_fill_manual(values = c('slategray3', 'sienna'), guide = 'none') +
  scale_linetype_manual(values = c(2, 1), guide = 'none') +
  coord_cartesian(ylim = c(0.29, 1)) +
  scale_y_continuous(breaks = c(0.3, 0.5, 0.7, 0.9)) +
  scale_x_continuous(breaks = c(-2.5, -1.25, 0, 1.25, 2.5)) +
  labs(y = 'Probability of Persistence', x = 'Increase in Air Temperature\n(z-transformed)')

bgrd01 =
  theme(axis.text = element_text(color="Black"),
        axis.title.x = element_text(face = "plain", size = 12, margin = margin(t = 10)), 
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(face = "plain", size = 12, margin = margin(r = 10)), 
        axis.text.y = element_text(size = 10),
        panel.background = element_rect(fill = "White"),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        axis.line.x = element_line(linetype = 'solid', color = 'black', size = 0.8),
        axis.line.y = element_line(linetype = 'solid', color = 'black', size = 0.8),
        legend.key = element_rect(fill = 'white')) 

o.temp.persist.fig <- o.persist.plot + bgrd01
png('revised.orn.temp.persist.png', width = 3.7, height = 5.28, units = 'in', res = 600)
print(o.temp.persist.fig)
dev.off()


#### fig 3b

### dummy dataset to predict values onto
M_wing_color <- c(rep('n', 1000), rep('y', 1000))
z.temp.change <- rep(mean(all.dat$z.temp.change), 2000)
z.all.burned <- rep(seq(min(all.dat$z.all.burned), max(all.dat$z.all.burned), length.out = 1000), 2)
burn.frame <- data.frame(M_wing_color, z.temp.change, z.all.burned)
burn.frame$pred.persist <- predict(new.a.mod00, newdata = burn.frame, re.form = NA, type = 'response')


## function for bootstrapping the 95% confidence intervals
myFunc3 <- function(mm) {
  predict(mm, newdata = burn.frame, re.form = NA, type = 'response')
}

wf.boot <- bootMer(new.a.mod00, FUN = myFunc3, use.u = TRUE, type = 'parametric', nsim = 500) # need to use use.u = TRUE, otherwise it spits out prediction intervals instead of confidence intervals (which we don't want)

pred.burn.boot <- t(apply(wf.boot$t, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)))
burn.frame$lci <- pred.burn.boot [,1] # lower bound of the confidence band (allegedly) - this looks sensible to me based on the model output
burn.frame$uci <- pred.burn.boot [,2] # upper bound of the confidence band (allegedly) - this looks sensible to me based on the model output
burn.frame$se <- apply(wf.boot$t, MARGIN = 2, FUN = sd) # allegedly this is the standard error


### make columns that allow us to plot the raw data onto the figures
all.dat$persist.burn.orn <- ifelse(all.dat$M_wing_color == 'y' & all.dat$Persistence == 1, all.dat$z.all.burned, NA)
all.dat$not.persist.burn.orn <- ifelse(all.dat$M_wing_color == 'y' & all.dat$Persistence == 0, all.dat$z.all.burned, NA)
all.dat$persist.burn.non.orn <- ifelse(all.dat$M_wing_color == 'n' & all.dat$Persistence == 1, all.dat$z.all.burned, NA)
all.dat$not.persist.burn.non.orn <- ifelse(all.dat$M_wing_color == 'n' & all.dat$Persistence == 0, all.dat$z.all.burned, NA)


### make figure for local persistence/extinction for ornamented vs non-ornamented species in response to wildfire activity
b.persist.plot <- 
  ggplot(data = burn.frame, aes(x = z.all.burned, group = M_wing_color)) +
  geom_rug(data = all.dat, aes(x = persist.burn.non.orn), color = 'slategray', sides = 't', outside = FALSE, length = unit(0.4, 'lines')) +
  geom_rug(data = all.dat, aes(x = not.persist.burn.non.orn), color = 'slategray', sides = 'b', outside = FALSE, length = unit(0.4, 'lines')) +
  geom_rug(data = all.dat, aes(x = persist.burn.orn), color = 'sienna', sides = 't', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
  geom_rug(data = all.dat, aes(x = not.persist.burn.orn), color = 'sienna', sides = 'b', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
  geom_ribbon(data = burn.frame, aes(x = z.all.burned, ymin = lci, ymax = uci, fill = M_wing_color), alpha = 0.32, size = 1) +
  geom_line(data = burn.frame, aes(x = z.all.burned, y = pred.persist, color = M_wing_color, linetype = M_wing_color), size = 1) +
  scale_color_manual(values = c('slategray4', 'sienna'), guide = 'none') +
  scale_fill_manual(values = c('slategray3', 'sienna'), guide = 'none') +
  scale_linetype_manual(values = c(2, 1), guide = 'none') +
  coord_cartesian(ylim = c(0.29, 1)) +
  scale_y_continuous(breaks = c(0.3, 0.5, 0.7, 0.9)) +
  scale_x_continuous(breaks = c(-1.50, -0.75, 0, 0.75, 1.50)) +
  labs(y = 'Probability of Persistence', x = 'Total Hectares Burned\n(z-transformed)')

bgrd01 =
  theme(axis.text = element_text(color="Black"),
        axis.title.x = element_text(face = "plain", size = 12, margin = margin(t = 10)), 
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(face = "plain", size = 12, margin = margin(r = 10)), 
        axis.text.y = element_text(size = 10),
        panel.background = element_rect(fill = "White"),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        axis.line.x = element_line(linetype = 'solid', color = 'black', size = 0.8),
        axis.line.y = element_line(linetype = 'solid', color = 'black', size = 0.8),
        legend.key = element_rect(fill = 'white')) 

o.burn.persist.fig <- b.persist.plot + bgrd01
png('revised.o.burn.persist.png', width = 3.7, height = 5.28, units = 'in', res = 600)
print(o.burn.persist.fig)
dev.off()


### fig S2


library(lme4)
library(lmerTest)
library(ggplot2)
library(phytools)
library(plyr)
library(metafor)
library(lsmeans)
library(visreg)
library(car)
library(nlme)
library(MuMIn)
library(ggplot2)
library(phylosignal)
library(phylobase)
library(viridis)
library(raster)
library(sp)
library(spgwr)
library(gstat)

## load in data
f.coords <- read.csv('FIDcoordinates_states.csv')
all.dat <- read.csv('nalley_moore_upd_ext.csv')
coords <- merge(all.dat, f.coords, by = 'FID') # this merges the lat/longs onto the broader persistence data


### Make a dataframe that just has the environmental variables so we can look at the spatial autocorrelation of the environmental factors
env.vars <- coords[,c('FID', 'Longitude', 'Latitude', 'z.temp.change', 'z.high.burned', 'z.all.burned')] 
colnames(env.vars) <- c('FID', 'x', 'y', 'warm', 'hb', 'ab') # rename columns so they're easier to work with

# Make dataframe into spatial dataframe
coordinates(env.vars) <- ~x+y # need to do this to tell R that this is a spatial dataframe

## Semi-variogram for high burn severity
burn.spat <- variogram(hb~1, data = env.vars)
plot(burn.spat, main = "Semi-Variogram for\nHigh Burn Severity") 

# High spatial autocorrelation, which is particularly expected in the US


## Semi-variogram for warming
warm.spat <- variogram(warm~1, data = env.vars)
plot(warm.spat, main = "Semi-Variogram for\nClimatic Warming") 

# High spatial autocorrelation, which is expected



### Make data frame of extinction rates for ornamented species in each grid cell
orn.sub <- subset(coords, M_wing_color == 'y') 
o.tab <- table(orn.sub$FID, orn.sub$Persistence) # table of ornamented species that persisted (1)/went extinct (0) in each grid cell
orn.dat <- data.frame(row.names(o.tab), o.tab[,1], o.tab[,2]) 
colnames(orn.dat) <-c('FID', 'e', 'p')

# Create a proportion of ornamented species that persisted in each grid cell
orn.dat$prop <- orn.dat$p/(orn.dat$e + orn.dat$p) 
o.coords <- merge(orn.dat, f.coords, by = 'FID') # now add lat/longs for each FID to this extinction rate data frame
colnames(o.coords) <- c('FID', 'e', 'p', 'prop', 'y', 'x', 'abbr', 'state') 

# Make dataframe into spatial dataframe
coordinates(o.coords) <- ~x+y 

## Semi-variogram for extinction rates of ornamented species
o.var01 <- variogram(prop~1, data = o.coords)
plot(o.var01, main = 'Semi-Variogram for\nOrnamented Species Extinction Rates') 

# very minor spatial autocorrelation


### Repeat process of creating spatial dataframe of extinction rates for non-ornamented species
no.sub <- subset(coords, M_wing_color == 'n')
no.tab <- table(no.sub$FID, no.sub$Persistence)
no.dat <- data.frame(row.names(no.tab), no.tab[,1], no.tab[,2])
colnames(no.dat) <-c('FID', 'e', 'p')

# Create a proportion of ornamented species that persisted in each grid cell
no.dat$prop <- no.dat$p/(no.dat$e + no.dat$p)
no.coords <- merge(no.dat, f.coords, by = 'FID')
colnames(no.coords) <- c('FID', 'e', 'p', 'prop', 'y', 'x', 'abbr', 'state') 

# Make dataframe into spatial dataframe
coordinates(no.coords) <- ~x+y

## Semi-variogram for extinction rates of non-ornamented species
no.var01 <- variogram(prop~1, data = no.coords)
plot(no.var01, main = 'Semi-Variogram for\nNon-Ornamented Species Extinction Rates') 

# minor spatial autocorrelation



### Make dataframes in a way that allows us to plot all of these semivariogram outputs on the same graph


# Pull columns 2 and 3 from the output of the ornamented species semi-variogram
o.var.tab <- o.var01[,c(2,3)] 
o.var.dat <- data.frame(o.var.tab) 

# scale the semivariance by the max semivariance of the dataset (because different variables have different max semivariances, this will allow us to plot all of the semivariograms in the same panel)
o.var.dat$scaled.var <- o.var.dat$gamma/max(o.var.dat$gamma) 


# Pull columns 2 and 3 from the output of the non-ornamented species semivariogram
no.var.tab <- no.var01[,c(2,3)]
no.var.dat <- data.frame(no.var.tab)
no.var.dat$scaled.var <- no.var.dat$gamma/max(no.var.dat$gamma)


# Pull columns 2 and 3 from the output of theburned area semivariogram
hb.var.tab <- burn.spat[,c(2,3)]
hb.var.dat <- data.frame(hb.var.tab)
hb.var.dat$scaled.var <- hb.var.dat$gamma/max(hb.var.dat$gamma)

# Pull columns 2 and 3 from the output of the warming semivariogram
warm.var.tab <- warm.spat[,c(2,3)]
warm.var.dat <- data.frame(warm.var.tab)
warm.var.dat$scaled.var <- warm.var.dat$gamma/max(warm.var.dat$gamma)



### Make figure to depict semivariance of species extinction rates and environmental variables
orn.semi.plot <- 
  ggplot() +
  geom_point(data = o.var.dat, aes(x = dist, y = scaled.var, color = "Ornamented", shape = "Ornamented")) +
  geom_point(data = no.var.dat, aes(x = dist, y = scaled.var, color = "Non-Ornamented", shape = "Non-Ornamented")) +
  geom_point(data = hb.var.dat, aes(x = dist, y = scaled.var, color = "High Burn Severity", shape = "High Burn Severity")) +
  geom_point(data = warm.var.dat, aes(x = dist, y = scaled.var, color = "Temperature Change", shape = "Temperature Change")) +
  labs(#title = "Semivariogram of\nWing Ornamentation and Environmental Factors", 
    x = 'Distance', y = 'Scaled Semivariance', 
    color = NULL, shape = NULL) +
  scale_color_manual(values = c('Ornamented' = 'black', 
                                'Non-Ornamented' = 'black', 
                                'High Burn Severity' = '#e66101', 
                                'Temperature Change' = '#5e3c99'),
                     labels = c("Ornamented", "Non-Ornamented", "High Burn Severity", "Temperature Change"),
                     breaks = c("Ornamented", "Non-Ornamented", "High Burn Severity", "Temperature Change")) +
  scale_shape_manual(values = c('Ornamented' = 15, 
                                'Non-Ornamented' = 0, 
                                'High Burn Severity' = 16,  
                                'Temperature Change' = 16),
                     labels = c("Ornamented", "Non-Ornamented", "High Burn Severity", "Temperature Change"),
                     breaks = c("Ornamented", "Non-Ornamented", "High Burn Severity", "Temperature Change")) 
bgrd01 =
  theme(#plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(color="Black"),
    axis.title.x = element_text(face = "plain", size = 12, margin = margin(t = 10)), 
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(face = "plain", size = 12, margin = margin(r = 10)), 
    axis.text.y = element_text(size = 10),
    panel.background = element_rect(fill = "White"),
    panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
    axis.line.x = element_line(linetype = 'solid', color = 'black', linewidth = 0.8),
    axis.line.y = element_line(linetype = 'solid', color = 'black', linewidth = 0.8),
    legend.key = element_rect(fill = 'white'),
    legend.title.align = 0.5,
    legend.position = c(1, 0.41),  # Slightly to the left from the very right edge
    legend.justification = c(1, 0.9), # Center the legend vertically, justify to the right
    legend.box.just = "bottom", # Center the title of the legend
    legend.spacing = unit(5, "cm"),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = .25), # Adds a border around the legend
    legend.margin = margin(1, 1, 1, 1), 
    legend.text = element_text(size = 10, margin = margin(0, 0.5, 0, -1)) #(B, R, T, L)
  )

orn.semi.fig <- orn.semi.plot + bgrd01
png('orn.semi.png', width = 6, height = 3, units = 'in', res = 600)
print(orn.semi.fig)
dev.off()
