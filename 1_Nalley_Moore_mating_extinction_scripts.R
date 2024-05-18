setwd('/Users/michaelmoore/Desktop/Working Directory')

library(car)
library(MASS)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)
library(phytools)
library(MCMCglmm)
library(visreg)


dat <- read.csv('mating.extinction.csv')
head(dat)


############ analyses in main text

### test if the local persistence/extinction of perchers vs fliers depends on local changes in temperature and local wildfire activity. First, fit overall model
r.mod01 <- glmer(Persistence ~ Repro_strategy + z.temp.change + z.burned + Repro_strategy:z.temp.change + Repro_strategy:z.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = dat, control = glmerControl('bobyqa'), na.action = na.fail)
summary(r.mod01)

# estimate slopes of relationship between local persistence/extinction and regional warming for perchers vs fliers
emtrends(r.mod01, specs=~'Repro_strategy', var = 'z.temp.change')
 # Repro_strategy z.temp.change.trend    SE  df asymp.LCL asymp.UCL
 # f                           -0.326 0.154 Inf    -0.628   -0.0239
 # p                           -0.352 0.135 Inf    -0.615   -0.0880

# both perchers and fliers are negatively affected by regional warming, but no real difference in magnitude of slopes


# estimate slopes of relationship between local persistence/extinction and wildfire activity for perchers vs fliers
emtrends(r.mod01, specs =~'Repro_strategy', var = 'z.burned')
 # Repro_strategy z.burned.trend    SE  df asymp.LCL asymp.UCL
 # f                      -0.255 0.176 Inf    -0.600    0.0906
 # p                      -0.282 0.152 Inf    -0.579    0.0155

# neither perchers nor fliers are negatively affected by regional warming.


### test if the local persistence/extinction of ornamented vs non-ornamented species depends on local changes in temperature and local wildfire activity. First, fit overall model
o.mod01 <- glmer(Persistence ~ M_wing_color + z.temp.change + z.burned + M_wing_color:z.temp.change + M_wing_color:z.burned + (1|FID) + (1|Genus_species) + (1|Genus) + (1|Family), family = 'binomial', data = dat, control = glmerControl('bobyqa'), na.action = na.fail)


# estimate slopes of relationship between local persistence/extinction and regional warming for ornamented vs non-ornamented species
emtrends(o.mod01, 'M_wing_color', var = 'z.temp.change')
 # M_wing_color z.temp.change.trend    SE  df asymp.LCL asymp.UCL
 # n                         -0.154 0.132 Inf    -0.412     0.105
 # y                         -0.499 0.131 Inf    -0.755    -0.242

## Increases in extinction risk with increases in temperature for ornamented species but not non-ornamented species (based on 95% CIs)


# estimate slopes of relationship between local persistence/extinction and wildfire activity for ornamented vs non-ornamented species
emtrends(o.mod01, specs =~'M_wing_color', var = 'z.burned')
 # M_wing_color z.burned.trend    SE  df asymp.LCL asymp.UCL
 # n                    -0.135 0.159 Inf    -0.445     0.176
 # y                    -0.382 0.141 Inf    -0.659    -0.105

## Increases in extinction risk with increases in wildfire activity for ornamented species but not non-ornamented species (based on 95% CIs)


################## analyses in supporting information


#### check if models with phylogenetic variance-covariance matrix are preferred over those where the shared evo history is modeled as nested random effects
drags <- read.tree('mating.extinction.phylo.tre')
summary(drags)


dat$animal <- dat$Genus_species
drags$node.label <- seq(1:59)



## fit models for ornamentation, compare DIC. Easier to do it in MCMCglmm, but it's a Bayesian framework, so fit flat non-informative priors. Reminder that DIC bounces around a bit because of the Bayesian element
# w/ phylogeny
prior.1 <- list(R = list(V=1, fix=1), G = list(G1 = list(V = 1, nu=0.0001), G2 = list(V = 1, nu=0.0001), G3 = list(V = 1, nu=0.0001)))

mod01 <- MCMCglmm(as.factor(Persistence) ~ M_wing_color + z.temp.change + z.burned + M_wing_color:z.temp.change + M_wing_color:z.burned, family = 'categorical', random =~Genus_species + FID + animal, prior = prior.1, pedigree = drags, node = 'tips', burnin = 50000, nitt = 200000, thin = 200, data = dat, pr = TRUE, verbose = FALSE)
plot(mod01)
summary(mod01) # DIC = 1126.615


# w/ nested random effects
prior.2 <- list(R = list(V=1, fix =1), G = list(G1 = list(V = 1, nu = 0.0001), G2 = list(V = 1, nu=0.0001), G3 = list(V = 1, nu=0.0001), G4 = list(V = 1, nu=0.0001)))

mod02 <- MCMCglmm(as.factor(Persistence) ~ M_wing_color + z.temp.change + z.burned + M_wing_color:z.temp.change + M_wing_color:z.burned, family = 'categorical', random =~Genus_species + FID + Genus + Family, prior = prior.2, burnin = 50000, nitt = 200000, thin = 200, data = dat, pr = TRUE, verbose = FALSE)
plot(mod02)
summary(mod02) # DIC = 1127.038. About the same as with phylogeny. Model with nested random effects has more parameters BUT makes more assumptions about the exact nature of the shared evo history between species. Given that the phylogenetic model isn't doing better, safer to just do nested random effects 


## fit phylo vs hierarchical models for territorial behavior (percher/flier). 
# w/ phylogeny

mod03 <- MCMCglmm(as.factor(Persistence) ~ Repro_strategy + z.temp.change + z.burned + Repro_strategy:z.temp.change + Repro_strategy:z.burned, family = 'categorical', random =~Genus_species + FID + animal, prior = prior.1, pedigree = drags, node = 'tips', burnin = 50000, nitt = 200000, thin = 200, data = dat, pr = TRUE, verbose = FALSE)
plot(mod03)
summary(mod03) # DIC = 1125.335


# w/ nested random effects

mod04 <- MCMCglmm(as.factor(Persistence) ~ Repro_strategy + z.temp.change + z.burned + Repro_strategy:z.temp.change + Repro_strategy:z.burned, family = 'categorical', random =~Genus_species + FID + Genus + Family, prior = prior.2, burnin = 50000, nitt = 200000, thin = 200, data = dat, pr = TRUE, verbose = FALSE)
plot(mod04)
summary(mod04) # DIC = 1124.652. About the same as with phylogeny. Model with nested random effects has more parameters BUT makes more assumptions about the exact nature of the shared evo history between species. Given that the phylogenetic model isn't doing better, safer to just do nested random effects 


############### figures

####### fig 2a

## make dummy dataset to predict values onto
Repro_strategy <- c(rep('p', 1000), rep('f', 1000))
z.temp.change <- rep(seq(min(dat$z.temp.change), max(dat$z.temp.change), length.out = 1000), 2)
z.burned <- rep(mean(dat$z.burned), 2000)
repro.temp.frame <- data.frame(Repro_strategy, z.temp.change, z.burned)
repro.temp.frame$pred.persist <- predict(r.mod01, newdata = repro.temp.frame, re.form = NA, type = 'response')

## function for bootstrapped 95% confidence intervals
myFunc4 <- function(mm) {
	 predict(mm, newdata = repro.temp.frame, re.form = NA, type = 'response')
}

repro.temp.boot <- bootMer(r.mod01, FUN = myFunc2, use.u = TRUE, type = 'parametric', nsim = 500) # need to use use.u = TRUE, otherwise it spits out prediction intervals instead of confidence intervals (which we don't want)

pred.repro.temp.boot <- t(apply(repro.temp.boot$t, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)))
repro.temp.frame$lci <- pred.repro.temp.boot [,1] # lower bound of the confidence band (allegedly) - this looks sensible to me based on the model output
repro.temp.frame$uci <- pred.repro.temp.boot [,2] # upper bound of the confidence band (allegedly) - this looks sensible to me based on the model output
repro.temp.frame$se <- apply(repro.temp.boot$t, MARGIN = 2, FUN = sd) # allegedly this is the standard error


## make columns that allow us to plot the raw data onto the figures
dat$persist.temp.p <- ifelse(dat$Repro_strategy == 'p' & dat$Persistence == 1, dat$z.temp.change, NA)
dat$not.persist.temp.p <- ifelse(dat$Repro_strategy == 'p' & dat$Persistence == 0, dat$z.temp.change, NA)
dat$persist.temp.f <- ifelse(dat$Repro_strategy == 'f' & dat$Persistence == 1, dat$z.temp.change, NA)
dat$not.persist.temp.f <- ifelse(dat$Repro_strategy == 'f' & dat$Persistence == 0, dat$z.temp.change, NA)

## make figure that depicts local extinction/persistence for perchers vs fliers in reponse to regional warming
repro.t.persist.plot <- 
ggplot(data = repro.temp.frame, aes(x = z.temp.change, group = Repro_strategy)) +
geom_rug(data = dat, aes(x = persist.temp.p), color = 'firebrick4', sides = 't', outside = FALSE, length = unit(0.4, 'lines')) +
geom_rug(data = dat, aes(x = not.persist.temp.p), color = 'firebrick4', sides = 'b', outside = FALSE, length = unit(0.4, 'lines')) +
geom_rug(data = dat, aes(x = persist.temp.f), color = 'cadetblue', sides = 't', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
geom_rug(data = dat, aes(x = not.persist.temp.f), color = 'cadetblue', sides = 'b', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
geom_ribbon(data = repro.temp.frame, aes(x = z.temp.change, ymin = lci, ymax = uci, fill = Repro_strategy), alpha = 0.32, size = 1) +
geom_line(data = repro.temp.frame, aes(x = z.temp.change, y = pred.persist, color = Repro_strategy), size = 1) +
scale_color_manual(values = c('cadetblue', 'firebrick4'), guide = 'none') +
scale_fill_manual(values = c('cadetblue', 'firebrick4'), guide = 'none') +
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

repro.t.persist.fig <- repro.t.persist.plot + bgrd01
png('repro.t.persist.png', width = 3.7, height = 5.28, units = 'in', res = 600)
print(repro.t.persist.fig)
dev.off()

####### fig 2b


## make dummy dataset to predict values onto
Repro_strategy <- c(rep('p', 1000), rep('f', 1000))
z.temp.change <- rep(mean(dat$z.temp.change), 2000)
z.burned <- rep(seq(min(dat$z.burned), max(dat$z.burned), length.out = 1000), 2)
repro.burn.frame <- data.frame(Repro_strategy, z.temp.change, z.burned)
repro.burn.frame$pred.persist <- predict(r.mod01, newdata = repro.burn.frame, re.form = NA, type = 'response')


## function for bootstrap 95% confidence intervals
myFunc5 <- function(mm) {
	 predict(mm, newdata = repro.burn.frame, re.form = NA, type = 'response')
}

wf.r.boot <- bootMer(r.mod01, FUN = myFunc5, use.u = TRUE, type = 'parametric', nsim = 500) # need to use use.u = TRUE, otherwise it spits out prediction intervals instead of confidence intervals (which we don't want)

pred.burn.r.boot <- t(apply(wf.r.boot$t, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)))
repro.burn.frame$lci <- pred.burn.r.boot [,1] # lower bound of the confidence band (allegedly) - this looks sensible to me based on the model output
repro.burn.frame$uci <- pred.burn.r.boot [,2] # upper bound of the confidence band (allegedly) - this looks sensible to me based on the model output
repro.burn.frame$se <- apply(wf.r.boot$t, MARGIN = 2, FUN = sd) # allegedly this is the standard error


## make column for plotting raw data onto figure
dat$persist.burn.p <- ifelse(dat$Repro_strategy == 'p' & dat$Persistence == 1, dat$z.burned, NA)
dat$not.persist.burn.p <- ifelse(dat$Repro_strategy == 'p' & dat$Persistence == 0, dat$z.burned, NA)
dat$persist.burn.f <- ifelse(dat$Repro_strategy == 'f' & dat$Persistence == 1, dat$z.burned, NA)
dat$not.persist.burn.f <- ifelse(dat$Repro_strategy == 'f' & dat$Persistence == 0, dat$z.burned, NA)

## make figure that depicts local extinction/persistence for perchers vs fliers in reponse to wildfire activity
repro.b.persist.plot <- 
ggplot(data = repro.burn.frame, aes(x = z.burned, group = Repro_strategy)) +
geom_rug(data = dat, aes(x = persist.burn.p), color = 'firebrick4', sides = 't', outside = FALSE, length = unit(0.4, 'lines')) +
geom_rug(data = dat, aes(x = not.persist.burn.p), color = 'firebrick4', sides = 'b', outside = FALSE, length = unit(0.4, 'lines')) +
geom_rug(data = dat, aes(x = persist.burn.f), color = 'cadetblue', sides = 't', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
geom_rug(data = dat, aes(x = not.persist.burn.f), color = 'cadetblue', sides = 'b', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
geom_ribbon(data = repro.burn.frame, aes(x = z.burned, ymin = lci, ymax = uci, fill = Repro_strategy), alpha = 0.32, size = 1) +
geom_line(data = repro.burn.frame, aes(x = z.burned, y = pred.persist, color = Repro_strategy), size = 1) +
scale_color_manual(values = c('cadetblue', 'firebrick4'), guide = 'none') +
scale_fill_manual(values = c('cadetblue', 'firebrick4'), guide = 'none') +
coord_cartesian(ylim = c(0.29, 1)) +
scale_y_continuous(breaks = c(0.3, 0.5, 0.7, 0.9)) +
scale_x_continuous(breaks = c(-0.75, 0, 0.75, 1.5, 2.25)) +
labs(y = 'Probability of Persistence', x = 'Hectares Burned\n(z-transformed)')

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

repro.b.persist.fig <- repro.b.persist.plot + bgrd01
png('repro.b.persist.png', width = 3.7, height = 5.28, units = 'in', res = 600)
print(repro.b.persist.fig)
dev.off()




####### ornament figures

#### fig 3a

## make dummy dataframe to predict values onto
M_wing_color <- c(rep('n', 1000), rep('y', 1000))
z.temp.change <- rep(seq(min(dat$z.temp.change), max(dat$z.temp.change), length.out = 1000), 2)
z.burned <- rep(mean(dat$z.burned), 2000)
col.frame <- data.frame(M_wing_color, z.temp.change, z.burned)
col.frame$pred.persist <- predict(o.mod01, newdata = col.frame, re.form = NA, type = 'response')

## function and code to estimate bootstrap confidence intervals
myFunc2 <- function(mm) {
	 predict(mm, newdata = col.frame, re.form = NA, type = 'response')
}

col.boot <- bootMer(f.mod01, FUN = myFunc2, use.u = TRUE, type = 'parametric', nsim = 500) # need to use use.u = TRUE, otherwise it spits out prediction intervals instead of confidence intervals (which we don't want)

pred.col.boot <- t(apply(col.boot$t, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)))
col.frame$lci <- pred.col.boot [,1] # lower bound of the confidence band (allegedly) - this looks sensible to me based on the model output
col.frame$uci <- pred.col.boot [,2] # upper bound of the confidence band (allegedly) - this looks sensible to me based on the model output
col.frame$se <- apply(col.boot$t, MARGIN = 2, FUN = sd) # allegedly this is the standard error


### make columns to put raw data on graph
dat$persist.temp.orn <- ifelse(dat$M_wing_color == 'y' & dat$Persistence == 1, dat$z.temp.change, NA)
dat$not.persist.temp.orn <- ifelse(dat$M_wing_color == 'y' & dat$Persistence == 0, dat$z.temp.change, NA)
dat$persist.temp.non.orn <- ifelse(dat$M_wing_color == 'n' & dat$Persistence == 1, dat$z.temp.change, NA)
dat$not.persist.temp.non.orn <- ifelse(dat$M_wing_color == 'n' & dat$Persistence == 0, dat$z.temp.change, NA)

### make figure for local persistence/extinction for ornamented vs non-ornamented species in response to regional warming (figure 3a)
o.persist.plot <- 
ggplot(data = col.frame, aes(x = z.temp.change, group = M_wing_color)) +
geom_rug(data = dat, aes(x = persist.temp.non.orn), color = 'slategray', sides = 't', outside = FALSE, length = unit(0.4, 'lines')) +
geom_rug(data = dat, aes(x = not.persist.temp.non.orn), color = 'slategray', sides = 'b', outside = FALSE, length = unit(0.4, 'lines')) +
geom_rug(data = dat, aes(x = persist.temp.orn), color = 'sienna', sides = 't', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
geom_rug(data = dat, aes(x = not.persist.temp.orn), color = 'sienna', sides = 'b', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
geom_ribbon(data = col.frame, aes(x = z.temp.change, ymin = lci, ymax = uci, fill = M_wing_color), alpha = 0.32, size = 1) +
geom_line(data = col.frame, aes(x = z.temp.change, y = pred.persist, color = M_wing_color), size = 1) +
scale_color_manual(values = c('slategray4', 'sienna'), guide = 'none') +
scale_fill_manual(values = c('slategray4', 'sienna'), guide = 'none') +
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
png('orn.temp.persist.png', width = 3.7, height = 5.28, units = 'in', res = 600)
print(o.temp.persist.fig)
dev.off()


#### fig 3b

### dummy dataset to predict values onto
M_wing_color <- c(rep('n', 1000), rep('y', 1000))
z.temp.change <- rep(mean(dat$z.temp.change), 2000)
z.burned <- rep(seq(min(dat$z.burned), max(dat$z.burned), length.out = 1000), 2)
burn.frame <- data.frame(M_wing_color, z.temp.change, z.burned)
burn.frame$pred.persist <- predict(o.mod01, newdata = burn.frame, re.form = NA, type = 'response')


## function for bootstrapping the confidence intervals
myFunc3 <- function(mm) {
	 predict(mm, newdata = burn.frame, re.form = NA, type = 'response')
}

wf.boot <- bootMer(f.mod01, FUN = myFunc3, use.u = TRUE, type = 'parametric', nsim = 500) # need to use use.u = TRUE, otherwise it spits out prediction intervals instead of confidence intervals (which we don't want)

pred.burn.boot <- t(apply(wf.boot$t, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)))
burn.frame$lci <- pred.burn.boot [,1] # lower bound of the confidence band (allegedly) - this looks sensible to me based on the model output
burn.frame$uci <- pred.burn.boot [,2] # upper bound of the confidence band (allegedly) - this looks sensible to me based on the model output
burn.frame$se <- apply(wf.boot$t, MARGIN = 2, FUN = sd) # allegedly this is the standard error


## make columns to plot raw data onto figure
dat$persist.burn.orn <- ifelse(dat$M_wing_color == 'y' & dat$Persistence == 1, dat$z.burned, NA)
dat$not.persist.burn.orn <- ifelse(dat$M_wing_color == 'y' & dat$Persistence == 0, dat$z.burned, NA)
dat$persist.burn.non.orn <- ifelse(dat$M_wing_color == 'n' & dat$Persistence == 1, dat$z.burned, NA)
dat$not.persist.burn.non.orn <- ifelse(dat$M_wing_color == 'n' & dat$Persistence == 0, dat$z.burned, NA)


## make figure that depicts relationship between local extinction/persistence for orn vs non-orn spp in response to wildfire activity
b.persist.plot <- 
ggplot(data = burn.frame, aes(x = z.burned, group = M_wing_color)) +
geom_rug(data = dat, aes(x = persist.burn.non.orn), color = 'slategray', sides = 't', outside = FALSE, length = unit(0.4, 'lines')) +
geom_rug(data = dat, aes(x = not.persist.burn.non.orn), color = 'slategray', sides = 'b', outside = FALSE, length = unit(0.4, 'lines')) +
geom_rug(data = dat, aes(x = persist.burn.orn), color = 'sienna', sides = 't', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
geom_rug(data = dat, aes(x = not.persist.burn.orn), color = 'sienna', sides = 'b', outside = FALSE, length = unit(0.4, 'lines'), alpha = 0.4) +
geom_ribbon(data = burn.frame, aes(x = z.burned, ymin = lci, ymax = uci, fill = M_wing_color), alpha = 0.32, size = 1) +
geom_line(data = burn.frame, aes(x = z.burned, y = pred.persist, color = M_wing_color), size = 1) +
scale_color_manual(values = c('slategray4', 'sienna'), guide = 'none') +
scale_fill_manual(values = c('slategray4', 'sienna'), guide = 'none') +
coord_cartesian(ylim = c(0.29, 1)) +
scale_y_continuous(breaks = c(0.3, 0.5, 0.7, 0.9)) +
scale_x_continuous(breaks = c(-0.75, 0, 0.75, 1.50, 2.25)) +
labs(y = 'Probability of Persistence', x = 'Hectares Burned\n(z-transformed)')

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
png('o.burn.persist.png', width = 3.7, height = 5.28, units = 'in', res = 600)
print(o.burn.persist.fig)
dev.off()
