library(tidyverse)
library(mgcv)
library(gamair)


#te() instead pf #s() effects + interactions


# Read in temperature data ------------------------------------------------

gtemp <- read_csv("hadcrutv4.csv")

plot(gtemp$year, gtemp$temperature)



# Fit GAM -----------------------------------------------------------------

m <- gam(temperature ~ s(year), data = gtemp, method = "REML")

summary(m)
plot(m, shade = T)


## set basis complexity (maximum wigglyness) k
# want to make this large enough to approximate the true function
# don't worry much about setting it too high -- penalty term (lambda) takes care of it
# REML chooses lambda for you
# bigger k does increase computational cost
# for regression predictors, usually doesn't need to be very high (10-15 max)
# for time series or spatial data, might need to be quite high
m <- gam(temperature ~ s(year, k = 10), data = gtemp, method = "REML")

# check whether k is high enough
gam.check(m)

# low p-value and effective degrees of freedom (EDF) close to k indicate k is too low
# increase k
m2 <- gam(temperature ~ s(year, k = 20), data = gtemp, method = "REML")

gam.check(m2)


# Bird example ------------------------------------------------------------

data(bird)

bird <- bird |> 
  mutate(crestlark = factor(crestlark), 
         linnet = factor(linnet), 
         e = x/1000,
         n = y/1000)

ggplot(bird, aes(x = e, y = n, color = crestlark)) + 
  geom_point() +
  coord_fixed()


## Fit model ----
crest <- gam(crestlark ~ s(e, n, k = 100), 
             data = bird, 
             family = binomial, 
             method = 'REML')

# s(e,n) fits 2D spatial smooth (isotropic thin plate spline--smoothness is the same in both dimensions) 
# k sets upper limit on model complexity (EDF)
# smoothness parameters estimated via REML

gam.check(crest)   # EDF is well below k, so k is high enough
summary(crest)

# plot predicted
crest$model |> 
  mutate(pred = predict(crest)) |> 
  ggplot(aes(x = e, y = n, color = pred)) + 
  geom_point() +
  coord_fixed()

# plot residuals
crest$model |> 
  mutate(resid = resid(crest)) |> 
  ggplot(aes(x = e, y = n, color = resid)) + 
  geom_point() +
  coord_fixed()

# plot binned residuals
bird |> 
  filter(!is.na(crestlark)) |> 
  mutate(pred = fitted(crest), 
         resid = resid(crest, type = "response"),
         bin = ntile(pred, 20)) |> 
  ggplot(aes(x = bin, y = resid)) + stat_summary()



# Cornucopia of smooths ---------------------------------------------------

# type of smooth controlled by bs argument (think basis)

# Low-rank thin plate spline (default): bs = 'tp'
# Cubic splines (faster for large datasets): bs = 'cr'
# Cyclic splines (e.g. daily or annual): bs = 'cc' or bs = 'cp'
# Random effects: bs = 're'
# Duchon splines (good for spatial): bs = 'ds'
# Spline on the sphere (e.g., global coordinates): bs = 'sos'
# Gaussian process: bs = 'gp'
# Adaptive splines (variable wigglyness): bs = 'ad'

## Fit bird model with Duchon spline----
crest2 <- gam(crestlark ~ s(e, n, bs = 'ds', k = 100, m = c(1, 0.5)),
             data = bird, 
             family = binomial, 
             method = 'REML')

gam.check(crest2)   # EDF is well below k, so k is high enough
summary(crest2)
plot(crest2)

# Bestiary of conditional distributions -----------------------------------

# binomal()           
# poisson() 
# Gamma()
# nb()        negative binomial (counts w/ overdispersion)
# tw()        Tweedie (continuous w/ zeroes)
# mvn()       multivariate normal
# multinom()  multinomial
# betar()     beta (continuous proportions)
# scat()      scale t (Gaussian w/ thick tails)
# gaulss()    Gaussian location-scale (fit mean and variance)
# ziplss()    zero-inflated (hurdle) Poisson (separate models for probability of zero and mean given non-zero)
# ziP()       zero-inflated (hurdle) Poisson (zero probability is function of Poisson mean parameter)
# twlss()     Tweedie location-scale
# cox.ph()    gamma location-scale
# ocat()      ordered categorical
# cox.ph()    Cox model (survival)
# inverse.gaussian



# Smooth interactions -----------------------------------------------------

# Two ways to fit:

#For two continuous variables

# 1. Bivariate (or higher order) thin plate splines
      # s(x, z, bs = 'tp')
      # isotropic: single smoothness parameter
      # use when x and z measured on same scale
# 2. Tensor product smooths
      # separate smoothness parameter for each variable
      # use for interactions when variables in different units
      # te(x, z)

# Three ways to build tensor products in mgcv
      # te(x, z) is most general form but not usable in gamm4::gamm4() or brms
      # t2(x, z) is alternative parameterization that does work in gamm4::gamm4() and brms
      # s(x) + s(z) + ti(x, z) decomposes interaction into main effects and interaction, ti(x, z)

# Two ways to fit interaction between smooth and factor variable

# 1. y ~ f + s(x, by = f)
      # entirely separate smooth for each level of the factor
      # each has its own smoothness parameter
      # smooths are centered, so need to include factor as a fixed effect
# 2. y ~ s(x, f, bs = 'fs')
      # smooth function for each level of the factor
      # share common smoothness parameter
      # includes group means, fully penalized
      # closer to random effects


# Random effects in mgcv

# s(f, bs = 're') will fit random intercept for f
# s(f, bs = 're') + s(f, x, bs = 're') fits (uncorrelated) random intercept and slope 

# for simple random effects without a ton of factor levels, use gam() or bam() function
  # can use all distribution families listed above, including those not in nlme or lme4
# for complex random effects or random effects with many levels, use
  # mgcv::gamm()   - fits using nlme::lme()
  # gamm4::gamm4() - fits using lmer() or glmer()
# for non-Gaussian models with complex random effects, use gamm4::gamm4()



# Testing model fit -------------------------------------------------------

# simulate data 
set.seed(2)
n <- 400
x1 <- rnorm(n)
x2 <- rnorm(n)
y_val <- 1 + 2*cos(pi*x1) + 2/(1+exp(-5*x2))
y_norm <- y_val + rnorm(n, 0, 0.5)
y_nbinom <- rnbinom(n, mu = exp(y_val), size = 10)
y_binom <- rbinom(n, 1, prob = plogis(y_val))

plot(x1, y_norm)
plot(x2, y_norm)
 
plot(x1, y_nbinom)
plot(x2, y_nbinom)


## Gaussian data ----

# fit Gaussian model with basis size = 4
norm_m1 <- gam(y_norm ~ s(x1, k = 4) + s(x2, k = 4), method = 'REML')
gam.check(norm_m1)  # k-index < 1 and p < 0.05 indicates k is too small 

# increase k of x1
norm_m2 <- gam(y_norm ~ s(x1, k = 12) + s(x2, k = 5), method = 'REML')
gam.check(norm_m2)  # now k needs to be bigger for x2

# increase k of x2
norm_m3 <- gam(y_norm ~ s(x1, k = 12) + s(x2, k = 12), method = 'REML')
gam.check(norm_m3, rep = 500)  # pretty good


#gam.check makes 4 plots
# 1) Quantile-quantile plots of residuals, if correct should follow a 1 to 1 line
# 2) Histogram of residuals -> want it to look Gaussian 
# 3) Residuals vs linear predictors -> want it to look even around 0
# 4) Response vs fitted values -> want a mostly 1 to 1 fit
# Note: gam.check uses deviance residuals by default

# can also fit gam to deviance residuals
# check if model finds any nonlinearity (edf > 1) and increase k for those variables

## Negative binomial data

# fit Poisson model
pois_model <- gam(y_nbinom ~ s(x1, k = 12) + s(x2, k = 12), family = poisson, method = 'REML')
gam.check(pois_model, rep = 500)   # clear problem with QQ plot

# fit negative binomial model
nb_model <- gam(y_nbinom ~ s(x1, k = 12) + s(x2, k = 12), family = nb, method = 'REML')
gam.check(nb_model, rep = 500)   # QQ plot looks better

# can also get diagnostic plots using gratia package
gratia::appraise(nb_model, method = "simulate")


# Model selection ---------------------------------------------------------

# default penalty only penalizes wiggly basis functions (doesn't penalize a constant or linear function), 
# so can't completely remove terms from the model

# two ways to penalize non-wiggly part of smooth (i.e., do model selection)
    #1. double penalty approach: gam(..., select = T)
          # tends to work best
          # applies penalty to all terms in model
          # doubles number of smoothness parameters to estimate
    #2. use special smoother with shrinkage
          # thin plate spline: s(..., bs = 'ts')
          # cubic spline: s(..., bs = 'cs')

# simulate data with 4 known function (1 flat) and 2 spurious covariates
set.seed(3)
n <- 200
dat <- gamSim(1, n = n, scale = 0.15, dist = 'poisson', verbose = F) |> 
  mutate(x4 = runif(n), 
         x5 = runif(n), 
         f4 = rep(0, n), 
         f5 = rep(0, n))

b <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3) + s(x4) + s(x5), 
         data = dat, family = poisson, method = 'REML', 
         select = T)

summary(b)   # model correctly drops the 3 variables with no effect
gratia::draw(b, scales = 'fixed')



# Confidence intervals ----------------------------------------------------

# 95% intervals plotted by plot.gam() are actually Bayesian credible intervals
# but they have good frequentist coverge properties across the function
# unless the models are over-smoothed... so don't over-smooth

# default confidence intervals fail if smooth is close to a linear function
  # use plot.gam(..., seWithMean = TRUE)
  # or gratia::draw(..., overall_uncertainty = TRUE)

# simulate data with linear function
x <- rnorm(200)
y <- rnorm(200, 0.5*x)

fit <- gam(y ~ s(x), method = 'REML')
plot(fit)                   # wrong confidence intervals
plot(fit, seWithMean = T)   # right confidence intervals


# p-values for smooths ----------------------------------------------------

# test of null hypothesis that the function is flat
# doesn't account for estimation of lambda (smoothing parameter), so p-values are biased low
# test statistic is a form of chi-square statistic, but with complicated degrees of freedom


#This whole section was above me
#Stuff about AIC <- No idea comparing models? 


# Atmospheric CO2 example (time series) -------------------------------------------------

data(co2s)
head(co2s)

ggplot(co2s, aes(x = c.month, y = co2)) + geom_line()


# fit model with single smooth
b <- gam(co2 ~ s(c.month, k=300, bs="cr"), data = co2s, method = 'REML')
summary(b)

# predict next 36 months
pd <- tibble(c.month = 1:nrow(co2s)+36) 
pred <- predict(b, pd, se = T)
pd <- pd |> 
  mutate(fit = pred$fit, 
         upr = fit + 2*pred$se.fit, 
         lwr = fit - 2*pred$se.fit)

ggplot(pd, aes(x = c.month, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.5) +
  geom_line() +
  geom_line(data = co2s, aes(y = co2), col = 'red')
#Model is not actually good for prediction, don't need it so wiggly so we do the next step

# fit model with separate smooth for overall trend and seasonal effect
b2 <- gam(co2 ~ s(c.month, k=200, bs="cr") + s(month, bs = 'cc'), 
          data = co2s, method = 'REML', 
          knots = list(month = c(0.5, 12.5)))
summary(b2)
#Cant set knots to 1 and 12 cuase then january is not diff then december

# predict next 36 months
pd <- tibble(c.month = 1:nrow(co2s)+36,
             month = rep(1:12, 50)[1:nrow(pd)]) 
pred <- predict(b2, pd, se = T)
pd <- pd |> 
  mutate(fit = pred$fit, 
         upr = fit + 2*pred$se.fit, 
         lwr = fit - 2*pred$se.fit)

ggplot(pd, aes(x = c.month, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.5) +
  geom_line() +
  geom_line(data = co2s, aes(y = co2), col = 'red')
#Much more realistic


# Water temperature example (spatial time series) -------------------------

galveston <- read_csv("gbtemp.csv") |> 
  rename_all(tolower) |> 
  mutate(datetime = mdy_hms(paste(date, time)), 
         station_id = factor(station_id),
         doy = yday(datetime), 
         tod = hour(datetime) + minute(datetime)/60)

knots <- list(doy = c(0.5, 366.5))  # 366.5 instead of 365.5 to account for leap years

# warning: take a long time to fit
m <- bam(measurement ~
           s(tod, bs = 'cc', k = 10) + 
           s(doy, bs = 'cc', k = 12) + 
           s(year, k = 30) + 
           s(latitude, longitude, bs = 'ds', k = 100, m = c(1, 0.5)) +
           ti(doy, year, bs = c('cc', 'tp')) + 
           ti(latitude, longitude, tod, d = c(2,1), bs = c('ds', 'tp'), k = c(20, 10)) +
           ti(latitude, longitude, doy, d = c(2,1), bs = c('ds', 'cc'), k = c(25, 12)) +
           ti(latitude, longitude, year, d = c(2,1), bs = c('ds', 'tp'), k = c(25, 15)),
         knots = knots,
         method = 'REML', 
         data = galveston)
         
