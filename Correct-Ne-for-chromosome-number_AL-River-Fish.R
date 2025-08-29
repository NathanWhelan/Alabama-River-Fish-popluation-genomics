# Process NeEstimator results:
# 9/27/2020
# Note "observed" = "overall" in NeEstimator output
rm(list=ls())
setwd("~/Fish/AL-River_migratory-fish/white-crappie/")

dat <- read.csv("NeEstimator_white-crappie_all.csv", stringsAsFactors = F)
dat

# STEP 1
# ------

# Start by recalculating r^2' and Ne since NeEstimator does not output many decimal places, and this can make a difference downstream:

dat$r2p <- dat$r2 - dat$exp_r2

# NOTE: Ne's with negative square root argument have that value set to zero
dat$Ne_new <- ifelse((0.308^2-2.08*dat$r2p)>0, (0.308+sqrt(0.308^2-2.08*dat$r2p))/(2*dat$r2p), (0.308+0)/(2*dat$r2p))

# ignore warnings...



# STEP 2
# ------

# Ne correction factor from Waples, Larson, Waples. 2016 Heredity
# y = 0.09775 + 0.21888*ln (Chr)
# where Chr = # of chromosome pairs (14 in plethodon, 13 in Rana pretiosa)

# Calculate CF then just divide Ne by CF

chromosomes <- 24

# correction factor:
CF <- 0.09775 + (0.21888*log(chromosomes))

dat$NeCorrected <- apply(dat["Ne_new"], 1, function(x) x/CF)



# STEP 3
# ------

# Adjust the jackknife CIs - calculated from the raw data:
# Use the Effective Degrees of Freedom output only in the TABULAR format!

# (1) quantile chi-square of individual comparisons at 0.025 and 0.975
dat$df025 <- apply(dat["EffDF"], 1, function(x) qchisq(p=0.025, df=x, lower.tail=FALSE))
dat$df975 <- apply(dat["EffDF"], 1, function(x) qchisq(p=0.975, df=x, lower.tail=FALSE))

# (2) calculate adjusted r2'
dat$adjR2p <- apply(dat["NeCorrected"], 1, function(x) (0.308/x)-0.52/x^2)

# (3) calculate adjusted overall r2 from #2
dat$adjobR2 <- dat$adjR2p + dat$exp_r2

# (4) apply adjusted overall r2 (#3) to recalculate CI r2
dat$r2025 <- dat$EffDF*dat$adjobR2/dat$df025
dat$r2975 <- dat$EffDF*dat$adjobR2/dat$df975

# (5) apply recalculated CI r2 to calculate r2'
dat$r2p025 <- dat$r2025-dat$exp_r2
dat$r2p975 <- dat$r2975-dat$exp_r2
  
# (6) finally, calc the CI Nes

# NOTE: low CIs with negative square root argument have that value set to zero
dat$lowCIcorrected <- ifelse((0.308^2-2.08*dat$r2p975)>0, (0.308+sqrt(0.308^2-2.08*dat$r2p975))/(2*dat$r2p975), (0.308+0)/(2*dat$r2p975))
dat$highCIcorrected  <- (0.308+sqrt(0.308^2-2.08*dat$r2p025))/(2*dat$r2p025)

write.csv(dat, "NeEst_white-crappie_all_corrected.csv", row.names = F)
# NEGATIVE high CIs are "infinity" estimates
