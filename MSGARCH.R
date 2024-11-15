# Loading package pacman and use its function p_load for loading all required packages
library(pacman)
p_load(tibbletime, timetk, tidyverse, MSGARCH, zoo, DEoptim, Rcpp)

# Import data set
Yen_Dollar <- read.csv("Yen_Dollar.csv") %>% 
  mutate(date = ymd(date))

# Convert Yen per Dollar object to a tibble
Yen_Dollar <- Yen_Dollar %>% 
  as_tibble()

head(Yen_Dollar)

# Calculate monthly log returns
Yen_Dollar <- Yen_Dollar %>% 
  mutate(return = (log(Yen_Dollar$Yen_per_Dollar) - log(lag(Yen_Dollar$Yen_per_Dollar)))) %>% 
  na.omit() %>% 
  select(-Yen_per_Dollar)

head(Yen_Dollar)

# Convert Yen-Dollar object to dataframe object
Yen_Dollar_df <- data.frame(Date = c(Yen_Dollar$date), return = c(Yen_Dollar$return))
head(Yen_Dollar_df)

# Convert Yen-Dollar dataframe object to zoo object
Yen_Dollar_z <- read.zoo(Yen_Dollar_df)
head(Yen_Dollar_z)            

# Set up a 3-regime Markov-switching specification with the help of variable K
# MS(3)-GARCH(1,1)- Student
msgarch <- CreateSpec(variance.spec = list(model = c("sGARCH")),
           distribution.spec = list(distribution = c("std")), 
           constraint.spec = list(regime.const = "nu"),
           switch.spec = list(do.mix = FALSE, K = 3))

# fit the model on data using the method which computes the Deviance Information Criterion (DIC) from a fit object created with the FitMCMC function
set.seed(123)

msgarch_dic <- FitMCMC(spec = msgarch, data = Yen_Dollar_z, ctr = list(nburn = 500L, nmcmc = 500L))

# Compute DIC
DIC(msgarch_dic)

#Check on the result of the DIC method
summary(msgarch_dic)

# Extract volatility of the model
msgarch_dic_vol <- Volatility(msgarch_dic)

# Make a plot of the volatility
plot(msgarch_dic_vol)

# Set up function for optimization from the stats package
f_custom_optim <- function(vPw, f_nll, spec, data, do.plm){
out <- stats::optim(vPw, f_nll, spec = spec, data = data,
  do.plm = do.plm, method = "Nelder-Mead")
return(out)
}

set.seed(123)

# Run optimization to the Yen-Dollar data set
msgarch_fit_optim <- FitML(msgarch, data = Yen_Dollar_z, ctr = list(OptimFUN = f_custom_optim))

# Check on the result of the optimization process
summary(msgarch_fit_optim)

# Check on coefficients of the msgarch model
coef(msgarch_fit_optim)

# Extract volatility of the model
msgarch_fit_optim_vol <- Volatility(msgarch_fit_optim)

# Make a plot of the volatility
plot(msgarch_fit_optim_vol)

# Set up Bayesian approach
msgarch_fit_mcmc <- FitMCMC(spec = msgarch, data = Yen_Dollar_z, ctr = list(nburn = 500L, nmcmc = 500L, nthin = 1L))

# Check on model result
summary(msgarch_fit_mcmc)

# Check on coefficients of Bayesian approach
summary(msgarch_fit_mcmc)

coef(msgarch_fit_mcmc)

# Extract volatility of the model
msgarch_fit_mcmc_vol <- Volatility(msgarch_fit_mcmc)

# Make a plot of the volatility
plot(msgarch_fit_mcmc_vol)

# Optimization of MSGARCH with DEoptim package setting up a function
f_DEoptim <- function(vPw, f_nll, spec, data, do.plm) {
  
  NP = 15 * length(vPw)
  
  mInitialPop = matrix(rep(vPw, NP) * rnorm(length(vPw) * NP) * 0.1, nrow = NP, byrow = TRUE)
  tmp <- DEoptim::DEoptim(f_nll, lower = rep(-10, length(vPw)), upper = rep(10, length(vPw)),
    control = DEoptim.control(initialpop = mInitialPop, NP = NP, trace = TRUE, itermax = 500),
    spec = spec, data = data, do.plm = do.plm)
  
  vPw_optim <- tmp$optim$bestmem
  dnllk <- tmp$optim$bestval
  names(vPw_optim) = names(vPw)
  
  out = list(par = vPw_optim, value = dnllk)
  
  return(out)
}

# Fit DEoptim to Yen-Dollar data set
msgarch_fit_DE = FitML(spec = msgarch, data = Yen_Dollar_z ,ctr = list(OptimFUN = f_DEoptim))

# Check for summary of the msgarch model
summary(msgarch_fit_DE)

coef(msgarch_fit_DE)

# Extract volatility of the model
msgarch_fit_DE_vol <- Volatility(msgarch_fit_DE)

# Make a plot of the volatility
plot(msgarch_fit_DE_vol)

# MSGARCH simulation
# simulation from ML fit
msgarch_ML_fit <- FitML(msgarch, data = Yen_Dollar_z, ctr = list(OptimFUN = f_custom_optim))

set.seed(1234)

msgarch_sim_ML <- simulate(object = msgarch_ML_fit, nsim = 1L, nahead = 1000L,
  nburn = 500L)

head(msgarch_sim_ML)

# Plot ML simulation 
plot(msgarch_sim_ML)

# simulation from MCMC fit
msgarch_fit_MCMC <- FitMCMC(spec = msgarch, data = Yen_Dollar_z)

set.seed(1234)

msgarch_sim_MCMC <- simulate(object = msgarch_fit_MCMC, nahead = 100L, nburn = 500L)

head(msgarch_sim_MCMC)

# Plot of the MCMC simulation
plot(msgarch_sim_MCMC)






