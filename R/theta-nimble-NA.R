library(tidyverse)
library(nimble)
library(viridis)

nimbleOptions(disallow_multivariate_argument_expressions = FALSE)


#-----------------------------------------------#
# Basic NIMBLE functions for population growth
#-----------------------------------------------#

predict_N <- nimbleFunction(
  run = function(
    N_current = double(),
    mu_r1 = double(),
    epsilon_r1 = double(),
    sigma_e = double(),
    sigma_d = double(),
    theta = double(),
    K = double()#,
    #beta_r1 = double(1),
    #env_cov = double(1)
  ) {

    r1 <- mu_r1 + epsilon_r1
    #r1 <- mu_r1 + sum(beta_r1 * env_cov) + epsilon_r1

    # K <- [some function of beta]

    s <- r1 - (0.5 * sigma_e^2) - (0.5 * sigma_d^2 / N_current)

    r <- s - mu_r1 * (((N_current^theta)-1) / ((K^theta)-1))

    N_next <- exp(log(N_current) + r)

    if(is.na(N_next)) stop('Predicted population size (N_next) is NA')

    returnType(double())
    return(N_next)
  }
)


predict_r <- nimbleFunction(
  run = function(
    N_current = double(),
    mu_r1 = double(),
    sigma_e2 = double(),
    sigma_d2 = double(),
    theta = double(),
    #beta = double(1),
    #env_cov = double(1),
    K = double()
  ) {

    s <- mu_r1 - (0.5 * sigma_e2) - (0.5 * sigma_d2 / N_current) #+ sum(beta * env_cov)

    pred_r <- s - mu_r1 * (((N_current^theta)-1) / ((K^theta)-1))

    if(is.na(pred_r)) stop('Predicted population growth rate (pred_r) is NA')

    returnType(double())
    return(pred_r)
  }
)


#-----------------------------#
# Multiple data simulation ####
#-----------------------------#

## Number of populations
pops <- 5

## Set parameter values
tmax <- 30
mu_r1 <- 0.5 # up to 1
sigma_e2 <- 0.01
sigma_e <- sqrt(sigma_e2)
sigma_d2 <- 0.2
sigma_d <- sqrt(sigma_d2)
K <- 100
theta <- 1.2

## Simulate random effects

## Prepare population vector and set initial population size
N <- matrix(NA, nrow = pops, ncol = tmax)
N[,1] <- 10

epsilon_r1 <- rep(NA, tmax-1)

## Use nimble function to predict population size over time
for(r in 1:nrow(N)) {
  for(t in 1:(tmax-1)){

    epsilon_r1[t] <- rnorm(1, 0, sqrt(sigma_e^2 + ((sigma_d^2) / N[r, t])))

    N[r, t+1] <- predict_N(N_current = N[r, t], mu_r1 = mu_r1, epsilon_r1 = epsilon_r1[t],
                           K = K, theta = theta, #beta.r1 = 0, EnvCov = 0,
                           sigma_e = sigma_e,
                           sigma_d = sigma_d)

  }
}

gamma <- mu_r1*theta/(1-K^(-theta))

obs_r <- matrix(NA, nrow = nrow(N), ncol = tmax-1)

for(r in 1:nrow(N)) {

  obs_r[r,] <- diff(log(N[r,]))

}

## log(N) and observation error
sigma_Y <- 0.1

log_Y <- matrix(NA, nrow = nrow(N), ncol = ncol(N))
for(r in 1:nrow(log_Y)) {

  log_Y[r,] <- rnorm(length(N[r,]), log(N[r,]) - 0.5 * sigma_Y^2, sigma_Y)

}

#---------------------#
# Plot simulated data
#---------------------#

par(mfrow = c(1, 2))

cl <- viridis(pops)

plot(N[1,], type = "l", xlab = "Time", ylab = "N", col = adjustcolor(cl[1], alpha.f = 0.6))

for(r in 2:nrow(N)) {

  lines(N[r,], type = "l", col = adjustcolor(cl[r], alpha.f = 0.6))

}

title(main = bquote(K~"="~.(K)*", "*bar(r)~"="~.(mu_r1)*", "*sigma[e]^2~"="~.(sigma_e2)*", "*sigma[d]^2~"="~.(sigma_d2)))
#legend("bottomright", col = c("black", "red"), lty = 1, legend = c("True N", "Obs N"), bty = "n", cex = 0.75)

# Plot density dependence
nn <- seq(0, max(N[1,]), length.out = 1000)
rr <- mu_r1 - (0.5 * sigma_e2) - mu_r1 * (((nn^theta)-1) / ((K^theta)-1))

plot(N[1, -tmax], diff(log(N[1,])), xlab = "N", ylab = "r", col = "black")
lines(nn, rr, col = "black")

title(main = bquote(theta~"="~.(theta)*", "*gamma~"="~.(round(gamma, 2))))


#-------------------------------#
# ... NIMBLE model code - B2 ####
#-------------------------------#

predict_r_mult_nimble <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  for (i in 1:pops){

    for (t in 1:(tmax-1)){

      pred_r[i, t] <- predict_r(N_current = N[i, t], mu_r1 = mu_r1[i],
                                sigma_e2 = sigma_e2[i], sigma_d2 = sigma_d2[i],
                                theta = theta[i], K = K[i])

    }

  }


  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  for (i in 1:pops){

    for (t in 1:(tmax-1)){

      obs_r[i, t] ~ dnorm(pred_r[i,t], var = var_r1[i, t])

      var_r1[i, t]  <- sigma_e2[i] + ((sigma_d2[i]) / N[i, t])

    }

  }


  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  for(i in 1:pops) {

    #mu_r1[i] ~ dlnorm(meanlog = 0, sdlog = 0.5)
    mu_r1[i] ~ dunif(-5, 5)
    sigma_e2[i] ~ dunif(0, 10)
    K[i] ~ dunif(1, max_K[i])
    theta[i] ~ dunif(-10, 10)

  }


  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops) {

    gamma[i] <- mu_r1[i] * theta[i] / (1 - K[i] ^ (-theta[i]))

  }


})

## Introduce NA into obs_r matrix
obs_r[,7] <- NA

## Function to sample initial values
sample_inits_b2 <- function(){

  # Setting initial values for missing data points
  init.obs_r <- matrix(NA, nrow = nrow(obs_r), ncol = ncol(obs_r))
  init.obs_r[which(is.na(obs_r))] <- 0

  list(
    mu_r1 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 1),
    theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops),
    obs_r = init.obs_r
  )

}

## Sample initial values
#inits_b2 <- list(sample_inits_b2())
inits_b2 <- list(sample_inits_b2(), sample_inits_b2(), sample_inits_b2())

#-----------------------------#
# NIMBLE model and MCMC setup
#-----------------------------#

## Set data and constants
input_data_b2 <- list(N = N, obs_r = obs_r)

input_constants_b2 <- list(tmax = tmax, max_K = rep(K,pops), sigma_d2 = rep(sigma_d2, pops))

## Set parameters to monitor
params_b2 <- c("K", "theta", "sigma_e2", "mu_r1", "gamma")

## Set MCMC parameters
niter <- 10
nburnin <- 0
nthin <- 1
nchains <- 3

#------------#
# Run model
#------------#

#start <- Sys.time()
mod_b2 <- nimbleMCMC(code = predict_r_mult_nimble,
                     constants = input_constants_b2,
                     data = input_data_b2,
                     inits = inits_b2,
                     monitors = params_b2,
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     nchains = nchains,
                     #setSeed = mySeed,
                     samplesAsCodaMCMC = TRUE)
#dur_b2 <- Sys.time() - start

#coda::gelman.diag(mod_b2)

