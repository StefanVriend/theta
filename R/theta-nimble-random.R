library(tidyverse)
library(nimble)
library(viridis)
# library(extrafont)
# extrafont::loadfonts(device = "win")

nimbleOptions(disallow_multivariate_argument_expressions = FALSE)

#------------------------------------------------#
# Simulate data from thetalogistic approximation #
#------------------------------------------------#

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

pops <- 1
tmax <- 75

mu_r1 <- 0.5
K <- 150
sigma_e2 <- 0.01
sigma_e <- sqrt(sigma_e2)
sigma_d2 <- 0.2
sigma_d <- sqrt(sigma_d2)
theta <- 1.2

N <- rep(NA, tmax)
N[1] <- 10

epsilon_r1 <- rep(NA, tmax-1)

for(t in 1:(tmax-1)){

  epsilon_r1[t] <- rnorm(1, 0, sqrt(sigma_e^2 + ((sigma_d^2) / N[t])))

  N[t+1] <- predict_N(N_current = N[t], mu_r1 = mu_r1, epsilon_r1 = epsilon_r1[t],
                      K = K, theta = theta,
                      sigma_e = sigma_e,
                      sigma_d = sigma_d)

}

N <- round(N)

obs_r <- diff(log(N))


#-----------------------------------------------#
# Simulate data from random-effects formulation #
#-----------------------------------------------#

N_sim <- matrix(rep(N, 1000), nrow = 1000, ncol = tmax, byrow = T)

pred_r <- matrix(NA, nrow = nrow(N_sim), ncol = tmax-1)

for(i in 1:nrow(N_sim)) {
  for(t in 1:(tmax-1)) {

    pred_r[i, t] <- mu_r1 - (mu_r1 * (((N_sim[i, t]^theta) - 1) / ((K^theta) - 1))) + rnorm(1, mean = 0, sd = sqrt(sigma_d2 / N_sim[i, t])) + rnorm(1, mean = 0, sd = sqrt(sigma_e2))

  }
}

plot(1:(tmax-1), pred_r[1,], type = "l", ylim = c(min(pred_r), max(pred_r))*1.1, col = adjustcolor("black", alpha.f = 0.1),
     ylab = "r", xlab = "Year")
for(i in 2:nrow(pred_r)) {

  lines(1:(tmax-1), pred_r[i,], type = "l", col = adjustcolor("black", alpha.f = 0.1))

}

obs_r2 <- rep(NA, tmax-1)
for(t in 1:(tmax-1)) {
  obs_r2[t] <- mu_r1 - (0.5 * sigma_e2) - (0.5 * sigma_d2 / N[t]) - (mu_r1 * (((N[t]^theta) - 1) / ((K^theta) - 1)))
}

lines(1:(tmax-1), obs_r2, type = "l", col = adjustcolor("red", alpha.f = 0.5))


#--------------#
# NIMBLE model #
#--------------#

predict_N_random <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  for(t in 1:(tmax-1)){

    log(N[t+1]) <- log(N[t]) + r0 * (1 - (N[t] / K)^theta) + eps_e[t] + eps_d[t]

    eps_e[t] ~ dnorm(0, var = sigma_e2)
    eps_d[t] ~ dnorm(0, var = sigma_d2[t] / N[t])

  }

  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  for(t in 1:tmax) {

    #obs_x[t] ~ dnorm(x[t], sd = sigma_obs)
    obs_N[t] ~ dpois(N[t])

  }

  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  initial_N ~ dunif(0, max_N)
  N[1] <- initial_N

  r0 ~ dunif(-5, 5)
  sigma_e2 ~ dunif(0, 10)
  K ~ dunif(1, max_K)
  theta ~ dunif(-10, 10)

  #sigma_obs ~ dunif(0, 100)

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  gamma <- r0 * theta

})

sample_inits <- function(){

  list(
    r0 = rnorm(1, 1, 0.5),
    sigma_e2  = runif(1, 0, 1),
    theta = rnorm(1, 2, 0.5),
    K = K,
    eps_e = rep(0, tmax-1),
    eps_d = rep(0, tmax-1),
    initial_N = N[1]#,
    #sigma_obs = runif(1, 0, 10)
  )

}

input_data <- list(obs_N = N)

input_constants <- list(tmax = tmax, max_K = K * 2, max_N = max(N) * 2, sigma_d2 = rep(sigma_d2, tmax))

inits <- list(sample_inits(), sample_inits(), sample_inits())

params <- c("K", "theta", "sigma_e2", "r0", "gamma")

# Set MCMC parameters
niter <- 10
nburnin <- 5
nthin <- 1
nchains <- 3

# Model
rod <- nimbleMCMC(code = predict_N_random,
                  constants = input_constants,
                  data = input_data,
                  inits = inits,
                  monitors = params,
                  niter = niter,
                  nburnin = nburnin,
                  thin = nthin,
                  nchains = nchains,
                  #setSeed = mySeed,
                  samplesAsCodaMCMC = TRUE)
