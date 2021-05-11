library(tidyverse)
library(nimble)
library(viridis)
library(extrafont)
extrafont::loadfonts()

nimbleOptions(disallow_multivariate_argument_expressions = FALSE)

seed <- 18
set.seed(seed)

#---------------------------------------------------#
# Simulate data from thetalogistic approximation ####
#---------------------------------------------------#

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

pops <- 5
tmax <- 50

mu_r1 <- 0.3
K <- 100
sigma_e2 <- 0.01
sigma_e <- sqrt(sigma_e2)
sigma_d2 <- 0.2
sigma_d <- sqrt(sigma_d2)
theta <- 1.2

r0 <- mu_r1 / (1 - K^(-theta))

gamma <- r0 * theta

# True N

N <- matrix(NA, nrow = pops, ncol = tmax)
N[,1] <- 10

epsilon_r1 <- matrix(NA, nrow = pops, ncol = tmax - 1)

for(i in 1:pops) {
  for(t in 1:(tmax-1)){

    epsilon_r1[i, t] <- rnorm(1, 0, sqrt(sigma_e^2 + ((sigma_d^2) / N[i, t])))

    N[i, t+1] <- predict_N(N_current = N[i, t], mu_r1 = mu_r1, epsilon_r1 = epsilon_r1[i, t],
                           K = K, theta = theta,
                           sigma_e = sigma_e,
                           sigma_d = sigma_d)

  }
}

N <- round(N)

obs_r <- matrix(NA, nrow = pops, ncol = tmax - 1)

for(i in 1:pops) {

  obs_r[i,] <- diff(log(N[i,]))

}

# Observed N
# Draw from Poisson model

obs_N <- matrix(NA, nrow = pops, ncol = tmax)

for(i in 1:pops){

    obs_N[i,] <- rpois(tmax, N[i,])

}

# obs_N <- N

# Plot N
cl <- viridis(pops)

plot(obs_N[1,], type = "l", xlab = "Time", ylab = "N", col = adjustcolor(cl[1], alpha.f = 0.6))

for(r in 2:nrow(N)) {

  lines(obs_N[r,], type = "l", col = adjustcolor(cl[r], alpha.f = 0.6))

}

#-----------------------------------------------#
# Simulate data from random-effects formulation #
#-----------------------------------------------#

# N_sim <- matrix(rep(N, 1000), nrow = 1000, ncol = tmax, byrow = T)
#
# pred_r <- matrix(NA, nrow = nrow(N_sim), ncol = tmax-1)
#
# for(i in 1:nrow(N_sim)) {
#   for(t in 1:(tmax-1)) {
#
#     pred_r[i, t] <- mu_r1 - (mu_r1 * (((N_sim[i, t]^theta) - 1) / ((K^theta) - 1))) + rnorm(1, mean = 0, sd = sqrt(sigma_d2 / N_sim[i, t])) + rnorm(1, mean = 0, sd = sqrt(sigma_e2)) - (0.5 * sigma_e2) - (0.5 * sigma_d2 / N_sim[i, t])
#
#   }
# }
#
# plot(1:(tmax-1), pred_r[1,], type = "l", ylim = c(min(pred_r), max(pred_r))*1.1, col = adjustcolor("black", alpha.f = 0.1),
#      ylab = "r", xlab = "Year")
# for(i in 2:nrow(pred_r)) {
#
#   lines(1:(tmax-1), pred_r[i,], type = "l", col = adjustcolor("black", alpha.f = 0.1))
#
# }
#
# obs_r2 <- rep(NA, tmax-1)
# for(t in 1:(tmax-1)) {
#   obs_r2[t] <- mu_r1 - (0.5 * sigma_e2) - (0.5 * sigma_d2 / N[t]) - (mu_r1 * (((N[t]^theta) - 1) / ((K^theta) - 1)))
# }
#
# lines(1:(tmax-1), obs_r2, type = "l", col = adjustcolor("red", alpha.f = 0.5))



#---------------#
# RE model r ####
#---------------#

predict_N_random <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  for(i in 1:pops){
    for(t in 1:(tmax-1)){

      log(N[i, t+1]) <- log(N[i, t]) + mu_r1[i] * (1 - (((N[i, t]^theta[i]) - 1) / ((K[i]^theta[i]) - 1))) + eps_e[i, t] + eps_d[i, t]

      eps_e[i, t] ~ dnorm(0, var = sigma_e2[i])
      eps_d[i, t] ~ dnorm(0, var = sigma_d2[i] / N[i, t])

    }
  }

  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  for(i in 1:pops) {
    for(t in 1:tmax) {

      #obs_N[i, t] ~ dnorm(N[i, t], sd = 0.00001)
      obs_N[i, t] ~ dpois(N[i, t])

    }
  }

  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  for(i in 1:pops) {

    initial_N[i] ~ dunif(0, max_N[i])
    N[i,1] <- initial_N[i]

    mu_r1[i] ~ dunif(-5, 5)
    sigma_e2[i] ~ dunif(0, 10)
    K[i] ~ dunif(1, max_K[i])
    theta[i] ~ dunif(-10, 10)

  }



  #sigma_obs ~ dunif(0, 100)

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops) {

    gamma[i] <- theta[i] * mu_r1[i] / (1 - K[i]^(-theta[i]))

  }


})

sample_inits1 <- function(){

  list(
    #r0 = rnorm(1, 1, 0.5),
    mu_r1 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 1),
    theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops),
    eps_e = matrix(0, pops, tmax-1),
    eps_d = matrix(0, pops, tmax-1),
    initial_N = obs_N[,1]#,
    #sigma_obs = runif(1, 0, 10)
  )

}

input_data1 <- list(obs_N = obs_N)

input_constants1 <- list(tmax = tmax, pops = pops, max_K = rep(K * 2, pops), max_N = apply(N, 1, max) * 2,
                         sigma_d2 = rep(sigma_d2, pops))

inits1 <- list(sample_inits1(), sample_inits1(), sample_inits1())

params1 <- c("K", "theta", "sigma_e2", "gamma", "mu_r1", "N")

# Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

# Model
start <- Sys.time()

seed1 <- 51
set.seed(seed1)

mod1 <- nimbleMCMC(code = predict_N_random,
                   constants = input_constants1,
                   data = input_data1,
                   inits = inits1,
                   monitors = params1,
                   niter = niter,
                   nburnin = nburnin,
                   thin = nthin,
                   nchains = nchains,
                   setSeed = seed1,
                   samplesAsCodaMCMC = TRUE)
dur1 <- Sys.time() - start


#---------------#
# RE model s ####
#---------------#

predict_N_random_unbiased <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  for(i in 1:pops){
    for(t in 1:(tmax-1)){

      #log(N[t+1]) <- log(N[t]) + r0 * (1 - (N[t] / K)^theta) + eps_e[t] + eps_d[t]
      log(N[i, t+1]) <- log(N[i, t]) + mu_r1[i] * (1 - (((N[i, t]^theta[i]) - 1) / ((K[i]^theta[i]) - 1))) - (sigma_e2[i] / 2) - (sigma_d2[i] / (2 * N[i, t])) + eps_e[i, t] + eps_d[i, t]

      eps_e[i, t] ~ dnorm(0, var = sigma_e2[i])
      eps_d[i, t] ~ dnorm(0, var = sigma_d2[i] / N[i, t])

    }
  }

  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  for(i in 1:pops) {
    for(t in 1:tmax) {

      #obs_N[i, t] ~ dnorm(N[i, t], sd = 0.00001)
      obs_N[i, t] ~ dpois(N[i, t])

    }
  }

  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  for(i in 1:pops) {

    initial_N[i] ~ dunif(0, max_N[i])
    N[i,1] <- initial_N[i]

    mu_r1[i] ~ dunif(-5, 5)
    sigma_e2[i] ~ dunif(0, 10)
    K[i] ~ dunif(1, max_K[i])
    theta[i] ~ dunif(-10, 10)

  }



  #sigma_obs ~ dunif(0, 100)

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops) {

    gamma[i] <- theta[i] * mu_r1[i] / (1 - K[i]^(-theta[i]))

  }


})

sample_inits2 <- function(){

  list(
    #r0 = rnorm(1, 1, 0.5),
    mu_r1 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 1),
    theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops),
    eps_e = matrix(0, pops, tmax-1),
    eps_d = matrix(0, pops, tmax-1),
    initial_N = obs_N[,1]#,
    #sigma_obs = runif(1, 0, 10)
  )

}

input_data2 <- list(obs_N = obs_N)

input_constants2 <- list(tmax = tmax, pops = pops, max_K = rep(K * 2, pops), max_N = apply(N, 1, max) * 2,
                         sigma_d2 = rep(sigma_d2, pops))

params2 <- c("K", "theta", "sigma_e2", "gamma", "mu_r1", "N")

# Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

# Model
# start <- Sys.time()
#
# seed2 <- 32
# set.seed(seed2)
#
# inits2 <- list(sample_inits2(), sample_inits2(), sample_inits2())
#
# mod2 <- nimbleMCMC(code = predict_N_random_unbiased,
#                    constants = input_constants2,
#                    data = input_data2,
#                    inits = inits2,
#                    monitors = params2,
#                    niter = niter,
#                    nburnin = nburnin,
#                    thin = nthin,
#                    nchains = nchains,
#                    setSeed = seed2,
#                    samplesAsCodaMCMC = TRUE)
#
# dur2 <- Sys.time() - start

# Three separate chains
start <- Sys.time()
nchains <- 1

seed2a <- 34
set.seed(seed2a)

inits2a <- list(sample_inits2())

mod2a <- nimbleMCMC(code = predict_N_random_unbiased,
                   constants = input_constants2,
                   data = input_data2,
                   inits = inits2a,
                   monitors = params2,
                   niter = niter,
                   nburnin = nburnin,
                   thin = nthin,
                   nchains = nchains,
                   setSeed = seed2a,
                   samplesAsCodaMCMC = TRUE)

seed2b <- 45
set.seed(seed2b)

inits2b <- list(sample_inits2())

# RUN 2B again

mod2b <- nimbleMCMC(code = predict_N_random_unbiased,
                    constants = input_constants2,
                    data = input_data2,
                    inits = inits2b,
                    monitors = params2,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    setSeed = seed2b,
                    samplesAsCodaMCMC = TRUE)

seed2c <- seed2b + 1
set.seed(seed2c)

inits2c <- list(sample_inits2())

mod2c <- nimbleMCMC(code = predict_N_random_unbiased,
                    constants = input_constants2,
                    data = input_data2,
                    inits = inits2c,
                    monitors = params2,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    setSeed = seed2c,
                    samplesAsCodaMCMC = TRUE)

dur2 <- Sys.time() - start

mod2 <- coda::mcmc.list(chain1 = mod2a, chain2 = mod2b, chain3 = mod2c)

# ggsave(filename = "post.pdf",
#        plot = gridExtra::marrangeGrob(purrr::map(colnames(as.matrix(mod2)),
#                                                  ~{bayesplot::mcmc_combo(mod2, pars = .x)}),
#                                       nrow=2, ncol=1),
#        width = 8, height = 8)


#------------------------#
# Approximation model ####
#------------------------#

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

      obs_r[i, t] ~ dnorm(pred_r[i, t], var = var_r1[i, t])

      var_r1[i, t]  <- sigma_e2[i] + sigma_d2[i] / N[i, t]

    }
  }



  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  for(i in 1:pops){

    mu_r1[i] ~ dunif(-5, 5)
    sigma_e2[i] ~ dunif(0, 10)
    K[i] ~ dunif(1, max_K[i])
    theta[i] ~ dunif(-10, 10)

  }


  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops){

    gamma[i] <- theta[i] * mu_r1[i] / (1 - K[i]^(-theta[i]))

  }


})

sample_inits3 <- function(){

  list(
    mu_r1 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 1),
    theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops)
  )

}

## Sample initial values
#inits_b3 <- list(sample_inits_b3())

#-----------------------------#
# NIMBLE model and MCMC setup
#-----------------------------#

## Set data and constants
input_data3 <- list(N = obs_N, obs_r = obs_r)

input_constants3 <- list(tmax = tmax, pops = pops, max_K = rep(K * 2, pops), sigma_d2 = rep(sigma_d2, pops))

## Set parameters to monitor
params3 <- c("K", "theta", "sigma_e2", "mu_r1", "gamma")

# Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

start <- Sys.time()

inits3 <- list(sample_inits3(), sample_inits3(), sample_inits3())
mod3 <- nimbleMCMC(code = predict_r_mult_nimble,
                   constants = input_constants3,
                   data = input_data3,
                   inits = inits3,
                   monitors = params3,
                   niter = niter,
                   nburnin = nburnin,
                   thin = nthin,
                   nchains = nchains,
                   #setSeed = mySeed,
                   samplesAsCodaMCMC = TRUE)
dur3 <- Sys.time() - start


#--------------#
# MLE model ####
#--------------#

mod4 <- purrr::map_dfc(.x = seq_len(pops),
                       .f = ~{

                         mod <- thetalogistic.est(N = obs_N[.x,],
                                                  year = 1:tmax,
                                                  sd = sigma_d2,
                                                  nboot = 1500)

                         tibble::tibble(
                           !!paste0("sigma_e2", "[", .x, "]") := mod$boot.se,
                           !!paste0("K", "[", .x, "]") := mod$boot.K,
                           !!paste0("mu_r1", "[", .x, "]") := mod$boot.r1,
                           !!paste0("gamma", "[", .x, "]") := mod$boot.gamma,
                           !!paste0("theta", "[", .x, "]") := NA
                         )

                       })


#---------------------#
# RE model r - hyp ####
#---------------------#

predict_N_random_hyp <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  for(i in 1:pops){
    for(t in 1:(tmax-1)){

      log(N[i, t+1]) <- log(N[i, t]) + mu_r1[i] * (1 - (((N[i, t]^theta[i]) - 1) / ((K[i]^theta[i]) - 1))) + eps_e[i, t] + eps_d[i, t]

      eps_e[i, t] ~ dnorm(0, var = sigma_e2[i])
      eps_d[i, t] ~ dnorm(0, var = sigma_d2[i] / N[i, t])

    }
  }

  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  for(i in 1:pops) {
    for(t in 1:tmax) {

      #obs_N[i, t] ~ dnorm(N[i, t], sd = 0.00001)
      obs_N[i, t] ~ dpois(N[i, t])

    }
  }

  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  for(i in 1:pops) {

    initial_N[i] ~ dunif(0, max_N[i])
    N[i,1] <- initial_N[i]

    # Priors for population parameters
    mu_r1[i] <- mean_mu_r1 + epsilon_mu_r1[i]
    epsilon_mu_r1[i] ~ dnorm(0, sd = sd_mu_r1)

    log(sigma_e2[i]) <- log(mean_sigma_e2) + epsilon_sigma_e2[i]
    epsilon_sigma_e2[i] ~ dnorm(0, sd = sd_sigma_e2)

    K[i] ~ dunif(1, max_K[i])
    theta[i] ~ dunif(-10, 10)

  }

  mean_mu_r1 ~ dunif(-5, 5)
  sd_mu_r1 ~ dunif(0, 10)

  mean_sigma_e2 ~ dunif(0, 10)
  sd_sigma_e2 ~ dunif(0, 10)

  #sigma_obs ~ dunif(0, 100)

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops) {

    gamma[i] <- theta[i] * mu_r1[i] / (1 - K[i]^(-theta[i]))

  }


})

sample_inits5 <- function(){

  list(
    #r0 = rnorm(1, 1, 0.5),
    mean_mu_r1 = rnorm(1, 1, 0.5),
    sd_mu_r1 = runif(1, 0, 1),
    epsilon_mu_r1 = rep(0, pops),
    mean_sigma_e2 = runif(1, 0, 1),
    sd_sigma_e2 = runif(1, 0, 1),
    epsilon_sigma_e2 = rep(0, pops),
    theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops),
    eps_e = matrix(0, pops, tmax-1),
    eps_d = matrix(0, pops, tmax-1),
    initial_N = obs_N[,1]#,
    #sigma_obs = runif(1, 0, 10)
  )

}

input_data5 <- list(obs_N = obs_N)

input_constants5 <- list(tmax = tmax, pops = pops, max_K = rep(K * 2, pops), max_N = apply(N, 1, max) * 2,
                         sigma_d2 = rep(sigma_d2, pops))

params5 <- c("K", "theta", "mean_sigma_e2", "sd_sigma_e2", "sigma_e2", "gamma", "mean_mu_r1", "sd_mu_r1", "mu_r1", "N")

# Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

# Model
start <- Sys.time()

seed5 <- 403
set.seed(seed5)

inits5 <- list(sample_inits5(), sample_inits5(), sample_inits5())

mod5 <- nimbleMCMC(code = predict_N_random_hyp,
                   constants = input_constants5,
                   data = input_data5,
                   inits = inits5,
                   monitors = params5,
                   niter = niter,
                   nburnin = nburnin,
                   thin = nthin,
                   nchains = nchains,
                   setSeed = seed5,
                   samplesAsCodaMCMC = TRUE)

dur5 <- Sys.time() - start


#---------------------#
# RE model s - hyp ####
#---------------------#

predict_N_random_unbiased_hyp <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  for(i in 1:pops){
    for(t in 1:(tmax-1)){

      #log(N[t+1]) <- log(N[t]) + r0 * (1 - (N[t] / K)^theta) + eps_e[t] + eps_d[t]
      log(N[i, t+1]) <- log(N[i, t]) + mu_r1[i] * (1 - (((N[i, t]^theta[i]) - 1) / ((K[i]^theta[i]) - 1))) - (sigma_e2[i] / 2) - (sigma_d2[i] / (2 * N[i, t])) + eps_e[i, t] + eps_d[i, t]

      eps_e[i, t] ~ dnorm(0, var = sigma_e2[i])
      eps_d[i, t] ~ dnorm(0, var = sigma_d2[i] / N[i, t])

    }
  }

  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  for(i in 1:pops) {
    for(t in 1:tmax) {

      #obs_N[i, t] ~ dnorm(N[i, t], sd = 0.00001)
      obs_N[i, t] ~ dpois(N[i, t])

    }
  }

  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  for(i in 1:pops) {

    initial_N[i] ~ dunif(0, max_N[i])
    N[i,1] <- initial_N[i]

    # Priors for population parameters
    mu_r1[i] <- mean_mu_r1 + epsilon_mu_r1[i]
    epsilon_mu_r1[i] ~ dnorm(0, sd = sd_mu_r1)

    log(sigma_e2[i]) <- log(mean_sigma_e2) + epsilon_sigma_e2[i]
    epsilon_sigma_e2[i] ~ dnorm(0, sd = sd_sigma_e2)

    K[i] ~ dunif(1, max_K[i])
    theta[i] ~ dunif(-10, 10)

  }

  mean_mu_r1 ~ dunif(-5, 5)
  sd_mu_r1 ~ dunif(0, 10)

  mean_sigma_e2 ~ dunif(0, 10)
  sd_sigma_e2 ~ dunif(0, 10)

  #sigma_obs ~ dunif(0, 100)

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops) {

    gamma[i] <- theta[i] * mu_r1[i] / (1 - K[i]^(-theta[i]))

  }


})

sample_inits6 <- function(){

  list(
    #r0 = rnorm(1, 1, 0.5),
    mean_mu_r1 = rnorm(1, 1, 0.5),
    sd_mu_r1 = runif(1, 0, 1),
    epsilon_mu_r1 = rep(0, pops),
    mean_sigma_e2 = runif(1, 0, 1),
    sd_sigma_e2 = runif(1, 0, 1),
    epsilon_sigma_e2 = rep(0, pops),
    theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops),
    eps_e = matrix(0, pops, tmax-1),
    eps_d = matrix(0, pops, tmax-1),
    initial_N = obs_N[,1]#,
    #sigma_obs = runif(1, 0, 10)
  )

}

input_data6 <- list(obs_N = obs_N)

input_constants6 <- list(tmax = tmax, pops = pops, max_K = rep(K * 2, pops), max_N = apply(N, 1, max) * 2,
                         sigma_d2 = rep(sigma_d2, pops))

params6 <- c("K", "theta", "mean_sigma_e2", "sd_sigma_e2", "sigma_e2", "gamma", "mean_mu_r1", "sd_mu_r1", "mu_r1", "N")

# Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

# Model
start <- Sys.time()

seed6 <- 78
set.seed(seed6)

inits6 <- list(sample_inits6(), sample_inits6(), sample_inits6())

mod6 <- nimbleMCMC(code = predict_N_random_unbiased_hyp,
                   constants = input_constants6,
                   data = input_data6,
                   inits = inits6,
                   monitors = params6,
                   niter = niter,
                   nburnin = nburnin,
                   thin = nthin,
                   nchains = nchains,
                   setSeed = seed6,
                   samplesAsCodaMCMC = TRUE)

dur6 <- Sys.time() - start

#--------------#
# X model s ####
#--------------#

predict_X <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  for(i in 1:pops){
    for(t in 1:(tmax-1)){

      N[i, t+1] <- exp(X[i, t+1])

      EX[i, t+1] <- X[i, t] + mu_r1[i] * (1 - (((N[i, t]^theta[i]) - 1) / ((K[i]^theta[i]) - 1))) - (sigma_e2[i] / 2) - (sigma_d2[i] / (2 * N[i, t]))

      X[i, t+1] ~ dnorm(EX[i, t+1], var = sigma_e2[i] + sigma_d2[i] / N[i, t])

      # eps now becomes a derived parameter
      # I keep it here instead of in "derived parameter section" because
      # we later may want to fit a distribution to epsilons (but how to do that
      # without making eps an parameter???)
      # eps[i,t+1] or eps[i, t] is a matter of preference / definition -
      # but important to be aware of / remember the indexing
      eps[i, t] <- X[i, t+1] - EX[i, t+1]
      demcomp[i, t] <- sigma_d2[i] / N[i, t]

      # eps_e[i, t] ~ dnorm(0, var = sigma_e2[i])
      # eps_d[i, t] ~ dnorm(0, var = sigma_d2[i] / N[i, t])

    }
  }

  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  for(i in 1:pops) {
    for(t in 1:tmax) {

      #obs_N[i, t] ~ dnorm(N[i, t], sd = 0.00001)
      obs_N[i, t] ~ dpois(N[i, t])

    }
  }

  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  for(i in 1:pops) {

    # initial_N[i] ~ dunif(0, max_N[i])
    # N[i,1] <- initial_N[i]

    X[i, 1] ~ dunif(0, max_X[i])
    N[i, 1] <- exp(X[i, 1])

    r0[i] ~ dunif(0, 5)
    sigma_e2[i] ~ dunif(0, 10)
    K[i] ~ dunif(10, max_K[i])
    theta[i] ~ dunif(-10, 10)

  }



  #sigma_obs ~ dunif(0, 100)

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops) {

    mu_r1[i] <- r0[i] / (1 - K[i]^(-theta[i]))
    gamma[i] <- theta[i] * mu_r1[i] / (1 - K[i]^(-theta[i]))

  }


})

sample_inits7 <- function(){

  list(
    r0 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 1),
    theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops),
    X = log(obs_N)
  )

}

input_data7 <- list(obs_N = obs_N)

input_constants7 <- list(tmax = tmax, pops = pops, max_K = rep(K * 2, pops), max_X = apply(obs_N, 1, max) * 2,
                         sigma_d2 = rep(sigma_d2, pops))

inits7 <- list(sample_inits7(), sample_inits7(), sample_inits7())

params7 <- c("K", "theta", "sigma_e2", "gamma", "mu_r1", "r0")

# Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

# Model
start <- Sys.time()

seed7 <- 191
set.seed(seed7)

mod7 <- nimbleMCMC(code = predict_X,
                   constants = input_constants7,
                   data = input_data7,
                   inits = inits7,
                   monitors = params7,
                   niter = niter,
                   nburnin = nburnin,
                   thin = nthin,
                   nchains = nchains,
                   setSeed = seed7,
                   samplesAsCodaMCMC = TRUE)
dur7 <- Sys.time() - start


#--------------------#
# X model s - eps ####
#--------------------#

predict_X_random <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  for(i in 1:pops){
    for(t in 1:(tmax-1)){

      N[i, t+1] <- exp(X[i, t+1])

      X[i, t+1] <- X[i, t] + mu_r1[i] * (1 - (((N[i, t]^theta[i]) - 1) / ((K[i]^theta[i]) - 1))) - (sigma_e2[i] / 2) - (sigma_d2[i] / (2 * N[i, t])) + eps[i, t]

      eps[i, t] ~ dnorm(0, var = sigma_e2[i] + sigma_d2[i] / N[i, t])

      # eps now becomes a derived parameter
      # I keep it here instead of in "derived parameter section" because
      # we later may want to fit a distribution to epsilons (but how to do that
      # without making eps an parameter???)
      # eps[i,t+1] or eps[i, t] is a matter of preference / definition -
      # but important to be aware of / remember the indexing
      #eps[i, t] <- X[i, t+1] - EX[i, t+1]
      demcomp[i, t] <- sigma_d2[i] / N[i, t]

      # eps_e[i, t] ~ dnorm(0, var = sigma_e2[i])
      # eps_d[i, t] ~ dnorm(0, var = sigma_d2[i] / N[i, t])

    }
  }

  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  for(i in 1:pops) {
    for(t in 1:tmax) {

      #obs_N[i, t] ~ dnorm(N[i, t], sd = 0.00001)
      obs_N[i, t] ~ dpois(N[i, t])

    }
  }

  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  for(i in 1:pops) {

    # initial_N[i] ~ dunif(0, max_N[i])
    # N[i,1] <- initial_N[i]

    X[i, 1] ~ dunif(0, max_X[i])
    N[i, 1] <- exp(X[i, 1])

    r0[i] ~ dunif(0, 5)
    sigma_e2[i] ~ dunif(0, 10)
    K[i] ~ dunif(10, max_K[i])
    theta[i] ~ dunif(-10, 10)

  }



  #sigma_obs ~ dunif(0, 100)

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops) {

    mu_r1[i] <- r0[i] / (1 - K[i]^(-theta[i]))
    gamma[i] <- theta[i] * mu_r1[i] / (1 - K[i]^(-theta[i]))

  }


})

sample_inits8 <- function(){

  list(
    r0 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 1),
    theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops),
    X = log(obs_N),
    eps = matrix(0, pops, tmax-1)
  )

}

input_data8 <- list(obs_N = obs_N)

input_constants8 <- list(tmax = tmax, pops = pops, max_K = rep(K * 2, pops), max_X = apply(obs_N, 1, max) * 2,
                         sigma_d2 = rep(sigma_d2, pops))

inits8 <- list(sample_inits8(), sample_inits8(), sample_inits8())

params8 <- c("K", "theta", "sigma_e2", "gamma", "mu_r1", "r0")

# Set MCMC parameters
niter <- 30000
nburnin <- 5000
nthin <- 10
nchains <- 3

# Model
start <- Sys.time()

seed8 <- 404
set.seed(seed8)

mod8 <- nimbleMCMC(code = predict_X_random,
                   constants = input_constants8,
                   data = input_data8,
                   inits = inits8,
                   monitors = params8,
                   niter = niter,
                   nburnin = nburnin,
                   thin = nthin,
                   nchains = nchains,
                   setSeed = seed8,
                   samplesAsCodaMCMC = TRUE)
dur8 <- Sys.time() - start


#-----------------#
# Plot outputs ####
#-----------------#

plot_theta_mult(filename = "Approx-RE-MLE-outputs2",
                model_list = list("Approx" = mod3, "MLE" = mod4, "RE-r" = mod1, "RE-s" = mod2,
                                  "RE-r-hyp" = mod5, "RE-s-hyp" = mod6))


plot_theta_mult(filename = "Approx-MLE-X-outputs",
                model_list = list("Approx" = mod3, "MLE" = mod4, "X-s" = mod7))
