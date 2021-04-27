# fixed theta = 1

library(tidyverse)
library(nimble)
library(viridis)
library(extrafont)
extrafont::loadfonts()

nimbleOptions(disallow_multivariate_argument_expressions = FALSE)

seed <- 189
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
theta <- 1

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

# obs_N <- matrix(NA, nrow = pops, ncol = tmax)
#
# for(i in 1:pops){
#
#     obs_N[i,] <- rpois(tmax, N[i,])
#
# }

obs_N <- N

# Plot N
cl <- viridis(pops)

plot(N[1,], type = "l", xlab = "Time", ylab = "N", col = adjustcolor(cl[1], alpha.f = 0.6))

for(r in 2:nrow(N)) {

  lines(N[r,], type = "l", col = adjustcolor(cl[r], alpha.f = 0.6))

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

predict_N_random_fixed <- nimbleCode({

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
    #theta[i] ~ dunif(-10, 10)

  }



  #sigma_obs ~ dunif(0, 100)

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops) {

    gamma[i] <- theta[i] * mu_r1[i] / (1 - K[i]^(-theta[i]))

  }


})

sample_inits10 <- function(){

  list(
    #r0 = rnorm(1, 1, 0.5),
    mu_r1 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 1),
    #theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops),
    eps_e = matrix(0, pops, tmax-1),
    eps_d = matrix(0, pops, tmax-1),
    initial_N = obs_N[,1]#,
    #sigma_obs = runif(1, 0, 10)
  )

}

input_data10 <- list(obs_N = obs_N)

input_constants10 <- list(tmax = tmax, pops = pops, max_K = rep(K * 2, pops), max_N = apply(N, 1, max) * 2,
                          sigma_d2 = rep(sigma_d2, pops), theta = rep(theta, pops))

inits10 <- list(sample_inits10(), sample_inits10(), sample_inits10())

params10 <- c("K", "sigma_e2", "gamma", "mu_r1", "N")

# Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

# Model
start <- Sys.time()
mod10 <- nimbleMCMC(code = predict_N_random_fixed,
                    constants = input_constants10,
                    data = input_data10,
                    inits = inits10,
                    monitors = params10,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    setSeed = seed,
                    samplesAsCodaMCMC = TRUE)
dur10 <- Sys.time() - start


#---------------#
# RE model s ####
#---------------#

predict_N_random_unbiased_fixed <- nimbleCode({

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
    #theta[i] ~ dunif(-10, 10)

  }



  #sigma_obs ~ dunif(0, 100)

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops) {

    gamma[i] <- theta[i] * mu_r1[i] / (1 - K[i]^(-theta[i]))

  }


})

sample_inits11 <- function(){

  list(
    #r0 = rnorm(1, 1, 0.5),
    mu_r1 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 1),
    #theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops),
    eps_e = matrix(0, pops, tmax-1),
    eps_d = matrix(0, pops, tmax-1),
    initial_N = obs_N[,1]#,
    #sigma_obs = runif(1, 0, 10)
  )

}

input_data11 <- list(obs_N = obs_N)

input_constants11 <- list(tmax = tmax, pops = pops, max_K = rep(K * 2, pops), max_N = apply(N, 1, max) * 2,
                          sigma_d2 = rep(sigma_d2, pops), theta = rep(theta, pops))

params11 <- c("K", "sigma_e2", "gamma", "mu_r1", "N")

# Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

# Model
start <- Sys.time()

seed11 <- 319
set.seed(seed11)

inits11 <- list(sample_inits11(), sample_inits11(), sample_inits11())

mod11 <- nimbleMCMC(code = predict_N_random_unbiased_fixed,
                    constants = input_constants11,
                    data = input_data11,
                    inits = inits11,
                    monitors = params11,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    setSeed = seed11,
                    samplesAsCodaMCMC = TRUE)

dur11 <- Sys.time() - start

# Three separate chains
start <- Sys.time()
nchains <- 1

seed6a <- 329
set.seed(seed6a)

inits6a <- list(sample_inits6())

mod6a <- nimbleMCMC(code = predict_N_random_unbiased_fixed,
                    constants = input_constants6,
                    data = input_data6,
                    inits = inits6a,
                    monitors = params6,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    setSeed = seed6a,
                    samplesAsCodaMCMC = TRUE)

seed6b <- 459
set.seed(seed6b)

inits6b <- list(sample_inits6())

mod6b <- nimbleMCMC(code = predict_N_random_unbiased_fixed,
                    constants = input_constants6,
                    data = input_data6,
                    inits = inits6b,
                    monitors = params6,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    setSeed = seed6b,
                    samplesAsCodaMCMC = TRUE)

seed6c <- seed6b + 1
set.seed(seed6c)

inits6c <- list(sample_inits6())

mod6c <- nimbleMCMC(code = predict_N_random_unbiased_fixed,
                    constants = input_constants6,
                    data = input_data6,
                    inits = inits6c,
                    monitors = params6,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    setSeed = seed6c,
                    samplesAsCodaMCMC = TRUE)

dur6 <- Sys.time() - start

mod6 <- coda::mcmc.list(chain1 = mod6a, chain2 = mod6b, chain3 = mod6c)

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

predict_r_mult_nimble_fixed <- nimbleCode({

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
    #theta[i] ~ dunif(-10, 10)

  }


  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops){

    gamma[i] <- theta[i] * mu_r1[i] / (1 - K[i]^(-theta[i]))

  }


})

sample_inits12 <- function(){

  list(
    mu_r1 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 1),
    #theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops)
  )

}

## Sample initial values
#inits_b3 <- list(sample_inits_b3())

#-----------------------------#
# NIMBLE model and MCMC setup
#-----------------------------#

## Set data and constants
input_data12 <- list(N = obs_N, obs_r = obs_r)

input_constants12 <- list(tmax = tmax, max_K = rep(K * 2, pops), sigma_d2 = rep(sigma_d2, pops), theta = rep(theta, pops))

## Set parameters to monitor
params12 <- c("K", "sigma_e2", "mu_r1", "gamma", "pred_r")

# Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

start <- Sys.time()

seed12 <- 559
set.seed(seed12)

inits12 <- list(sample_inits12(), sample_inits12(), sample_inits12())

mod12 <- nimbleMCMC(code = predict_r_mult_nimble_fixed,
                    constants = input_constants12,
                    data = input_data12,
                    inits = inits12,
                    monitors = params12,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    setSeed = seed12,
                    samplesAsCodaMCMC = TRUE)
dur12 <- Sys.time() - start


#--------------#
# MLE model ####
#--------------#

mod13 <- purrr::map_dfc(.x = seq_len(pops),
                        .f = ~{

                          mod <- thetalogistic.est(N = obs_N[.x,],
                                                   year = 1:tmax,
                                                   sd = sigma_d2,
                                                   theta = theta,
                                                   nboot = 1500)

                          tibble::tibble(
                            !!paste0("sigma_e2", "[", .x, "]") := mod$boot.se,
                            !!paste0("K", "[", .x, "]") := mod$boot.K,
                            !!paste0("mu_r1", "[", .x, "]") := mod$boot.r1,
                            !!paste0("gamma", "[", .x, "]") := mod$boot.gamma
                          )

                        })

#---------------------#
# RE model r - hyp ####
#---------------------#

predict_N_random_fixed_hyp <- nimbleCode({

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
    #theta[i] ~ dunif(-10, 10)

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

sample_inits14 <- function(){

  list(
    #r0 = rnorm(1, 1, 0.5),
    mean_mu_r1 = rnorm(1, 1, 0.5),
    sd_mu_r1 = runif(1, 0, 1),
    epsilon_mu_r1 = rep(0, pops),
    mean_sigma_e2 = runif(1, 0, 1),
    sd_sigma_e2 = runif(1, 0, 1),
    epsilon_sigma_e2 = rep(0, pops),
    #theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops),
    eps_e = matrix(0, pops, tmax-1),
    eps_d = matrix(0, pops, tmax-1),
    initial_N = obs_N[,1]#,
    #sigma_obs = runif(1, 0, 10)
  )

}

input_data14 <- list(obs_N = obs_N)

input_constants14 <- list(tmax = tmax, pops = pops, max_K = rep(K * 2, pops), max_N = apply(N, 1, max) * 2,
                          sigma_d2 = rep(sigma_d2, pops), theta = rep(theta, pops))

params14 <- c("K", "mean_sigma_e2", "sd_sigma_e2", "sigma_e2", "gamma", "mean_mu_r1", "sd_mu_r1", "mu_r1", "N")

# Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

# Model
start <- Sys.time()

seed14 <- 774
set.seed(seed14)

inits14 <- list(sample_inits14(), sample_inits14(), sample_inits14())

mod14 <- nimbleMCMC(code = predict_N_random_fixed_hyp,
                    constants = input_constants14,
                    data = input_data14,
                    inits = inits14,
                    monitors = params14,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    setSeed = seed14,
                    samplesAsCodaMCMC = TRUE)
dur14 <- Sys.time() - start


#---------------------#
# RE model s - hyp ####
#---------------------#

predict_N_random_unbiased_fixed_hyp <- nimbleCode({

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
    #theta[i] ~ dunif(-10, 10)

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

sample_inits15 <- function(){

  list(
    #r0 = rnorm(1, 1, 0.5),
    mean_mu_r1 = rnorm(1, 1, 0.5),
    sd_mu_r1 = runif(1, 0, 1),
    epsilon_mu_r1 = rep(0, pops),
    mean_sigma_e2 = runif(1, 0, 1),
    sd_sigma_e2 = runif(1, 0, 1),
    epsilon_sigma_e2 = rep(0, pops),
    #theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops),
    eps_e = matrix(0, pops, tmax-1),
    eps_d = matrix(0, pops, tmax-1),
    initial_N = obs_N[,1]#,
    #sigma_obs = runif(1, 0, 10)
  )

}

input_data15 <- list(obs_N = obs_N)

input_constants15 <- list(tmax = tmax, pops = pops, max_K = rep(K * 2, pops), max_N = apply(N, 1, max) * 2,
                          sigma_d2 = rep(sigma_d2, pops), theta = rep(theta, pops))

params15 <- c("K", "mean_sigma_e2", "sd_sigma_e2", "sigma_e2", "gamma", "mean_mu_r1", "sd_mu_r1", "mu_r1", "N")

# Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

# Model
start <- Sys.time()

seed15 <- 931
set.seed(seed15)

inits15 <- list(sample_inits15(), sample_inits15(), sample_inits15())

mod15 <- nimbleMCMC(code = predict_N_random_unbiased_fixed_hyp,
                    constants = input_constants15,
                    data = input_data15,
                    inits = inits15,
                    monitors = params15,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    setSeed = seed15,
                    samplesAsCodaMCMC = TRUE)

dur15 <- Sys.time() - start

#-----------------#
# Plot outputs ####
#-----------------#

plot_theta_mult(filename = "Approx-RE-MLE-fixed-theta",
                model_list = list("Approx" = mod12, "MLE" = mod13, "RE-r" = mod10, "RE-s" = mod11,
                                  "RE-r-hyp" = mod14, "RE-s-hyp" = mod15),
                pars = c("sigma_e2", "K", "mu_r1", "gamma"), true = c(sigma_e2, K, mu_r1, gamma))
