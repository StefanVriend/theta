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

pops <- 4
tmax <- 75

mu_r1 <- 0.5
K <- 150
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

# Plot N
cl <- viridis(pops)

plot(N[1,], type = "l", xlab = "Time", ylab = "N", col = adjustcolor(cl[1], alpha.f = 0.6))

for(r in 2:nrow(N)) {

  lines(N[r,], type = "l", col = adjustcolor(cl[r], alpha.f = 0.6))

}

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

  for(i in 1:pops){
    for(t in 1:(tmax-1)){

      #log(N[t+1]) <- log(N[t]) + r0 * (1 - (N[t] / K)^theta) + eps_e[t] + eps_d[t]
      log(N[i, t+1]) <- log(N[i, t]) + mu_r1[i] * (1 - (((N[i, t]^theta[i]) - 1) / ((K[i]^theta[i]) - 1))) + eps_e[i, t] + eps_d[i, t]

      eps_e[i, t] ~ dnorm(0, var = sigma_e2[i])
      eps_d[i, t] ~ dnorm(0, var = sigma_d2[i, t] / N[i, t])

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

sample_inits <- function(){

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

input_data <- list(obs_N = obs_N)

input_constants <- list(tmax = tmax, max_K = rep(K * 2, pops), max_N = apply(N, 1, max) * 2, sigma_d2 = matrix(sigma_d2, pops, tmax))

inits <- list(sample_inits(), sample_inits(), sample_inits())

params <- c("K", "theta", "sigma_e2", "gamma", "mu_r1", "N")

# Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

# Model
start <- Sys.time()
rod4 <- nimbleMCMC(code = predict_N_random,
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
dr <- Sys.time() - start



## Approximation model for comparison
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
  #-----

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

sample_inits2 <- function(){

  list(
    mu_r1 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 1),
    theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops)
  )

}

## Sample initial values
#inits_b2 <- list(sample_inits_b2())
inits2 <- list(sample_inits2(), sample_inits2(), sample_inits2())

#-----------------------------#
# NIMBLE model and MCMC setup
#-----------------------------#

## Set data and constants
input_data2 <- list(N = obs_N, obs_r = obs_r)

input_constants2 <- list(tmax = tmax, max_K = rep(K * 2, pops), sigma_d2 = rep(sigma_d2, pops))

## Set parameters to monitor
params2 <- c("K", "theta", "sigma_e2", "mu_r1", "gamma")

start <- Sys.time()
mod4 <- nimbleMCMC(code = predict_r_mult_nimble,
                   constants = input_constants2,
                   data = input_data2,
                   inits = inits2,
                   monitors = params2,
                   niter = niter,
                   nburnin = nburnin,
                   thin = nthin,
                   nchains = nchains,
                   #setSeed = mySeed,
                   samplesAsCodaMCMC = TRUE)
dm <- Sys.time() - start



## Compare random-effects model with approximation model
compare_mods <- function(par, true, models = c("Approx", "Rand-Poisson")) {

  data <- purrr::map2_dfr(.x = list(mod4, rod4),
                          .y = models,
                          .f = ~{

                            tibble(
                              y = as.matrix(.x)[, par],
                              group = as.character(.y),
                            )

                          }) %>%
    dplyr::mutate(group = forcats::fct_relevel(group, models))

  mean_data <- data %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(median = quantile(y, probs = 0.5, na.rm = TRUE),
                     low = quantile(y, probs = 0.025, na.rm = TRUE),
                     high = quantile(y, probs = 0.975, na.rm = TRUE),
                     .groups = "drop")

  ggplot(data = data, aes(y = group, x = y)) +
    geom_vline(xintercept = true, linetype = "dashed") +
    ggridges::geom_density_ridges(aes(fill = group), scale = 0.9, alpha = 0.5) +
    geom_segment(aes(y = group, yend = group, x = low, xend = high, color = group), data = mean_data, size = 1) +
    geom_point(aes(y = group, x = median, color = group), data = mean_data, size = 3) +
    labs(x = par, y = "") +
    scale_color_manual(values = c("#164850", "#883041", "#54611a")) +
    scale_fill_manual(values = c("#319eaf", "#cf7789", "#c5d86d")) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none")

}

# pl <- purrr::map2(.x = c("sigma_e2", "mu_r1", "K", "theta", "gamma"),
#                   .y = c(sigma_e2, mu_r1, K, theta, gamma),
#                   .f = ~{compare_mods(.x, .y)})
#
# cowplot::plot_grid(plotlist = pl)
# cowplot::save_plot(plot = ggplot2::last_plot(), filename = "rand-approx.png", nrow = 2, ncol = 3, base_asp = 1)

pl1 <- cowplot::plot_grid(plotlist = {purrr::map(paste0("sigma_e2[", 1:pops, "]"), ~compare_mods(.x, sigma_e2))}, nrow = 1, align = "h")
pl2 <- cowplot::plot_grid(plotlist = {purrr::map(paste0("mu_r1[", 1:pops, "]"), ~compare_mods(.x, mu_r1))}, nrow = 1, align = "h")
pl3 <- cowplot::plot_grid(plotlist = {purrr::map(paste0("K[", 1:pops, "]"), ~compare_mods(.x, K))}, nrow = 1, align = "h")
pl4 <- cowplot::plot_grid(plotlist = {purrr::map(paste0("theta[", 1:pops, "]"), ~compare_mods(.x, theta))}, nrow = 1, align = "h")
pl5 <- cowplot::plot_grid(plotlist = {purrr::map(paste0("gamma[", 1:pops, "]"), ~compare_mods(.x, gamma))}, nrow = 1, align = "h")

cowplot::plot_grid(pl1, pl2, pl3, pl4, pl5, nrow = 5, align = "hv")
cowplot::save_plot(plot = ggplot2::last_plot(), filename = "rand-approx4.png", nrow = 5, ncol = pops, base_asp = 1)
