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
tmax <- 40
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


# Remove a set of observations from N
# for(r in 1:nrow(N)){
#
#   idxNA <- sample(2:(tmax-1), 1)
#   N[r,idxNA] <- NA
#
# }

# Calculate obs_r
obs_r <- matrix(NA, nrow = nrow(N), ncol = tmax-1)

for(r in 1:nrow(N)) {

  obs_r[r,] <- diff(log(N[r,]))

}

# Save N & obs_r
N1 <- N
obs_r1 <- obs_r

## log(N) and observation error
# sigma_Y <- 0.1
#
# log_Y <- matrix(NA, nrow = nrow(N), ncol = ncol(N))
# for(r in 1:nrow(log_Y)) {
#
#   log_Y[r,] <- rnorm(length(N[r,]), log(N[r,]) - 0.5 * sigma_Y^2, sigma_Y)
#
# }

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
nn <- seq(0, max(N[1,], na.rm=TRUE), length.out = 1000)
rr <- mu_r1 - (0.5 * sigma_e2) - mu_r1 * (((nn^theta)-1) / ((K^theta)-1))

plot(N[1, -tmax], diff(log(N[1,])), xlab = "N", ylab = "r", col = "black")
lines(nn, rr, col = "black")

title(main = bquote(theta~"="~.(theta)*", "*gamma~"="~.(round(gamma, 2))))


#-------------------------------#
# ... NIMBLE model code - B2 ####
#-------------------------------#

predict_r_mult_na_nimble <- nimbleCode({

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

    # Prior for missing values in N
    for(t in 1:(tmax-1)){
      N[i, t] ~ dunif(0, max_K[i]*1.5)
    }
    # NOTE: A variety of distributions may be suitable here.
    #       Alternatives include:
    #       - discrete uniform (with more or less informative boundaries)
    #       - truncated normal using mean and sd of whole time-series
    #       - Poisson with expected value = mean/median of N[t-1] and N[t+1]
    #       - etc.

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


## Function to sample initial values
sample_inits_na <- function(){

  # Setting initial values for missing data points
  init.N <- matrix(NA, nrow = nrow(N), ncol = ncol(N))
  init.N[which(is.na(N))] <- mean(N, na.rm = T)

  init.obs_r <- matrix(NA, nrow = nrow(obs_r), ncol = ncol(obs_r))
  init.obs_r[which(is.na(obs_r))] <- 0

  list(
    mu_r1 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 1),
    theta = rnorm(pops, 2, 0.5),
    K = rep(K, pops),
    N = init.N,
    obs_r = init.obs_r
  )

}

## Sample initial values
inits_na <- list(sample_inits_na(), sample_inits_na(), sample_inits_na())

#-----------------------------#
# NIMBLE model and MCMC setup
#-----------------------------#

## Set data and constants
input_data_na <- list(N = N, obs_r = obs_r)

input_constants_na <- list(tmax = tmax, max_K = rep(K, pops)*2, sigma_d2 = rep(sigma_d2, pops))

## Set parameters to monitor
params_na <- c("K", "theta", "sigma_e2", "mu_r1", "gamma")

## Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

#------------#
# Run model
#------------#


## Run 1 ####
## 0 NA     #

#start <- Sys.time()
mod_na1 <- nimbleMCMC(code = predict_r_mult_na_nimble,
                     constants = input_constants_na,
                     data = input_data_na,
                     inits = inits_na,
                     monitors = params_na,
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     nchains = nchains,
                     #setSeed = mySeed,
                     samplesAsCodaMCMC = TRUE)
#dur_b2 <- Sys.time() - start

#coda::gelman.diag(mod_b2)


##      Run 2       ####
## 1 NA per population #

# Reset N & obs_r
N <- N1
obs_r <- obs_r1

# Remove a set of observations from N
for(r in 1:nrow(N)){

  idxNA <- sample(2:(tmax-1), 1)
  N[r,idxNA] <- NA

}

# Calculate obs_r
obs_r <- matrix(NA, nrow = nrow(N), ncol = tmax-1)

for(r in 1:nrow(N)) {

  obs_r[r,] <- diff(log(N[r,]))

}

# Save N & obs_r
N2 <- N
obs_r2 <- obs_r

## Sample initial values
inits_na <- list(sample_inits_na(), sample_inits_na(), sample_inits_na())

## Input data
input_data_na <- list(N = N, obs_r = obs_r)

## Run model
mod_na2 <- nimbleMCMC(code = predict_r_mult_na_nimble,
                     constants = input_constants_na,
                     data = input_data_na,
                     inits = inits_na,
                     monitors = params_na,
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     nchains = nchains,
                     #setSeed = mySeed,
                     samplesAsCodaMCMC = TRUE)


##      Run 3       ####
## 2 NA per population #

# Reset N & obs_r
N <- N1
obs_r <- obs_r1

# Remove a set of observations from N
for(r in 1:nrow(N)){

  idxNA <- sample(2:(tmax-1), 2)
  N[r,idxNA] <- NA

}

# Calculate obs_r
obs_r <- matrix(NA, nrow = nrow(N), ncol = tmax-1)

for(r in 1:nrow(N)) {

  obs_r[r,] <- diff(log(N[r,]))

}

# Save N & obs_r
N3 <- N
obs_r3 <- obs_r

## Sample initial values
inits_na <- list(sample_inits_na(), sample_inits_na(), sample_inits_na())

## Input data
input_data_na <- list(N = N, obs_r = obs_r)

## Run model
mod_na3 <- nimbleMCMC(code = predict_r_mult_na_nimble,
                      constants = input_constants_na,
                      data = input_data_na,
                      inits = inits_na,
                      monitors = params_na,
                      niter = niter,
                      nburnin = nburnin,
                      thin = nthin,
                      nchains = nchains,
                      #setSeed = mySeed,
                      samplesAsCodaMCMC = TRUE)


##      Run 4       ####
## 5 NA per population #

# Reset N & obs_r
N <- N1
obs_r <- obs_r1

# Remove a set of observations from N
for(r in 1:nrow(N)){

  idxNA <- sample(2:(tmax-1), 5)
  N[r,idxNA] <- NA

}

# Calculate obs_r
obs_r <- matrix(NA, nrow = nrow(N), ncol = tmax-1)

for(r in 1:nrow(N)) {

  obs_r[r,] <- diff(log(N[r,]))

}

# Save N & obs_r
N4 <- N
obs_r4 <- obs_r

## Sample initial values
inits_na <- list(sample_inits_na(), sample_inits_na(), sample_inits_na())

## Input data
input_data_na <- list(N = N, obs_r = obs_r)

## Run model
mod_na4 <- nimbleMCMC(code = predict_r_mult_na_nimble,
                      constants = input_constants_na,
                      data = input_data_na,
                      inits = inits_na,
                      monitors = params_na,
                      niter = niter,
                      nburnin = nburnin,
                      thin = nthin,
                      nchains = nchains,
                      #setSeed = mySeed,
                      samplesAsCodaMCMC = TRUE)


##      Run 5       ####
## 10 NA per population #

# Reset N & obs_r
N <- N1
obs_r <- obs_r1

# Remove a set of observations from N
for(r in 1:nrow(N)){

  idxNA <- sample(2:(tmax-1), 10)
  N[r,idxNA] <- NA

}

# Calculate obs_r
obs_r <- matrix(NA, nrow = nrow(N), ncol = tmax-1)

for(r in 1:nrow(N)) {

  obs_r[r,] <- diff(log(N[r,]))

}

# Save N & obs_r
N5 <- N
obs_r5 <- obs_r

## Sample initial values
inits_na <- list(sample_inits_na(), sample_inits_na(), sample_inits_na())

## Input data
input_data_na <- list(N = N, obs_r = obs_r)

## Run model
mod_na5 <- nimbleMCMC(code = predict_r_mult_na_nimble,
                      constants = input_constants_na,
                      data = input_data_na,
                      inits = inits_na,
                      monitors = params_na,
                      niter = niter,
                      nburnin = nburnin,
                      thin = nthin,
                      nchains = nchains,
                      #setSeed = mySeed,
                      samplesAsCodaMCMC = TRUE)


##      Run 6       ####
## 5 NA per population in cluster #

# Reset N & obs_r
N <- N1
obs_r <- obs_r1

# Remove a set of observations from N
for(r in 1:nrow(N)){

  idxNA <- seq(sample(2:(tmax-5), 1), length.out = 5)
  N[r,idxNA] <- NA

}

# Calculate obs_r
obs_r <- matrix(NA, nrow = nrow(N), ncol = tmax-1)

for(r in 1:nrow(N)) {

  obs_r[r,] <- diff(log(N[r,]))

}

# Save N & obs_r
N6 <- N
obs_r6 <- obs_r

## Sample initial values
inits_na <- list(sample_inits_na(), sample_inits_na(), sample_inits_na())

## Input data
input_data_na <- list(N = N, obs_r = obs_r)

## Run model
mod_na6 <- nimbleMCMC(code = predict_r_mult_na_nimble,
                      constants = input_constants_na,
                      data = input_data_na,
                      inits = inits_na,
                      monitors = params_na,
                      niter = niter,
                      nburnin = nburnin,
                      thin = nthin,
                      nchains = nchains,
                      #setSeed = mySeed,
                      samplesAsCodaMCMC = TRUE)


##      Run 7       ####
## 10 NA per population in cluster #

# Reset N & obs_r
N <- N1
obs_r <- obs_r1

# Remove a set of observations from N
for(r in 1:nrow(N)){

  idxNA <- seq(sample(2:(tmax-10), 1), length.out = 10)
  N[r,idxNA] <- NA

}

# Calculate obs_r
obs_r <- matrix(NA, nrow = nrow(N), ncol = tmax-1)

for(r in 1:nrow(N)) {

  obs_r[r,] <- diff(log(N[r,]))

}

# Save N & obs_r
N7 <- N
obs_r7 <- obs_r

## Sample initial values
inits_na <- list(sample_inits_na(), sample_inits_na(), sample_inits_na())

## Input data
input_data_na <- list(N = N, obs_r = obs_r)

## Run model
mod_na7 <- nimbleMCMC(code = predict_r_mult_na_nimble,
                      constants = input_constants_na,
                      data = input_data_na,
                      inits = inits_na,
                      monitors = params_na,
                      niter = niter,
                      nburnin = nburnin,
                      thin = nthin,
                      nchains = nchains,
                      #setSeed = mySeed,
                      samplesAsCodaMCMC = TRUE)


##      Run 8       ####
## 15 NA per population in cluster #

# Reset N & obs_r
N <- N1
obs_r <- obs_r1

# Remove a set of observations from N
for(r in 1:nrow(N)){

  idxNA <- seq(sample(2:(tmax-15), 1), length.out = 15)
  N[r,idxNA] <- NA

}

# Calculate obs_r
obs_r <- matrix(NA, nrow = nrow(N), ncol = tmax-1)

for(r in 1:nrow(N)) {

  obs_r[r,] <- diff(log(N[r,]))

}

# Save N & obs_r
N8 <- N
obs_r8 <- obs_r

## Sample initial values
inits_na <- list(sample_inits_na(), sample_inits_na(), sample_inits_na())

## Input data
input_data_na <- list(N = N, obs_r = obs_r)

## Run model
mod_na8 <- nimbleMCMC(code = predict_r_mult_na_nimble,
                      constants = input_constants_na,
                      data = input_data_na,
                      inits = inits_na,
                      monitors = params_na,
                      niter = niter,
                      nburnin = nburnin,
                      thin = nthin,
                      nchains = nchains,
                      #setSeed = mySeed,
                      samplesAsCodaMCMC = TRUE)


# Compare runs
compare_runs <- function(par, true, xlim) {

  data <- purrr::map2_dfr(.x = list(mod_na1, mod_na2, mod_na3, mod_na4,
                                    mod_na5, mod_na6, mod_na7, mod_na8),
                          .y = c(0, 1, 2, 5, 10, "CL-5", "CL-10", "CL-15"),
                          .f = ~{

                            #data <- as.matrix(.x)[, par]

                            tibble(
                              # mode = max_den(data),
                              # mean = mean(data),
                              # median = median(data),
                              # low = quantile(data, probs = 0.025),
                              # high = quantile(data, probs = 0.975),
                              y = as.matrix(.x)[, par],
                              group = as.character(.y),
                            )

                          }) %>%
    dplyr::mutate(group = forcats::fct_relevel(group, c("0", "1", "2", "5", "10", "CL-5", "CL-10"))) %>%
    dplyr::mutate(NA_type = dplyr::case_when(grepl("CL", group) ~ "Clustered",
                                             TRUE ~ "Random"))

  mean_data <- data %>%
    dplyr::group_by(NA_type, group) %>%
    dplyr::summarise(median = quantile(y, probs = 0.5),
                     low = quantile(y, probs = 0.025),
                     high = quantile(y, probs = 0.975),
                     .groups = "drop")

  ggplot(data = data, aes(y = group, x = y)) +
    geom_vline(xintercept = true, linetype = "dashed") +
    ggridges::geom_density_ridges_gradient(aes(fill = NA_type), scale = 0.9) +
    #ggridges::stat_density_ridges(quantile_lines = TRUE, scale = 0.95) +
    geom_segment(aes(y = group, yend = group, x = low, xend = high, color = NA_type), data = mean_data, size = 1) +
    geom_point(aes(y = group, x = median, color = NA_type), data = mean_data, size = 3) +
    labs(y = "# NAs", x = par) +
    coord_cartesian(xlim = xlim) +
    scale_color_manual(values = c("#123a40", "#883041")) +
    scale_fill_manual(values = c("#319eaf", "#cf7789")) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none")

}

# Sigma e2
pl1 <- purrr::map(paste0("sigma_e2[", 1:pops, "]"), ~compare_runs(.x, sigma_e2, c(0, 0.05)))

# Theta
pl2 <- purrr::map(paste0("theta[", 1:pops, "]"), ~compare_runs(.x, theta, c(-4, 4)))

# Mu_r1
pl3 <- purrr::map(paste0("mu_r1[", 1:pops, "]"), ~compare_runs(.x, mu_r1, c(0, 3)))

# K
pl4 <- purrr::map(paste0("K[", 1:pops, "]"), ~compare_runs(.x, K, c(75, 125)))


cowplot::save_plot(filename = "pars_na2.pdf", plot = cowplot::plot_grid(plotlist = c(pl1, pl2, pl3, pl4),
                                                                       nrow = 4, ncol = 5),
                   nrow = 4, ncol = 5, base_asp = 1.2)
