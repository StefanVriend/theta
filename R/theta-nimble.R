library(tidyverse)
library(nimble)
library(viridis)
library(extrafont)
extrafont::loadfonts(device = "win")

nimbleOptions(disallow_multivariate_argument_expressions = FALSE)

mySeed <- 74
set.seed(mySeed)

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

#--------------------#
# Data simulation ####
#--------------------#

## NOTE: Simulate data with population sizes both very low and fluctuating around K

## Set parameter values
tmax <- 100
mu_r1 <- 0.5 # up to 1
sigma_e2 <- 0.01
sigma_e <- sqrt(sigma_e2)
sigma_d2 <- 0.2
sigma_d <- sqrt(sigma_d2)
K <- 100
theta <- 1.2

## Simulate random effects

## Prepare population vector and set initial population size
N <- rep(NA, tmax)
N[1] <- 10

epsilon_r1 <- rep(NA, tmax-1)

## Use nimble function to predict population size over time
for(t in 1:(tmax-1)){

  epsilon_r1[t] <- rnorm(1, 0, sqrt(sigma_e^2 + ((sigma_d^2) / N[t])))

  N[t+1] <- predict_N(N_current = N[t], mu_r1 = mu_r1, epsilon_r1 = epsilon_r1[t],
                      K = K, theta = theta, #beta.r1 = 0, EnvCov = 0,
                      sigma_e = sigma_e,
                      sigma_d = sigma_d)

}

gamma <- mu_r1*theta/(1-K^(-theta))

## log(N) and observation error
sigma_Y <- 0.1
log_Y <- rnorm(length(N), log(N) - 0.5 * sigma_Y^2, sigma_Y)

#---------------------#
# Plot simulated data
#---------------------#

#png(here::here("inst", "images", "sim-data.png"), width = 7, height = 4, units = "in", res = 350)
par(mfrow = c(1, 2))
plot(N, type = "l", xlab = "Time")
#lines(exp(log_Y), type = "l", col = "red")
title(main = bquote(K~"="~.(K)*", "*bar(r)~"="~.(mu_r1)*", "*sigma[e]^2~"="~.(sigma_e2)*", "*sigma[d]^2~"="~.(sigma_d2)))
legend("bottomright", col = c("black", "red"), lty = 1, legend = c("True N", "Obs N"), bty = "n", cex = 0.75)

# Plot density dependence
nn <- seq(0, max(N), length.out = 1000)
rr <- mu_r1 - (0.5 * sigma_e2) - mu_r1 * (((nn^theta)-1) / ((K^theta)-1)) # no sigma_d, no epsilon_r1

plot(N[-tmax], diff(log(N)), xlab = "N", ylab = "r")
lines(nn, rr)
title(main = bquote(theta~"="~.(theta)*", "*gamma~"="~.(round(gamma, 2))))
#dev.off()


#-----------------------------#
# MODEL A ####
# Model for population growth
# Predicting N
#-----------------------------#


#--------------------------#
# ... NIMBLE model code ####
#--------------------------#

predict_N_nimble <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  ## Initial population size
  N[1] <- initial_N

  ## Population growth
  for(t in 1:(tmax-1)){

    N[t+1] <- predict_N(N_current = N[t],
                        mu_r1 = mu_r1, epsilon_r1 = epsilon_r1[t],
                        theta = theta, K = K,
                        #env_cov = 0, beta_r1 = 0,
                        sigma_e = sigma_e, sigma_d = sigma_d)


  }

  # NOTE: If env_cov and beta_r1 are set to 0, have to pass that 0 as a vector (look up how to do again)
  # NOTE: When fitting a model without REs, we can simply set epsilon_r1 = 0 in the function call

  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  ## Binomial observation model
  #for(t in 1:tmax){
  #  N_obs[t] ~ dbin(p_obs[t], round(N[t]))
  # NOTE: The binomial observation model does not work yet.
  #       At the moment, it's "too rigid" this way because the predict_N
  #       function already does not allow for stochastic outcomes (so whenever
  #       predict_N returns an N that does not satify the condition round(N) = N_obs
  #       then there will be an infinite log probability.

  ## Alternative: Poisson observation model
  for(t in 1:(tmax-1)){

    N_obs[t] ~ dpois(p_obs[t] * N[t])

  }

  ## Alternative: Normal observation model
  #for(t in 1:(tmax-1)){
  #
  #  N_obs[t] ~ dnorm(p_obs[t] * N[t], sd = sigma_obs)
  #
  #}

  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  ## Random effects
  for(t in 1:tmax){

    epsilon_r1[t] ~ dnorm(0, var = var_r1[t])
    var_r1[t] <- sigma_e^2 + ((sigma_d^2) / N[t])

  }

  ## Priors
  mu_r1 ~ dlnorm(meanlog = 0, sdlog = 0.5)
  sigma_e ~ dunif(0, 1)
  K ~ dunif(1, max_K)
  #theta ~ dunif(-2, 10)
  theta ~ dlogis(location = 3, scale = 1)
  initial_N ~ dunif(1, max_N1)

  # for(i in 1:NoCovariates){
  #
  #   beta_r1[i] ~ dunif(-5, 5)
  #
  # }

  #sigma_obs ~ dunif(0, 50)

  gamma <- mu_r1*theta/(1-K^(-theta))


})

#-------------------#
# Initial values
#-------------------#

## Function to sample initial values
sample_inits_a <- function(){

  list(
    mu_r1 = rnorm(1, 1, 0.5),
    sigma_e  = runif(1, 0, 1),
    theta = rnorm(1, 2, 0.5),
    epsilon_r1 = runif(tmax, -0.5, 0.5),
    K = K,
    initial_N = N[1]#,
    #sigma_obs = runif(1, 0, 5)
  )

}

## Sample initial values
#inits <- list(sample_inits())
inits_a <- list(sample_inits_a(), sample_inits_a(), sample_inits_a())

#-----------------------------#
# NIMBLE model and MCMC setup
#-----------------------------#

## Set data and constants
input_data_a <- list(N_obs = round(N), p_obs = rep(1, tmax))

input_constants_a <- list(tmax = tmax, max_N1 = N[1]*2, max_K = K*2, sigma_d = sigma_d)

## Set parameters to monitor
params_a <- c("K", "theta", "sigma_e", "mu_r1", "gamma", "epsilon_r1", "N")

## Set MCMC parameters
niter <- 35000
nburnin <- 25000
nthin <- 20
nchains <- 3

#------------#
# Run model
#------------#

start <- Sys.time()
mod_a <- nimbleMCMC(code = predict_N_nimble,
                    constants = input_constants_a,
                    data = input_data_a,
                    inits = inits_a,
                    monitors = params_a,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    #setSeed = mySeed,
                    samplesAsCodaMCMC = TRUE)
dur_a <- Sys.time() - start

#--------------#
# Plot results
#--------------#

# plot(mod_a, ask = T)
#
# coda::gelman.diag(mod_a)
#
# l <- purrr::map2(.x = c("K", "mu_r1", "sigma_e", "theta"),
#                  .y = c(K, mu_r1, sigma_e, theta),
#            .f = ~{
#
#              #bayesplot::mcmc_combo(mod_a, pars = .x)
#              p <- bayesplot::mcmc_dens(mod_a, pars = .x) +
#                geom_vline(xintercept = .y, color = "#03396c")
#              q <- bayesplot::mcmc_trace(mod_a, pars = .x)
#
#              cowplot::plot_grid(p, q, align = "h", nrow = 1, ncol = 2)
#
#            })
#
# #ggsave("test.pdf", gridExtra::arrangeGrob(grobs = l))
#
# # Create multi-page plots
# pdf("est.pdf", width = 6, height = 4)
# invisible(lapply(l[1:4], print))
# dev.off()
#
# ggsave("N.pdf", gridExtra::arrangeGrob(grobs = l[5:104]),
#        width = 45, height = 25, dpi = 450)

#-----------------------------------------------------------#
# Extra: Model assembly step (for interactive exploration)
#-----------------------------------------------------------#

nimble_model <- nimbleModel(code = predict_N_nimble,
                            constants = input_constants,
                            data = input_data,
                            inits = inits[[1]])

nimble_model$check()
nimble_model$expandNodeNames("N")
nimble_model$getDependencies("N")
nimble_model$initializeInfo()

#----------------------------#
# ... Comparison with MLE ####
#----------------------------#

# Maximum likelihood estimation procedure for theta-logistic model
mle_a <- thetalogistic.est(N, year = seq_len(tmax), sd = sigma_d2, nboot = nrow(as.matrix(mod_a)))
thetaplot(mle_a)

# Posterior distributions from NIMBLE model
mod_a_est <- tibble::tibble(
  x = rep(c("K", "mu_r1", "sigma_e", "gamma"),
          each = nrow(as.matrix(mod_a))),
  val = c(as.matrix(mod_a)[, "K"], as.matrix(mod_a)[, "mu_r1"], as.matrix(mod_a)[, "sigma_e"], as.matrix(mod_a)[, "gamma"]),
  type = "nimble",
  chain = rep(rep(seq_len(nchains), each = (niter-nburnin)/nthin, 4))
)

# Bootstraps from MLE model
mle_a_est <- tibble::tibble(
  x = c(rep(c("K", "mu_r1", "sigma_e", "gamma"), each = nrow(as.matrix(mod_a)))),
  val = c(mle_a$boot.K, mle_a$boot.r1, sqrt(mle_a$boot.se), mle_a$boot.gamma),
  type = "mle"
)

# True values for parameters
true_est <- tibble::tibble(
  x = c("K", "mu_r1", "sigma_e", "gamma"),
  val = c(K, mu_r1, sigma_e, gamma),
  type = "true"
)

# Plot
bind_rows(mod_a_est, mle_a_est) %>%
  ggplot(aes(x = val, color = type)) +
  geom_segment(data = true_est, mapping = aes(x = val, xend = val,
                                              y = -Inf, yend = Inf),
               size = 0.75) +
  geom_density(data = mod_a_est, mapping = aes(x = val, group = chain),
               color = colorspace::lighten("#9ed1cc"), size = 0.5) +
  geom_density(size = 0.75) +
  facet_wrap(~x, scales = "free") +
  theme_classic(base_family = "Roboto") +
  labs(x = "", y = "Density") +
  scale_color_manual(values = c("#006D77", "#53b5ab", "#E29578"),
                     name = "") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks())

# ggsave(here::here("inst", "images", "parameter-comparison.pdf"), device = cairo_pdf, width = 6, height = 4)
# ggsave(here::here("inst", "images", "parameter-comparison.png"), width = 6, height = 4, dpi = 350)



#-----------------------------#
# MODEL B ####
# Model for population growth
# Predicting r
# No observation error
#-----------------------------#

#--------------------------#
# ... NIMBLE model code ####
#--------------------------#

predict_r_nimble <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  for (t in 1:(tmax-1)){

    s[t] <- mu_r1 - (0.5 * sigma_e2) - (0.5 * sigma_d2 / N[t])

    pred_r[t] <- s[t] - mu_r1 * (((N[t]^theta)-1) / ((K^theta)-1))

  }

  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  for (t in 1:(tmax-1)){

    obs_r[t] ~ dnorm(pred_r[t], var = var_r1[t])

    var_r1[t]  <- sigma_e2 + ((sigma_d2) / N[t])

  }


  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  mu_r1 ~ dunif(-5, 5)
  #mu_r1 ~ dlnorm(meanlog = 0, sdlog = 0.5)
  sigma_e2 ~ dunif(0, 10)
  K ~ dunif(1, max_K)
  theta ~ dunif(-10, 10)
  #theta ~ dlogis(location = 3, scale = 1)


  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  gamma <- mu_r1*theta/(1-K^(-theta))

})

#-------------------#
# Initial values
#-------------------#

## Function to sample initial values
sample_inits_b <- function(){

  list(
    mu_r1 = rnorm(1, 1, 0.5),
    sigma_e2  = runif(1, 0, 10),
    theta = rnorm(1, 2, 0.5),
    K = K
  )

}

## Sample initial values
#inits <- list(sample_inits())
inits_b <- list(sample_inits_b(), sample_inits_b(), sample_inits_b())

#-----------------------------#
# NIMBLE model and MCMC setup
#-----------------------------#

## Set data and constants
input_data_b <- list(N = round(N), obs_r = diff(log(N)))

input_constants_b <- list(tmax = tmax, max_K = K*2, sigma_d2 = sigma_d2)

## Set parameters to monitor
params_b <- c("K", "theta", "sigma_e2", "mu_r1", "gamma")

## Set MCMC parameters
niter <- 35000
nburnin <- 25000
nthin <- 20
nchains <- 3

#------------#
# Run model
#------------#
start <- Sys.time()
mod_b <- nimbleMCMC(code = predict_r_nimble,
                    constants = input_constants_b,
                    data = input_data_b,
                    inits = inits_b,
                    monitors = params_b,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    #setSeed = mySeed,
                    samplesAsCodaMCMC = TRUE)
dur_b <- Sys.time() - start

coda::gelman.diag(mod_b)

bayesplot::mcmc_combo(as.list(mod_b), pars = "theta")
bayesplot::mcmc_acf(as.list(mod_b))

# nim_b <- nimbleModel(code = predict_r_nimble, data = input_data_b,
#                      constants = input_constants_b, inits = inits_b[[1]])
#
# nim_conf <- configureMCMC(nim_b, monitors = params_b)

#----------------------------#
# ... Comparison with MLE ####
#----------------------------#

# Maximum likelihood estimation procedure for theta-logistic model
mle_b <- thetalogistic.est(N, year = seq_len(tmax), sd = sigma_d2, nboot = nrow(as.matrix(mod_b)))
thetaplot(mle_b)

# Posterior distributions from NIMBLE model
mod_b_est <- tibble::tibble(
  x = rep(c("K", "mu_r1", "sigma_e2", "gamma"),
          each = nrow(as.matrix(mod_b))),
  val = c(as.matrix(mod_b)[, "K"], as.matrix(mod_b)[, "mu_r1"], as.matrix(mod_b)[, "sigma_e2"], as.matrix(mod_b)[, "gamma"]),
  type = "nimble",
  chain = rep(rep(seq_len(nchains), each = (niter-nburnin)/nthin, 4))
)

# Bootstraps from MLE model
mle_b_est <- tibble::tibble(
  x = c(rep(c("K", "mu_r1", "sigma_e2", "gamma"), each = nrow(as.matrix(mod_b)))),
  val = c(mle_b$boot.K, mle_b$boot.r1, mle_b$boot.se, mle_b$boot.gamma),
  type = "mle"
)

# True values for parameters
true_est <- tibble::tibble(
  x = c("K", "mu_r1", "sigma_e2", "gamma"),
  val = c(K, mu_r1, sigma_e2, gamma),
  type = "true"
)

# Plot
bind_rows(mod_b_est, mle_b_est) %>%
  ggplot(aes(x = val, color = type)) +
  geom_segment(data = true_est, mapping = aes(x = val, xend = val,
                                              y = -Inf, yend = Inf),
               size = 0.75) +
  geom_density(data = mod_b_est, mapping = aes(x = val, group = chain),
               color = colorspace::lighten("#9ed1cc"), size = 0.5) +
  geom_density(size = 0.75) +
  facet_wrap(~x, scales = "free") +
  theme_classic(base_family = "Roboto") +
  labs(x = "", y = "Density") +
  scale_color_manual(values = c("#006D77", "#53b5ab", "#E29578"),
                     name = "") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks())

# ggsave(here::here("inst", "images", "parameter-comparison.pdf"), device = cairo_pdf, width = 6, height = 4)
# ggsave(here::here("inst", "images", "parameter-comparison.png"), width = 6, height = 4, dpi = 350)


#-----------------------------#
# MODEL C ####
# Model for population growth
# Predicting r
# Including observation error
#-----------------------------#

#--------------------------#
# ... NIMBLE model code ####
#--------------------------#

predict_r_nimble_obs <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  # Initialize N
  #logN[1] <- initial_N
  N[1] <- initial_N

  for (t in 1:(tmax-1)){

    s[t] <- mu_r1 - (0.5 * sigma_e2) - (0.5 * sigma_d2 / N[t])

    pred_r[t] <- s[t] - mu_r1 * (((N[t]^theta)-1) / ((K^theta)-1))

    var_r1[t] <- sigma_e2 + ((sigma_d2) / N[t])

    log(N[t+1]) ~ dnorm(log(N[t]) + pred_r[t], sd = sqrt(var_r1[t]))

  }

  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  for (t in 1:(tmax-1)){

    #logN_obs[t] ~ dnorm(p_obs[t] * (log(N[t]) - 0.5 * sigma_obs^2), sd = sigma_obs)
    #N_obs[t] ~ dnorm(p_obs[t] * N[t], sd = sigma_obs)

    # NB: lognormal observation error can be realistically expected in many ecological settings
    # see (e.g. Dennis et al. 2006 Ecological Monographs)

    N_obs[t] ~ dpois(N[t] * p_obs[t])
  }

  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  mu_r1 ~ dunif(-5, 5)
  #mu_r1 ~ dlnorm(meanlog = 0, sdlog = 0.5)
  sigma_e2 ~ dunif(0, 10)
  K ~ dunif(1, max_K)
  theta ~ dunif(-10, 10)
  #theta ~ dlogis(location = 3, scale = 1)

  #sigma_obs ~ dunif(0, 1)
  initial_N ~ dunif(1, max_N1)

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  gamma <- mu_r1*theta/(1-K^(-theta))

})

#-------------------#
# Initial values
#-------------------#

## Function to sample initial values
sample_inits_c <- function(){

  list(
    mu_r1 = runif(1, 0.1, 0.5),
    sigma_e2 = runif(1, 0, 0.1),
    theta = runif(1, 1.2, 2.5),
    K = K#,
    #sigma_obs = runif(1, 0, 0.1)
  )

}

## Sample initial values
inits_c <- list(sample_inits_c(), sample_inits_c(), sample_inits_c())

#-----------------------------#
# NIMBLE model and MCMC setup
#-----------------------------#

## Set data and constants
#input_data_c <- list(log_Y = log_Y)
#input_data_c <- list(logN_obs = log_Y, p_obs = rep(1, tmax))
input_data_c <- list(N_obs = round(exp(log_Y)), p_obs = rep(1, tmax))

input_constants_c <- list(tmax = tmax, max_K = K*2,
                          max_N1 = max(N),
                          sigma_d2 = sigma_d2)

## Set parameters to monitor
#params_c <- c("K", "theta", "sigma_e2", "mu_r1", "gamma", "sigma_obs")
params_c <- c("K", "theta", "sigma_e2", "mu_r1", "gamma")

## Set MCMC parameters
niter <- 95000
nburnin <- 75000
nthin <- 50
nchains <- 3

#------------#
# Run model
#------------#

mod_c <- nimbleMCMC(code = predict_r_nimble_obs,
                    constants = input_constants_c,
                    data = input_data_c,
                    inits = inits_c,
                    monitors = params_c,
                    niter = niter,
                    nburnin = nburnin,
                    thin = nthin,
                    nchains = nchains,
                    #setSeed = mySeed,
                    samplesAsCodaMCMC = TRUE)

coda::gelman.diag(mod_c)

bayesplot::mcmc_combo(as.list(mod_c), pars = "theta")
bayesplot::mcmc_acf(as.list(mod_c))

# nim_c <- nimbleModel(predict_r_nimble_obs, constants = input_constants_c,
#                      data = input_data_c, inits = inits_c[[1]]) # works
#
# nim_comp <- compileNimble(nim_c) # works
#
# nim_conf <- configureMCMC(nim_c, monitors = params_c, print = TRUE)
#
# nim_build <- buildMCMC(nim_conf)


#----------------------------#
# ... Comparison with MLE ####
#----------------------------#

# Maximum likelihood estimation procedure for theta-logistic model
mle_c <- thetalogistic.est(N, year = seq_len(tmax), sd = sigma_d2, nboot = nrow(as.matrix(mod_c)))
thetaplot(mle_c)

# Posterior distributions from NIMBLE model
# mod_c_est <- tibble::tibble(
#   x = rep(c("K", "mu_r1", "sigma_e2", "gamma", "sigma_obs"),
#           each = nrow(as.matrix(mod_c))),
#   val = c(as.matrix(mod_c)[, "K"], as.matrix(mod_c)[, "mu_r1"], as.matrix(mod_c)[, "sigma_e2"], as.matrix(mod_c)[, "gamma"], as.matrix(mod_c)[, "sigma_obs"]),
#   type = "nimble",
#   chain = rep(rep(seq_len(nchains), each = (niter-nburnin)/nthin, 5))
# )

mod_c_est <- tibble::tibble(
  x = rep(c("K", "mu_r1", "sigma_e2", "gamma"),
          each = nrow(as.matrix(mod_c))),
  val = c(as.matrix(mod_c)[, "K"], as.matrix(mod_c)[, "mu_r1"], as.matrix(mod_c)[, "sigma_e2"], as.matrix(mod_c)[, "gamma"]),
  type = "nimble",
  chain = rep(rep(seq_len(nchains), each = (niter-nburnin)/nthin, 4))
)

# Bootstraps from MLE model
mle_c_est <- tibble::tibble(
  x = c(rep(c("K", "mu_r1", "sigma_e2", "gamma"), each = nrow(as.matrix(mod_c)))),
  val = c(mle_c$boot.K, mle_c$boot.r1, mle_c$boot.se, mle_c$boot.gamma),
  type = "mle"
)

# True values for parameters
# true_est <- tibble::tibble(
#   x = c("K", "mu_r1", "sigma_e2", "gamma", "sigma_obs"),
#   val = c(K, mu_r1, sigma_e2, gamma, sigma_Y),
#   type = "true"
# )

true_est <- tibble::tibble(
  x = c("K", "mu_r1", "sigma_e2", "gamma"),
  val = c(K, mu_r1, sigma_e2, gamma),
  type = "true"
)

# Plot
bind_rows(mod_c_est, mle_c_est) %>%
  ggplot(aes(x = val, color = type)) +
  geom_segment(data = true_est, mapping = aes(x = val, xend = val,
                                              y = -Inf, yend = Inf),
               size = 0.75) +
  geom_density(data = mod_c_est, mapping = aes(x = val, group = chain),
               color = colorspace::lighten("#9ed1cc"), size = 0.5) +
  geom_density(size = 0.75) +
  facet_wrap(~x, scales = "free") +
  theme_classic(base_family = "Roboto") +
  labs(x = "", y = "Density") +
  scale_color_manual(values = c("#006D77", "#53b5ab", "#E29578"),
                     name = "") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks())

# ggsave(here::here("inst", "images", "parameter-comparison.pdf"), device = cairo_pdf, width = 6, height = 4)
# ggsave(here::here("inst", "images", "parameter-comparison.png"), width = 6, height = 4, dpi = 350)

#-----------------------------#
# Multiple data simulation ####
#-----------------------------#

## NOTE: Simulate data with population sizes both very low and fluctuating around K

## Number of populations
pops <- 25

## Set parameter values
tmax <- 100
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

#-----------------------------#
# MODEL B2 ####
# Model for population growth
# Predicting r
# No observation error
# Multiple populations
#-----------------------------#

#--------------------------#
# ... NIMBLE model code ####
#--------------------------#

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


## Function to sample initial values
sample_inits_b2 <- function(){

  list(
    mu_r1 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 1),
    theta = rnorm(pops, 2, 0.5),
    K = K
  )

}

## Sample initial values
#inits_b2 <- list(sample_inits_b2())
inits_b2 <- list(sample_inits_b2(), sample_inits_b2(), sample_inits_b2())

#-----------------------------#
# NIMBLE model and MCMC setup
#-----------------------------#

## Set data and constants
input_data_b2 <- list(N = round(N), obs_r = obs_r)

input_constants_b2 <- list(tmax = tmax, max_K = K*2, sigma_d2 = rep(sigma_d2, pops))

## Set parameters to monitor
params_b2 <- c("K", "theta", "sigma_e2", "mu_r1", "gamma")

## Set MCMC parameters
niter <- 300000
nburnin <- 250000
nthin <- 100
nchains <- 3

#------------#
# Run model
#------------#

start <- Sys.time()
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
dur_b2 <- Sys.time() - start

coda::gelman.diag(mod_b2)

#----------------------------#
# ... Comparison with MLE ####
#----------------------------#

mle_b2 <- purrr::map2(.x = asplit(N, 1),
                      .y = 1:pops,
                      .f = ~{

                        cat("\nMLE theta-logistic estimation (", .y, "/", pops, ")\n", sep = "")
                        thetalogistic.est(.x, year = seq_len(tmax), sd = sigma_d2, nboot = nrow(as.matrix(mod_b2)))

                      })

# True values
true_est <- tibble::tibble(
  x = c("K", "mu_r1", "sigma_e2", "gamma", "sigma_obs"),
  val = c(K, mu_r1, sigma_e2, gamma, sigma_Y),
  type = "true"
)

# Recoding table
cd <- tibble::tibble(
  Par = c("K", "mu_r1", "sigma_e2", "gamma"),
  Arg = c("boot.K", "boot.r1", "boot.se", "boot.gamma"),
  Nr = c(16, 17, 19, 18)
)

# Plot NIMBLE outputs
plot_mult_nim_est <- function(par, nim, mle) {

  xmin <- min(c(min(as.matrix(nim)[,grepl(par, colnames(as.matrix(nim)))]), min(unlist(map(mle, cd[cd$Par == par,]$Nr)))))*0.95
  xmax <- max(c(max(as.matrix(nim)[,grepl(par, colnames(as.matrix(nim)))]), max(unlist(map(mle, cd[cd$Par == par,]$Nr)))))*1.05

  data <- purrr::map_dfr(.x = 1:pops,
                 .f =  ~{

                   tibble::tibble(
                     x = par,
                     val = as.matrix(nim)[, paste0(par, "[", .x, "]")],
                     repl = as.character(.x),
                     chain = rep(seq_len(nchains), each = (niter-nburnin)/nthin)
                   )

                 })

  ggplot(data, aes(x = val)) +
    geom_vline(xintercept = true_est[true_est$x == par, ]$val, linetype = "dashed") +
    geom_line(mapping = aes(color = repl, group = repl), stat = "density", alpha = 0.5) +
    #geom_density(mapping = aes(color = repl, group = repl), alpha = 0.5) +
    theme_classic(base_family = "Roboto") +
    labs(x = par, y = "Density") +
    coord_cartesian(xlim = c(xmin, xmax)) +
    #scale_color_manual(values = cl) +
    #ggthemes::scale_color_tableau(palette = "Tableau 20") +
    scale_color_hue() +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    theme(legend.position = "none")

}

multi_nim_b <- cowplot::plot_grid(plotlist = purrr::map(c("K", "mu_r1", "sigma_e2", "gamma"), ~{plot_mult_nim_est(.x, mod_b2, mle_b2)}), align = "v")

# Plot MLE outputs
plot_mult_mle_est <- function(par, nim, mle) {

  xmin <- min(c(min(as.matrix(nim)[,grepl(par, colnames(as.matrix(nim)))]), min(unlist(map(mle, cd[cd$Par == par,]$Nr)))))*0.95
  xmax <- max(c(max(as.matrix(nim)[,grepl(par, colnames(as.matrix(nim)))]), max(unlist(map(mle, cd[cd$Par == par,]$Nr)))))*1.05

  data <- purrr::map2_dfr(.x = mle,
                         .y = 1:length(mle),
                         .f =  ~{

                           tibble::tibble(
                             x = par,
                             val = .x[[cd[cd$Par == par,]$Arg]],
                             repl = as.character(.y)
                           )

                         })

  ggplot(data, aes(x = val)) +
    geom_vline(xintercept = true_est[true_est$x == par, ]$val) +
    geom_line(mapping = aes(color = repl, group = repl), stat = "density", alpha = 0.5) +
    #geom_density(mapping = aes(color = repl, group = repl), ) +
    theme_classic(base_family = "Roboto") +
    labs(x = par, y = "Density") +
    coord_cartesian(xlim = c(xmin, xmax)) +
    #scale_color_manual(values = dl) +
    #ggthemes::scale_color_tableau(palette = "Tableau 20") +
    scale_color_hue() +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    theme(legend.position = "none")

}

multi_mle_b <- cowplot::plot_grid(plotlist = purrr::map(c("K", "mu_r1", "sigma_e2", "gamma"), ~{plot_mult_mle_est(.x, mod_b2, mle_b2)}), align = "v")

# Combine in one plot
cowplot::plot_grid(multi_nim_b, multi_mle_b, nrow = 1, labels = list("NIM", "MLE"))
ggsave(here::here("inst", "images", "multi-outputs.png"), width = 12, height = 7, units = "in")


#-----------------------------#
# MODEL C2 ####
# Model for population growth
# Predicting r
# Including observation error
# Multiple populations
#-----------------------------#

#--------------------------#
# ... NIMBLE model code ####
#--------------------------#

predict_r_mult_nimble_obs <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  # Initialize N

  for (i in 1:pops){

    N[i, 1] <- initial_N[i]

    for (t in 1:(tmax-1)){

      s[i, t] <- mu_r1[i] - (0.5 * sigma_e2[i]) - (0.5 * sigma_d2[i] / N[i, t])

      pred_r[i, t] <- s[i, t] - mu_r1[i] * (((N[i, t] ^ theta[i]) - 1) / ((K[i] ^ theta[i]) - 1))

      var_r1[i, t] <- sigma_e2[i] + ((sigma_d2[i]) / N[i, t])

      log(N[i, t+1]) ~ dnorm(log(N[i, t]) + pred_r[i, t], sd = sqrt(var_r1[i, t]))

    }

  }

  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  for (i in 1:pops){

    for (t in 1:(tmax-1)){

      #logN_obs[t] ~ dnorm(p_obs[t] * (log(N[t]) - 0.5 * sigma_obs^2), sd = sigma_obs)
      #N_obs[i, t] ~ dnorm(p_obs[i, t] * N[i, t], sd = sigma_obs[i])

      # NB: lognormal observation error can be realistically expected in many ecological settings
      # see (e.g. Dennis et al. 2006 Ecological Monographs)

      N_obs[i, t] ~ dpois(N[i, t] * p_obs[i, t])
    }

  }

  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  for (i in 1:pops){

    mu_r1[i] ~ dunif(-5, 5)
    #mu_r1 ~ dlnorm(meanlog = 0, sdlog = 0.5)
    sigma_e2[i] ~ dunif(0, 10)
    K[i] ~ dunif(1, max_K)
    theta[i] ~ dunif(-10, 10)
    #theta ~ dlogis(location = 3, scale = 1)

    #sigma_obs[i] ~ dunif(0, 10)
    initial_N[i] ~ dunif(1, max_N1)

  }

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for (i in 1:pops){

    gamma[i] <- mu_r1[i] * theta[i] / (1 - K[i] ^ (-theta[i]))

  }

})

#-------------------#
# Initial values
#-------------------#

## Function to sample initial values
sample_inits_c2 <- function(){

  list(
    mu_r1 = runif(pops, 0.1, 0.5),
    sigma_e2 = runif(pops, 0, 0.1),
    theta = runif(pops, 1.2, 2.5),
    K = rep(K, pops)#,
    #sigma_obs = runif(pops, 0, 0.1)
  )

}

## Sample initial values
inits_c2 <- list(sample_inits_c2(), sample_inits_c2(), sample_inits_c2())

input_data_c2 <- list(N_obs = round(exp(log_Y)), p_obs = matrix(1, nrow = pops, ncol = tmax))

input_constants_c2 <- list(tmax = tmax, max_K = K*2,
                           max_N1 = max(N),
                           sigma_d2 = rep(sigma_d2, pops))

## Set parameters to monitor
#params_c2 <- c("K", "theta", "sigma_e2", "mu_r1", "gamma", "sigma_obs")
params_c2 <- c("K", "theta", "sigma_e2", "mu_r1", "gamma")

## Set MCMC parameters
niter <- 105000
nburnin <- 85000
nthin <- 50
nchains <- 3

#------------#
# Run model
#------------#

mod_c2 <- nimbleMCMC(code = predict_r_mult_nimble_obs,
                     constants = input_constants_c2,
                     data = input_data_c2,
                     inits = inits_c2,
                     monitors = params_c2,
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     nchains = nchains,
                     #setSeed = mySeed,
                     samplesAsCodaMCMC = TRUE)

coda::gelman.diag(mod_c2)

#----------------------------#
# ... Comparison with MLE ####
#----------------------------#

mle_c2 <- purrr::map2(.x = asplit(N, 1),
                      .y = 1:pops,
                      .f = ~{

                        cat("\nMLE theta-logistic estimation (", .y, "/", pops, ")\n", sep = "")
                        thetalogistic.est(.x, year = seq_len(tmax), sd = sigma_d2, nboot = nrow(as.matrix(mod_c2)))

                      })

true_est <- tibble::tibble(
  x = c("K", "mu_r1", "sigma_e2", "gamma"),
  val = c(K, mu_r1, sigma_e2, gamma),
  type = "true"
)

multi_nim_c <- cowplot::plot_grid(plotlist = purrr::map(c("K", "mu_r1", "sigma_e2", "gamma"), ~{plot_mult_nim_est(.x, mod_c2, mle_b2)}), align = "h", nrow = 1)
multi_mle_c <- cowplot::plot_grid(plotlist = purrr::map(c("K", "mu_r1", "sigma_e2", "gamma"), ~{plot_mult_mle_est(.x, mod_c2, mle_b2)}), align = "h", nrow = 1)
cowplot::plot_grid(multi_nim_c, multi_mle_c, nrow = 2, labels = list("NIM", "MLE"))
ggsave(here::here("inst", "images", "multi-outputs_c.png"), width = 12, height = 7, units = "in")



#-----------------------------#
# Multiple data simulation ####
#-----------------------------#

## NOTE: Simulate data with population sizes both very low and fluctuating around K

## Number of populations
pops <- 100

## Set parameter values
tmax <- 60

mean_mu_r1 <- 0.5
sd_mu_r1 <- 0.2
epsilon_mu_r1 <- rnorm(pops, mean = 0, sd = sd_mu_r1)
mu_r1 <- mean_mu_r1 + epsilon_mu_r1

mean_sigma_e2 <- 0.001
sd_sigma_e2 <- 0.5
epsilon_sigma_e2 <- rnorm(pops, mean = 0, sd = sd_sigma_e2)
sigma_e2 <- exp(log(mean_sigma_e2) + epsilon_sigma_e2)
sigma_e <- sqrt(sigma_e2)

sigma_d2 <- rep(0.4, pops)
sigma_d <- sqrt(sigma_d2)

K <- round(runif(pops, 200, 400))

theta <- rnorm(pops, 1.2, 0.5)

## Simulate random effects

## Prepare population vector and set initial population size
N <- matrix(NA, nrow = pops, ncol = tmax)
N[,1] <- 10

epsilon_r1 <- rep(NA, tmax-1)

## Use nimble function to predict population size over time
for(r in 1:nrow(N)) {
  for(t in 1:(tmax-1)){

    epsilon_r1[t] <- rnorm(1, 0, sqrt(sigma_e[r]^2 + ((sigma_d[r]^2) / N[r, t])))

    N[r, t+1] <- predict_N(N_current = N[r, t], mu_r1 = mu_r1[r], epsilon_r1 = epsilon_r1[t],
                           K = K[r], theta = theta[r], #beta.r1 = 0, EnvCov = 0,
                           sigma_e = sigma_e[r],
                           sigma_d = sigma_d[r])

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

plot(N[1,], type = "l", xlab = "Time", ylab = "N", col = adjustcolor(cl[1], alpha.f = 0.6), ylim = c(0, max(N)))

for(r in 2:nrow(N)) {

  lines(N[r,], type = "l", col = adjustcolor(cl[r], alpha.f = 0.6))

}

#title(main = bquote(K~"="~.(K)*", "*bar(r)~"="~.(mu_r1)*", "*sigma[e]^2~"="~.(sigma_e2)*", "*sigma[d]^2~"="~.(sigma_d2)))
#legend("bottomright", col = c("black", "red"), lty = 1, legend = c("True N", "Obs N"), bty = "n", cex = 0.75)

# Plot density dependence

nn <- rr <- matrix(NA, nrow = pops, ncol = 1000)

for(r in 1:nrow(N)) {

  nn[r,] <- seq(0, max(N[r,]), length.out = 1000)
  rr[r,] <- mu_r1[r] - (0.5 * sigma_e2[r]) - mu_r1[r] * (((nn[r,]^theta[r])-1) / ((K[r]^theta[r])-1))

}

plot(N[1, -tmax] / max(N[1,]), diff(log(N[1,])), xlab = "N", ylab = "r", col = adjustcolor(cl[1], alpha.f = 0.3),
     xlim = c(0, 1), ylim = c(-3, 2))

for(r in 2:nrow(nn)) {

  points(N[r, -tmax] / max(N[r,]), diff(log(N[r,])), col = adjustcolor(cl[r], alpha.f = 0.3))

}

lines(nn[1,] / max(N[1,]), rr[1,], col = adjustcolor(cl[1], alpha.f = 0.6))

for(r in 2:nrow(nn)) {

  lines(nn[r,] / max(N[r,]), rr[r,], col = adjustcolor(cl[r], alpha.f = 0.6))

}

#title(main = bquote(theta~"="~.(theta)*", "*gamma~"="~.(round(gamma, 2))))

#-----------------------------#
# MODEL B3 ####
# Model for population growth
# Predicting r
# No observation error
# Multiple populations
# Hierarchical structure
#-----------------------------#


#--------------------------#
# ... NIMBLE model code ####
#--------------------------#

predict_r_mh_nimble <- nimbleCode({

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

    mu_r1[i] <- mean_mu_r1 + epsilon_mu_r1[i]
    epsilon_mu_r1[i] ~ dnorm(0, sd = sd_mu_r1)

    log(sigma_e2[i]) <- log(mean_sigma_e2) + epsilon_sigma_e2[i]
    epsilon_sigma_e2[i] ~ dnorm(0, sd = sd_sigma_e2)

    theta[i] ~ dunif(-10, 10)
    #theta[i] <- mean_theta[species[i]]

    K[i] ~ dunif(1, max_K[i])

  }

  mean_mu_r1 ~ dunif(-10, 10)
  sd_mu_r1 ~ dunif(0, 10)

  mean_sigma_e2 ~ dunif(0, 10)
  sd_sigma_e2 ~ dunif(0, 10)

  # for(i in 1:sps) {
  #
  #   mean_theta[i] ~ dunif(-10, 10)
  #
  # }

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops) {

    gamma[i] <- mu_r1[i] * theta[i] / (1 - K[i] ^ (-theta[i]))

  }


})

## Function to sample initial values
sample_inits_b3 <- function(){

  list(
    mean_mu_r1 = rnorm(1, 0.5, 0.25),
    sd_mu_r1 = runif(1, 0, 1),
    mean_sigma_e2 = runif(1, 0, 1),
    sd_sigma_e2 = runif(1, 0, 1),
    epsilon_sigma_e2 = rep(0, pops),
    epsilon_mu_r1 = rep(0, pops),
    theta = runif(pops, -2, 5),
    K = K
  )

}

## Sample initial values
#inits_b3 <- list(sample_inits_b3())
inits_b3 <- list(sample_inits_b3(), sample_inits_b3(), sample_inits_b3())

#-----------------------------#
# NIMBLE model and MCMC setup
#-----------------------------#

## Set data and constants
input_data_b3 <- list(N = round(N), obs_r = obs_r)

input_constants_b3 <- list(tmax = tmax, max_K = K*2, sigma_d2 = sigma_d2, pops = pops)

## Set parameters to monitor
params_b3 <- c("K", "theta", "mean_sigma_e2", "sd_sigma_e2", "mean_mu_r1", "sd_mu_r1", "gamma", "mu_r1", "sigma_e2")
#params_b3 <- c("K", "theta", "mean_sigma_e2", "sd_sigma_e2", "gamma", "mu_r1", "sigma_e2")

## Set MCMC parameters
niter <- 300000
nburnin <- 250000
nthin <- 100
nchains <- 3

#------------#
# Run model
#------------#

start <- Sys.time()
mod_b3 <- nimbleMCMC(code = predict_r_mh_nimble,
                     constants = input_constants_b3,
                     data = input_data_b3,
                     inits = inits_b3,
                     monitors = params_b3,
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     nchains = nchains,
                     #setSeed = mySeed,
                     samplesAsCodaMCMC = TRUE)
dur_b3 <- Sys.time() - start

coda::gelman.diag(mod_b3)

# Compare estimates for mu_r1 between B2 & B3
# cowplot::plot_grid(
#
#   plot_mult_nim_est("mu_r1", mod_b2, mle_b2),
#   {purrr::map_dfr(.x = 1:pops,
#                   .f =  ~{
#
#                     tibble::tibble(
#                       x = "mu_r1",
#                       val = as.matrix(mod_b3)[, paste0("mu_r1", "[", .x, "]")],
#                       repl = as.character(.x),
#                       chain = rep(seq_len(nchains), each = (niter-nburnin)/nthin)
#                     )
#
#                   }) %>%
#       ggplot(aes(x = val)) +
#       geom_vline(xintercept = mean(mu_r1), linetype = "dashed") +
#       geom_line(data = {tibble::tibble(
#         x = "mu_r1",
#         val = as.matrix(mod_b3)[, "mean_mu_r1"],
#         chain = rep(seq_len(nchains), each = (niter-nburnin)/nthin)
#       )}, mapping = aes(x = val), stat = "density", color = "gray30") +
#       geom_line(mapping = aes(color = repl, group = repl), stat = "density", alpha = 0.5) +
#       #geom_density(mapping = aes(color = repl, group = repl), alpha = 0.5) +
#       theme_classic(base_family = "Roboto") +
#       labs(x = "mu_r1", y = "Density") +
#       #coord_cartesian(xlim = c(0, 3)) +
#       #scale_color_manual(values = cl) +
#       #ggthemes::scale_color_tableau(palette = "Tableau 20") +
#       scale_color_grey() +
#       scale_x_continuous(breaks = scales::pretty_breaks()) +
#       scale_y_continuous(breaks = scales::pretty_breaks()) +
#       theme(legend.position = "none")},
#   nrow = 1, align = "h"
# ) %>%
#   cowplot::save_plot(plot = ., filename = here::here("inst", "images", "mu_r1.png"), nrow = 1, ncol = 2, base_asp = 1.2)

pl <- purrr::map(.x = 1:pops,
                 .f =  ~{

                   tibble::tibble(
                     x = "mu_r1",
                     val = c(as.matrix(mod_b2)[, paste0("mu_r1", "[", .x, "]")],
                             as.matrix(mod_b3)[, paste0("mu_r1", "[", .x, "]")],
                             as.matrix(mod_b3)[, "mean_mu_r1"]),
                   pop = as.character(.x),
                   mod = rep(c("indep", "hyper", "mean"), each = nrow(as.matrix(mod_b2)))) %>%
                     ggplot(aes(x = val)) +
                     geom_vline(xintercept = mu_r1[.x], linetype = "dashed") +
                     geom_line(mapping = aes(color = mod), stat = "density", alpha = 0.5) +
                     scale_color_manual(values =  c("#00587a", "#f2af29", "#ad343e")) +
                     facet_wrap(~pop) +
                     theme_classic() +
                     guides(color = guide_legend(override.aes = list(alpha = 1)))

                 })

ggsave(filename = "plots-mu-r1.pdf",
       plot = gridExtra::marrangeGrob(pl, nrow=3, ncol=3),
       width = 15, height = 15,
       dpi = 450)

pl <- purrr::map(.x = 1:pops,
                 .f =  ~{

                   tibble::tibble(
                     x = "sigma_e2",
                     val = c(as.matrix(mod_b2)[, paste0("sigma_e2", "[", .x, "]")],
                             as.matrix(mod_b3)[, paste0("sigma_e2", "[", .x, "]")],
                             as.matrix(mod_b3)[, "mean_sigma_e2"]),
                     pop = as.character(.x),
                     mod = rep(c("indep", "hyper", "mean"), each = nrow(as.matrix(mod_b2)))) %>%
                     ggplot(aes(x = val)) +
                     geom_vline(xintercept = sigma_e2[.x], linetype = "dashed") +
                     geom_line(mapping = aes(color = mod), stat = "density", alpha = 0.5) +
                     scale_color_manual(values =  c("#00587a", "#f2af29", "#ad343e")) +
                     facet_wrap(~pop) +
                     theme_classic() +
                     guides(color = guide_legend(override.aes = list(alpha = 1)))

                 })

ggsave(filename = "plots-sigma-e2.pdf",
       plot = gridExtra::marrangeGrob(pl, nrow=3, ncol=3),
       width = 15, height = 15,
       dpi = 450)



#-----------------------------#
# Multiple data simulation ####
#-----------------------------#

## Number of populations and replicates
pops <- 10
repl <- 10

pop_vec <- rep(1:pops, each = repl)

## Set parameter values
tmax <- 60

mean_mu_r1 <- 0.5
sd_mu_r1 <- 0.2
epsilon_mu_r1 <- rnorm(pops, mean = 0, sd = sd_mu_r1)
mu_r1 <- mean_mu_r1 + epsilon_mu_r1

mean_sigma_e2 <- 0.001
sd_sigma_e2 <- 0.5
epsilon_sigma_e2 <- rnorm(pops, mean = 0, sd = sd_sigma_e2)
sigma_e2 <- exp(log(mean_sigma_e2) + epsilon_sigma_e2)
sigma_e <- sqrt(sigma_e2)

sigma_d2 <- rep(0.4, pops)
sigma_d <- sqrt(sigma_d2)

K <- round(runif(pops, 200, 400))

theta <- rnorm(pops, 1.2, 0.5)

## Simulate random effects

## Prepare population vector and set initial population size
N <- matrix(NA, nrow = pops*repl, ncol = tmax)
N[,1] <- 10

epsilon_r1 <- rep(NA, tmax-1)

## Use nimble function to predict population size over time
for(r in 1:nrow(N)) {
  for(t in 1:(tmax-1)){

    epsilon_r1[t] <- rnorm(1, 0, sqrt(sigma_e[pop_vec[r]]^2 + ((sigma_d[pop_vec[r]]^2) / N[r, t])))

    N[r, t+1] <- predict_N(N_current = N[r, t], mu_r1 = mu_r1[pop_vec[r]], epsilon_r1 = epsilon_r1[t],
                           K = K[pop_vec[r]], theta = theta[pop_vec[r]], #beta.r1 = 0, EnvCov = 0,
                           sigma_e = sigma_e[pop_vec[r]],
                           sigma_d = sigma_d[pop_vec[r]])

  }
}

gamma <- mu_r1*theta/(1-K^(-theta))

obs_r <- matrix(NA, nrow = nrow(N), ncol = tmax-1)

for(r in 1:nrow(N)) {

  obs_r[r,] <- diff(log(N[r,]))

}


#------------------------------------------------------#
# MODEL B4 ####
# Model for population growth
# Predicting r
# No observation error
# Multiple populations, each with multiple replicates
# Hierarchical structure
#------------------------------------------------------#

#--------------------------#
# ... NIMBLE model code ####
#--------------------------#

predict_r_mhr_nimble <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  for (j in 1:(repl*pops)){

    for (t in 1:(tmax-1)){

      pred_r[j, t] <- predict_r(N_current = N[j, t], mu_r1 = mu_r1[pop[j]],
                                sigma_e2 = sigma_e2[pop[j]], sigma_d2 = sigma_d2[pop[j]],
                                theta = theta[pop[j]], K = K[pop[j]])

    }

  }

  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#

  for (j in 1:(repl*pops)){

    for (t in 1:(tmax-1)){

      obs_r[j, t] ~ dnorm(pred_r[j,t], var = var_r1[j, t])

      var_r1[j, t]  <- sigma_e2[pop[j]] + ((sigma_d2[pop[j]]) / N[j, t])

    }

  }


  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  for(i in 1:pops) {

    mu_r1[i] <- mean_mu_r1 + epsilon_mu_r1[i]
    epsilon_mu_r1[i] ~ dnorm(0, sd = sd_mu_r1)

    log(sigma_e2[i]) <- log(mean_sigma_e2) + epsilon_sigma_e2[i]
    epsilon_sigma_e2[i] ~ dnorm(0, sd = sd_sigma_e2)

    theta[i] ~ dunif(-10, 10)

    K[i] ~ dunif(1, max_K[i])

  }

  mean_mu_r1 ~ dunif(-10, 10)
  sd_mu_r1 ~ dunif(0, 10)

  mean_sigma_e2 ~ dunif(0, 10)
  sd_sigma_e2 ~ dunif(0, 10)

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops) {

    gamma[i] <- mu_r1[i] * theta[i] / (1 - K[i] ^ (-theta[i]))

  }


})

## Function to sample initial values
sample_inits_b4 <- function(){

  list(
    mean_mu_r1 = rnorm(1, 0.5, 0.25),
    sd_mu_r1 = runif(1, 0, 1),
    epsilon_mu_r1 = rep(0, pops),
    mean_sigma_e2 = runif(1, 0, 1),
    sd_sigma_e2 = runif(1, 0, 1),
    epsilon_sigma_e2 = rep(0, pops),
    theta = runif(pops, -2, 5),
    K = K
  )

}

## Sample initial values
#inits_b4 <- list(sample_inits_b4())
inits_b4 <- list(sample_inits_b4(), sample_inits_b4(), sample_inits_b4())

#-----------------------------#
# NIMBLE model and MCMC setup
#-----------------------------#

## Set data and constants
input_data_b4 <- list(N = round(N), obs_r = obs_r)

input_constants_b4 <- list(tmax = tmax, max_K = K*2, sigma_d2 = sigma_d2, pop = pop_vec, pops = pops, repl = repl)

## Set parameters to monitor
params_b4 <- c("K", "theta", "mean_sigma_e2", "sd_sigma_e2", "mean_mu_r1", "sd_mu_r1", "gamma", "mu_r1", "sigma_e2")

## Set MCMC parameters
niter <- 200000
nburnin <- 150000
nthin <- 100
nchains <- 3

#------------#
# Run model
#------------#

start <- Sys.time()
mod_b4 <- nimbleMCMC(code = predict_r_mhr_nimble,
                     constants = input_constants_b4,
                     data = input_data_b4,
                     inits = inits_b4,
                     monitors = params_b4,
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     nchains = nchains,
                     #setSeed = mySeed,
                     samplesAsCodaMCMC = TRUE)
dur_b4 <- Sys.time() - start

