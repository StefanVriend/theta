library(tidyverse)
library(nimble)
library(extrafont)
extrafont::loadfonts(device = "win")

nimbleOptions(disallow_multivariate_argument_expressions = FALSE)

mySeed <- 74
set.seed(mySeed)

#-----------------------------------------------#
# Basic NIMBLE function for population growth
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

predict_N_c <- compileNimble(predict_N)

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
# MODEL A. ####
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
niter <- 10000
nburnin <- 5000
nthin <- 10
nchains <- 3

#------------#
# Run model
#------------#

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
# MODEL B. ####
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

  mu_r1 ~ dlnorm(meanlog = 0, sdlog = 0.5)
  sigma_e2 ~ dunif(0, 1)
  K ~ dunif(1, max_K)
  theta ~ dunif(-2, 10)
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
    sigma_e2  = runif(1, 0, 1),
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
# MODEL C. ####
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
    N_obs[t] ~ dnorm(p_obs[t] * N[t], sd = sigma_obs)

    # NB: lognormal observation error can be realistically expected in many ecological settings
    # see (e.g. Dennis et al. 2006 Ecological Monographs)

    #N_obs[t] ~ dpois(N[t] * p_obs[t])
  }

  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#

  mu_r1 ~ dunif(0, 1)
  #mu_r1 ~ dlnorm(meanlog = 0, sdlog = 0.5)
  sigma_e2 ~ dunif(0, 1)
  K ~ dunif(1, max_K)
  theta ~ dunif(-2, 10)
  #theta ~ dlogis(location = 3, scale = 1)

  sigma_obs ~ dunif(0, 1)
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
    K = K,
    sigma_obs = runif(1, 0, 0.1)
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
input_data_c <- list(N_obs = exp(log_Y), p_obs = rep(1, tmax))

input_constants_c <- list(tmax = tmax, max_K = K*2,
                          max_N1 = max(N),
                          sigma_d2 = sigma_d2)

## Set parameters to monitor
params_c <- c("K", "theta", "sigma_e2", "mu_r1", "gamma", "sigma_obs")
#params_c <- c("K", "theta", "sigma_e2", "mu_r1", "gamma")

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
mod_c_est <- tibble::tibble(
  x = rep(c("K", "mu_r1", "sigma_e2", "gamma", "sigma_obs"),
          each = nrow(as.matrix(mod_c))),
  val = c(as.matrix(mod_c)[, "K"], as.matrix(mod_c)[, "mu_r1"], as.matrix(mod_c)[, "sigma_e2"], as.matrix(mod_c)[, "gamma"], as.matrix(mod_c)[, "sigma_obs"]),
  type = "nimble",
  chain = rep(rep(seq_len(nchains), each = (niter-nburnin)/nthin, 5))
)

# Bootstraps from MLE model
mle_c_est <- tibble::tibble(
  x = c(rep(c("K", "mu_r1", "sigma_e2", "gamma"), each = nrow(as.matrix(mod_c)))),
  val = c(mle_c$boot.K, mle_c$boot.r1, mle_c$boot.se, mle_c$boot.gamma),
  type = "mle"
)

# True values for parameters
true_est <- tibble::tibble(
  x = c("K", "mu_r1", "sigma_e2", "gamma", "sigma_obs"),
  val = c(K, mu_r1, sigma_e2, gamma, sigma_Y),
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


#-------------------------#
# Multiple comparisons ####
#-------------------------#

pars <- expand.grid(mu_r1 = c(0.4, 0.7),
                    sigma_e = c(sqrt(0.005), sqrt(0.015), sqrt(0.00001)),
                    theta = c(0.5, 1, 2),
                    K = c(200, 500))

# Simulate data
sim_data <- purrr::pmap(.l = pars,
                        .f = ~{

                          tmax <- 150
                          sigma_d <- sqrt(0.01)

                          # Prepare population vector and set initial population size
                          N <- rep(NA, tmax)
                          N[1] <- 10

                          epsilon_r1 <- rep(NA, tmax-1)

                          # Use nimble function to predict population size over time
                          for(t in 1:(tmax-1)){

                            epsilon_r1[t] <- rnorm(1, 0, sqrt(..2^2 + ((sigma_d^2) / N[t])))

                            N[t+1] <- predict_N(N_current = N[t],
                                                mu_r1 = ..1,
                                                epsilon_r1 = epsilon_r1[t],
                                                K = ..4,
                                                theta = ..3, #beta.r1 = 0, EnvCov = 0,
                                                sigma_e = ..2,
                                                sigma_d = sigma_d)

                          }

                          return(list(N = N,
                                      mu_r1 = ..1,
                                      epsilon_r1 = epsilon_r1,
                                      sigma_e = ..2,
                                      sigma_d = sigma_d,
                                      theta = ..3,
                                      K = ..4,
                                      tmax = tmax))

                        })

# Plot simulated data
plot_sim_data <- function(data) {

  p1 <- tibble::tibble(N = data$N,
                 x = seq_len(data$tmax)) %>%
    ggplot(aes(x = x, y = N)) +
    geom_line() +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    theme_classic() +
    labs(title = bquote(K~"="~.(data$K)*", "*bar(r)~"="~.(data$mu_r1)*", "*sigma[e]^2~"="~.(data$sigma_e^2)*", "*sigma[d]^2~"="~.(data$sigma_d^2)),
         x = "Time") +
    theme(axis.text = element_text(colour = "black"))

  point_data <- tibble::tibble(
    N = data$N[-data$tmax],
    r = diff(log(data$N))
  )

  line_data <- tibble::tibble(
    N = seq(0, max(data$N), length.out = 1000),
    r = data$mu_r1 - (0.5 * data$sigma_e^2) - data$mu_r1 * (((N^data$theta)-1) / ((data$K^data$theta)-1))
  )

  p2 <- ggplot(data = point_data, aes(x = N, y = r)) +
    geom_hline(yintercept = 0, colour = "gray70", linetype = "dashed") +
    geom_point(shape = 21) +
    geom_line(data = line_data, aes(x = N, y = r)) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    theme_classic() +
    labs(title = bquote(theta~"="~.(data$theta))) +
    theme(axis.text = element_text(colour = "black"))

  cowplot::plot_grid(p1, p2, ncol = 2, align = "h")

}

# Run both models
# run_models <- function(data) {
#
#   #-----------#
#   # A. NIMBLE #
#   #-----------#
#
#   ## Function to sample initial values
#   sample_inits <- function(){
#
#     list(
#       mu_r1 = rnorm(1, 1, 0.5),
#       sigma_e  = runif(1, 0, 1),
#       theta = rnorm(1, data$theta, 0.5),
#       epsilon_r1 = runif(data$tmax, -0.5, 0.5),
#       K = data$K,
#       initial_N = data$N[1]
#     )
#
#   }
#
#
#   ## Sample initial values
#   inits <- list(sample_inits(), sample_inits(), sample_inits())
#
#   ## Set data and constants
#   input_data <- list(N_obs = round(data$N), p_obs = rep(1, data$tmax))
#
#   input_constants <- list(tmax = data$tmax, max_N1 = data$N[1]*2,
#                           max_K = data$K*2, sigma_d = data$sigma_d)
#
#   ## Set parameters to monitor
#   params <- c("K", "theta", "sigma_e", "mu_r1", "gamma", "epsilon_r1", "N")
#
#   ## Set MCMC parameters
#   niter <- 35000
#   nburnin <- 30000
#   nthin <- 10
#   nchains <- 3
#
#   #-----------#
#   # Run model #
#   #-----------#
#
#   model_output <- nimbleMCMC(code = predict_N_nimble,
#                              constants = input_constants,
#                              data = input_data,
#                              inits = inits,
#                              monitors = params,
#                              niter = niter,
#                              nburnin = nburnin,
#                              thin = nthin,
#                              nchains = nchains,
#                              samplesAsCodaMCMC = TRUE)
#
#   #--------#
#   # B. MLE #
#   #--------#
#
#   tl_mle <- thetalogistic.est(data$N, year = 1:data$tmax, sd = data$sigma_d,
#                               nboot = nrow(as.matrix(model_output)))
#   #thetaplot(tl_mle)
#
#   return(list(nim_mod = model_output,
#               mle_mod = tl_mle))
#
# }
#
# # Plot model estimates against true values
# plot_mod_outputs <- function(data, mod) {
#
#   # Extract model outputs
#   model_output <- mod$nim_mod
#   tl_mle <- mod$mle_mod
#
#   # Posterior distributions from NIMBLE model
#   nim_est <- tibble::tibble(
#     x = rep(c("K", "mu_r1", "sigma_e", "gamma"),
#             each = nrow(as.matrix(model_output))),
#     val = c(as.matrix(model_output)[, "K"], as.matrix(model_output)[, "mu_r1"], as.matrix(model_output)[, "sigma_e"], as.matrix(model_output)[, "gamma"]),
#     type = "nimble",
#     chain = rep(rep(1:3, each = 250), 4)
#   )
#
#   # Bootstraps from MLE model
#   mle_est <- tibble::tibble(
#     x = c(rep(c("K", "mu_r1", "sigma_e", "gamma"), each = nrow(as.matrix(model_output)))),
#     val = c(tl_mle$boot.K, tl_mle$boot.r1, sqrt(tl_mle$boot.se), tl_mle$boot.gamma),
#     type = "mle"
#   )
#
#   # True values for parameters
#   data$gamma <- data$mu_r1*data$theta/(1-data$K^(-data$theta))
#
#   true_est <- tibble::tibble(
#     x = c("K", "mu_r1", "sigma_e", "gamma"),
#     val = c(data$K, data$mu_r1, data$sigma_e, data$gamma),
#     type = "true"
#   )
#
#   # Plot
#   bind_rows(nim_est, mle_est) %>%
#     ggplot(aes(x = val, color = type)) +
#     geom_segment(data = true_est, mapping = aes(x = val, xend = val,
#                                                 y = -Inf, yend = Inf),
#                  size = 0.75) +
#     geom_density(data = nim_est, mapping = aes(x = val, group = chain),
#                  color = colorspace::lighten("#9ed1cc"), size = 0.5) +
#     geom_density(size = 0.75) +
#     facet_wrap(~x, scales = "free") +
#     theme_classic(base_family = "Roboto") +
#     labs(x = "", y = "Density") +
#     scale_color_manual(values = c("#006D77", "#53b5ab", "#E29578"),
#                        name = "") +
#     scale_x_continuous(breaks = scales::pretty_breaks()) +
#     scale_y_continuous(breaks = scales::pretty_breaks())
#
# }
#
# sim_plot_list <- purrr::map(.x = sim_data,
#                             .f = ~{plot_sim_data(.x)})
#
# mod_list <- purrr::map2(.x = sim_data,
#                         .y = 1:length(sim_data),
#                         .f = ~{
#
#                           cat("Running models on simulated data (", .y, "/16)...\n", sep = "")
#                           run_models(.x)
#
#                         })
#
# sim_plot_list[[1]]
# mod <- run_models(sim_data[[1]])
# coda::gelman.diag(mod$nim_mod)
# #bayesplot::mcmc_dens_overlay(mod$nim_mod, pars = "theta")
# plot_mod_outputs(sim_data[[1]], mod)
#
# #purrr::map(mod_list, ~{coda::gelman.diag(.x$nim_mod)$mpsrf})
#
# save(pars, sim_data, mod_list, predict_N, predict_N_nimble, plot_sim_data, run_models, plot_mod_outputs, thetalogistic.est, thetaplot, file = "theta-comp.rda")
