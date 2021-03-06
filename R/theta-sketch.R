library(nimble)
nimbleOptions(disallow_multivariate_argument_expressions = FALSE)

mySeed <- 74
set.seed(mySeed)

#------------------------------------------------#
# Basic Nimble function for population growth ####
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

#--------------------#
# Data simulation ####
#--------------------#

## NOTE: Simulate data with population sizes both very low and fluctuating around K

## Set parameter values
tmax <- 30
mu_r1 <- 0.7 # up to 1
sigma_e <- sqrt(0.01)
sigma_d <- sqrt(0.4)
K <- 20
theta <- 2

## Simulate random effects

## Prepare population vector and set initial population size
N <- rep(NA, tmax)
N[1] <- 2

epsilon_r1 <- rep(NA, tmax-1)

## Use nimble function to predict population size over time
for(t in 1:(tmax-1)){

  epsilon_r1[t] <- rnorm(1, 0, sqrt(sigma_e^2 + ((sigma_d^2) / N[t])))

  N[t+1] <- predict_N(N_current = N[t], mu_r1 = mu_r1, epsilon_r1 = epsilon_r1[t],
                      K = K, theta = theta, #beta.r1 = 0, EnvCov = 0,
                      sigma_e = sigma_e,
                      sigma_d = sigma_d)

}
ts.plot(N)
plot(diff(log(N))~N[-length(N)])
sigmaY <- 0.2
logY <- rnorm(length(N), log(N) - 0.5*sigmaY^2, sigmaY)
ts.plot(cbind(N,exp(logY)), col=1:2)

#----------------------#
# Nimble model code ####
#----------------------#

predict_N_nimble <- nimbleCode({

  for (t in 1:(tmax-1)){
    #obs.r[t] ~ dnorm(pred.r[t], var = var_r1[t])
    # state process
    s[t] <- mu_r1 - (0.5 * sigma_e^2) - (0.5 * sigma_d^2 / N[t])
    pred.r[t] <- s[t] - mu_r1 * (((N[t]^theta)-1) / ((K^theta)-1))
    var_r1[t]  <- sigma_e^2 + ((sigma_d^2) / N[t])
    eps[t] ~ dnorm(0, var = var_r1[t])
    logN[t+1] <- logN[t] + pred.r[t] + eps[t]
    N[t+1] <- exp(logN[t+1])
    # observation process
    logY[t] ~ dnorm(logN[t] - 0.5 * sigmaY^2, sd = sigmaY)
    # Y[t] ~ dpois(N[t])
  }

  ## Priors
  mu_r1 ~ dlnorm(meanlog = 0, sdlog = 0.5)
  sigma_e ~ dunif(0, 1)
  K ~ dunif(1, max_K)
  #theta ~ dunif(-1, 10)
  theta ~ dlogis(location = 3, scale = 1)
  sigmaY ~ dunif(0, 1)
  N[1] ~ dunif(1, max_K)
  logN[1] <- log(N[1])
})


predict_N_nimble_WO_obserror <- nimbleCode({

  for (t in 1:(tmax-1)){
    obs.r[t] ~ dnorm(pred.r[t], var = var_r1[t])
    s[t] <- mu_r1 - (0.5 * sigma_e^2) - (0.5 * sigma_d^2 / N[t])
    pred.r[t] <- s[t] - mu_r1 * (((N[t]^theta)-1) / ((K^theta)-1))
    var_r1[t]  <- sigma_e^2 + ((sigma_d^2) / N[t])
  }

  ## Priors
  mu_r1 ~ dlnorm(meanlog = 0, sdlog = 0.5)
  sigma_e ~ dunif(0, 1)
  K ~ dunif(1, max_K)
  #theta ~ dunif(-1, 10)
  theta ~ dlogis(location = 3, scale = 1)
})



#-------------------#
# Initial values ####
#-------------------#


#--------------------------------#
# Nimble model and MCMC setup ####
#--------------------------------#

## Set data and constants
input_data <- list(logY=logY)

input_constants <- list(tmax = tmax, max_K = K*2, sigma_d = sigma_d)

## Set parameters to monitor
params <- c("K", "theta", "sigma_e", "mu_r1", "sigmaY")

## Set MCMC parameters
niter <- 10000
nburnin <- 5000
nthin <- 1
nchains <- 1

#niter <- 2
#nburnin <- 0
#nthin <- 1
#nchains <- 1

#--------------#
# Run model ####
#--------------#

mu_r1 ~ dlnorm(meanlog = 0, sdlog = 0.5)
sigma_e ~ dunif(0, 1)
K ~ dunif(1, max_K)
#theta ~ dunif(-1, 10)
theta ~ dlogis(location = 3, scale = 1)
sigmaY ~ dunif(0, 1)
N[1] ~ dunif(1, max_K)
logN[1] <- log(N[1])

## Function to sample initial values
sample_inits <- function(){

  list(
    mu_r1 = rnorm(1, 1, 0.5),
    sigma_e  = runif(1, 0, 0.2),
    theta = rnorm(1, 2, 0.5),
    K = K,
    pred.r = diff(logY),
    s = diff(logY),
    var_r1 = 0.001,
    sigmaY=0.1,
    eps = rnorm(length(logY)-1, 0, 0.0001)
    #eps = runif(length(logY)-1, 0, 0.0001)
    #sigma_obs = runif(1, 0, 5)
  )

}

## Sample initial values
#inits <- list(sample_inits())
inits <- list(sample_inits())

mod <- nimbleModel(predict_N_nimble, constants = input_constants,
                   data = input_data, inits=inits[[1]])
configureMCMC(mod)
modR <- buildMCMC(mod, niter =1000, nburnin= 1000)
modC <- compileNimble(mod)
tr <- runMCMC(modC)

mod$initializeInfo()


Rmodel <- nimbleModel(code)
Rmodel$setData(list(x = c(2, 5, 3, 4, 1, 0, 1, 3, 5, 3)))
Rmcmc <- buildMCMC(Rmodel)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
inits <- function() list(mu = rnorm(1,0,1), sigma = runif(1,0,10))
samplesList <- runMCMC(Cmcmc, niter = 10000, nchains = 3, inits = inits)


model_output <- nimbleMCMC(code = predict_N_nimble,
                           constants = input_constants,
                           data = input_data,
                           inits = inits,
                           monitors = params,
                           niter = niter,
                           nburnin = nburnin,
                           thin = nthin,
                           nchains = nchains,
                           setSeed = mySeed,
                           samplesAsCodaMCMC = TRUE)

#-----------------#
# Plot results ####
#-----------------#

plot(model_output, ask = T)

coda::gelman.diag(model_output)

l <- purrr::map2(.x = c("K", "mu_r1", "sigma_e", "theta"),
                 .y = c(K, mu_r1, sigma_e, theta),
                 .f = ~{

                   #bayesplot::mcmc_combo(model_output, pars = .x)
                   p <- bayesplot::mcmc_dens(model_output, pars = .x) +
                     geom_vline(xintercept = .y, color = "#03396c")
                   q <- bayesplot::mcmc_trace(model_output, pars = .x)

                   cowplot::plot_grid(p, q, align = "h", nrow = 1, ncol = 2)

                 })

#ggsave("test.pdf", gridExtra::arrangeGrob(grobs = l))

# Create multi-page plots
pdf("est.pdf", width = 6, height = 4)
invisible(lapply(l[1:4], print))
dev.off()

ggsave("N.pdf", gridExtra::arrangeGrob(grobs = l[5:104]),
       width = 45, height = 25, dpi = 450)

#-------------------------------------------------------------#
# Extra: Model assembly step (for interactive exploration) ####
#-------------------------------------------------------------#

nimble_model <- nimbleModel(code = predict_N_nimble,
                            constants = input_constants,
                            data = input_data,
                            inits = inits[[1]])

nimble_model$check()
nimble_model$expandNodeNames("N")
nimble_model$getDependencies("N")
nimble_model$initializeInfo()

#---------------#
# Comparison ####
#---------------#

# Maximum likelihood estimation procedure for theta-logistic model
tl_mle <- thetalogistic.est(N, year = 1:35, sd = sigma_d, nboot = nrow(as.matrix(model_output)))
thetaplot(tl_mle)

# Posterior distributions from NIMBLE model
nim_est <- tibble::tibble(
  x = rep(c("K", "mu_r1", "sigma_e", "theta"), each = nrow(as.matrix(model_output))),
  val = c(as.matrix(model_output)[, "K"], as.matrix(model_output)[, "mu_r1"], as.matrix(model_output)[, "sigma_e"], as.matrix(model_output)[, "theta"]),
  type = "nimble",
  chain = rep(rep(1:3, each = 500), 4)
)

# Bootstraps from MLE model
mle_est <- tibble::tibble(
  x = c(rep(c("K", "mu_r1", "sigma_e"), each = nrow(as.matrix(model_output)))),
  val = c(tl_mle$boot.K, tl_mle$boot.r1, sqrt(tl_mle$boot.se)),
  type = "mle"
)

# Estimate for theta from MLE model
mle_theta <- tibble::tibble(
  x = "theta",
  val = tl_mle$theta,
  type = "mle"
)

# True values for parameters
true_est <- tibble::tibble(
  x = c("K", "mu_r1", "sigma_e", "theta"),
  val = c(K, mu_r1, sigma_e, theta),
  type = "true"
)

# Plot
bind_rows(nim_est, mle_est) %>%
  ggplot(aes(x = val, color = type)) +
  geom_density(data = nim_est, mapping = aes(x = val, group = chain),
               color = colorspace::lighten("#9ed1cc"), size = 0.5) +
  geom_density(size = 0.75) +
  geom_segment(data = true_est, mapping = aes(x = val, xend = val,
                                              y = -Inf, yend = Inf),
               size = 0.75) +
  geom_segment(data = mle_theta, mapping = aes(x = val, xend = val,
                                               y = -Inf, yend = Inf),
               size = 0.75) +
  facet_wrap(~x, scales = "free") +
  theme_classic(base_family = "Roboto") +
  labs(x = "", y = "Density") +
  scale_color_manual(values = c("#006D77", "#53b5ab", "#E29578"),
                     name = "") +
  scale_x_continuous(breaks = scales::pretty_breaks())
