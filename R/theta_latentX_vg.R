
#---------------------------------------------------#
# Simulate data from thetalogistic approximation ####
#---------------------------------------------------#
library(nimble)
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


simData <- function(){
  pops <<- 3
  tmax <<- 50

  # simulate from prior distribution

  #initial_N[i] ~ dunif(0, max_N[i])
  #N[i,1] <- initial_N[i]
  #mu_r1[i] ~ dunif(0, 2)
  #sigma_e2[i] ~ dunif(0, 0.2)
  #K[i] ~ dunif(1, 200)
  #theta[i] ~ dunif(-1, 3)

  r0 <<- runif(pops, 0.3, 1)
  K <<- runif(pops, 50, 200)
  sigma_e2 <<- runif(pops, 0.01, 0.02)
  sigma_e <<- sqrt(sigma_e2)
  sigma_d2 <<- 0.5
  sigma_d <<- sqrt(sigma_d2)
  theta <<- runif(pops, 0.1, 2.2)

  mu_r1 <<- r0 / (1 - K^(-theta))
  gamma <<- r0 * theta

  # True N

  N <- matrix(NA, nrow = pops, ncol = tmax)
  max_N <<- K*1.2
  N[,1] <- runif(pops, 10, 20)
  #N[,1] <- K

  epsilon_r1 <- matrix(NA, nrow = pops, ncol = tmax - 1)

  for(i in 1:pops) {
    for(t in 1:(tmax-1)){

      epsilon_r1[i, t] <- rnorm(1, 0, sqrt(sigma_e[i]^2 + ((sigma_d^2) / N[i, t])))

      N[i, t+1] <- predict_N(N_current = N[i, t], mu_r1 = mu_r1[i], epsilon_r1 = epsilon_r1[i, t],
                             K = K[i], theta = theta[i],
                             sigma_e = sigma_e[i],
                             sigma_d = sigma_d)

    }
  }


  # Observed N
  # Draw from Poisson model
  obs_N <<- matrix(rpois(nrow(N)*ncol(N),N), nrow(N), ncol(N))
  return(NULL)

}


TRY <- TRUE
gamma <- NA
while(TRY){
  simTry <- try(simData(), silent = TRUE)
  if (class(simTry)!="try-error"){
    TRY <- FALSE
  }
}

str(data_list)
str(parameters)

ts.plot(t(obs_N))



#-----------------------------#
# X latent variable model     #
#-----------------------------#

predict_N_random_unbiased <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  for(i in 1:pops){
    for(t in 1:(tmax-1)){

      #log(N[t+1]) <- log(N[t]) + r0 * (1 - (N[t] / K)^theta) + eps_e[t] + eps_d[t]
      N[i, t+1] <- exp(X[i, t+1])
      E_X[i, t+1] <- X[i, t] + mu_r1[i] * (1 - (((N[i, t]^theta[i]) - 1) / ((K[i]^theta[i]) - 1))) - (sigma_e2[i] / 2) - (sigma_d2[i] / (2 * N[i, t]))
      VAR_X[i, t+1] <- sigma_e2[i] + sigma_d2[i] / N[i, t]
      X[i, t+1] ~ dnorm(E_X[i, t+1], var = VAR_X[i, t+1])

      # eps now becomes a derived parameter
      # I keep it here instead of in "derived parameter section" because
      # we later may want to fit a distribution to epsilons (but how to do that
      # without making eps an parameter???)
      # eps[i,t+1] or eps[i, t] is a matter of preference / definition -
      # but important to be aware of / remember the indexing
      eps[i, t] <- X[i, t+1] - E_X[i, t+1]
      # we may want to save the demographic component of the variance later:
      demcomp[i, t] <- sigma_d2[i] / N[i, t]
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

    X[i, 1] ~ dunif(0, maxK[i])
    N[i, 1] <- exp(X[i, 1])
    r0[i] ~ dunif(0.01, 2)
    sigma_e2[i] ~ dunif(0.0001, 0.5)
    K[i] ~ dunif(10, maxK[1])
    theta[i] ~ dunif(0.1, 10)
    mu_r1[i] <- r0[i] / (1 - K[i]^(-theta[i]))
    gamma[i] <- theta[i] * mu_r1[i] / (1 - K[i]^(-theta[i]))

  }

})

sample_inits2 <- function(){

  list(
    #r0 = rnorm(1, 1, 0.5),
    r0 = runif(pops, 0.2, 0.8),
    sigma_e2  = runif(pops, 0, 0.1),
    theta = rnorm(pops, 2, 0.5),
    K = runif(pops, apply(obs_N, 1, quantile, prob=1/3), apply(obs_N, 1, quantile, prob=4/5)),
    X = log(obs_N)
  )

}

input_data2 <- list(obs_N = obs_N)

input_constants2 <- list(tmax = tmax, sigma_d2 = rep(sigma_d2, pops), pops = pops, maxK = apply(obs_N,1,max) * 1.5)

params2 <- c("K", "theta", "sigma_e2", "gamma", "mu_r1", "r0")

# Set MCMC parameters
niter <- 30000
nburnin <- 5000
nthin <- 10
nchains <- 3

Model
start <- Sys.time()


inits2 <- list(sample_inits2(), sample_inits2(), sample_inits2())


mod2 <- nimbleMCMC(code = predict_N_random_unbiased,
                   constants = input_constants2,
                   data = input_data2,
                   inits = inits2,
                   monitors = params2,
                   niter = niter,
                   nburnin = nburnin,
                   thin = nthin,
                   nchains = nchains,
                   samplesAsCodaMCMC = TRUE)

dur2 <- Sys.time() - start


library(coda)

gelman.diag(mod2)

round(summary(mod2)$statistics,2)

traceplot(mod2)





