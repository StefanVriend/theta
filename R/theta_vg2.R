setwd("C:/Users/vidargr/OneDrive - NTNU/Prosjekt/theta")

library(tidyverse)
library(nimble)
library(viridis)
library(coda)
library(bayesplot)
# library(extrafont)
# extrafont::loadfonts(device = "win")

nimbleOptions(disallow_multivariate_argument_expressions = FALSE)

#seed <- 18
#set.seed(seed)



#---------------#
# RE model s ####
#---------------#

predict_N_random_unbiased <- nimbleCode({

  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#

  for(i in 1:pops){
    for(t in 1:(tmax-1)){

      log(N[i, t+1]) <- log(N[i, t]) + mu_r1[i] * (1 - (((N[i, t]^theta[i]) - 1) / ((K[i]^theta[i]) - 1))) -
        (sigma_e2[i] / 2) - (sigma_d2[i] / (2 * N[i, t])) + eps_ed[i, t]

      eps_ed[i, t] ~ dnorm(0, var = sigma_e2[i] + sigma_d2[i] / N[i, t])
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

    mu_r1[i] ~ dunif(0, 2)
    sigma_e2[i] ~ dunif(0, 0.2)
    K[i] ~ dunif(1, 200)
    theta[i] ~ dunif(-1, 3)

  }



  #sigma_obs ~ dunif(0, 100)

  #--------------------#
  # DERIVED PARAMETERS #
  #--------------------#

  for(i in 1:pops) {

    gamma[i] <- theta[i] * mu_r1[i] / (1 - K[i]^(-theta[i]))

  }


})





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

  mu_r1 <<- runif(pops, 0, 2)
  K <<- runif(pops, 1, 200)
  sigma_e2 <<- runif(pops, 0, 0.2)
  sigma_e <<- sqrt(sigma_e2)
  sigma_d2 <<- 0.5
  sigma_d <<- sqrt(sigma_d2)
  theta <<- runif(pops, -1, 3)

  r0 <<- mu_r1 / (1 - K^(-theta))
  gamma <<- r0 * theta

  # True N

  N <- matrix(NA, nrow = pops, ncol = tmax)
  max_N <<- K*1.2
  N[,1] <- runif(pops, 0, max_N)

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


sample_inits2 <- function(){

  list(
    #r0 = rnorm(1, 1, 0.5),
    mu_r1 = rnorm(pops, 1, 0.5),
    sigma_e2  = runif(pops, 0, 0.2),
    theta = rnorm(pops, 2, 0.5),
    K = K,
    eps_ed = matrix(0, pops, tmax-1),
    initial_N = apply(cbind(obs_N[,1],max_N-0.01), 1, min)#,
    #sigma_obs = runif(1, 0, 10)
  )

}


write_q_vals <- function(){
  # include only pops that show approximate convergence
  tr <- matrix(gelman.diag(mod2)$psrf[,2], byrow=TRUE, ncol=pops)
  print(tr)
  sel <- apply(tr, 2, function(x) all(x<1.1))

  bind <- rbind(mod2$chain1, mod2$chain2,  mod2$chain3)

  get_q <- function(par, simVals){
    indx <- grep(par, colnames(bind))
    mat <- bind[,indx]

    for (j in 1:length(sel)){
      if (sel[j]){
        q <- NA
        q <- sum(mat[,j]<simVals[j])/nrow(mat)
        cat(q, "\n", file=paste0("q_", par, ".out"), append = TRUE)
      }
    }

  }

  get_q("K", K)
  get_q("theta", theta)
  get_q("sigma_e2", sigma_e2)
  get_q("gamma", gamma)
  get_q("mu_r1", mu_r1)

}


params2 <- c("K", "theta", "sigma_e2", "gamma", "mu_r1")

# Set MCMC parameters
niter <- 300000
nburnin <- 150000
nthin <- 75
nchains <- 3


for (i in 1:1000){
  gamma <- NA
  simTry <- try(simData(), silent = TRUE)
  if (class(simTry)!="try-error"){
    input_data2 <- list(obs_N = obs_N)

    input_constants2 <- list(tmax = tmax, max_N = max_N, sigma_d2 = rep(sigma_d2, pops), pops = pops)

    mod2 <- try(nimbleMCMC(code = predict_N_random_unbiased,
                           constants = input_constants2,
                           data = input_data2,
                           inits = sample_inits2(),
                           monitors = params2,
                           niter = niter,
                           nburnin = nburnin,
                           thin = nthin,
                           nchains = nchains,
                           #setSeed = seed2,
                           samplesAsCodaMCMC = TRUE), silent = TRUE)
    if (class(mod2)!="try-error"){
      try(write_q_vals(), silent = TRUE)
    }

  }
}




theta.hat <- unique(scan("q_theta.out"))
K.hat <- unique(scan("q_K.out"))
mu_r1.hat <- unique(scan("q_mu_r1.out"))
gamma.hat <- unique(scan("q_gamma.out"))
sigma_e2.hat <- unique(scan("q_sigma_e2.out"))


hist(theta.hat)
hist(K.hat)
hist(mu_r1.hat)
hist(gamma.hat)
hist(sigma_e2.hat)

K.hat[K.hat>0.999999] <- 0.99999
theta.hat[theta.hat>0.999999] <- 0.99999

shapiro.test(qnorm(theta.hat))
shapiro.test(qnorm(K.hat))
shapiro.test(qnorm(mu_r1.hat))
shapiro.test(qnorm(gamma.hat))
shapiro.test(qnorm(sigma_e2.hat))





# F****K NIMBLE!!!
Rmodel <- nimbleModel(code = predict_N_random_unbiased,
                      constants = input_constants2, data = input_data2,
                      inits = sample_inits2())
#Rmodel <- nimbleModel(code = predict_N_random_unbiased,data = input_data2)



MCMCsamples <- as.matrix(Cmodel$Rmcmc$mvSamples)


Cmcmc <- compileNimble(Rmcmc, project="Rmodel")




