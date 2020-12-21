library(nimble)
nimbleOptions(disallow_multivariate_argument_expressions = FALSE)

mySeed <- 0
set.seed(mySeed)

# Basic Nimble function for population growth #
#---------------------------------------------#

N.predict <- nimbleFunction(
  run = function(
    N_current = double(),
    Mu.r1 = double(),
    
    epsilon.r1 = double(),
    sigma_e = double(),
    sigma_d = double(),
    theta = double(),
    K = double(),
    beta.r1 = double(1),
    EnvCov = double(1)
    ) {
    
    r1 <- Mu.r1 + sum(beta.r1*EnvCov) + epsilon.r1 # epsilon.r1 ~ dnorm(0, sd = sigma_e)
    
    # K <- [some function of beta] 
    
    s <- r1 - (0.5 * sigma_e^2) - (0.5 * sigma_d^2 / N_current)
    
    r <- s - Mu.r1 * ((N_current^theta)-1 / (K^theta)-1) 
    
    N_next <- exp(log(N_current) + r)

    if(is.na(N_next)) stop('Predicted population size (N_next) is NA')

    returnType(double())
    return(N_next)
  }
)


# Data simulation # - NEEDS UPDATING
#-----------------#

## NOTE: Simulate data with population sizes both very low and fluctuating around K

## Set parameter values 
Tmax <- 100
Mu.r1 <- 0.3 # up to 1
sigma.e <- sqrt(0.01)
sigma.d <- sqrt(0.5)

## Simulate random effects
epsilon.r <- rnorm(Tmax, 0, sigma.r)

## Prepare population vector and set initial population size
N <- rep(NA, Tmax)
N[1] <- 50

## Use nimble function to predict population size over time
for(t in 1:(Tmax-1)){
  N[t+1] <- N.predict(N_current = N[t], Mu.r = Mu.r, epsilon.r = epsilon.r[t])
}

plot(N, type = 'l')


# Nimble model code #
#-------------------#

N.predict.nimble <- nimbleCode({
  
  #-----------------------------------#
  # PROCESS MODEL (POPULATION GROWTH) #
  #-----------------------------------#
  
  ## Initial population size
  N[1] <- initialN
  
  ## Population growth
  for(t in 1:(Tmax-1)){
    N[t+1] <- N.predict(N_current = N[t], Mu.r = Mu.r, epsilon.r = epsilon.r[t], EnvCov = EnvCov[1:NoCovariates,t], beta.r1 = , ...)
    # If EnvCov and beta.r1 are set to 0, have to pass that 0 as a vector (look up how to do again)
  }
  
  # NOTE: When fitting a model without REs, we can simply set epsilon.r = 0 in the function call
  
  #-------------------#
  # OBSERVATION MODEL #
  #-------------------#
  
  ## Binomial observation model
  #for(t in 1:Tmax){
  #  N_obs[t] ~ dbin(pObs[t], round(N[t]))
  # NOTE: The binomial observation model does not work yet.
  #       At the moment, it's "too rigid" this way because the N.predict
  #       function already does not allow for stochastic outcomes (so whenever
  #       N.predict returns an N that does not satify the condition round(N) = N_obs
  #       then there will be an infinite log probability. 
  
  ## Alternative: Poisson observation model
  for(t in 1:(Tmax-1)){
    N_obs[t] ~ dpois(pObs[t]*N[t])
  }
  
  ## Alternative: Normal observation model
  #for(t in 1:(Tmax-1)){
  #  N_obs[t] ~ dnorm(pObs[t]*N[t], sd = sigma.obs)
  #}
  
  #------------------------#
  # PRIORS AND CONSTRAINTS #
  #------------------------#
  
  ## Random effects
  for(t in 1:Tmax){
    epsilon.r1[t] ~ dnorm(0, var = sigma.e^2 + ((sigma.d^2)/N))
  }
  
  ## Priors
  Mu.r1 ~ dunif(-50, 50)
  sigma.e ~ dunif(0, 10)
  # sigma.d (--> added as data)
  K ~ dunif(1, maxK)
  theta ~ dunif(-1, 20)
  
  for(i in 1:NoCovariates){ # Check if dimension call works in BUGS?
    beta.r1[i] ~ dunif(-5, 5)
  }
  
  initialN ~ dunif(1, maxN1)
  
  #sigma.obs ~ dunif(0, 50)
  
})
  

# Initial values #
#----------------#

## Function to sample initial values
initsFun <- function(){
  list(
    Mu.r = runif(1, -0.01, 0.01),
    sigma.r = runif(1, 0, 0.5),
    epsilon.r = rep(0, Tmax),
    initialN = N[1]#,
    #sigma.obs = runif(1, 0, 5)
  )
}

## Sample initial values
#Inits <- list(initsFun())
Inits <- list(initsFun(), initsFun(), initsFun())


# Nimble model and MCMC setup #
#-----------------------------#

## Set data and constants
input.data <- list(N_obs = round(N), pObs = rep(1, Tmax))

input.constants <- list(Tmax = Tmax, maxN1 = 200)

## Set parameters to monitor
params <- c('Mu.r', 'sigma.r', 'N')

## Set MCMC parameters
niter <- 50000
nburnin <- 10000
nthin <- 10

#niter <- 2
#nburnin <- 0
#nthin <- 1

#nchains <- 1
nchains <- 3

# Run model #
#-----------#
ModelRun <- nimbleMCMC(code = N.predict.nimble, 
                       constants = input.constants, 
                       data = input.data, 
                       inits = Inits, 
                       monitors = params, 
                       niter = niter, 
                       nburnin = nburnin, 
                       nchains = nchains, 
                       setSeed = mySeed, 
                       samplesAsCodaMCMC = TRUE)


# Plot results #
#--------------#
plot(ModelRun, ask = T)


# Extra: Model assembly step (for interactive exploration) #
#----------------------------------------------------------#

NimModel <- nimbleModel(code = N.predict.nimble, 
                        constants = input.constants, 
                        data = input.data, 
                        inits = Inits[[1]])

NimModel$check()
NimModel$expandNodeNames('N')
NimModel$getDependencies('N')
NimModel$initializeInfo()
