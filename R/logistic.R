# Simulate data
simulate_data <- function(init_X, sigma_e2, sigma_d2, r, beta, tmax, plot = TRUE, seed = NULL) {
  # init_X: initial value for X
  # sigma_e2: desired value for sigma_e2
  # sigma_d2: desired value for sigma_d2
  # r: desired value for r
  # beta: desired value for beta
  # tmax: time series length
  # plot: plot time series: TRUE or FALSE
  # seed: set seed to RNG

  if(!is.null(seed)) {

    set.seed(seed)

  }

  X <- rep(NA, tmax)
  X[1] <- init_X

  # Simulate Xs
  for(i in 1:(tmax-1)) {

    N <- exp(X[i])

    E_X <- X[i] + r - (sigma_e2 / 2) - (sigma_d2 / (2 * N)) + beta * N

    Var_X <- sigma_e2 + (sigma_d2 / N)

    X[i+1] <- E_X + rnorm(1, 0, sqrt(Var_X))

  }

  # Plot time series
  plot(1:tmax, X, type = "l")

  return(list(X = X, sigma_e2 = sigma_e2, sigma_d2 = sigma_d2, r = r, beta = beta, tmax = tmax))

}

# Run logistic or density-independent population model
run_population_model <- function(X, sigma_d2, n_boot = NULL, n_time = length(X), hessian = FALSE) {
  # X: observed log population sizes
  # sigma_d2: demographic variance
  # dd: form of density dependence:
  # --- "logistic": logistic form of density dependence
  # --- "independent": density-independent
  # n_boot: number of replicates for simulation/bootstrap
  # n_time: number of time steps to simulate data for
  # hessian: return Hessian matrix (TRUE or FALSE)

  par_names <- c("sigma_e2", "r", "beta", "K") # Parameters of interest
  boot <- NULL

  #-----------#
  # FUNCTIONS #
  #-----------#

  # Log-likelihood function of density-dependent/-independent population model
  log_likelihood <- function(par, X, sigma_d2, dd = "logistic") {
    # par: model parameters. Number is dependent on the form of density dependence
    # X: observed log population sizes
    # sigma_d2: demographic variance
    # dd: form of density dependence:
    # --- "logistic": logistic form of density dependence
    # --- "independent": density-independent

    sigma_e2 <- exp(par[1]) # environmental variance
    r <- par[2] # deterministic growth rate

    if(dd == "logistic") {

      beta <- par[3] # "intraspecific competition coefficient"

    } else if(dd == "independent") {

      beta <- 0

    }

    N <- exp(X) # Population size on arithmetic scale

    E_X <- X + r - (sigma_e2 / 2) - (sigma_d2 / (2 * N)) + beta * N # Expectation

    Var_X <- sigma_e2 + (sigma_d2 / N) # Variance

    -sum(dnorm(x = X[-1], mean = E_X[-length(X)], sd = sqrt(Var_X[-length(X)]), log = TRUE), na.rm = TRUE) # Log-likelihood

  }

  # Estimate parameters via optimization
  estimate_pars <- function(X, sigma_d2, dd, start, hessian) {
    # X: observed log population sizes
    # sigma_d2: demographic variance
    # dd: form of density dependence:
    # --- "logistic": logistic form of density dependence
    # --- "independent": density-independent
    # start: initial values for parameters to be optimized over; e.g. output of `set_starting_values()`
    # hessian: return Hessian matrix (TRUE or FALSE)

    optim(par = start, fn = log_likelihood, X = X, sigma_d2 = sigma_d2, dd = dd,
          hessian = hessian)

  }

  # Set starting values for parameters
  set_starting_values <- function(X, dd = "logistic") {
    # X: observed log population sizes
    # dd: form of density dependence:
    # --- "logistic": logistic form of density dependence
    # --- "independent": density-independent

    # Fit linear regression of per capita growth rate versus population size
    # log(N[t+1] / N[t]) versus N
    fit <- lm(diff(X) ~ exp(X[-length(X)]))

    sigma_e2_hat <- var(residuals(fit)) # Residuals of lm
    r_hat <- coef(fit)[1] # Intercept of lm

    if(dd == "logistic") {

      beta_hat <- coef(fit)[2] # Slope of lm

      # Starting values
      start <- c(sigma_e2_hat, r_hat, beta_hat)
      names(start) <- c("log_sigma_e2", "r", "beta")

    } else if(dd == "independent") {

      # Starting values
      start <- c(sigma_e2_hat, r_hat)
      names(start) <- c("log_sigma_e2", "r")

    }

    return(start)

  }

  # Calculate information criteria
  calculate_criteria <- function(X, est) {
    # X: observed log population sizes
    # est: output of `estimate_pars()`

    n_pars <- length(est$par) # Number of parameters
    n_obs <- length(X[!is.na(X)]) # Number of observations

    # NB: optim minimizes, so est$value is a minimum - not a maximum

    AIC <- (2 * est$value) + (2 * n_pars)
    AICc <- AIC + (((2 * n_pars) * (n_pars + 1)) / (n_obs - n_pars - 1))
    BIC <- (2 * est$value) + (n_pars * log(n_obs))

    return(list(AIC = AIC,
                AICc = AICc,
                BIC = BIC))

  }

  # Simulate time series based on estimated parameters
  simulate_ts <- function(X, sigma_d2, est, dd = "logistic", init_X = NULL, n_time) {
    # X: observed log population sizes
    # sigma_d2: demographic variance
    # est: output of `estimate_pars()`
    # dd: form of density dependence:
    # --- "logistic": logistic form of density dependence
    # --- "independent": density-independent
    # init_X: initial value for X. By default, it takes the first value in X.
    #         Can be set to any value (e.g. to K when calculating the stationary distribution)
    # n_time: number of time steps to simulate for

    # Retrieve parameter estimates
    sigma_e2 <- exp(est$par[1])
    r <- est$par[2]

    if(dd == "logistic") {

      beta <- est$par[3]

    } else if(dd == "independent") {

      beta <- 0

    }

    first <- min(which(!is.na(X))) # Get first non-NA observation
    X_sim <- rep(NA, n_time)

    # Set initial value of X_sim
    if(is.null(init_X)) {


      X_sim[first] <- X[first] # First value in observations (X)

    } else if(!is.null(init_X))  {

      X_sim[first] <- init_X

    }


    while(any(is.na(X_sim)) | any(X_sim == 0)){

      for (i in first:(n_time-1)){

        N <- exp(X_sim[i]) # Population size on arithmetic scale

        E_X <- X_sim[i] + r - (sigma_e2 / 2) - (sigma_d2 / (2 * N)) + beta * N # Expectation

        Var_X <- sigma_e2 + (sigma_d2 / N) # Variance

        X_sim[i+1] <- E_X + rnorm(1, 0, sqrt(Var_X))

        if(X_sim[i+1] < 0) X_sim[i+1] <- 0

      }

    }

    X_sim[is.na(X)] <- NA # NA in observations -> NA in simulations

    return(X_sim)

  }

  # Bootstrap uncertainty
  bootstrap <- function(X, sigma_d2, est, dd = "logistic", n_boot, n_time, hessian) {
    # X: observed log population sizes
    # sigma_d2: demographic variance
    # est: output of `estimate_pars()`
    # dd: form of density dependence:
    # --- "logistic": logistic form of density dependence
    # --- "independent": density-independent
    # n_boot: number of replicates for simulation/bootstrap
    # n_time: number of time steps to simulate data for
    # hessian: return Hessian matrix (TRUE or FALSE)

    # Simulate data
    sim_data <- matrix(NA, nrow = n_boot, ncol = n_time)

    pb_sim <- progress::progress_bar$new(total = n_boot, clear = FALSE,
                                         format = "Simulating data [:bar] :percent eta: :eta")

    for(i in 1:n_boot) {

      pb_sim$tick() # Update progress bar

      sim_data[i,] <- simulate_ts(X = X, sigma_d2 = sigma_d2, est = est, dd = dd, n_time = n_time) # Simulate time series

    }

    # Bootstrap
    boot <- matrix(NA, nrow = n_boot, ncol = 4)
    colnames(boot) <- par_names

    pb_boot <- progress::progress_bar$new(total = n_boot, clear = FALSE,
                                          format = "Bootstrapping [:bar] :percent eta: :eta")

    for(i in 1:n_boot) {

      pb_boot$tick() # Update progress bar

      X_boot <- sim_data[i,] # Retrieve set of simulated data

      est_boot <- estimate_pars(X_boot, sigma_d2 = sigma_d2, dd = dd,
                                start = set_starting_values(X_boot, dd = dd), hessian = hessian) # Estimate parameters

      # Save parameters
      boot[i, 1] <- exp(est_boot$par[1]) # sigma_e2
      boot[i, 2] <- est_boot$par[2] # r

      if(length(est_boot$par) > 2) {

        boot[i, 3] <- est_boot$par[3] # beta
        boot[i, 4] <- -est_boot$par[2] / est_boot$par[3] # carrying capacity K

      }

    }

    return(list(boot = boot,
                sim = sim_data,
                n_boot = n_boot,
                n_time = n_time))

  }

  #------------#
  # RUN MODELS #
  #------------#

  # Estimate parameters
  estimates <- estimate_pars(X = X, sigma_d2 = sigma_d2, dd = "logistic",
                             start = set_starting_values(X, dd = "logistic"),
                             hessian = hessian)

  # Bootstrap uncertainty
  if(!is.null(n_boot)) {

    boot <- bootstrap(X, sigma_d2, estimates, dd = "logistic", n_boot = n_boot, n_time = n_time, hessian = hessian)

  }

  #----------------#
  # RETURN OUTPUTS #
  #----------------#

  sigma_e2 <- exp(estimates$par[1]) # Estimate for sigma_e2
  names(sigma_e2) <- "sigma_e2" # Set name to sigma_e2 (instead of log_sigma_e2)
  r <- estimates$par[2] # Estimate for r
  beta <- estimates$par[3] # Estimate for beta
  K <- -r / beta
  LL <- estimates$value # Log-likelihood
  N <- exp(X) # Population size on the arithmetic scale
  pred <- c(NA, X + r - (sigma_e2 / 2) - (sigma_d2 / (2 * N)) + beta * N) # Predicted values
  res <- c(X, NA) - pred # Residuals
  dem_comp <- sigma_d2 / N # Demographic component

  # Adjust indexing
  pred <- pred[-length(pred)]
  res <- res[-1]

  ic <- calculate_criteria(X, estimates) # Calculate information criteria

  return(list(estimates = estimates, sigma_e2 = sigma_e2, r = r, beta = beta, K = K,
              X = X, N = N, sigma_d2 = sigma_d2,
              pred = pred, res = res, dem_comp = dem_comp, boot = boot, LL = LL,
              AIC = ic$AIC, AICc = ic$AICc, BIC = ic$BIC))

}

#--------------------------#
# NOTE: UNDER CONSTRUCTION #
#--------------------------#

# Stationary distribution
stat_dist_logis <- function(mod, init_X, n_rep, n_time) {
  # mod: output from `run_population_model()`
  # init_X: initial value for X.
  # n_rep: number of replicates
  # n_time: number of time steps to simulate for

  if(missing(init_X)) {

    init_X <- log(-mod$r / mod$beta) # Initial value of simulation at K

  }

  # Output from `estimate_pars()`
  est <- mod$estimates

  # Simulate data
  sim_data <- matrix(NA, nrow = n_rep, ncol = n_time)

  pb_rep <- progress::progress_bar$new(total = n_rep,
                                       format = "Simulating data [:bar] :percent eta: :eta")

  for(i in 1:n_rep) {

    pb_rep$tick() # Update progress bar

    log_sim <- simulate_ts(X = mod$X, init_X = init_X, sigma_d2 = mod$sigma_d2,
                           est = est, dd = "logistic", n_time = n_time) # Simulate time series

    sim_data[i,] <- exp(log_sim)

  }

  # Mean and variance over simulated time series (first 10% of TS excluded)
  mean_stat_dist <- apply(sim_data[, -seq_len(n_time / 10)], 1, mean)
  var_stat_dist <- apply(sim_data[, -seq_len(n_time / 10)], 1, var)

  return(list(sim = sim_data[, -seq_len(n_time / 10)],
              mean_stat_dist = mean_stat_dist,
              var_stat_dist = var_stat_dist))

}

# Quasi-stationary distribution
quasi_stat_dist_logis <- function(mod, from = 1, to, n, lower = 0, upper) {
  # mod: output from `run_population_model()`
  # from: lower bound population size (Default: 1)
  # to: upper bound population size (Default: 2K)
  # n: number of population sizes
  # lower: lower limit of integration in Green function (Default: 0)
  # upper: upper limit of integration in Green function (Default: K)

  Green_function <- function(x, lower, upper, r, K, sigma_e2, sigma_d2) {
    # x: population size (arithmetic scale)
    # lower: lower limit of integration in Green function
    # upper: upper limit of integration in Green function

    m <- function(n, r, K) {r * n-(r/K) * n^2} # Infinitesimal mean
    v <- function(n, sigma_e2, sigma_d2) {sigma_e2 * n^2 + sigma_d2 * n} # Infinitesimal variance

    sv <- function(x) {m(x, r, K) / v(x, sigma_e2, sigma_d2)}
    mx <- function(x) {1/(v(x, sigma_e2 = sigma_e2, sigma_d2 = sigma_d2) * s(x))}

    s <- Vectorize(function(y) {exp(-2 * integrate(sv, lower = lower, upper = y, subdivisions = 100)$value)})
    S <- function(x) {integrate(s, lower = lower, upper = x)$value}

    ifelse(x<=upper, 2*mx(x)*S(x), 2*mx(x)*S(upper))

  }

  # Retrieve estimates from `run_population_model()`
  K <- -mod$r / mod$beta
  r <- mod$r
  sigma_e2 <- mod$sigma_e2
  sigma_d2 <- mod$sigma_d2

  # Vector of population sizes
  if(missing(to)) {

    to <- K * 2

  }

  vec <- seq(from = from, to = to, length.out = n)

  # Vector of probabilities
  res <- NA

  lower <- ifelse(sigma_d2 <= 0.01, 1, 0)

  if(missing(upper)) {

    upper <- K

  }

  for(i in 1:length(vec)) {

    res[i] <- Green_function(vec[i], lower, upper, r, K, sigma_e2, sigma_d2)

  }

  res <- res / (sum(res)) / diff(vec[1:2])

  tibble::tibble(
    N = vec,
    Pr = res
  )

}

# Test run ####
library(tidyverse)
#extrafont::loadfonts(device = "win")

# Source "original" script
source(here::here("R", "Thetalogistic_BlueTit_VL.R"))

# Simulate data
data <- simulate_data(log(10), sigma_e2 = 0.02, sigma_d2 = 0.3, r = 0.6, beta = -0.005, tmax = 50, seed = 18)

# Run logistic models
logis <- run_population_model(data$X, data$sigma_d2, n_boot = 1000)
logis2 <- logismodfit(x = exp(data$X), sd = data$sigma_d2, nboot = 1000) # From Thetalogistic_BlueTit_VL.R

# Compare my code and Vidar's code
mod_outputs <- tibble::tibble(
  val = c(logis$boot$boot[, "r"], logis$boot$boot[, "sigma_e2"], logis$boot$boot[, "beta"], logis$boot$boot[, "K"],
          logis2$boot[, "r"] + (logis2$boot[, "sig2"] * 0.5), logis2$boot[, "sig2"],
          -logis2$boot[, "r"] / logis2$boot[, "K"], logis2$boot[, "K"]),
  par = rep(rep(c("r", "sigma_e2", "beta", "K"), each = 1000), 2),
  mod = rep(c("New", "Original"), each = 4000)
)

pl <- purrr::map(.x = c("r", "sigma_e2", "beta", "K"),
                 .f = ~{

                   # Filter data for selected parameter
                   plot_data <- mod_outputs %>%
                     dplyr::filter(par == .x)

                   # True value for selected parameter
                   true <- data[[.x]]

                   if(.x == "K") {

                     true <- -data$r / data$beta

                   }

                   # Plot
                   ggplot(plot_data, aes(x = val, y = mod, fill = mod, color = mod)) +
                     geom_vline(xintercept = true, linetype = "dashed") +
                     ggridges::geom_density_ridges(scale = 0.85, alpha = 0.5) +
                     theme_classic(base_family = "Signika") +
                     labs(x = .x, y = "Model") +
                     scale_color_manual(values = c("#2c91a0", "#de6e4b")) +
                     scale_fill_manual(values = c("#2c91a0", "#de6e4b")) +
                     theme(axis.text = element_text(size = 11, color = "black"),
                           axis.title = element_text(size = 12, color = "black"),
                           plot.margin = margin(l = 5, b = 5, t = 10, r = 15),
                           legend.position = "none")

           })

cowplot::plot_grid(plotlist = pl) %>%
  cowplot::save_plot(filename = here::here("inst", "images", "Logis-outputs.pdf"),
                     device = cairo_pdf,
                     plot = .,
                     nrow = 2,
                     ncol = 2,
                     base_asp = 1.4)


#--------------------------#
# NOTE: UNDER CONSTRUCTION #
#--------------------------#

# Stationary distribution via simulation
stat_dist <- stat_dist_logis(logis, n_rep = 100, n_time = 10000)
qs <- quasi_stat_dist_logis(logis, n = 100)

# Compare my code and original code
# Input needed for original version of the stationary distribution
x <- list(s = logis$r - logis$sigma_e2 * 0.5,
          sig2 = logis$sigma_e2,
          sigd = logis$sigma_d2,
          K = -logis$r / logis$beta)

stasj.logismod(x) # From Thetalogistic_BlueTit_VL.R

hist(stat_dist$sim, prob = TRUE, ylim = c(0, max(qs$Pr) * 1.1), xlim = c(min(qs$N), 252),
     main = "Title", nclass = 30)
lines(qs$N, qs$Pr, col = "red")
