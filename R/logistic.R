# Simulate data
simulate_data <- function(init_X, sigma_e2, sigma_d2, r, beta, tmax, nn = 0, plot = TRUE, seed = NULL) {
  # init_X: initial value for X
  # sigma_e2: desired value for sigma_e2
  # sigma_d2: desired value for sigma_d2
  # r: desired value for r
  # beta: desired value for beta
  # tmax: time series length
  # nn: number of NAs
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

  if(nn > 0) {

    na_years <- sample(2:(tmax-1), nn)

    X[na_years] <- NA

  }

  # Plot time series
  plot(1:tmax, X, type = "l", xlab = "Time")

  return(list(X = X, sigma_e2 = sigma_e2, sigma_d2 = sigma_d2, r = r, beta = beta, tmax = tmax))

}

# Run logistic or density-independent population model
run_population_model <- function(X, sigma_d2, dd = "logistic",
                                 n_boot = NULL, n_time = length(X), hessian = FALSE) {
  # X: observed log population sizes
  # sigma_d2: demographic variance
  # dd: form of density dependence:
  # --- "logistic": logistic form of density dependence
  # --- "independent": density-independent
  # n_boot: number of replicates for bootstrap
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

    } else if(dd == "independent") {

      # Starting values
      start <- c(sigma_e2_hat, r_hat)

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
  simulate_ts <- function(X, sigma_d2, est, dd = "logistic", init_X = NULL, n_rep, n_time) {
    # X: observed log population sizes
    # sigma_d2: demographic variance
    # est: output of `estimate_pars()`
    # dd: form of density dependence:
    # --- "logistic": logistic form of density dependence
    # --- "independent": density-independent
    # init_X: initial value for X. By default, it takes the first value in X.
    #         Can be set to any value (e.g. to K when calculating the stationary distribution)
    # n_rep = number of replications
    # n_time: number of time steps to simulate for in each replicate

    # Retrieve parameter estimates
    sigma_e2 <- exp(est$par[1])
    r <- est$par[2]

    if(dd == "logistic") {

      beta <- est$par[3]

    } else if(dd == "independent") {

      beta <- 0

    }

    first <- min(which(!is.na(X))) # Get first non-NA observation
    X_sim <- matrix(NA, nrow = n_rep, ncol = n_time)

    # Set initial value of X_sim and XX
    if(is.null(init_X)) {

      X_sim[,first] <- XX <- X[first] # First value in observations (X)

    } else if(!is.null(init_X))  {

      X_sim[,first] <- XX <- init_X

    }

    # Progess bar
    #pb_sim <- progress::progress_bar$new(total = n_time, clear = FALSE,
    #                                     format = "Simulating [:bar] :percent eta: :eta")

    # Random variation
    # NB: Calculate random variation for all replicates and time steps outside loop for shorter run time
    eps <- matrix(rnorm(n_rep * n_time, mean = 0, sd = 1), nrow = n_rep, ncol = n_time)

    try <- 0 # Try counter

    # Simulate
    while(any(is.na(X_sim)) | any(X_sim == 0)){

      # If try counter reaches 1000, stop simulation
      try <- try + 1
      if (try > 1000) stop("Unable to simulate time series.")

      for (t in first:(n_time-1)){

        #pb_sim$tick()

        E_X <- XX + r - (sigma_e2 / 2) - (sigma_d2 / (2 * exp(XX))) + beta * exp(XX) # Expectation

        Var_X <- sigma_e2 + (sigma_d2 / exp(XX)) # Variance

        XX <- E_X + sqrt(Var_X) * eps[,t] # X next time step

        XX <- ifelse(XX <= 1, NA, XX) # If time series is extinct (log(1) == 0), set to NA

        X_sim[,t+1] <- XX # Store calculations for t+1 in X_sim matrix

      }

    }

    X_sim[,is.na(X)] <- NA # NA in observations -> NA in simulations

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
    # n_boot: number of bootstrap replicates
    # n_time: number of time steps to simulate data for
    # hessian: return Hessian matrix (TRUE or FALSE)

    # Simulate data
    sim_data <- matrix(NA, nrow = n_boot, ncol = n_time)
    sim_data <- simulate_ts(X = X, sigma_d2 = sigma_d2, est = est, dd = dd, n_rep = n_boot, n_time = n_time)

    # Bootstrap
    boot <- matrix(NA, nrow = n_boot, ncol = 4)
    colnames(boot) <- par_names

    # Progess bar
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
  estimates <- estimate_pars(X = X, sigma_d2 = sigma_d2, dd = dd,
                             start = set_starting_values(X, dd = dd),
                             hessian = hessian)

  # Bootstrap uncertainty
  if(!is.null(n_boot)) {

    boot <- bootstrap(X, sigma_d2, estimates, dd = dd, n_boot = n_boot, n_time = n_time, hessian = hessian)

  }

  #----------------#
  # RETURN OUTPUTS #
  #----------------#

  sigma_e2 <- exp(estimates$par[1]) # Estimate for sigma_e2

  r <- estimates$par[2] # Estimate for r

  if(dd == "logistic") {

    beta <- estimates$par[3] # Estimate for beta
    K <- -r / beta

  } else if (dd == "independent") {

    beta <- 0
    K <- NULL

  }

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
              X = X, N = N, sigma_d2 = sigma_d2, dd = dd,
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

  # Simulate time series
  sim_data <- matrix(NA, nrow = n_rep, ncol = n_time)
  log_sim <- simulate_ts(X = mod$X, init_X = init_X, sigma_d2 = mod$sigma_d2,
                         est = est, dd = "logistic", n_rep = n_rep, n_time = n_time)
  sim_data <- exp(log_sim)

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
data <- simulate_data(log(10), sigma_e2 = 0.07, sigma_d2 = 0.4, r = 0.6, beta = -0.01, tmax = 50)

# Run logistic models
logis <- run_population_model(data$X, data$sigma_d2, n_boot = 10000)
logis2 <- logismodfit(x = exp(data$X), sd = data$sigma_d2, nboot = 10000) # From Thetalogistic_BlueTit_VL.R

# Compare my code and original code
mod_outputs <- tibble::tibble(
  val = c(logis$boot$boot[, "r"], logis$boot$boot[, "sigma_e2"], logis$boot$boot[, "beta"], logis$boot$boot[, "K"],
          logis2$boot[, "r"] + (logis2$boot[, "sig2"] * 0.5), logis2$boot[, "sig2"],
          -logis2$boot[, "r"] / logis2$boot[, "K"], logis2$boot[, "K"]),
  par = rep(rep(c("r", "sigma_e2", "beta", "K"), each = 10000), 2),
  mod = rep(c("New", "Original"), each = 40000)
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
stat_dist <- stat_dist_logis(logis, n_rep = 1000, n_time = 10000)
qs <- quasi_stat_dist_logis(logis, n = 100)

# Compare my code and original code
# Input needed for original version of the stationary distribution
x <- list(s = logis$r - logis$sigma_e2 * 0.5,
          sig2 = logis$sigma_e2,
          sigd = logis$sigma_d2,
          K = logis$K)

stasj.logismod(x) # From Thetalogistic_BlueTit_VL.R

hist(stat_dist$sim, prob = TRUE, ylim = c(0, max(qs$Pr) * 1.1), xlim = c(min(qs$N), 252),
     main = "New", nclass = 30)
lines(qs$N, qs$Pr, col = "red")


# Multiple test runs ####
rep <- 15
purrr::walk(.x = 1:rep,
            .f = ~{

              cat(crayon::yellow(paste("Model", .x, sep = " ")))

              data <- simulate_data(log(10), sigma_e2 = 0.05, sigma_d2 = 0.4, r = 0.6, beta = -0.01, tmax = 50)

              # Run logistic models
              logis <- run_population_model(data$X, data$sigma_d2, n_boot = 2000)
              logis2 <- logismodfit(x = exp(data$X), sd = data$sigma_d2, nboot = 2000) # From Thetalogistic_BlueTit_VL.R

              # Compare my code and Vidar's code
              mod_outputs <- tibble::tibble(
                val = c(logis$boot$boot[, "r"], logis$boot$boot[, "sigma_e2"], logis$boot$boot[, "beta"], logis$boot$boot[, "K"],
                        logis2$boot[, "r"] + (logis2$boot[, "sig2"] * 0.5), logis2$boot[, "sig2"],
                        -(logis2$boot[, "r"] + (logis2$boot[, "sig2"] * 0.5)) / logis2$boot[, "K"], logis2$boot[, "K"]),
                par = rep(rep(c("r", "sigma_e2", "beta", "K"), each = 2000), 2),
                mod = rep(c("New", "Original"), each = 8000)
              )

              # LL estimates
              ll_outputs <- tibble::tibble(
                val = c(logis$r, logis$sigma_e2, logis$beta, logis$K,
                        logis2$s + logis2$sig2 * 0.5, logis2$sig2, -(logis2$s + (logis2$sig2 * 0.5)) / logis2$K, logis2$K),
                par = rep(c("r", "sigma_e2", "beta", "K"), 2),
                mod = rep(c("New", "Original"), each = 4),
                id = rep(1:2, each = 4)
              )

              pl <- purrr::map(.x = c("r", "sigma_e2", "beta", "K"),
                               .f = ~{

                                 # Filter data for selected parameter
                                 plot_data <- mod_outputs %>%
                                   dplyr::filter(par == .x)

                                 plot_ll <- ll_outputs %>%
                                   dplyr::filter(par == .x)

                                 # True value for selected parameter
                                 true <- data[[.x]]

                                 if(.x == "K") {

                                   true <- -data$r / data$beta

                                 }

                                 if(.x == "K") {

                                   xlim <- c(20, 100)

                                 } else if(.x == "r") {

                                   xlim <- c(0.2, 1.5)

                                 } else if(.x == "beta") {

                                   xlim <- c(-0.03, 0)

                                 } else if(.x == "sigma_e2") {

                                   xlim <- c(0, 0.2)

                                 }

                                 # Plot
                                 ggplot(plot_data, aes(x = val, y = mod, fill = mod, color = mod)) +
                                   geom_segment(data = plot_ll, mapping = aes(x = val, xend = val, y = mod, yend = id + 0.8,
                                                                              color = mod)) +
                                   geom_vline(xintercept = true, linetype = "dashed") +
                                   ggridges::geom_density_ridges(scale = 0.8, alpha = 0.5) +
                                   theme_classic(base_family = "Signika") +
                                   labs(x = .x, y = "Model") +
                                   coord_cartesian(xlim = xlim) +
                                   scale_color_manual(values = c("#2c91a0", "#de6e4b")) +
                                   scale_fill_manual(values = c("#2c91a0", "#de6e4b")) +
                                   theme(axis.text = element_text(size = 11, color = "black"),
                                         axis.title = element_text(size = 12, color = "black"),
                                         plot.margin = margin(l = 5, b = 5, t = 10, r = 15),
                                         legend.position = "none")

                               })

              cowplot::plot_grid(plotlist = pl) %>%
                cowplot::save_plot(filename = here::here("inst", "images", paste0("Logis-outputs", "-", .x, ".pdf")),
                                   device = cairo_pdf,
                                   plot = .,
                                   nrow = 2,
                                   ncol = 2,
                                   base_asp = 1.4)


            })

pdftools::pdf_combine(here::here("inst", "images", paste0("Logis-outputs", "-", seq_len(rep), ".pdf")),
                      output = here::here("inst", "images", paste0("Logis-outputs-all", ".pdf")))

unlink(here::here("inst", "images", paste0("Logis-outputs", "-", seq_len(rep), ".pdf")))


# NA test run ####
data1 <- simulate_data(log(10), sigma_e2 = 0.03, sigma_d2 = 0.5, r = 0.6, beta = -0.01, tmax = 50, nn = 1, seed = 112)
data2 <- simulate_data(log(10), sigma_e2 = 0.03, sigma_d2 = 0.5, r = 0.6, beta = -0.01, tmax = 50, nn = 2, seed = 112)
data5 <- simulate_data(log(10), sigma_e2 = 0.03, sigma_d2 = 0.5, r = 0.6, beta = -0.01, tmax = 50, nn = 5, seed = 112)
data10 <- simulate_data(log(10), sigma_e2 = 0.03, sigma_d2 = 0.5, r = 0.6, beta = -0.01, tmax = 50, nn = 10, seed = 112)
data15 <- simulate_data(log(10), sigma_e2 = 0.03, sigma_d2 = 0.5, r = 0.6, beta = -0.01, tmax = 50, nn = 15, seed = 112)

logis1 <- run_population_model(data1$X, data1$sigma_d2, n_boot = 1000)
logis2 <- run_population_model(data2$X, data2$sigma_d2, n_boot = 1000)
logis5 <- run_population_model(data5$X, data5$sigma_d2, n_boot = 1000)
logis10 <- run_population_model(data10$X, data10$sigma_d2, n_boot = 1000)
logis15 <- run_population_model(data15$X, data15$sigma_d2, n_boot = 1000)

mod_outputs_na <- tibble::tibble(
  val = c(logis1$boot$boot[, "r"], logis1$boot$boot[, "sigma_e2"], logis1$boot$boot[, "beta"], logis1$boot$boot[, "K"],
          logis2$boot$boot[, "r"], logis2$boot$boot[, "sigma_e2"], logis2$boot$boot[, "beta"], logis2$boot$boot[, "K"],
          logis5$boot$boot[, "r"], logis5$boot$boot[, "sigma_e2"], logis5$boot$boot[, "beta"], logis5$boot$boot[, "K"],
          logis10$boot$boot[, "r"], logis10$boot$boot[, "sigma_e2"],
          logis10$boot$boot[, "beta"], logis10$boot$boot[, "K"],
          logis15$boot$boot[, "r"], logis15$boot$boot[, "sigma_e2"],
          logis15$boot$boot[, "beta"], logis15$boot$boot[, "K"]),
  par = rep(rep(c("r", "sigma_e2", "beta", "K"), each = 1000), 5),
  mod = rep(c("1", "2", "5", "10", "15"), each = 4000)
) %>%
  dplyr::mutate(mod = forcats::fct_relevel(mod, c("1", "2", "5", "10", "15")))

pl_na <- purrr::map(.x = c("r", "sigma_e2", "beta", "K"),
                    .f = ~{

                      # Filter data for selected parameter
                      plot_data <- mod_outputs_na %>%
                        dplyr::filter(par == .x)

                      # True value for selected parameter
                      true <- data1[[.x]]

                      if(.x == "K") {

                        true <- -data1$r / data1$beta

                      }

                      # Plot
                      ggplot(plot_data, aes(x = val, y = mod, fill = mod, color = mod)) +
                        geom_vline(xintercept = true, linetype = "dashed") +
                        ggridges::geom_density_ridges(scale = 0.85, alpha = 0.5) +
                        theme_classic(base_family = "Titillium Web") +
                        labs(x = .x, y = "Number of NAs") +
                        scale_color_manual(values = c("#053c5e", "#559cad", "#856084", "#f39c6b", "#562423")) +
                        scale_fill_manual(values = c("#053c5e", "#559cad", "#856084", "#f39c6b", "#562423")) +
                        theme(axis.text = element_text(size = 11, color = "black"),
                              axis.title = element_text(size = 12, color = "black"),
                              plot.margin = margin(l = 5, b = 5, t = 10, r = 15),
                              legend.position = "none")

                    })

cowplot::plot_grid(plotlist = pl_na) %>%
  cowplot::save_plot(filename = here::here("inst", "images", "Logis-outputs-NA.pdf"),
                     device = cairo_pdf,
                     plot = .,
                     nrow = 2,
                     ncol = 2,
                     base_asp = 1.4)
