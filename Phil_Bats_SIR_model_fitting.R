### Example code #1 ####
########################

# Required libraries
library(deSolve)
library(ggplot2)

# SIR model differential equations
sir_model <- function(time, state, parameters) {
  S <- state[1]
  I <- state[2]
  R <- state[3]
  
  beta <- parameters[1]
  gamma <- parameters[2]
  
  dS <- -beta * S * I
  dI <- beta * S * I - gamma * I
  dR <- gamma * I
  
  list(c(dS, dI, dR))
}

# Solve SIR model using the lsoda solver (adaptive step size)
solve_sir <- function(time, beta, gamma, S0, I0, R0) {
  initial_state <- c(S = S0, I = I0, R = R0)
  parameters <- c(beta, gamma)
  
  out <- try(ode(y = initial_state, times = time, func = sir_model, parms = parameters, 
                 method = "lsoda", atol = 1e-8, rtol = 1e-8))  # Tighter tolerances
  
  # If the solver fails, return NA
  if (inherits(out, "try-error") || any(is.na(out))) {
    return(NA)
  } else {
    return(out)
  }
}

# Objective function to minimize: Sum of Squared Errors (SSE) between observed and model-predicted R values
objective_function <- function(params, time, R_data) {
  beta <- params[1]
  gamma <- params[2]
  S0 <- params[3]
  I0 <- params[4]
  
  # Solve the SIR model with given parameters
  sir_output <- solve_sir(time, beta, gamma, S0, I0, R0 = 0)
  
  # If the solver fails or returns NA, return a large error value
  if (is.na(sir_output[1]) || any(is.na(sir_output))) {
    return(1e10)  # Large penalty for invalid solutions
  }
  
  R_predicted <- sir_output[, "R"]
  
  # Calculate SSE (sum of squared errors)
  SSE <- sum((R_data - R_predicted)^2)
  
  return(SSE)
}

# Simulated time series over 100 days
time_data <- seq(0, 100, by = 1)  # Simulate for 100 days

# Simulated observed recovered data (plausible epidemic dynamics)
R_data <- c(0, 0.002, 0.005, 0.012, 0.03, 0.08, 0.15, 0.28, 0.42, 0.55, 0.68, 
            0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.93, 0.94, 0.95, 0.95, 0.96, 
            0.97, 0.97, 0.97, 0.98, 0.98, 0.98, 0.99, 0.99, 0.99, 1.0)

# Ensure length of R_data matches time_data
R_data <- c(R_data, rep(1, length(time_data) - length(R_data)))

# Initial guesses for parameters
initial_params <- c(beta = 0.4, gamma = 0.1, S0 = 0.99, I0 = 0.01)

# Optimize the parameters using optim
fit <- optim(par = initial_params, fn = objective_function, time = time_data, R_data = R_data,
             method = "L-BFGS-B", lower = c(0, 0, 0, 0), upper = c(2, 1, 1, 1))

# Extract fitted parameters
beta_est <- fit$par[1]
gamma_est <- fit$par[2]
S0_est <- fit$par[3]
I0_est <- fit$par[4]

# Solve SIR model using the fitted parameters
SIR_fitted <- solve_sir(time_data, beta_est, gamma_est, S0_est, I0_est, R0 = 0)

# Check for NA values in the solver output
if (any(is.na(SIR_fitted))) {
  stop("SIR model solver failed with the fitted parameters.")
}

# Extract each state
S_fitted <- SIR_fitted[, "S.S0"]
I_fitted <- SIR_fitted[, "I.I0"]
R_fitted <- SIR_fitted[, "R"]

# Combine data into a single data frame for plotting
data <- data.frame(
  time = time_data,
  S = S_fitted,
  I = I_fitted,
  R = R_fitted,
  R_observed = R_data  # For comparison of observed vs fitted R
)

# Plot all three states S, I, and R
ggplot(data, aes(x = time)) +
  geom_line(aes(y = S, color = "S"), size = 1) +
  geom_line(aes(y = I, color = "I"), size = 1) +
  geom_line(aes(y = R, color = "R (fitted)"), size = 1) +
  geom_point(aes(y = R_observed, color = "R (observed)"), size = 3, shape = 16) +
  labs(title = "SIR Model: Susceptible, Infected, and Recovered States (100 Days)",
       x = "Time (days)", y = "Proportion of Population") +
  scale_color_manual(name = "State",
                     values = c("S" = "blue", "I" = "red", "R (fitted)" = "green", "R (observed)" = "black")) +
  theme_minimal()

### Example code #2 ####
########################

# Required libraries
library(deSolve)
library(ggplot2)

# SIR model differential equations with births and deaths
sir_model_with_births_deaths <- function(time, state, parameters) {
  S <- state[1]
  I <- state[2]
  R <- state[3]
  
  beta <- parameters[1]
  gamma <- parameters[2]
  lambda <- parameters[3]  # Birth rate
  mu <- parameters[4]      # Death rate
  
  dS <- lambda - beta * S * I - mu * S
  dI <- beta * S * I - gamma * I - mu * I
  dR <- gamma * I - mu * R
  
  list(c(dS, dI, dR))
}

# Solve SIR model with births and deaths
solve_sir_with_births_deaths <- function(time, beta, gamma, lambda, mu, S0, I0, R0) {
  initial_state <- c(S = S0, I = I0, R = R0)
  parameters <- c(beta, gamma, lambda, mu)
  
  out <- try(ode(y = initial_state, times = time, func = sir_model_with_births_deaths, parms = parameters, 
                 method = "lsoda", atol = 1e-8, rtol = 1e-8))  # Tighter tolerances
  
  # If the solver fails, return NA
  if (inherits(out, "try-error") || any(is.na(out))) {
    return(NA)
  } else {
    return(out)
  }
}

# Objective function to minimize for model with births and deaths
objective_function_with_births_deaths <- function(params, time, R_data) {
  beta <- params[1]
  gamma <- params[2]
  lambda <- params[3]
  mu <- params[4]
  S0 <- params[5]
  I0 <- params[6]
  
  # Solve the SIR model with given parameters
  sir_output <- solve_sir_with_births_deaths(time, beta, gamma, lambda, mu, S0, I0, R0 = 0)
  
  # If the solver fails or returns NA, return a large error value
  if (is.na(sir_output[1]) || any(is.na(sir_output))) {
    return(1e10)  # Large penalty for invalid solutions
  }
  
  R_predicted <- sir_output[, "R"]
  
  # Calculate SSE (sum of squared errors)
  SSE <- sum((R_data - R_predicted)^2)
  
  return(SSE)
}

# Simulated time series over 100 days
time_data <- seq(0, 100, by = 1)  # Simulate for 100 days

# Simulated observed recovered data (plausible epidemic dynamics)
R_data <- c(0, 0.002, 0.005, 0.012, 0.03, 0.08, 0.15, 0.28, 0.42, 0.55, 0.68, 
            0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.93, 0.94, 0.95, 0.95, 0.96, 
            0.97, 0.97, 0.97, 0.98, 0.98, 0.98, 0.99, 0.99, 0.99, 1.0)

# Ensure length of R_data matches time_data
R_data <- c(R_data, rep(1, length(time_data) - length(R_data)))

# Initial guesses for parameters for the basic SIR model
initial_params_sir <- c(beta = 0.4, gamma = 0.1, S0 = 0.99, I0 = 0.01)

# Optimize the parameters for the basic SIR model
fit_sir <- optim(par = initial_params_sir, fn = objective_function, time = time_data, R_data = R_data,
                 method = "L-BFGS-B", lower = c(0, 0, 0, 0), upper = c(2, 1, 1, 1))

# Extract fitted parameters for the basic SIR model
beta_est_sir <- fit_sir$par[1]
gamma_est_sir <- fit_sir$par[2]
S0_est_sir <- fit_sir$par[3]
I0_est_sir <- fit_sir$par[4]

# Solve SIR model using the fitted parameters
SIR_fitted <- solve_sir(time_data, beta_est_sir, gamma_est_sir, S0_est_sir, I0_est_sir, R0 = 0)

# Check for NA values in the solver output for basic SIR model
if (any(is.na(SIR_fitted))) {
  stop("SIR model solver failed with the fitted parameters.")
}

# Extract each state for basic SIR model
S_fitted_sir <- SIR_fitted[, "S.S0"]
I_fitted_sir <- SIR_fitted[, "I.I0"]
R_fitted_sir <- SIR_fitted[, "R"]

# Initial guesses for parameters for the extended SIR model with births and deaths
initial_params_sir_bdb <- c(beta = 0.4, gamma = 0.1, lambda = 0.02, mu = 0.01, S0 = 0.99, I0 = 0.01)

# Optimize the parameters for the extended SIR model with births and deaths
fit_sir_bdb <- optim(par = initial_params_sir_bdb, fn = objective_function_with_births_deaths, 
                     time = time_data, R_data = R_data, method = "L-BFGS-B", 
                     lower = c(0, 0, 0, 0, 0, 0), upper = c(2, 1, 1, 1, 1, 1))

# Extract fitted parameters for the extended SIR model
beta_est_bdb <- fit_sir_bdb$par[1]
gamma_est_bdb <- fit_sir_bdb$par[2]
lambda_est_bdb <- fit_sir_bdb$par[3]
mu_est_bdb <- fit_sir_bdb$par[4]
S0_est_bdb <- fit_sir_bdb$par[5]
I0_est_bdb <- fit_sir_bdb$par[6]

# Solve the extended SIR model using the fitted parameters
SIR_fitted_bdb <- solve_sir_with_births_deaths(time_data, beta_est_bdb, gamma_est_bdb, lambda_est_bdb, mu_est_bdb, S0_est_bdb, I0_est_bdb, R0 = 0)

# Check for NA values in the solver output for extended SIR model
if (any(is.na(SIR_fitted_bdb))) {
  stop("Extended SIR model solver failed with the fitted parameters.")
}

# Extract each state for extended SIR model
S_fitted_bdb <- SIR_fitted_bdb[, "S.S0"]
I_fitted_bdb <- SIR_fitted_bdb[, "I.I0"]
R_fitted_bdb <- SIR_fitted_bdb[, "R"]

# Calculate AIC for both models
n <- length(R_data)
aic_sir <- n * log(sum((R_data - R_fitted_sir)^2) / n) + 2 * 4  # 4 parameters for basic SIR
aic_sir_bdb <- n * log(sum((R_data - R_fitted_bdb)^2) / n) + 2 * 6  # 6 parameters for extended SIR

# Print AIC values
cat("AIC for basic SIR model: ", aic_sir, "\n")
cat("AIC for extended SIR model with births and deaths: ", aic_sir_bdb, "\n")

# Combine data into a single data frame for plotting
data <- data.frame(
  time = time_data,
  S_S0 = S_fitted_sir,
  I_I0 = I_fitted_sir,
  R_fitted = R_fitted_sir,
  R_observed = R_data,
  S_S0_bdb = S_fitted_bdb,
  I_I0_bdb = I_fitted_bdb,
  R_fitted_bdb = R_fitted_bdb
)

# Plot all three states S, I, and R for both models
ggplot(data, aes(x = time)) +
  geom_line(aes(y = S_S0, color = "S (basic)"), size = 1, linetype = "dashed") +
  geom_line(aes(y = I_I0, color = "I (basic)"), size = 1, linetype = "dashed") +
  geom_line(aes(y = R_fitted, color = "R (fitted basic)"), size = 1, linetype = "dashed") +
  geom_point(aes(y = R_observed, color = "R (observed)"), size = 3, shape = 16) +
  geom_line(aes(y = S_S0_bdb, color = "S (extended)"), size = 1) +
  geom_line(aes(y = I_I0_bdb, color = "I (extended)"), size = 1) +
  geom_line(aes(y = R_fitted_bdb, color = "R (fitted extended)"), size = 1) +
  labs(title = "SIR Model Comparison: Basic vs. Extended (with Births and Deaths)",
       x = "Time (days)", y = "Proportion of Population") +
  scale_color_manual(name = "State",
                     values = c("S (basic)" = "blue", "I (basic)" = "red", "R (fitted basic)" = "green", 
                                "R (observed)" = "black", "S (extended)" = "cyan", 
                                "I (extended)" = "orange", "R (fitted extended)" = "purple")) +
  theme_minimal()

### Example code #3 ####
########################

# SIR model differential equations with time-varying births and deaths
sir_model_with_time_varying_births <- function(time, state, parameters) {
  S <- state[1]
  I <- state[2]
  R <- state[3]
  
  beta <- parameters[1]
  gamma <- parameters[2]
  lambda0 <- parameters[3]  # Baseline birth rate
  lambda1 <- parameters[4]  # Amplitude of the birth rate variation
  omega <- parameters[5]    # Frequency of the birth rate variation
  phi <- parameters[6]      # Phase shift of the birth rate variation
  mu <- parameters[7]       # Death rate
  
  # Time-varying birth rate modeled with a sinusoidal function
  lambda_t <- lambda0 + lambda1 * sin(omega * time + phi)
  
  # SIR model equations with time-varying birth rate and deaths
  dS <- lambda_t - beta * S * I - mu * S
  dI <- beta * S * I - gamma * I - mu * I
  dR <- gamma * I - mu * R
  
  list(c(dS, dI, dR))
}

# Solve SIR model with time-varying births and deaths
solve_sir_with_time_varying_births <- function(time, beta, gamma, lambda0, lambda1, omega, phi, mu, S0, I0, R0) {
  initial_state <- c(S = S0, I = I0, R = R0)
  parameters <- c(beta, gamma, lambda0, lambda1, omega, phi, mu)
  
  out <- try(ode(y = initial_state, times = time, func = sir_model_with_time_varying_births, parms = parameters, 
                 method = "lsoda", atol = 1e-8, rtol = 1e-8))  # Tighter tolerances
  
  # If the solver fails, return NA
  if (inherits(out, "try-error") || any(is.na(out))) {
    return(NA)
  } else {
    return(out)
  }
}

# Objective function to minimize for model with time-varying births and deaths
objective_function_with_time_varying_births <- function(params, time, R_data) {
  beta <- params[1]
  gamma <- params[2]
  lambda0 <- params[3]
  lambda1 <- params[4]
  omega <- params[5]
  phi <- params[6]
  mu <- params[7]
  S0 <- params[8]
  I0 <- params[9]
  
  # Solve the SIR model with given parameters
  sir_output <- solve_sir_with_time_varying_births(time, beta, gamma, lambda0, lambda1, omega, phi, mu, S0, I0, R0 = 0)
  
  # If the solver fails or returns NA, return a large error value
  if (is.na(sir_output[1]) || any(is.na(sir_output))) {
    return(1e10)  # Large penalty for invalid solutions
  }
  
  R_predicted <- sir_output[, "R"]
  
  # Calculate SSE (sum of squared errors)
  SSE <- sum((R_data - R_predicted)^2)
  
  return(SSE)
}

# Simulated time series over 100 days
time_data <- seq(0, 100, by = 1)  # Simulate for 100 days

# Simulated observed recovered data (plausible epidemic dynamics)
R_data <- c(0, 0.002, 0.005, 0.012, 0.03, 0.08, 0.15, 0.28, 0.42, 0.55, 0.68, 
            0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.93, 0.94, 0.95, 0.95, 0.96, 
            0.97, 0.97, 0.97, 0.98, 0.98, 0.98, 0.99, 0.99, 0.99, 1.0)

# Ensure length of R_data matches time_data
R_data <- c(R_data, rep(1, length(time_data) - length(R_data)))

# Initial guesses for parameters
initial_params_sir_time_varying <- c(beta = 0.4, gamma = 0.1, lambda0 = 0.02, lambda1 = 0.01, 
                                     omega = 2 * pi / 365, phi = 0, mu = 0.01, S0 = 0.99, I0 = 0.01)

# Optimize the parameters for the time-varying births model
fit_sir_time_varying <- optim(par = initial_params_sir_time_varying, fn = objective_function_with_time_varying_births, 
                              time = time_data, R_data = R_data, method = "L-BFGS-B", 
                              lower = c(0, 0, 0, 0, 0, -pi, 0, 0, 0), 
                              upper = c(2, 1, 1, 1, 2 * pi, pi, 1, 1, 1))

# Extract fitted parameters for the time-varying births model
beta_est_tv <- fit_sir_time_varying$par[1]
gamma_est_tv <- fit_sir_time_varying$par[2]
lambda0_est_tv <- fit_sir_time_varying$par[3]
lambda1_est_tv <- fit_sir_time_varying$par[4]
omega_est_tv <- fit_sir_time_varying$par[5]
phi_est_tv <- fit_sir_time_varying$par[6]
mu_est_tv <- fit_sir_time_varying$par[7]
S0_est_tv <- fit_sir_time_varying$par[8]
I0_est_tv <- fit_sir_time_varying$par[9]

# Solve the time-varying births model using the fitted parameters
SIR_fitted_tv <- solve_sir_with_time_varying_births(time_data, beta_est_tv, gamma_est_tv, lambda0_est_tv, lambda1_est_tv, omega_est_tv, phi_est_tv, mu_est_tv, S0_est_tv, I0_est_tv, R0 = 0)

# Check for NA values in the solver output for time-varying births model
if (any(is.na(SIR_fitted_tv))) {
  stop("Time-varying births model solver failed with the fitted parameters.")
}

# Extract each state for time-varying births model
S_fitted_tv <- SIR_fitted_tv[, "S.S0"]
I_fitted_tv <- SIR_fitted_tv[, "I.I0"]
R_fitted_tv <- SIR_fitted_tv[, "R"]

# Combine data into a single data frame for plotting
data_tv <- data.frame(
  time = time_data,
  S_fitted_tv = S_fitted_tv,
  I_fitted_tv = I_fitted_tv,
  R_fitted_tv = R_fitted_tv,
  R_observed = R_data
)

# Plot the results for the time-varying births model
ggplot(data_tv, aes(x = time)) +
  geom_line(aes(y = S_fitted_tv, color = "S (time-varying)"), size = 1) +
  geom_line(aes(y = I_fitted_tv, color = "I (time-varying)"), size = 1) +
  geom_line(aes(y = R_fitted_tv, color = "R (fitted time-varying)"), size = 1) +
  geom_point(aes(y = R_observed, color = "R (observed)"), size = 3, shape = 16) +
  labs(title = "SIR Model with Time-Varying Births",
       x = "Time (days)", y = "Proportion of Population") +
  scale_color_manual(name = "State",
                     values = c("S (time-varying)" = "blue", "I (time-varying)" = "red", 
                                "R (fitted time-varying)" = "green", "R (observed)" = "black")) +
  theme_minimal()

### Example code #3 ####
########################

# SIR model differential equations with time-varying births and deaths
sir_model_with_time_varying_births <- function(time, state, parameters) {
  S <- state[1]
  I <- state[2]
  R <- state[3]
  
  beta <- parameters[1]
  gamma <- parameters[2]
  lambda0 <- parameters[3]  # Baseline birth rate
  lambda1 <- parameters[4]  # Amplitude of the birth rate variation
  omega <- parameters[5]    # Frequency of the birth rate variation
  phi <- parameters[6]      # Phase shift of the birth rate variation
  mu <- parameters[7]       # Death rate
  
  # Time-varying birth rate modeled with a sinusoidal function
  lambda_t <- lambda0 + lambda1 * sin(omega * time + phi)
  
  # SIR model equations with time-varying birth rate and deaths
  dS <- lambda_t - beta * S * I - mu * S
  dI <- beta * S * I - gamma * I - mu * I
  dR <- gamma * I - mu * R
  
  list(c(dS, dI, dR))
}

# Solve SIR model with time-varying births and deaths
solve_sir_with_time_varying_births <- function(time, beta, gamma, lambda0, lambda1, omega, phi, mu, S0, I0, R0) {
  initial_state <- c(S = S0, I = I0, R = R0)
  parameters <- c(beta, gamma, lambda0, lambda1, omega, phi, mu)
  
  out <- try(ode(y = initial_state, times = time, func = sir_model_with_time_varying_births, parms = parameters, 
                 method = "lsoda", atol = 1e-8, rtol = 1e-8))  # Tighter tolerances
  
  # If the solver fails, return NA
  if (inherits(out, "try-error") || any(is.na(out))) {
    return(NA)
  } else {
    return(out)
  }
}

# Objective function to minimize for model with time-varying births and deaths
objective_function_with_time_varying_births <- function(params, time, R_data) {
  beta <- params[1]
  gamma <- params[2]
  lambda0 <- params[3]
  lambda1 <- params[4]
  omega <- params[5]
  phi <- params[6]
  mu <- params[7]
  S0 <- params[8]
  I0 <- params[9]
  
  # Solve the SIR model with given parameters
  sir_output <- solve_sir_with_time_varying_births(time, beta, gamma, lambda0, lambda1, omega, phi, mu, S0, I0, R0 = 0)
  
  # If the solver fails or returns NA, return a large error value
  if (is.na(sir_output[1]) || any(is.na(sir_output))) {
    return(1e10)  # Large penalty for invalid solutions
  }
  
  R_predicted <- sir_output[, "R"]
  
  # Calculate SSE (sum of squared errors)
  SSE <- sum((R_data - R_predicted)^2)
  
  return(SSE)
}

# Simulated time series over 100 days
time_data <- seq(0, 100, by = 1)  # Simulate for 100 days

# Simulated observed recovered data (plausible epidemic dynamics)
R_data <- c(0, 0.002, 0.005, 0.012, 0.03, 0.08, 0.15, 0.28, 0.42, 0.55, 0.68, 
            0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.93, 0.94, 0.95, 0.95, 0.96, 
            0.97, 0.97, 0.97, 0.98, 0.98, 0.98, 0.99, 0.99, 0.99, 1.0)

# Ensure length of R_data matches time_data
R_data <- c(R_data, rep(1, length(time_data) - length(R_data)))

# Initial guesses for parameters
initial_params_sir_time_varying <- c(beta = 0.4, gamma = 0.1, lambda0 = 0.02, lambda1 = 0.01, 
                                     omega = 2 * pi / 365, phi = 0, mu = 0.01, S0 = 0.99, I0 = 0.01)

# Optimize the parameters for the time-varying births model
fit_sir_time_varying <- optim(par = initial_params_sir_time_varying, fn = objective_function_with_time_varying_births, 
                              time = time_data, R_data = R_data, method = "L-BFGS-B", 
                              lower = c(0, 0, 0, 0, 0, -pi, 0, 0, 0), 
                              upper = c(2, 1, 1, 1, 2 * pi, pi, 1, 1, 1))

# Extract fitted parameters for the time-varying births model
beta_est_tv <- fit_sir_time_varying$par[1]
gamma_est_tv <- fit_sir_time_varying$par[2]
lambda0_est_tv <- fit_sir_time_varying$par[3]
lambda1_est_tv <- fit_sir_time_varying$par[4]
omega_est_tv <- fit_sir_time_varying$par[5]
phi_est_tv <- fit_sir_time_varying$par[6]
mu_est_tv <- fit_sir_time_varying$par[7]
S0_est_tv <- fit_sir_time_varying$par[8]
I0_est_tv <- fit_sir_time_varying$par[9]

# Solve the time-varying births model using the fitted parameters
SIR_fitted_tv <- solve_sir_with_time_varying_births(time_data, beta_est_tv, gamma_est_tv, lambda0_est_tv, lambda1_est_tv, omega_est_tv, phi_est_tv, mu_est_tv, S0_est_tv, I0_est_tv, R0 = 0)

# Check for NA values in the solver output for time-varying births model
if (any(is.na(SIR_fitted_tv))) {
  stop("Time-varying births model solver failed with the fitted parameters.")
}

# Extract each state for time-varying births model
S_fitted_tv <- SIR_fitted_tv[, "S.S0"]
I_fitted_tv <- SIR_fitted_tv[, "I.I0"]
R_fitted_tv <- SIR_fitted_tv[, "R"]

# Combine data into a single data frame for plotting
data_tv <- data.frame(
  time = time_data,
  S_fitted_tv = S_fitted_tv,
  I_fitted_tv = I_fitted_tv,
  R_fitted_tv = R_fitted_tv,
  R_observed = R_data
)

# Plot the results for the time-varying births model
ggplot(data_tv, aes(x = time)) +
  geom_line(aes(y = S_fitted_tv, color = "S (time-varying)"), size = 1) +
  geom_line(aes(y = I_fitted_tv, color = "I (time-varying)"), size = 1) +
  geom_line(aes(y = R_fitted_tv, color = "R (fitted time-varying)"), size = 1) +
  geom_point(aes(y = R_observed, color = "R (observed)"), size = 3, shape = 16) +
  labs(title = "SIR Model with Time-Varying Births",
       x = "Time (days)", y = "Proportion of Population") +
  scale_color_manual(name = "State",
                     values = c("S (time-varying)" = "blue", "I (time-varying)" = "red", 
                                "R (fitted time-varying)" = "green", "R (observed)" = "black")) +
  theme_minimal()

### Example code #3 version 2 ####
##################################

# SIR model differential equations with time-varying births and deaths
sir_model_with_time_varying_births <- function(time, state, parameters) {
  S <- state[1]
  I <- state[2]
  R <- state[3]
  
  beta <- parameters[1]
  gamma <- parameters[2]
  lambda_0 <- parameters[3]   # Average birth rate
  lambda_amp <- parameters[4]  # Amplitude of birth rate oscillations
  T <- parameters[5]           # Period of oscillations
  mu <- parameters[6]          # Death rate
  
  # Time-varying birth rate using a sine function
  lambda_t <- lambda_0 + lambda_amp * sin(2 * pi * time / T)
  
  dS <- lambda_t - beta * S * I - mu * S
  dI <- beta * S * I - gamma * I - mu * I
  dR <- gamma * I - mu * R
  
  list(c(dS, dI, dR))
}

# Solve SIR model with time-varying births
solve_sir_with_time_varying_births <- function(time, beta, gamma, lambda_0, lambda_amp, T, mu, S0, I0, R0) {
  initial_state <- c(S = S0, I = I0, R = R0)
  parameters <- c(beta, gamma, lambda_0, lambda_amp, T, mu)
  
  out <- try(ode(y = initial_state, times = time, func = sir_model_with_time_varying_births, parms = parameters, 
                 method = "lsoda", atol = 1e-8, rtol = 1e-8))  # Tighter tolerances
  
  if (inherits(out, "try-error") || any(is.na(out))) {
    return(NA)
  } else {
    return(out)
  }
}

# Objective function to minimize for the time-varying births model
objective_function_with_time_varying_births <- function(params, time, R_data) {
  beta <- params[1]
  gamma <- params[2]
  lambda_0 <- params[3]
  lambda_amp <- params[4]
  T <- params[5]
  mu <- params[6]
  S0 <- params[7]
  I0 <- params[8]
  
  # Solve the SIR model with the given parameters
  sir_output <- solve_sir_with_time_varying_births(time, beta, gamma, lambda_0, lambda_amp, T, mu, S0, I0, R0 = 0)
  
  if (is.na(sir_output[1]) || any(is.na(sir_output))) {
    return(1e10)  # Large penalty for invalid solutions
  }
  
  R_predicted <- sir_output[, "R"]
  
  # Calculate SSE (sum of squared errors)
  SSE <- sum((R_data - R_predicted)^2)
  
  return(SSE)
}

# Simulated time series over 100 days
time_data <- seq(0, 100, by = 1)  # Simulate for 100 days

# Simulated observed recovered data (plausible epidemic dynamics)
R_data <- c(0, 0.002, 0.005, 0.012, 0.03, 0.08, 0.15, 0.28, 0.42, 0.55, 0.68, 
            0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.93, 0.94, 0.95, 0.95, 0.96, 
            0.97, 0.97, 0.97, 0.98, 0.98, 0.98, 0.99, 0.99, 0.99, 1.0)

# Ensure length of R_data matches time_data
R_data <- c(R_data, rep(1, length(time_data) - length(R_data)))

# Initial guesses for parameters for the SIR model with time-varying births
initial_params_sir_tvb <- c(beta = 0.4, gamma = 0.1, lambda_0 = 0.02, lambda_amp = 0.01, T = 365, mu = 0.01, S0 = 0.99, I0 = 0.01)

# Optimize the parameters for the time-varying births model
fit_sir_tvb <- optim(par = initial_params_sir_tvb, fn = objective_function_with_time_varying_births, 
                     time = time_data, R_data = R_data, method = "L-BFGS-B", 
                     lower = c(0, 0, 0, 0, 10, 0, 0, 0), upper = c(2, 1, 1, 1, 365, 1, 1, 1))

# Extract fitted parameters for the time-varying births model
beta_est_tvb <- fit_sir_tvb$par[1]
gamma_est_tvb <- fit_sir_tvb$par[2]
lambda_0_est_tvb <- fit_sir_tvb$par[3]
lambda_amp_est_tvb <- fit_sir_tvb$par[4]
T_est_tvb <- fit_sir_tvb$par[5]
mu_est_tvb <- fit_sir_tvb$par[6]
S0_est_tvb <- fit_sir_tvb$par[7]
I0_est_tvb <- fit_sir_tvb$par[8]

# Solve the model with the fitted parameters
SIR_fitted_tvb <- solve_sir_with_time_varying_births(time_data, beta_est_tvb, gamma_est_tvb, lambda_0_est_tvb, lambda_amp_est_tvb, T_est_tvb, mu_est_tvb, S0_est_tvb, I0_est_tvb, R0 = 0)

# Check for NA values in the solver output for the time-varying births model
if (any(is.na(SIR_fitted_tvb))) {
  stop("Time-varying births model solver failed with the fitted parameters.")
}

# Extract each state for the time-varying births model
S_fitted_tvb <- SIR_fitted_tvb[, "S.S0"]
I_fitted_tvb <- SIR_fitted_tvb[, "I.I0"]
R_fitted_tvb <- SIR_fitted_tvb[, "R"]

# Combine data into a single data frame for plotting
data <- data.frame(
  time = time_data,
  R_observed = R_data,
  S_fitted_tvb = S_fitted_tvb,
  I_fitted_tvb = I_fitted_tvb,
  R_fitted_tvb = R_fitted_tvb
)

# Plot the model output and the observed R data
ggplot(data, aes(x = time)) +
  geom_line(aes(y = S_fitted_tvb, color = "S (time-varying births)"), size = 1) +
  geom_line(aes(y = I_fitted_tvb, color = "I (time-varying births)"), size = 1) +
  geom_line(aes(y = R_fitted_tvb, color = "R (fitted time-varying births)"), size = 1) +
  geom_point(aes(y = R_observed, color = "R (observed)"), size = 3, shape = 16) +
  labs(title = "SIR Model with Time-Varying Births", x = "Time (days)", y = "Proportion of Population") +
  scale_color_manual(name = "State",
                     values = c("S (time-varying births)" = "blue", "I (time-varying births)" = "red", 
                                "R (fitted time-varying births)" = "green", "R (observed)" = "black")) +
  theme_minimal()

