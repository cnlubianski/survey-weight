# ------------------------------------------------------------------------------
# Title: Simulation Utility Functions
# General functions to be sourced in by the simulation scripts. These include
# functions to generate the population data and for sampling from the population.
# ------------------------------------------------------------------------------

# Function to generate data from Wang et al.'s (2023) Study 1
generate_data_study1 <- function(N, sigma, alpha, delta) {
  X <- runif(N, 0, 1)
  u <- runif(N, 0, 1)
  epsilon <- rnorm(N, 0, sd = sigma)
  Y <- 1 + X + epsilon
  w <- alpha * Y + 0.3 * X + delta * u
  data.frame(y = Y, x = X, w = w)
}

# Function to generate data from Wang et al.'s (2023) Study 2
generate_data_study2 = function(N, sigma, alpha) {
  X <- runif(N, 0, 1)
  u <- runif(N, 0, 1)
  epsilon <- rnorm(N, 0, sd = sigma)
  
  Y <- 1 + X + epsilon
  w <- alpha * (Y - 1.5 * alpha)^2 + 0.3 * X - 0.3 * X^2 + u
  data <- data.frame(y = Y, x = X, w)
  return(data)
}

# Function to generate data from Wang et al.'s (2023) Study 3
generate_data_study3 = function(N, alpha, psi) {
  X <- rnorm(N, 0, sd = sqrt(0.5))
  epsilon <- rnorm(N, 0, sd = sqrt(0.5))
  z <- rnorm(N, 0, sd = sqrt(0.5))
  beta <- 2 - alpha
  
  eta <- function(x) {
    ifelse(x < 0.2, 0.025, 
           ifelse(0.2 <= x & x <= 1.2, 0.475 * (x - 0.2) + 0.025, 0.5))
  }
  Y <- 0.5 + X + epsilon
  w <- alpha * eta(X) + beta * eta(psi * epsilon + (1 - psi) * z)
  data <- data.frame(y = Y, x = X, w)
  return(data)
}

# Function to generate the data from Wang et al.'s (2023) Study 2 with different error distributions
generate_data_study2_perm <- function(N, sigma = 0.2, delta = 1, alpha, distribution) {
  X <- runif(N, 0, 1)
  u <- runif(N, 0, 1)
  sigmas <- abs(rnorm(N, sigma, sigma / 3))
  
  epsilon <- NULL
  if (distribution == "Normal") {
    epsilon <- rnorm(N, 0, sd = sigmas)
  } 
  if (distribution == "Uniform") {
    epsilon <- runif(N, min = -sqrt(3 * sigmas^2), max = sqrt(3 * sigmas^2))
  } 
  if (distribution == "Gamma") {
    generated <- rgamma(N, shape = 10, scale = sqrt(sigmas^2 / 10))
    epsilon <- generated - mean(generated)
  } 
  if (distribution == "t") {
    scale <- sqrt((sigmas^2 * 3) / 5)
    epsilon <- scale * rt(N, df = 5)
  }
  
  Y <- 1 + X + epsilon
  w <- alpha * Y + 0.3 * X + 1 * u
  data <- data.frame(y = Y, x = X, w)
  return(data)
}

# Function to retrieve sample from population using Brewer sampling given n sample size
generate_sample_brewer <- function(data, w, n) {
  pik <- sampling::inclusionprobabilities(w, n)
  choosen <- sampling::UPbrewer(pik)
  samp <- data[1:nrow(data) * choosen, ] |>
    mutate(w = 1 / pik[1:nrow(data) * choosen])
  samp
}

# Function to retrieve sample from population using Poisson sampling given n expected sample size
generate_sample_poisson <- function(data, w, n, rescale = FALSE) {
  choosen <- as.numeric(runif(length(w), 0, (1 / n) * sum(w)) < w)
  samp <- cbind(data, choosen) %>%
    filter(choosen == 1) %>%
    select(-choosen) %>%
    mutate(w = 1 / w) # Redefine from pi to weights w
  if (rescale == TRUE) mutate(samp, w = w / sum(w))
  return(samp)
}

# The following two functions are versions above but they do not redefine weights
generate_sample_poisson_not_presumed <- function(data, w, n, rescale = FALSE) {
  choosen <- as.numeric(runif(length(w), 0, (1 / n) * sum(w)) < w)
  samp <- cbind(data, choosen) %>%
    filter(choosen == 1) %>%
    select(-choosen)
  if (rescale == TRUE) mutate(samp, w = w / sum(w))
  return(samp)
}

generate_sample_brewer_not_presumed <- function(data, w, n, rescale = FALSE) {
  pik <- inclusionprobabilities(w, n)
  choosen <- UPbrewer(pik)
  samp <- data[1:nrow(data) * choosen,]
  if (rescale == TRUE) mutate(samp, w = w / sum(w))
  return(samp)
}
