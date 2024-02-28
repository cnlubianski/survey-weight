# Install packages -------------------------------------------------------------
if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if (!require(rpms)) install.packages("rpms"); library(rpms)
if (!require(sampling)) install.packages("sampling"); library(sampling)
if (!require(survey)) install.packages("survey"); library(survey)


# Sampling Function - Wang et al (2023) Study 1 for testing --------------------
generate_data_study1 = function(N, sigma, alpha, delta) {
  X <- runif(N, 0, 1)
  u <- runif(N, 0, 1)
  epsilon <- rnorm(N, 0, sd = sigma)
  
  Y <- 1 + X + epsilon
  w <- alpha * Y + 0.3 * X + delta * u
  data = data.frame(y = Y, x = X, w)
  return(data)
}

generate_sample_brewer = function(data, w, n, rescale = FALSE) {
  pik = inclusionprobabilities(w, n)
  choosen = UPbrewer(pik)
  samp = data[1:nrow(data) * choosen,] %>%
    mutate(w = 1 / pik[1:nrow(data) * choosen])
  if (rescale == TRUE) mutate(samp, w = w / sum(w))
  return(samp)
}

# Set cases and case grid
N <- 3000
n <- c(100, 200)
sigma <- c(0.1, 0.2)
alpha <- c(0, 0.2, 0.4, 0.6)
delta <- c(1, 1.5)

cases <- expand_grid(N, n, sigma, delta, alpha)

# If we want to sample outside of the function for testing purposes
pop = generate_data_study1(N = cases$N[m],
                           sigma = cases$sigma[m],
                           alpha = cases$alpha[m],
                           delta = cases$delta[m])

# Keeping all units to test out the test code
samp = generate_sample_brewer(data = pop,
                              w = pop$w, 
                              n = cases$n[m], #nrow(pop) / 2,
                              rescale = FALSE)


# Difference in Coefficients Test ----------------------------------------------

## Hausman-Pfeffermann DC Test =================================================

# TO-DO: incorporate categorical variables
HP_DC_test = function(data, y, x, wts) {
  # Unweighted Regression
  X = cbind(1, x)
  betas_u = solve(t(X) %*% X) %*% t(X) %*% y
  
  # Weighted Regression
  W =  diag(x = wts, nrow = length(wts), ncol = length(wts))
  betas_w = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
  
  # Calculate Test Statistic
  
  # Calculate sigma^2 estimate from OLS under null
  y_hat = X %*% betas_u
  residuals = y - y_hat
  SSE = sum(residuals^2)
  sigma_sq_hat = SSE / (length(y) - length(betas_u) - 1)
  
  # A (note that in paper, H = diag(w) which is our W)
  A = (solve(t(X) %*% W %*% X) %*% t(X) %*% W) - (solve(t(X) %*% X) %*% t(X))
  
  # Calculate variance
  V_hat = sigma_sq_hat * A %*% t(A)
  
  # Chi Squared Test Statistic
  Chi_statistic = t(betas_u - betas_w) %*% solve(V_hat) %*% (betas_u - betas_w)
  p_value = pchisq(q = Chi_statistic, df = length(betas_u), lower.tail = FALSE)
  return(p_value)
}

# HP_DC_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)




# Weight Association Tests -----------------------------------------------------

## DuMouchel-Duncan WA Test ====================================================

DD_WA_test = function(data, y, x, wts) {
  # Create matrices
  X = cbind(1, x)
  W = diag(x = wts, ncol = length(wts), nrow = length(wts))
  X_tilde = W %*% X 
  X_comb = cbind(X, X_tilde)
  
  # Full model
  betas_comb = solve(t(X_comb) %*% X_comb) %*% t(X_comb) %*% y
  y_hat_full = X_comb %*% betas_comb
  RSS_full = sum((y - y_hat_full)^2)
  
  # Reduced model
  X_reduced = X_comb[, c(1,2)]
  betas_reduced = solve(t(X_reduced) %*% X_reduced) %*% t(X_reduced) %*% y
  y_hat_reduced = X_reduced %*% betas_reduced
  RSS_reduced = sum((y - y_hat_reduced)^2)
  
  # F-test
  F_statistic = (RSS_reduced - RSS_full) / (length(betas_comb) - length(betas_reduced)) /
    (RSS_full / (length(y) - ncol(X_comb)))
  p_value = 1 - pf(F_statistic, df1 = ncol(X_tilde), df2 = length(y) - ncol(X_comb))
  return(p_value)
}

# DD_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)

## Pfeffermann-Sverchkov WA Test 1 - Correlation Testing =======================

PS1_WA_test = function(data, y, x, wts) {
  W = wts
  
  # Estimate residuals with unweighted regression
  X = cbind(1, x)
  betas_u = solve(t(X) %*% X) %*% t(X) %*% y
  
  y_hat = X %*% betas_u
  residuals = y - y_hat
  residuals_diag = diag(x = as.vector(residuals), ncol = length(residuals), nrow = length(residuals))
  
  # Create Design matrix
  E = residuals
  E_sq = residuals^2
  X_tilde = residuals_diag %*% x
  
  X_design = cbind(1, x, E, E_sq, X_tilde)
  betas = solve(t(X_design) %*% X_design) %*% t(X_design) %*% W
  
  # Compute F-statistic
  W_hat = X_design %*% betas
  RSS = sum((W - W_hat)^2)
  TSS = sum((W - mean(W))^2)
  F_statistic = ((TSS - RSS) / (length(betas) - ncol(X))) / 
    (RSS / (length(W) - length(betas) - 1))
  p_value <- 1 - pf(F_statistic,
                    df1 = length(betas) - ncol(X),
                    df2 = length(W) - length(betas) - 1)
  return(p_value)
}

# PS1_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)

# For quadratic X terms
PS1q_WA_test = function(data, y, x, wts) { #TO-DO: Figure out circumstances when p != 0
  W = wts
  
  # Estimate residuals with unweighted regression
  X = cbind(1, x)
  betas_u = solve(t(X) %*% X) %*% t(X) %*% y
  
  y_hat = X %*% betas_u
  residuals = y - y_hat
  residuals_diag = diag(x = as.vector(residuals), ncol = length(residuals), nrow = length(residuals))
  
  # Create Design matrix
  X_sq = x^2
  E = residuals
  E_sq = residuals^2
  X_tilde = residuals_diag %*% x
  
  X_design = cbind(1, x, X_sq, E, E_sq, X_tilde)
  betas = solve(t(X_design) %*% X_design) %*% t(X_design) %*% W
  
  # Compute F-statistic
  W_hat = X_design %*% betas
  RSS = sum((W - W_hat)^2)
  TSS = sum((W - mean(W))^2)
  F_statistic = ((TSS - RSS) / (length(betas) - ncol(cbind(X, X_sq)))) / 
    (RSS / (length(W) - length(betas) - 1))
  p_value <- 1 - pf(F_statistic,
                    df1 = length(betas) - ncol(cbind(X, X_sq)),
                    df2 = length(W) - length(betas) - 1)
  return(p_value)
}

# PS1q_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)

## Pfeffermann-Sverchkov WA Test 2 - Weight Informative ========================

PS2_WA_test = function(data, y, x, wts) {
  W = wts
  
  # Full model: Regressing W on both X and Y with quadratic terms 
  XY_full = cbind(1, x, y, y^2)
  betas_full = solve(t(XY_full) %*% XY_full) %*% t(XY_full) %*% W
  w_hat_full = XY_full %*% betas_full
  RSS_full = sum((W - w_hat_full)^2)
  
  # Reduced model: Regressing W only on X and X^2
  XY_reduced = cbind(1, x)
  betas_reduced = solve(t(XY_reduced) %*% XY_reduced) %*% t(XY_reduced) %*% W
  w_hat_reduced = XY_reduced %*% betas_reduced
  RSS_reduced = sum((W - w_hat_reduced)^2)
  
  # F-test
  F_statistic = ((RSS_reduced - RSS_full) / length(betas_reduced[-1,])) /
    (RSS_full / (length(W) - length(betas_full[-1,])))
  p_value = 1 - pf(F_statistic, df1 = length(betas_reduced[-1,]),
                   df2 = length(W) - length(betas_full[-1,]))
  return(p_value)
}

# PS2_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)

PS2q_WA_test = function(data, y, x, wts) {
  W = wts
  
  # Full model: Regressing W on both X and Y with quadratic terms 
  XY_full = cbind(1, x, x^2, y, y^2)
  betas_full = solve(t(XY_full) %*% XY_full) %*% t(XY_full) %*% W
  w_hat_full = XY_full %*% betas_full
  RSS_full = sum((W - w_hat_full)^2)
  
  # Reduced model: Regressing W only on X and X^2
  XY_reduced = cbind(1, x, x^2)
  betas_reduced = solve(t(XY_reduced) %*% XY_reduced) %*% t(XY_reduced) %*% W
  w_hat_reduced = XY_reduced %*% betas_reduced
  RSS_reduced = sum((W - w_hat_reduced)^2)
  
  # F-test
  F_statistic = ((RSS_reduced - RSS_full) / length(betas_reduced[-1,])) /
    (RSS_full / (length(W) - length(betas_full[-1,])))
  p_value = 1 - pf(F_statistic, df1 = length(betas_reduced[-1,]),
                   df2 = length(W) - length(betas_full[-1,]))
  return(p_value)
}

# PS2q_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)

## Wu-Fuller WA Test ===========================================================

WF_WA_test = function(data, y, x, wts) {
  W = wts
  
  # Auxiliary Regression
  X_design = cbind(1, x, x^2)
  etas = solve(t(X_design) %*% X_design) %*% t(X_design) %*% W
  W_hat = X_design %*% etas
  q = W / W_hat
  
  # Full Regression
  X_tilde = diag(as.vector(q)) %*% x
  XQ_design = cbind(1, x, X_tilde)
  betas_full = solve(t(XQ_design) %*% XQ_design) %*% t(XQ_design) %*% y
  y_hat_full = XQ_design %*% betas_full
  RSS_full = sum((y - y_hat_full)^2)
  
  # Reduced Regression
  X = cbind(1, x)
  betas_reduced = solve(t(X) %*% X) %*% t(X) %*% y
  y_hat_reduced = X %*% betas_reduced
  RSS_reduced = sum((y - y_hat_reduced)^2)
  
  # F-test - Testing whether beta coefficient for X_tilde is stat sig
  F_statistic = ((RSS_reduced - RSS_full) / length(betas_reduced[-1,])) /
    (RSS_full / (length(y) - length(betas_full)))
  p_value = 1 - pf(F_statistic, df1 = length(betas_reduced[-1,]), df2 = length(y) - length(betas_full))
  return(p_value)
}

# WF_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)


# Other ------------------------------------------------------------------------

## Pfeffermann-Sverchkov Other Test - Estimating Equations =====================

PS3_test = function(data, y, x, wts) {
  W = wts
  
  # Auxiliary Regression
  X_design = cbind(1, x, x^2)
  etas = solve(t(X_design) %*% X_design) %*% t(X_design) %*% W
  W_hat = X_design %*% etas
  q = W / W_hat
  
  # Estimating Regression
  X = cbind(1, x)
  betas = solve(t(X) %*% X) %*% t(X) %*% y
  
  # Unweighted estimating function
  delta = rep(NA, length(y))
  
  for (i in 1:length(y)) {
    delta[i] = x[i] * (y[i] - X[i,] %*% betas)
  }
  
  R_x = delta - q * delta
  
  # Hotelling Test Statistic
  F_statistic = ((length(y) - length(betas[-1])) / length(betas[-1])) *
    t(mean(R_x)) %*% solve(cov(R_x)) %*% mean(R_x)
  p_value = 1 - pf(F_statistic, df1 = length(betas[-1]), df2 = length(y) - length(betas[-1]))
  return(p_value)
}

# PS3_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)

## Pfeffermann-Nathan Predictive Power Test ===================================

# Fix to generalize the parameters x,y,wts from any dataframe

PN_test = function(data, y, x, wts, est_split) {
  est_index = sample(1:nrow(data), floor(est_split * length(x)))
  x_est = x[est_index]
  x_val = x[-est_index]
  y_est = y[est_index]
  y_val = y[-est_index]
  wts_est = wts[est_index]
  wts_val = wts[-est_index]
  
  est_samp = data[est_index,]
  val_samp = data[-est_index,]
  
  # Unweighted Regression --------------
  X = cbind(1, x_est)
  betas_u = solve(t(X) %*% X) %*% t(X) %*% y_est
  y_val_u = cbind(1, x_val) %*% betas_u
  v_u = y_val - y_val_u
  
  # Weighted Regression ----------------
  W =  diag(x = wts_est)
  betas_w = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y_est
  y_val_w = cbind(1, x_val) %*% betas_w
  v_w = y_val - y_val_w
  
  # Standard Z-test
  D = v_u^2 - v_w^2
  Z = (mean(D) - 0) / (sd(D) / sqrt(length(D)))
  p_value = as.numeric(2 * (1 - pnorm(abs(Z))))
  return(p_value)
}

# PN_test(data = samp, y = samp$y, x = samp$x, wts = samp$w, est_split = 0.5)



## Breidt Likelihood-Ratio Test ================================================

LR_test = function(data, y, x, wts) {
  X = cbind(1, x)
  
  # Define the log-likelihood function
  logLike <- function(params, y, X, wts) {
    beta <- params[1:ncol(X)]
    sigma2 <- exp(params[ncol(X)+1])  # we use exp to ensure sigma^2 is positive
    mu <- X %*% beta
    n <- length(y)
    logLik <- -0.5 * log(2 * pi * sigma2) * sum(wts) - 0.5 * sum(wts * (y - mu)^2) / sigma2
    return(-logLik)  # We negate because optim() minimizes by default
  }
  
  initial <- c(rep(0, ncol(X)), log(var(y)))
  
  # theta_u estimates
  res <- optim(initial, logLike, y = y, X = X, wts = rep(sum(wts) / length(wts), length(wts)))
  beta_u <- res$par[1:ncol(X)]
  sigma_sq_u <- exp(res$par[ncol(X)+1])
  theta_u = c(beta_u, sigma_sq_u)
  loglike_u = res$value
  
  # theta-w estimates
  res <- optim(initial, logLike, y = y, X = X, wts = wts)
  beta_w <- res$par[1:ncol(X)]
  sigma_sq_w <- exp(res$par[ncol(X)+1])
  theta_w = c(beta_w, sigma_sq_w)
  loglike_w = res$value
  
  # True parameter
  X = cbind(1, x)
  beta_true = solve(t(X) %*% X) %*% t(X) %*% y
  sigma_sq = mean((y - X %*% beta_true)^2)
  theta_true = c(beta_true, sigma_sq)
  
  # Test Statistic
  null_weights = rep(sum(wts) / length(wts), length(wts))
  l_u_theta_u = 0.5 * log(2 * pi * sigma_sq_u) * sum(null_weights) -
    0.5 * sum(null_weights * (y - X %*% beta_u)^2) / sigma_sq_u
  l_u_theta_w = 0.5 * log(2 * pi * sigma_sq_w) * sum(null_weights) - 
    0.5 * sum(null_weights * (y - X %*% beta_w)^2) / sigma_sq_w
  
  test_statistic = 2 * (l_u_theta_u - l_u_theta_w)
  
  # Asymptotic Test Distribution (Using Fisher Information matrix)
  J_u = diag(c(sum(t(X) %*% X / sigma_sq), # with respect to X, not X_i
               1 / (2 * (sigma_sq)^2)))
  J_w = diag(c(sum(t(X) %*% diag(wts) %*% X / sigma_sq),
               sum(wts) / (2 * length(wts) * (sigma_sq)^2)))
  K_w = diag(c(sum(t(X) %*% diag(wts^2) %*% X / sigma_sq),
               sum(wts^2) / (2 * length(wts) * (sigma_sq)^2)))
  
  Gamma = solve(J_w) %*% K_w %*% solve(J_w) - solve(J_u)
  
  stat_matrix = t(sqrt(Gamma)) %*% J_u %*% sqrt(Gamma)
  eigenvalues = eigen(stat_matrix, only.values = TRUE)$values
  
  df = 1:length(eigenvalues)
  chi_square_rv = sapply(df, function(d) rchisq(10000, df = d))
  asymptotic_comb = (1 / chi_square_rv) %*% eigenvalues
  ecdf_LR = ecdf(asymptotic_comb)
  p_value = 1 - ecdf_LR(test_statistic)
  return(p_value)
}

# LR_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
## Confidence Intervals Overlap Test ===========================================

confint_test = function(data, y, x, wts, alpha = 0.05) {
  X = cbind(1, x)
  
  # Unweighted beta CI
  betas_u = solve(t(X) %*% X) %*% t(X) %*% y
  residuals_u = y -  X %*% betas_u
  MSE_u = sum(residuals_u^2) / (nrow(X) - ncol(X))
  vcov_matrix = MSE_u * solve(t(X) %*% X)
  se_u = sqrt(diag(vcov_matrix))
  t_value = qt(1 - alpha / 2, df = nrow(X) - ncol(X))
  
  CI_u = cbind(betas_u - t_value * se_u, betas_u + t_value * se_u)
  beta_1_u = CI_u[2,]
  
  # Weighted beta CI
  W = diag(wts)
  betas_w = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
  residuals_w = y -  X %*% betas_w
  MSE_w = sum(residuals_w^2) / (nrow(X) - ncol(X))
  vcov_matrix = MSE_w * solve(t(X) %*% X)
  se_w = sqrt(diag(vcov_matrix))
  t_value = qt(1 - alpha / 2, df = nrow(X) - ncol(X))
  
  CI_w = cbind(betas_w - t_value * se_w, betas_w + t_value * se_w)
  beta_1_w = CI_w[2,]
  
  # Determine overlap
  overlap = !(beta_1_w[2] < beta_1_u[1] | beta_1_u[2] < beta_1_w[1])
  return(as.numeric(overlap)) # if overlap, returns 1 so fail to reject the null
}

# confint_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)

# Permutation Proposal Test (TO-DO) --------------------------------------------

weight_perm_test <- function(y, x, w, B) {
  stat_stor = rep(NA, B)
  data = data.frame(x = x, y = y, w = w)
  # Permutating
  for (b in 1:B) {
    shuffled_w = sample(w)
    b_data = data.frame(y = y, x = x, w = shuffled_w, id = length(w))
    
    # Unweighted Regression
    unweighted = lm(y ~ x, data = b_data)
    beta_u = unweighted$coefficients[2] # faster to manually calculate beta_1
    #beta_u = sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
    
    # Weighted Regression
    design = svydesign(id = ~1, weights = ~w, fpc = ~rep(N, nrow(b_data)),
                       data = b_data)
    weighted = svyglm(y ~ x, design = design)
    beta_w = weighted$coefficients[2]
    
    stat_stor[b] = beta_w - beta_u
  }
  
  # Estimating test statistic
  act_betau = lm(y ~ x, data = data)$coefficients[2]
  
  act_design = svydesign(id = ~1, weights = ~w, fpc = ~rep(N, nrow(data)),
                         data = data)
  act_betaw = svyglm(y ~ x, design = act_design)$coefficients[2]
  est_stat = act_betaw - act_betau
  
  # Calculating p-value
  dist = ecdf(stat_stor)
  p = dist(est_stat)
  pvalue = 2 * min(p, 1 - p)
  
  # return list for distribution and p-value
  return(pvalue)
}

B = 100
stat_stor = rep(NA, B)

data = samp

for (b in 1:B) {
  b_data = mutate(data, w = sample(w))
  
  # Unweighted Regression
  X = cbind(1, b_data$x)
  Y = b_data$y
  betas_u = solve(t(X) %*% X) %*% t(X) %*% Y
  
  # Weighted Regression
  W = diag(x = b_data$w)
  betas_w = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
  
  stat_stor[b] = betas_w[2] - betas_u[2]
}

# Calculating Test Statistic
# Unweighted Regression
X = cbind(1, data$x)
Y = data$y
betas_u = solve(t(X) %*% X) %*% t(X) %*% Y

# Weighted Regression
W = diag(x = data$w)
betas_w = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y

est_stat = betas_w[2] - betas_u[2]

# Calculating p-value
dist = ecdf(stat_stor)
p = dist(est_stat)
p_value = 2 * min(p, 1 - p)
p_value


# hist(stat_stor)
p_value <- sum(abs(stat_stor) >= abs(est_stat)) / length(stat_stor)
p_value


hist(stat_stor)
