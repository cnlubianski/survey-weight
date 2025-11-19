## Set-up

if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if (!require(rpms)) install.packages("rpms"); library(rpms)
if (!require(sampling)) install.packages("sampling"); library(sampling)
if (!require(survey)) install.packages("survey"); library(survey)


## Permutation 1: DC Test using HP variation

perm_HP <- function(y, x, wts, B = 1000, replacement = FALSE) {
  stat_stor = rep(NA, B)
  X = cbind(1, x)
  
  for (b in 1:B) {
    wts_b = sample(wts, replace = replacement)
    W = diag(wts_b)
    betas_u = solve(t(X) %*% X) %*% t(X) %*% y
    betas_w = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
    stat_stor[b] = betas_w[2] - betas_u[2]
  }
  
  # Actual estimates
  W = diag(wts)
  betas_u = solve(t(X) %*% X) %*% t(X) %*% y
  betas_w = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
  est_stat = betas_w[2] - betas_u[2]
  
  # Calculating p-value
  p_o = mean(ifelse(stat_stor >= est_stat, 1, 0))
  p_value = 2 * min(p_o, 1 - p_o)
  return(p_value)
}

#perm_HP(y = samp$y, x = samp$x, wts = samp$w, B = 1000)


## Permutation 2: WA Test using DD variation

perm_DD <- function(y, x, wts, B = 1000, replacement = FALSE) {
  stat_stor = rep(NA, B)
  X = cbind(1, x)
  
  for (b in 1:B) {
    wts_b = sample(wts, replace = replacement)
    
    W = diag(wts_b)
    X_tilde = W %*% X
    X_comb = cbind(X, X_tilde)
    betas_comb = solve(t(X_comb) %*% X_comb) %*% t(X_comb) %*% y
    stat_stor[b] = betas_comb[4]
  }
  
  W = diag(wts)
  X_tilde = W %*% X
  X_comb = cbind(X, X_tilde)
  betas_comb = solve(t(X_comb) %*% X_comb) %*% t(X_comb) %*% y
  est_stat = betas_comb[4]
  
  p_o = mean(ifelse(stat_stor >= est_stat, 1, 0))
  p_value = 2 * min(p_o, 1 - p_o)
  return(p_value)
}

#perm_DD(y = samp$y, x = samp$x, wts = samp$w, B = 1000)


## Permutation 3: WA Test using correlations

perm_PS1 <- function(y, x, wts, B = 1000, replacement = FALSE) {
  stat_stor = rep(NA, B)
  X = cbind(1, x)
  
  for (b in 1:B) {
    wts_b = sample(wts, replace = replacement)
    betas_u = solve(t(X) %*% X) %*% t(X) %*% y
    residuals = y - X %*% betas_u
    stat_stor[b] = cor(residuals, wts_b)
  }
  
  betas_u = solve(t(X) %*% X) %*% t(X) %*% y
  residuals = y - X %*% betas_u
  est_stat = cor(residuals, wts)
  
  p_o = mean(ifelse(stat_stor >= rep(est_stat, B), 1, 0))
  p_value = 2 * min(p_o, 1 - p_o)
  return(p_value)
}

#perm_PS1(y = samp$y, x = samp$x, wts = samp$w, B = 10000)


## Permutation 4 - RMSE

RMSE = function(y, yhat){
  SSE = sum((y - yhat)^2)
  return(sqrt(SSE / length(y)))  
}

perm_RMSE <- function(y, x, wts, B = 1000, replacement = FALSE) {
  stat_stor = rep(NA, B)
  X = cbind(1, x)
  
  for (b in 1:B) {
    wts_b = sample(wts, replace = replacement)
    W = diag(wts_b)
    betas_u = solve(t(X) %*% X) %*% t(X) %*% y
    betas_w = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
    stat_stor[b] = RMSE(y, X %*% betas_u) - RMSE(y, X %*% betas_w)
  }
  
  # Actual estimates
  W = diag(wts)
  betas_u = solve(t(X) %*% X) %*% t(X) %*% y
  betas_w = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
  est_stat = RMSE(y, X %*% betas_u) - RMSE(y, X %*% betas_w)
  
  # Calculating p-value
  p_o = mean(ifelse(stat_stor >= est_stat, 1, 0))
  p_value = 2 * min(p_o, 1 - p_o)
  return(p_value)
}

#perm_RMSE(y = samp$y, x = samp$x, wts = samp$w, B = 100)


## Permutation 5 - PN variation

perm_PN <- function(data, y, x, wts, B = 1000, est_split = 0.7, replacement = FALSE) {
  stat_stor = rep(NA, B)
  index = sample(1:nrow(data), floor(est_split * length(x)))
  
  for (b in 1:B) {
    index_b = index
    #index_b = sample(1:nrow(data), floor(est_split * length(x)))
    wts_b = sample(wts[index_b], replace = replacement)
    
    # Unweighted Regression
    X = cbind(1, x[index_b])
    betas_u = solve(t(X) %*% X) %*% t(X) %*% y[index_b]
    y_val_u = cbind(1, x[-index_b]) %*% betas_u
    v_u = y[-index_b] - y_val_u
    
    # Weighted Regression 
    W = diag(wts_b)
    betas_w = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y[index_b]
    y_val_w = cbind(1, x[-index_b]) %*% betas_w
    v_w = y[-index_b] - y_val_w
    
    stat_stor[b] = mean(v_u^2 - v_w^2)
  }
  X = cbind(1, x[index])
  betas_u = solve(t(X) %*% X) %*% t(X) %*% y[index]
  y_val_u = cbind(1, x[-index]) %*% betas_u
  v_u = y[-index] - y_val_u
  
  # Weighted Regression 
  W = diag(wts[index])
  betas_w = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y[index]
  y_val_w = cbind(1, x[-index]) %*% betas_w
  v_w = y[-index] - y_val_w
  
  # Standard Z-test
  est_stat = mean(v_u^2 - v_w^2)
  p_o = mean(ifelse(stat_stor >= rep(est_stat, B), 1, 0))
  p_value = 2 * min(p_o, 1 - p_o)
  return(p_value)
}

#perm_PN(data = samp, y = samp$y, x = samp$x, wts = samp$w, B = 1000)


## Permutation 6 - PS2

perm_PS2 <- function(y, x, wts, B = 1000, replacement = FALSE) {
  stat_stor = rep(NA, B)
  XY_design = cbind(1, x, y)
  
  for (b in 1:B) {
    wts_b = sample(wts, replace = replacement)
    betas = solve(t(XY_design) %*% XY_design) %*% t(XY_design) %*% wts_b
    stat_stor[b] = betas[3]
  }
  
  betas = solve(t(XY_design) %*% XY_design) %*% t(XY_design) %*% wts
  est_stat = betas[3]
  p_o = mean(ifelse(stat_stor >= rep(est_stat, B), 1, 0))
  p_value = 2 * min(p_o, 1 - p_o)
  return(p_value)
}

#perm_PS2(y = samp$FINCBTAX, x = samp$TOTEXPCQ, wts = samp$wts, B = 1000)


## Permutation 7 - GREG

perm_greg <- function(y, x_s, x_u, wts, B = 1000, replacement = FALSE) {
  stat_stor = rep(NA, B)
  X = cbind(1, x_s)
  
  for (b in 1:B) {
    wts_b = sample(wts, replace = replacement)
    betas = solve(t(X) %*% diag(wts_b) %*% X) %*% t(X) %*% diag(wts_b) %*% y
    stat_stor[b] = sum(y * wts_b) + t(c(length(y), sum(x_u) - sum(x_s * wts_b))) %*% betas
  }
  
  W = diag(wts)
  betas = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
  est_stat = sum(y * wts) + t(c(length(y), sum(x_u) - sum(x_s * wts))) %*% betas
  p_o = mean(ifelse(stat_stor >= rep(est_stat, B), 1, 0))
  p_value = 2 * min(p_o, 1 - p_o)
  return(p_value)
}

#perm_greg(y = samp$y, x_s = samp$x, x_u = pop$x, wts = samp$w, B = 10000)