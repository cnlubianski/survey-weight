# ------------------------------------------------------------------------------
# Title: Sampling Utility Functions
# General functions to be sourced in by the simulation scripts. These include
# functions to generate different inclusion probabilities with regards to the 
# different sampling methods. Functions return a list of inclusion probabilities
# and selected observations.
# ------------------------------------------------------------------------------

# Function to select samples using quantiles of some variable X and probability per group
grouping <- function(x, n, strata_prob) {
  N_h <- length(strata_prob) # Number of stratum
  quantiles <- quantile(x, probs = seq(0, 1, 1 / N_h), na.rm = TRUE)
  groups <- cut(x, breaks = quantiles, labels = FALSE, include.lowest = TRUE)
  pik <- inclusionprobastrata(strata = groups, nh = ceiling(n * strata_prob))
  selected <- UPbrewer(pik)
  return(list(pik = pik, sampled = selected))
}

# Function to select samples using Probability Proportional to Size (PPS)
pps <- function(x, n, noise_sd) {
  noise <- rnorm(length(x), 0, noise_sd)
  x_noisy <- x * (1 + noise)
  pik <- inclusionprobabilities(x_noisy, n)
  selected <- UPbrewer(pik)
  return(list(pik = pik, sampled = selected))
}

# Function to select samples using Stratification (STSRS) by factor levels
stratify <- function(stratum, nh) {
  pik <- inclusionprobastrata(as.numeric(stratum), nh)
  sampled <- rep(0, length(stratum))
  sampled[unlist(sapply(levels(as.factor(stratum)), function(s) {
    indices <- which(stratum == s)
    selected <- sample(indices, size = nh[as.numeric(s)], replace = FALSE)
    return(selected)
  }))] <- 1
  return(list(pik = pik, sampled = sampled))
}

# Function to select samples using Two-stage Clustering (CL-SRS)
clustering <- function(clusters, n, m) {
  selected_clusters <- sample(x = unique(clusters), size = n, replace = FALSE)
  
  sampled <- rep(0, length(clusters))
  for (c in selected_clusters) {
    indices <- which(clusters == c)
    selected <- sample(indices, size = m, replace = FALSE)
    sampled[selected] <- 1
  }
  
  clust_pik <- data.frame(clust = clusters) %>%
    group_by(clust) %>%
    summarize(pik = n / length(unique(clusters)) * m / n())
  
  pik_list <- left_join(data.frame(clusters = clusters),
                       mutate(clust_pik, clusters = clust),
                       by = "clusters")
  
  return(list(pik = pik_list$pik, sampled = sampled))
}

# Function to select samples using Three stage sampling method 
# by Two-stage Clustering then Stratified (CL-ST-SRS)
three_stage_clust_strat <- function(clusters, stratum, n, sample_size) {
  selected_clusters <- sample(x = unique(clusters), size = n, replace = FALSE)
  clust_indices <- which(clusters %in% selected_clusters)
  stratum_in_clust <- paste0(clusters[clust_indices], "_", stratum[clust_indices])
  threestage <- data.frame(clusters = clusters, stratum = stratum,
                          index = 1:length(clusters))
  
  m <- round(sample_size / length(unique(stratum_in_clust)), 0)
  sampled <- threestage %>%
    filter(clusters %in% selected_clusters) %>%
    group_by(clusters, stratum) %>%
    sample_n(size = m, replace = FALSE) %>%
    ungroup()
  
  selected <- rep(0, length(clusters))
  selected[sampled$index] <- 1
  
  stratum_pik <- data.frame(clust = clusters, strat = stratum) %>%
    group_by(clust, strat) %>%
    summarize(pik = (n / length(unique(clusters))) * (m / n()))
  
  pik_list <- left_join(threestage, stratum_pik,
                       by = c("clusters" = "clust", "stratum" = "strat")) %>%
    mutate(pik = pik * sum(selected) / sum(pik)) # pik does not always equal n, so scale
  
  return(list(pik = pik_list$pik, sampled = selected))
}
