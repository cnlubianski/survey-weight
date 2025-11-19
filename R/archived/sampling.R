# Packages
if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if (!require(rpms)) install.packages("rpms"); library(rpms)
if (!require(sampling)) install.packages("sampling"); library(sampling)
if (!require(survey)) install.packages("survey"); library(survey)

# Data Wrangling
ce = rpms::CE %>%
  filter(TOTEXPCQ > 0, FINCBTAX > 10, SALARYX > 0, !is.na(REGION), 
         FAM_SIZE %in% factor(1:10), ROOMSQ %in% factor(1:11), NO_EARNR %in% factor(1:4)) %>%
  mutate(TOTEXPCQ = log(TOTEXPCQ), FINCBTAX = log(FINCBTAX)) %>% #,
  #REG_POP = paste0(REGION, "_", POPSIZE)) %>%
  select(-c(QINTRVMO, PSU, INCNONWK, IRAX, LIQUIDX, STOCKX, STUDNTX,
            FOOTWRCQ, TOBACCCQ, TOTXEST, VEHQL, EARNER)) %>%
  group_by(REGION, MARITAL) %>%
  filter(n() >= 70) %>%
  ungroup()

### Grouping
grouping <- function(x, n, strata_prob) {
  N_h = length(strata_prob) # Number of stratum
  quantiles <- quantile(x, probs = seq(0, 1, 1 / N_h), na.rm = TRUE)
  groups <- cut(x, breaks = quantiles, labels = FALSE, include.lowest = TRUE)
  pik = inclusionprobastrata(strata = groups,
                             nh = ceiling(n * strata_prob))
  selected = UPbrewer(pik)
  return(list(pik = pik, sampled = selected))
}

#sample_size = 100
#strata_prob = c(0.1, 0.15, 0.25, 0.5)
#grouping(x = ce$FINCBTAX, n = sample_size, strata_prob = strata_prob)


### Probability Proportional to Size (PPS)
pps <- function(x, n, noise_sd) {
  noise = rnorm(length(x), 0, noise_sd)
  x_noisy = x * (1 + noise)
  pik = inclusionprobabilities(x_noisy, n)
  selected = UPbrewer(pik)
  return(list(pik = pik, sampled = selected))
}
#sampling = pps(x = exp(ce$TOTEXPCQ), n = 50, noise_sd = 0.025)


### Stratification (STSRS) 
stratify <- function(stratum, nh) {
  pik = inclusionprobastrata(as.numeric(stratum), nh)
  sampled <- rep(0, length(stratum))
  sampled[unlist(sapply(levels(as.factor(stratum)), function(s) {
    indices <- which(stratum == s)
    selected <- sample(indices, size = nh[as.numeric(s)], replace = FALSE)
    return(selected)
  }))] <- 1
  return(list(pik = pik, sampled = sampled))
}

#stratify(stratum = ce$REGION, nh = c(150, 150, 500, 200))


### Two-stage Clustering (CL-SRS)
clustering <- function(clusters, n, m) {
  selected_clusters <- sample(x = unique(clusters), size = n, replace = FALSE)
  
  sampled <- rep(0, length(clusters))
  for (c in selected_clusters) {
    indices <- which(clusters == c)
    selected <- sample(indices, size = m, replace = FALSE)
    sampled[selected] <- 1
  }
  
  clust_pik = data.frame(clust = clusters) %>%
    group_by(clust) %>%
    summarize(pik = n / length(unique(clusters)) * m / n())
  
  pik_list = left_join(data.frame(clusters = clusters),
                       mutate(clust_pik, clusters = clust),
                       by = "clusters")
  
  return(list(pik = pik_list$pik, sampled = sampled))
}

#clustering(clusters = ce$INCOMEY, n = 2, m = 500 / 2)


### Two-stage Clustering - Stratified (CL-ST-SRS)
three_stage_clust_strat <- function(clusters, stratum, n, sample_size) {
  selected_clusters = sample(x = unique(clusters), size = n, replace = FALSE)
  clust_indices = which(clusters %in% selected_clusters)
  stratum_in_clust = paste0(clusters[clust_indices], "_", stratum[clust_indices])
  threestage = data.frame(clusters = clusters, stratum = stratum,
                          index = 1:length(clusters))
  
  m = round(sample_size / length(unique(stratum_in_clust)), 0)
  sampled = threestage %>%
    filter(clusters %in% selected_clusters) %>%
    group_by(clusters, stratum) %>%
    sample_n(size = m, replace = FALSE) %>%
    ungroup()
  
  selected = rep(0, length(clusters))
  selected[sampled$index] <- 1
  
  stratum_pik = data.frame(clust = clusters, strat = stratum) %>%
    group_by(clust, strat) %>%
    summarize(pik = (n / length(unique(clusters))) * (m / n()))
  
  pik_list = left_join(threestage, stratum_pik,
                       by = c("clusters" = "clust", "stratum" = "strat")) %>%
    mutate(pik = pik * sum(selected) / sum(pik)) # pik does not always equal n, so scale
  
  return(list(pik = pik_list$pik, sampled = selected))
}

#three_stage_clust_strat(clusters = ce$REGION, stratum = ce$MARITAL, n = 3, sample_size = 1000)


