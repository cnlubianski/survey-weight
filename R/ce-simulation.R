# ------------------------------------------------------------------------------
# Title: Consumer Expenditure Sampling Design Simulation
# Description: This script runs the simulation using Consumer Expenditure data
# to test the performance of the survey weight diagnostic tests with different 
# sampling designs.
# ------------------------------------------------------------------------------

# Preamble ---------------------------------------------------------------------
# Load in libraries
library(dplyr)
library(sampling)
library(survey)
library(tidyr)
library(future.apply) # for parallelization
library(svytest)
library(rpms)

# Set working directory and define script/storage paths
working_dir <- getwd() # Make sure you are in root directory
sim_script_dir <- file.path(working_dir, "R")
sim_storage_dir <- file.path(working_dir, "output")

# Define simulation attributes
set.seed(1337)
sim_rep <- 10000 # How many times a specific case is run
critical_value <- 0.05

# Define cases and construct case grid
n <- c(50, 100, 250, 500, 1000)
methods <- c("grouping", "pps", "stratify", "cluster", "twostage")
cases <- tidyr:expand_grid(methods, n)

# Generate a subset of the CE data for the simulation
ce <- rpms::CE %>%
  filter(TOTEXPCQ > 0, FINCBTAX > 10, SALARYX > 0, !is.na(REGION), 
         FAM_SIZE %in% factor(1:10), ROOMSQ %in% factor(1:11), NO_EARNR %in% factor(1:4)) %>%
  mutate(TOTEXPCQ = log(TOTEXPCQ), FINCBTAX = log(FINCBTAX)) %>%
  select(-c(QINTRVMO, PSU, INCNONWK, IRAX, LIQUIDX, STOCKX, STUDNTX,
            FOOTWRCQ, TOBACCCQ, TOTXEST, VEHQL, EARNER)) %>%
  group_by(REGION, MARITAL) %>%
  filter(n() >= 70) %>%
  ungroup()


# One replication for a given case
run_one_rep <- function(method, n) {
  # Choose sampling scheme
  sampling <- switch(method,
    "grouping"  = grouping(x = ce$TOTEXPCQ, n = n, strata_prob = c(0.15, 0.2, 0.25, 0.4)),
    "pps"       = pps(x = ce$TOTEXPCQ, n = n, noise_sd = 0.025),
    "stratify"  = stratify(stratum = ce$NO_EARNR, nh = round(n * c(0.4, 0.35, 0.15, 0.10))),
    "cluster"   = clustering(clusters = ce$INCOMEY, n = 3, m = n / 3),
    "twostage"  = {
      n_I <- ifelse(n %in% c(50, 100), 2, 3)
      three_stage_clust_strat(clusters = ce$REGION, stratum = ce$MARITAL,
                              n = n_I, sample_size = n)
    },
    stop("Unknown method")
  )
  
  # Construct sample
  samp <- ce[1:nrow(ce) * sampling$sample, ] %>%
    mutate(wts = 1 / sampling$pik[1:nrow(ce) * sampling$sample])
  
  # Construct svyglm model object for diagnostic functions
  design <- svydesign(id = ~1, weights = ~wts, data = samp)
  svyglm_model <- svyglm(FINCBTAX ~ TOTEXPCQ, design = design)
  
  # Run tests
  results <- list(
    DD   = wa_test(model = svyglm_model, type = "DD")$p.value,
    HP   = diff_in_coef_test(model = svyglm_model, var_equal = TRUE)$p.value,
    PS1  = wa_test(model = svyglm_model, type = "PS1q")$p.value,
    PS1q  = wa_test(model = svyglm_model, type = "PS1q")$p.value,
    PS2  = wa_test(model = svyglm_model, type = "PS2")$p.value,
    PS2q  = wa_test(model = svyglm_model, type = "PS2q")$p.value,
    PS3  = estim_eq_test(model = svyglm_model, q_method = "linear")$p.value,
    WF   = wa_test(model = svyglm_model, type = "WF")$p.value
  )
  
  return(results)
}

# Simulation driver
run_simulation <- function(cases, B = 1000) {
  columns <- c("case", "iteration", "DD", "HP", "PS1", "PS1q", "PS2", "PS2q", "PS3", "WF")
  storage <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(storage) <- columns
  
  for (case in seq_len(nrow(cases))) {
    case_storage <- data.frame(iteration = seq_len(B))
    
    for (b in seq_len(B)) {
      results <- run_one_rep(cases$methods[case], cases$n[case])
      case_storage[b, names(results)] <- unlist(results)
    }
    
    storage <- rbind(storage, cbind(case = case, case_storage))
    print(case)
  }
  
  return(storage)
}

# Run simulation
sim_results <- run_simulation(cases, B = sim_rep)

# Summarize rejection rates
reject_table <- sim_results %>%
  mutate(across(c(DD, PS1, PS1q, PS2, PS2q, WF, HP, PS3),
                ~ ifelse(. <= critical_value, 1, 0)),
         case = rep(seq_len(nrow(cases)), each = sim_rep)) %>%
  group_by(case) %>%
  summarize(across(-row, mean)) %>%
  left_join(cases %>% mutate(case = row_number()), by = "case") %>%
  select(n, sigma, delta, alpha, DD, HP, PS1, PS1q, PS2, PS2q, PS3)

# Save the results
saveRDS(reject_table, file.path(sim_storage_dir, "ce_reject_table.rds"))
write.csv(reject_table, file.path(sim_storage_dir, "ce_reject_table.csv"), row.names = FALSE)