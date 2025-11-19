# ------------------------------------------------------------------------------
# Title: Permutation Tests Simulation 2
# Description: This script runs the second simulation for the permutation tests 
# chapter which tests the performance of the two permutation tests with regards 
# to different error distributions.
# ------------------------------------------------------------------------------

# Preamble ---------------------------------------------------------------------
# Load in libraries
library(dplyr)
library(sampling)
library(survey)
library(tidyr)
library(future.apply) # for parallelization
library(svytest)

# Set working directory and define script/storage paths
working_dir <- getwd() # Make sure you are in root directory
sim_script_dir <- file.path(working_dir, "R")
sim_storage_dir <- file.path(working_dir, "output")

# Define simulation attributes
set.seed(1337)
sim_rep <- 1000 # How many times a specific case is run
func_rep <- 1000 # how many times a permutation statistic is computed within a function run
critical_value <- 0.05

# Define cases and construct case grid
N <- 3000
n <- c(50, 100)
sigma <- 0.2
delta <- 1
alpha <- c(0, 0.2, 0.4, 0.6)
distribution <- c("Normal", "Uniform", "Gamma", "t")
cases <- tidyr::expand_grid(N, n, distribution, alpha)

# General simulation functions -------------------------------------------------
# Retrieve simulation utility functions
source(file.path(sim_script_dir, "simulation-utilities.R"))

# Helper to safely extract p-values
safe_pval <- function(expr) {
  tryCatch(
    expr,
    error = function(e) NA,
    warning = function(w) {
      # suppress warnings but still return result if possible
      invokeRestart("muffleWarning")
    }
  )
}

# Function to run one simulation replication
run_one_rep <- function(N, n, alpha, dist) {
  # Generate population and sample
  pop <- generate_data_study2_perm(N, sigma, delta, alpha, distribution = dist)
  samp <- generate_sample_brewer(pop, w = pop$w, n = n)
  
  # Construct svyglm model object for diagnostic functions
  design <- svydesign(id = ~1, weights = ~w, data = samp)
  svyglm_model <- svyglm(y ~ x, design = design)
  
  # Call the package functions and return the p-values (safe)
  out <- list(
    DD        = safe_pval(wa_test(model = svyglm_model, type = "DD")$p.value),
    PS1       = safe_pval(wa_test(model = svyglm_model, type = "PS1")$p.value),
    PS2       = safe_pval(wa_test(model = svyglm_model, type = "PS2")$p.value),
    WF        = safe_pval(wa_test(model = svyglm_model, type = "WF")$p.value),
    HP        = safe_pval(diff_in_coef_test(model = svyglm_model, var_equal = TRUE)$p.value),
    PS3       = safe_pval(estim_eq_test(model = svyglm_model, q_method = "linear")$p.value),
    pred_mean = safe_pval(perm_test(model = svyglm_model, stat = "pred_mean", 
                                    B = func_rep, engine = "R")$p.value),
    coef_mahal= safe_pval(perm_test(model = svyglm_model, stat = "coef_mahal", 
                                    B = func_rep, engine = "R")$p.value)
  )
  
  return(out)
}

# Function to perform the simulation
simulate_study <- function(cases, B = sim_rep) {
  # Initiate a storage object
  results <- list()
  
  # Iterate through each case
  for (case in 1:nrow(cases)) {
    param <- cases[case, ]
    message("Running case ", case, " of ", nrow(cases))
    
    # Parallel replications
    reps <- future_lapply(seq_len(B), function(b) {
      run_one_rep(param$N, param$n, param$alpha, param$distribution)
    }, future.seed = TRUE)
    
    case_df <- tibble(row = seq_len(B), data = reps) |> unnest_wider(data)
    results[[case]] <- case_df
  }
  
  results <- bind_rows(results)
  return(results)
}

# Run simulation, wrangle results into table, then save table ------------------
res <- simulate_study(cases, B = sim_rep)

# Summarize rejection rates
reject_table <- res %>%
  mutate(across(c(DD, PS1, PS2, WF, HP, PS3, pred_mean, coef_mahal),
                ~ ifelse(. <= critical_value, 1, 0)),
         case = rep(seq_len(nrow(cases)), each = sim_rep)) %>%
  group_by(case) %>%
  summarize(across(-row, mean)) %>%
  left_join(cases %>% mutate(case = row_number()), by = "case") %>%
  select(distribution, n, alpha, DD, HP, PS1, PS2, PS3, pred_mean, coef_mahal) |> 
  arrange(distribution, n, alpha)

# Save the results
saveRDS(reject_table, file.path(sim_storage_dir, "perm_reject_table2.rds"))
write.csv(reject_table, file.path(sim_storage_dir, "perm_reject_table2.csv"), row.names = FALSE)
