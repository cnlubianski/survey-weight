# ------------------------------------------------------------------------------
# Title: Replication of Wang et al. Study 2
# Description: This script runs the second simulation for the replication 
# chapter which tests the performance of the diagnostic tests with regards 
# to different levels of informative weighting and quadratic informative weighting.
# ------------------------------------------------------------------------------

# Set working directory and define script/storage paths
working_dir <- getwd() # Make sure you are in root directory
sim_script_dir <- file.path(working_dir, "R")
sim_storage_dir <- file.path(working_dir, "output")

# Define simulation attributes
set.seed(1337)
sim_rep <- 1000 # How many times a specific case is run
critical_value <- 0.05

# Define cases and construct case grid
N <- 3000
n <- c(100, 200)
sigma <- c(0.1)
alpha <- c(0, 0.5, 1.0, 1.5)
cases <- tidyr::expand_grid(N, n, sigma, alpha)

# General simulation functions -------------------------------------------------
# Retrieve simulation utility functions
source(file.path(sim_script_dir, "simulation-utilities.R"))

# One replication for a given case
run_one_rep <- function(N, n, sigma, alpha) {
  # Generate population and sample
  pop <- generate_data_study2(N = N, sigma = sigma, alpha = alpha)
  samp <- generate_sample_brewer(pop, w = pop$w, n = n)
  
  # Construct svyglm model object for diagnostic functions
  design <- svydesign(id = ~1, weights = ~w, data = samp)
  svyglm_model <- svyglm(y ~ x, design = design)
  
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
  results <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(results) <- columns
  
  for (case in seq_len(nrow(cases))) {
    case_storage <- data.frame(iteration = seq_len(B))
    
    for (b in seq_len(B)) {
      res <- run_one_rep(
        N     = cases$N[case],
        n     = cases$n[case],
        sigma = cases$sigma[case],
        alpha = cases$alpha[case]
      )
      case_storage[b, names(res)] <- unlist(res)
    }
    
    results <- rbind(results, cbind(case = case, case_storage))
    print(case)
  }
  
  return(results)
}

# Run simulation, wrangle results into table, then save table ------------------
sim_results <- run_simulation(cases, B = sim_rep)

# Summarize rejection rates
reject_table <- sim_results %>%
  mutate(across(c(DD, PS1, PS1q, PS2, PS2q, WF, HP, PS3),
                ~ ifelse(. <= critical_value, 1, 0)),
         case = rep(seq_len(nrow(cases)), each = sim_rep)) %>%
  group_by(case) %>%
  summarize(across(-row, mean)) %>%
  left_join(cases %>% mutate(case = row_number()), by = "case") %>%
  select(n, sigma, delta, alpha, DD, PS1, PS1q, PS2, PS2q, WF, HP, PS3)

# Save the results
saveRDS(reject_table, file.path(sim_storage_dir, "replication_reject_table2.rds"))
write.csv(reject_table, file.path(sim_storage_dir, "replication_reject_table2.csv"), row.names = FALSE)
