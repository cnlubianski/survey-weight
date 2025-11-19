# Simulation: Wang et al.

set.seed(51483464)
B = 1000

N <- 3000
n <- c(25, 50, 75, 100)
sigma <- 0.2 #c(0.1, 0.2)
alpha <- c(0, 0.2, 0.4, 0.6)
delta <- 1 # c(1.5, 1)
cases <- expand_grid(N, n, sigma, delta, alpha)

columns = c("case", "iteration", "perm_HP", "HP", "perm_DD", "DD",
            "perm_PS1", "PS1", "p_PS2", "PS2", "greg") ###
mini_results = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(mini_results) = columns

for (case in 1:nrow(cases)) {
  p_HP = HP = p_DD = DD = p_PS1 = PS1 = p_PS2 = PS2 = greg = rep(NA, B) ###
  case_storage = data.frame(iteration = seq_len(B), p_HP, HP, p_DD, DD, p_PS1, PS1, p_PS2, PS2, greg) ###
  for (b in 1:B) {
    pop = generate_data_study1_perm(N = cases$N[case], 
                                    sigma = cases$sigma[case],
                                    alpha = cases$alpha[case],
                                    delta = cases$delta[case])
    samp = generate_sample_brewer(pop, w = pop$w, n = cases$n[case])
    
    case_storage$p_HP[b] = perm_HP(y = samp$y, x = samp$x, wts = samp$w, B = 1000, replacement = FALSE)
    case_storage$HP[b] = HP_DC_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$p_DD[b] = perm_DD(y = samp$y, x = samp$x, wts = samp$w, B = 1000, replacement = FALSE)
    case_storage$DD[b] = DD_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$p_PS1[b] = perm_PS1(y = samp$y, x = samp$x, wts = samp$w, B = 1000, replacement = FALSE)
    case_storage$PS1[b] = PS1_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$p_PS2[b] = perm_PS2(y = samp$y, x = samp$x, wts = samp$w, B = 1000, replacement = FALSE)
    case_storage$PS2[b] = PS2_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$greg[b] = perm_greg(y = samp$y, x_s = samp$x, x_u = pop$x, wts = samp$w, B = 1000, replacement = FALSE)
  }
  
  mini_results = rbind(mini_results, cbind(case, case_storage))
  print(case)
}

write.csv(mini_results, "perm_sim.csv")

reject_mini = mini_results %>% 
  mutate(p_HP = case_when(p_HP <= 0.05 ~ 1, TRUE ~ 0),
         HP = case_when(HP <= 0.05 ~ 1, TRUE ~ 0),
         p_DD = case_when(p_DD <= 0.05 ~ 1, TRUE ~ 0),
         DD = case_when(DD <= 0.05 ~ 1, TRUE ~ 0),
         p_PS1 = case_when(p_PS1 <= 0.05 ~ 1, TRUE ~ 0),
         PS1 = case_when(PS1 <= 0.05 ~ 1, TRUE ~ 0),
         p_PS2 = case_when(p_PS2 <= 0.05 ~ 1, TRUE ~ 0),
         PS2 = case_when(PS2 <= 0.05 ~ 1, TRUE ~ 0),
         greg = case_when(greg <= 0.05 ~ 1, TRUE ~ 0)) %>%
  select(-iteration) %>%
  group_by(case) %>%
  summarize(across(everything(), mean)) %>%
  mutate(p_HP = format(round(p_HP * 100, 1), nsmall = 1),
         HP = format(round(HP * 100, 1), nsmall = 1),
         p_DD = format(round(p_DD * 100, 1), nsmall = 1),
         DD = format(round(DD * 100, 1), nsmall = 1),
         p_PS1 = format(round(p_PS1 * 100, 1), nsmall = 1),
         PS1 = format(round(PS1 * 100, 1), nsmall = 1),
         p_PS2 = format(round(p_PS2 * 100, 1), nsmall = 1),
         PS2 = format(round(PS2 * 100, 1), nsmall = 1),
         greg = format(round(greg * 100, 1), nsmall = 1))

reject_mini_table = cbind(cases, reject_mini) %>% select(-c(N, case, delta))
reject_mini_table

write.csv(reject_mini_table, "perm_table.csv")


# Simulation: Error Distribution

set.seed(51483464)
B = 1000

N <- 3000
distribution <- c("Normal", "Uniform", "Gamma", "t")
n <- c(50, 100)
alpha <- c(0, 0.2, 0.4, 0.6)
cases <- expand_grid(N, distribution, n, alpha)

columns = c("case", "iteration", "perm_HP", "HP", "perm_DD", "DD",
            "perm_PS1", "PS1", "p_PS2", "PS2", "greg") ###
mini_results = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(mini_results) = columns

for (case in 1:nrow(cases)) {
  p_HP = HP = p_DD = DD = p_PS1 = PS1 = p_PS2 = PS2 = greg = rep(NA, B) ###
  case_storage = data.frame(iteration = seq_len(B), p_HP, HP, p_DD, DD, p_PS1, PS1, p_PS2, PS2, greg) ###
  for (b in 1:B) {
    pop = generate_data_study2_perm(N = cases$N[case], 
                                    distribution = cases$distribution[case],
                                    alpha = cases$alpha[case])
    samp = generate_sample_brewer(pop, w = pop$w, n = cases$n[case])
    
    case_storage$p_HP[b] = perm_HP(y = samp$y, x = samp$x, wts = samp$w,
                                   B = 1000, replacement = FALSE)
    case_storage$HP[b] = HP_DC_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$p_DD[b] = perm_DD(y = samp$y, x = samp$x, wts = samp$w,
                                   B = 1000, replacement = FALSE)
    case_storage$DD[b] = DD_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$p_PS1[b] = perm_PS1(y = samp$y, x = samp$x, wts = samp$w,
                                     B = 1000, replacement = FALSE)
    case_storage$PS1[b] = PS1_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$p_PS2[b] = perm_PS2(y = samp$y, x = samp$x, wts = samp$w,
                                     B = 1000, replacement = FALSE)
    case_storage$PS2[b] = PS2_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$greg[b] = perm_greg(y = samp$y, x_s = samp$x, x_u = pop$x, wts = samp$w,
                                     B = 1000, replacement = FALSE)
  }
  
  mini_results = rbind(mini_results, cbind(case, case_storage))
  print(case)
}

write.csv(mini_results, "perm_sim_error.csv")

reject_mini = mini_results %>% ###
  mutate(p_HP = case_when(p_HP <= 0.05 ~ 1, TRUE ~ 0),
         HP = case_when(HP <= 0.05 ~ 1, TRUE ~ 0),
         p_DD = case_when(p_DD <= 0.05 ~ 1, TRUE ~ 0),
         DD = case_when(DD <= 0.05 ~ 1, TRUE ~ 0),
         p_PS1 = case_when(p_PS1 <= 0.05 ~ 1, TRUE ~ 0),
         PS1 = case_when(PS1 <= 0.05 ~ 1, TRUE ~ 0),
         p_PS2 = case_when(p_PS2 <= 0.05 ~ 1, TRUE ~ 0),
         PS2 = case_when(PS2 <= 0.05 ~ 1, TRUE ~ 0),
         greg = case_when(greg <= 0.05 ~ 1, TRUE ~ 0)) %>%
  select(-iteration) %>%
  group_by(case) %>%
  summarize(across(everything(), mean)) %>%
  mutate(p_HP = format(round(p_HP * 100, 1), nsmall = 1),
         HP = format(round(HP * 100, 1), nsmall = 1),
         p_DD = format(round(p_DD * 100, 1), nsmall = 1),
         DD = format(round(DD * 100, 1), nsmall = 1),
         p_PS1 = format(round(p_PS1 * 100, 1), nsmall = 1),
         PS1 = format(round(PS1 * 100, 1), nsmall = 1),
         p_PS2 = format(round(p_PS2 * 100, 1), nsmall = 1),
         PS2 = format(round(PS2 * 100, 1), nsmall = 1),
         greg = format(round(greg * 100, 1), nsmall = 1))

reject_mini_table = cbind(cases, reject_mini) %>% select(-c(N, case))
reject_mini_table

write.csv(reject_mini_table, "perm_table_error.csv")