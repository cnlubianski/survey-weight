set.seed(51483464)
B = 5000

n <- c(50, 100, 250, 500, 1000)
methods <- c("grouping", "pps", "stratify", "cluster", "twostage")
cases <- expand_grid(methods, n)

columns = c("case", "iteration", "DD", "PN", "HP", "PS1", "PS1q", "PS2", "PS2q", 
            "PS3", "WF", "LR")
storage = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(storage) = columns

for (case in 1:nrow(cases)) {
  DD = PN = HP = PS1 = PS1q = PS2 = PS2q = PS3 = WF = LR = rep(NA, B)
  case_storage = data.frame(iteration = seq_len(B), DD, PN, HP, PS1, PS1q,
                            PS2, PS2q, PS3, WF, LR)
  
  for (b in 1:B) {
    
    sampling = NULL
    if (cases$methods[case] == "grouping") { # Bad
      sampling = grouping(x = ce$TOTEXPCQ, n = cases$n[case], strata_prob = c(0.15, 0.2, 0.25, 0.4))
    } 
    if (cases$methods[case] == "pps") { # Decent - converges only with some tests
      sampling = pps(x = ce$TOTEXPCQ, n = cases$n[case], noise_sd = 0.025)
    } 
    if (cases$methods[case] == "stratify") { # Good
      sampling = stratify(stratum = ce$NO_EARNR, nh = round(cases$n[case] * c(0.4, 0.35, 0.15, 0.10)))
    }
    if (cases$methods[case] == "cluster") { # Good
      #n = ifelse(cases$n[case] %in% c(50, 100), 2, 3) # Increase m per n
      sampling = clustering(clusters = ce$INCOMEY, n = 3, m = cases$n[case] / 3)
    }
    if (cases$methods[case] == "twostage") { # Very good!
      n = ifelse(cases$n[case] %in% c(50, 100), 2, 3) # Increase n_I per M
      sampling = three_stage_clust_strat(clusters = ce$REGION, stratum = ce$MARITAL,
                                         n = n, sample_size = cases$n[case])
    }
    samp = ce[1:nrow(ce) * sampling$sample,] %>%
      mutate(wts = 1 / sampling$pik[1:nrow(ce) * sampling$sample])
    
    case_storage$HP[b] = HP_DC_test(data = samp, y = samp$FINCBTAX, x = samp$TOTEXPCQ, wts = samp$wts)
    case_storage$DD[b] = DD_WA_test(data = samp, y = samp$FINCBTAX, x = samp$TOTEXPCQ, wts = samp$wts)
    case_storage$PS1[b] = PS1_WA_test(data = samp, y = samp$FINCBTAX, x = samp$TOTEXPCQ, wts = samp$wts)
    case_storage$PS2[b] = PS2_WA_test(data = samp, y = samp$FINCBTAX, x = samp$TOTEXPCQ, wts = samp$wts)
    case_storage$PS3[b] = PS3_test(data = samp, y = samp$FINCBTAX, x = samp$TOTEXPCQ, wts = samp$wts)
    case_storage$WF[b] = WF_WA_test(data = samp, y = samp$FINCBTAX, x = samp$TOTEXPCQ, wts = samp$wts)
    case_storage$PN[b] = PN_test(data = samp, y = samp$FINCBTAX, x = samp$TOTEXPCQ, wts = samp$wts)
    case_storage$PS1q[b] = PS1q_WA_test(data = samp, y = samp$FINCBTAX, x = samp$TOTEXPCQ, wts = samp$wts)
    case_storage$PS2q[b] = PS2q_WA_test(data = samp, y = samp$FINCBTAX, x = samp$TOTEXPCQ, wts = samp$wts)
    case_storage$LR[b] = LR_test(data = samp, y = samp$FINCBTAX, x = samp$TOTEXPCQ, wts = samp$wts)
    
  }
  storage = rbind(storage, cbind(case, case_storage))
  print(case)
}
write.csv(storage, "ce_results.csv")


test_results = read.csv("ce_results.csv")

test_reject = test_results %>%
  mutate(HP = case_when(HP <= 0.05 ~ 1, TRUE ~ 0),
         DD = case_when(DD <= 0.05 ~ 1, TRUE ~ 0),
         PS1 = case_when(PS1 <= 0.05 ~ 1, TRUE ~ 0),
         PS1q = case_when(PS1q <= 0.05 ~ 1, TRUE ~ 0),
         PS2 = case_when(PS2 <= 0.05 ~ 1, TRUE ~ 0),
         PS2q = case_when(PS2q <= 0.05 ~ 1, TRUE ~ 0),
         PS3 = case_when(PS3 <= 0.05 ~ 1, TRUE ~ 0),
         WF = case_when(WF <= 0.05 ~ 1, TRUE ~ 0),
         LR = case_when(LR <= 0.05 ~ 1, TRUE ~ 0),
         PN = case_when(PN <= 0.05 ~ 1, TRUE ~ 0)) %>%
  select(-iteration) %>%
  group_by(case) %>%
  summarize(across(everything(), mean)) %>%
  mutate(HP = format(round(HP * 100, 1), nsmall = 1),
         DD = format(round(DD * 100, 1), nsmall = 1),
         PS1 = format(round(PS1 * 100, 1), nsmall = 1),
         PS1q = format(round(PS1q * 100, 1), nsmall = 1),
         PS2 = format(round(PS2 * 100, 1), nsmall = 1),
         PS2q = format(round(PS2q * 100, 1), nsmall = 1),
         PS3 = format(round(PS3 * 100, 1), nsmall = 1),
         WF = format(round(WF * 100, 1), nsmall = 1),
         LR = format(round(LR * 100, 1), nsmall = 1),
         PN = format(round(PN * 100, 1), nsmall = 1))

test_reject_table = cbind(cases, test_reject) %>% select(-c(X, case, PN, LR))
test_reject_table

write.csv(test_reject_table, "ce_reject_table.csv")