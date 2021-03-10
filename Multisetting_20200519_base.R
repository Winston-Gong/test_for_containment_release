###############################################
# Simulation starts here
###############################################

#simulate Scenario 1: do nothing
all_new_daily <- matrix(nrow = num_rep, ncol = len_sim)
Rt <- matrix(nrow = num_rep, ncol = len_sim)
for(rp in 1:num_rep){
  nw_para <- load_sample_nw(setting,v_id,s_id, para)
  est1 <- nw_para[[1]]
  para <- nw_para[[2]]
  
  print(paste("Running Setting=",setting,"V=",v_id,"S=",s_id,", Scenario 1, Rep: ",rp))
  # C is the contact matrix, which is traced for 7 days. If i,j element is 7, it means the latest contact of i,j is today,
  # if it is 0, then there is no contact between i,j in the past 7 days,
  C <- matrix(0, para$pop_sz,para$pop_sz)
  
  # I is a vector indicating when an individual has been infected
  # NA means susceptible, 0 means infected date, 1 indicates the person has been infected 1 day
  I <- matrix(NA, para$pop_sz,1)
  # Z is a vector indicating if the infected is detectable
  Z <- matrix(F, para$pop_sz,1)
  ############################################
  ############################################
  trace_inf_n <- matrix(0, para$pop_sz,1)
  
  # onset just save the incubation time 
  O <- matrix(NA, para$pop_sz,1)
  
  # Q is the vector indicating when an individual has been quarantine, NA means not quarantine, 1 means first day
  Q <- matrix(NA, para$pop_sz,1)
  
  C_lst <- list()
  O_lst <- list()
  I_lst <- list()
  # initial case
  # E_0 <- floor(runif(1,min = 1, max = 3))
  init_idx <- sample(1:para$pop_sz, para$E_0)
  I[init_idx] <- 0
  O[init_idx] <- incub(para$E_0)
  Z[init_idx] <- F # initial case is never detectable
  for(t in 1:len_sim){
    C_lst[[t]] <- C
    C <- C_update(C, Q, est1, para)
    lst <- I_O_update(I, Q, C, O, Z, trace_inf_n, para)
    I <- lst[[1]]
    O <- lst[[2]]
    Z <- lst[[3]]
    trace_inf_n <- lst[[4]]
    I_lst[[t]] <- I
    O_lst[[t]] <- O
  }
  # Rt
  rt <- rep(NA, len_sim)
  for(t in 1:len_sim){
    i <- I_lst[[t]]
    idx <- which(i == 2) # from the second day they are infectious
    if(length(idx)) rt[t] <- mean(trace_inf_n[idx])
  }
  # plot the daily new case
  new_daily <- rep(NA, len_sim)
  new_daily[1] <- para$E_0
  for(t in 2:len_sim){
    new_daily[t] <- sum(is.na(I_lst[[t-1]]))  - sum(is.na(I_lst[[t]]))
  }
  all_new_daily[rp,] <- new_daily
  Rt[rp,] <- rt
}

longData <- melt(C_lst[[8]])
longData<-longData[longData$value!=0,]

# ggplot(longData, aes(x = Var2, y = Var1)) + 
#   geom_raster(aes(fill=value)) + 
#   scale_fill_gradient(low="grey90", high="red") +
#   labs(x="letters", y="LETTERS", title="Matrix") +
#   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
#                      axis.text.y=element_text(size=9),
#                      plot.title=element_text(size=11))

all_Q_new_daily <- matrix(nrow = num_rep, ncol = len_sim)
cap_Q <- matrix(nrow = num_rep, ncol = len_sim) # store the number of max quarantine
Rt_Q <- matrix(nrow = num_rep, ncol = len_sim)
PCR_Q <-matrix(nrow = num_rep, ncol = len_sim)

# Simulate Scenario 2: symptom & PCR
for(rp in 1:num_rep){
  print(paste("Running Setting=",setting,"V=",v_id,"S=",s_id,", Scenario 2, Rep: ",rp))
  nw_para <- load_sample_nw(setting,v_id,s_id, para)
  est1 <- nw_para[[1]]
  para <- nw_para[[2]]
  # C is the contact matrix, which is traced for 7 days. If i,j element is 7, it means the latest contact of i,j is today,
  # if it is 0, then there is no contact between i,j in the past 7 days,
  C <- matrix(0, para$pop_sz,para$pop_sz)
  
  # I is a vector indicating when an individual has been infected
  # NA means susceptible, 0 means infected date
  I <- matrix(NA, para$pop_sz,1)
  # Z is a vector indicating if the infected is detectable
  Z <- matrix(F, para$pop_sz,1)
  ############################################
  ############################################
  trace_inf_n <- matrix(0, para$pop_sz,1)
  
  # onset just save the incubation time 
  O <- matrix(NA, para$pop_sz,1)
  
  # Q is the vector indicating when an individual has been quarantine, NA means not quarantine, 1 means first day
  Q <- matrix(NA, para$pop_sz,1)
  RDT_n <- 0 # total number of RDT used
  PCR_n <- 0 # total number of PCR used
  rdt <- matrix(0, len_sim)
  pcr <- matrix(0, len_sim)
  
  C_lst <- list()
  O_lst <- list()
  I_lst <- list()
  Q_lst <- list()
  
  # initial case
  # para$E_0 <- floor(runif(1, min = 1, max = 3))
  init_idx <- sample(1:para$pop_sz, para$E_0)
  I[init_idx] <- 0
  O[init_idx] <- incub(para$E_0)
  Z[init_idx] <- F
  for(t in 1:len_sim){
    C_lst[[t]] <- C
    C <- C_update(C, Q, est1, para)
    lst <- I_O_update(I, Q, C, O, Z, trace_inf_n, para)
    I <- lst[[1]]
    O <- lst[[2]]
    Z <- lst[[3]]
    trace_inf_n <- lst[[4]]
    I_lst[[t]] <- I
    O_lst[[t]] <- O
    lst1 <- Q_update(I, Q, C, O, Z, RDT_n, PCR_n, para, flg_ab = FALSE, flg_multi_ab = FALSE, flg_PCR = TRUE, PCR_need_yesterday)
    Q <- lst1[[1]]
    RDT_n <-lst1[[2]]
    PCR_n <- lst1[[3]]
    PCR_need_yesterday <- lst1[[4]]
    Q_lst[[t]] <- Q
    rdt[t] <- RDT_n
    pcr[t] <- PCR_n
  }
  rt <- rep(NA, len_sim)
  for(t in 1:len_sim){
    i <- I_lst[[t]]
    idx <- which(i == 2) # from the second day they are infectious
    if(length(idx)) rt[t] <- mean(trace_inf_n[idx])
  }
  # plot the daily new case
  Q_new_daily <- rep(NA, len_sim)
  cap <- rep(NA, len_sim) # store the Q numbers
  Q_new_daily[1] <- para$E_0
  cap[1] <- sum(!is.na(Q_lst[[1]]) & Q_lst[[t]] < 14)
  for(t in 2:len_sim){
    Q_new_daily[t] <- sum(is.na(I_lst[[t-1]]))  - sum(is.na(I_lst[[t]]))
    cap[t] <- sum(!is.na(Q_lst[[t]]) & Q_lst[[t]] < 14)
  }
  all_Q_new_daily[rp,] <- Q_new_daily
  cap_Q[rp,] <- cap
  Rt_Q[rp,] <- rt
  PCR_Q[rp,] <- pcr
}

# Simulate Scenario 3: symptom & pcr & ab test
Rt_ab <- matrix(nrow = num_rep, ncol = len_sim)
all_ab_new_daily <- matrix(nrow = num_rep, ncol = len_sim)
cap_ab <- matrix(nrow = num_rep, ncol = len_sim) # store the number of max quarantine
# only save for this senario, where all three methods are combined
iso_num <- matrix(nrow = num_rep, ncol = len_sim) # store the number of isoloated
RDT_ab <-matrix(nrow = num_rep, ncol = len_sim)
PCR_ab <-matrix(nrow = num_rep, ncol = len_sim)
for(rp in 1:num_rep){
  print(paste("Running Setting=",setting,"V=",v_id,"S=",s_id,", Scenario 3, Rep: ",rp))
  
  nw_para <- load_sample_nw(setting,v_id,s_id, para)
  est1 <- nw_para[[1]]
  para <- nw_para[[2]]
  
  # C is the contact matrix, which is traced for 7 days. If i,j element is 7, it means the latest contact of i,j is today,
  # if it is 0, then there is no contact between i,j in the past 7 days,
  C <- matrix(0, para$pop_sz,para$pop_sz)
  
  # I is a vector indicating when an individual has been infected
  # NA means susceptible, 0 means infected date
  I <- matrix(NA, para$pop_sz,1)
  trace_inf_n <- matrix(0, para$pop_sz,1)
  # Z is a vector indicating if the infected is detectable
  Z <- matrix(F, para$pop_sz,1)
  
  # onset just save the incubation time 
  O <- matrix(NA, para$pop_sz,1)
  
  # Q is the vector indicating when an individual has been quarantine, NA means not quarantine, 1 means first day
  Q <- matrix(NA, para$pop_sz,1)
  RDT_n <- 0 # total number of RDT used
  PCR_n <- 0 # total number of PCR used
  rdt <- matrix(0, len_sim)
  pcr <- matrix(0, len_sim)
  
  C_lst <- list()
  O_lst <- list()
  I_lst <- list()
  Q_lst <- list()
  
  # initial case
  # para$E_0 <- floor(runif(1,min = 1, max = 3))
  init_idx <- sample(1:para$pop_sz, para$E_0)
  I[init_idx] <- 0
  O[init_idx] <- incub(para$E_0)
  Z[init_idx] <- F
  for(t in 1:len_sim){
    C_lst[[t]] <- C
    C <- C_update(C, Q, est1, para)
    lst <- I_O_update(I, Q, C, O, Z, trace_inf_n, para)
    I <- lst[[1]]
    O <- lst[[2]]
    Z <- lst[[3]]
    trace_inf_n <- lst[[4]]
    I_lst[[t]] <- I
    O_lst[[t]] <- O
    lst1 <- Q_update(I, Q, C, O, Z, RDT_n, PCR_n, para, flg_ab = TRUE, flg_multi_ab = FALSE, flg_PCR=TRUE, PCR_need_yesterday)
    Q <- lst1[[1]]
    RDT_n <-lst1[[2]]
    PCR_n <-lst1[[3]]
    PCR_need_yesterday <- lst1[[4]]
    Q_lst[[t]] <- Q
    rdt[t] <- RDT_n
    pcr[t] <- PCR_n
  }
  rt <- rep(NA, len_sim)
  for(t in 1:len_sim){
    i <- I_lst[[t]]
    idx <- which(i == 2) # from the second day they are infectious
    if(length(idx)) rt[t] <- mean(trace_inf_n[idx])
  }
  # plot the daily new case
  ab_new_daily <- rep(NA, len_sim)
  cap <- rep(NA, len_sim) # store the Q numbers
  iso <- rep(NA, len_sim) # store the isolated number 
  ab_new_daily[1] <- para$E_0
  cap[1] <- sum(!is.na(Q_lst[[1]]) & Q_lst[[t]] < 14)
  for(t in 2:len_sim){
    ab_new_daily[t] <- sum(is.na(I_lst[[t-1]]))  - sum(is.na(I_lst[[t]]))
    cap[t] <- sum(!is.na(Q_lst[[t]]) & Q_lst[[t]] < 14)
    iso[t] <- sum(!is.na(Q_lst[[t]]) & (!is.na(I_lst[[t]])) & (I_lst[[t]] <= 5))
  }
  all_ab_new_daily[rp,] <- ab_new_daily
  cap_ab[rp,] <- cap
  iso_num[rp,] <- iso
  Rt_ab[rp,] <- rt
  RDT_ab[rp,] <- rdt
  PCR_ab[rp,] <- pcr
}

# Simulate Scenario 4: symptom & pcr & antigen test
Rt_antig <- matrix(nrow = num_rep, ncol = len_sim)
all_antig_new_daily <- matrix(nrow = num_rep, ncol = len_sim)
cap_antig <- matrix(nrow = num_rep, ncol = len_sim) # store the number of max quarantine
RDT_antig <-matrix(nrow = num_rep, ncol = len_sim)
PCR_antig <-matrix(nrow = num_rep, ncol = len_sim)
###########################################
# reset parameter for the antigen case
para$ab_test_rate  <- para$pcr_test_rate # consent for antigen 
para$ab_rate <- function(x, t_onset) (1-para$samplefailure_pcr) # failure rate of sample for antigen 
# para$abiso # same as ab
para$sensitivity_ab <- para$sensitivity_antig # sensitivity is 0.8
para$tracing_cc_ab <- para$tracing_cc_onset
# para$cc_success_ab # using the same as ab
# para$qrate_ab # using the same as ab
para$delay_ab <- para$delay_symptom # we do the test as soon as the symptom is discovered
para$onsetiso <- 0 # in this case, we don't quarantine anyone based on symptom
###########################################

for(rp in 1:num_rep){
  print(paste("Running Setting=",setting,"V=",v_id,"S=",s_id,", Scenario 4, Rep: ",rp))
  
  nw_para <- load_sample_nw(setting,v_id,s_id, para)
  est1 <- nw_para[[1]]
  para <- nw_para[[2]]
  
  # C is the contact matrix, which is traced for 7 days. If i,j element is 7, it means the latest contact of i,j is today,
  # if it is 0, then there is no contact between i,j in the past 7 days,
  C <- matrix(0, para$pop_sz,para$pop_sz)
  
  # I is a vector indicating when an individual has been infected
  # NA means susceptible, 0 means infected date
  I <- matrix(NA, para$pop_sz,1)
  trace_inf_n <- matrix(0, para$pop_sz,1)
  # Z is a vector indicating if the infected is detectable
  Z <- matrix(F, para$pop_sz,1)
  
  # onset just save the incubation time 
  O <- matrix(NA, para$pop_sz,1)
  
  # Q is the vector indicating when an individual has been quarantine, NA means not quarantine, 1 means first day
  Q <- matrix(NA, para$pop_sz,1)
  RDT_n <- 0 # total number of RDT used
  PCR_n <- 0 # total number of PCR used
  rdt <- matrix(0, len_sim)
  pcr <- matrix(0, len_sim)
  
  C_lst <- list()
  O_lst <- list()
  I_lst <- list()
  Q_lst <- list()
  
  # initial case
  # para$E_0 <- floor(runif(1,min = 1, max = 3))
  init_idx <- sample(1:para$pop_sz, para$E_0)
  I[init_idx] <- 0
  O[init_idx] <- incub(para$E_0)
  Z[init_idx] <- F
  for(t in 1:len_sim){
    C_lst[[t]] <- C
    C <- C_update(C, Q, est1, para)
    lst <- I_O_update(I, Q, C, O, Z, trace_inf_n, para)
    I <- lst[[1]]
    O <- lst[[2]]
    Z <- lst[[3]]
    trace_inf_n <- lst[[4]]
    I_lst[[t]] <- I
    O_lst[[t]] <- O
    lst1 <- Q_update(I, Q, C, O, Z, RDT_n, PCR_n, para, flg_ab = TRUE, flg_multi_ab = FALSE, flg_PCR=TRUE, PCR_need_yesterday)
    Q <- lst1[[1]]
    RDT_n <-lst1[[2]]
    PCR_n <-lst1[[3]]
    PCR_need_yesterday <- lst1[[4]]
    Q_lst[[t]] <- Q
    rdt[t] <- RDT_n
    pcr[t] <- PCR_n
  }
  rt <- rep(NA, len_sim)
  for(t in 1:len_sim){
    i <- I_lst[[t]]
    idx <- which(i == 2) # from the second day they are infectious
    if(length(idx)) rt[t] <- mean(trace_inf_n[idx])
  }
  # plot the daily new case
  ab_new_daily <- rep(NA, len_sim)
  cap <- rep(NA, len_sim) # store the Q numbers
  ab_new_daily[1] <- para$E_0
  cap[1] <- sum(!is.na(Q_lst[[1]]) & Q_lst[[t]] < 14)
  for(t in 2:len_sim){
    ab_new_daily[t] <- sum(is.na(I_lst[[t-1]]))  - sum(is.na(I_lst[[t]]))
    cap[t] <- sum(!is.na(Q_lst[[t]]) & Q_lst[[t]] < 14)
  }
  all_antig_new_daily[rp,] <- ab_new_daily
  cap_antig[rp,] <- cap
  Rt_antig[rp,] <- rt
  RDT_antig[rp,] <- rdt
  PCR_antig[rp,] <- pcr
}


# # Simulate Scenario 4: symptom & pcr & multiple ab tests
# Rt_multiab <- matrix(nrow = num_rep, ncol = len_sim)
# all_multiab_new_daily <- matrix(nrow = num_rep, ncol = len_sim)
# cap_multiab <- matrix(nrow = num_rep, ncol = len_sim) # store the number of
# RDT_multiab <-matrix(nrow = num_rep, ncol = len_sim)
# PCR_multiab <-matrix(nrow = num_rep, ncol = len_sim)
# for(rp in 1:num_rep){
#     print(paste("Running Setting=",setting,"V=",v_id,"S=",s_id,", Scenario 4, Rep: ",rp))
#     # C is the contact matrix, which is traced for 7 days. If i,j element is 7, it means the latest contact of i,j is today,
#     # if it is 0, then there is no contact between i,j in the past 7 days,
#     C <- matrix(0, para$pop_sz,para$pop_sz)
#     
#     # I is a vector indicating when an individual has been infected
#     # NA means susceptible, 0 means infected date
#     I <- matrix(NA, para$pop_sz,1)
#     trace_inf_n <- matrix(0, para$pop_sz,1)
#     # Z is a vector indicating if the infected is detectable
#     Z <- matrix(F, para$pop_sz,1)
#     
#     # onset just save the incubation time
#     O <- matrix(NA, para$pop_sz,1)
#     
#     # Q is the vector indicating when an individual has been quarantine, NA means not quarantine, 1 means first day
#     Q <- matrix(NA, para$pop_sz,1)
#     RDT_n <- 0 # total number of RDT used
#     PCR_n <- 0 # total number of PCR used
#     rdt <- matrix(0, len_sim)
#     pcr <- matrix(0, len_sim)
#     
#     C_lst <- list()
#     O_lst <- list()
#     I_lst <- list()
#     Q_lst <- list()
#     multiab_new_daily <- rep(NA, len_sim)
#     cap <- rep(NA, len_sim) # store the Q numbers
#     rt <- rep(NA, len_sim)
#     
#  if (senario4_yn) {    
#     # initial case
#     # para$E_0 <- floor(runif(1,min = 1, max = 3))
#     init_idx <- sample(1:para$pop_sz, para$E_0)
#     I[init_idx] <- 0
#     O[init_idx] <- incub(para$E_0)
#     Z[init_idx] <- F
#     for(t in 1:len_sim){
#       C_lst[[t]] <- C
#       C <- C_update(C, Q, est1, para)
#       lst <- I_O_update(I, Q, C, O, Z, trace_inf_n, para)
#       I <- lst[[1]]
#       O <- lst[[2]]
#       Z <- lst[[3]]
#       trace_inf_n <- lst[[4]]
#       I_lst[[t]] <- I
#       O_lst[[t]] <- O
#       lst1 <- Q_update(I, Q, C, O, Z, RDT_n, PCR_n, para, flg_ab = FALSE, flg_multi_ab = TRUE, flg_PCR = TRUE, PCR_need_yesterday)
#       Q <- lst1[[1]]
#       RDT_n <-lst1[[2]]
#       PCR_n <-lst1[[3]]
#       PCR_need_yesterday <- lst1[[4]]
#       Q_lst[[t]] <- Q
#       rdt[t] <- RDT_n
#       pcr[t] <- PCR_n
#     }
#     for(t in 1:len_sim){
#       i <- I_lst[[t]]
#       idx <- which(i == 2) # from the second day they are infectious
#       if(length(idx)) rt[t] <- mean(trace_inf_n[idx])
#     }
#     # plot the daily new case
#     multiab_new_daily[1] <- para$E_0
#     cap[1] <- sum(!is.na(Q_lst[[1]]) & Q_lst[[t]] < 14)
#     for(t in 2:len_sim){
#       multiab_new_daily[t] <- sum(is.na(I_lst[[t-1]]))  - sum(is.na(I_lst[[t]]))
#       cap[t] <- sum(!is.na(Q_lst[[t]]) & Q_lst[[t]] < 14)
#     }
#   }
#   all_multiab_new_daily[rp,] <- multiab_new_daily
#   cap_multiab[rp,] <- cap
#   Rt_multiab[rp,] <- rt
#   RDT_multiab[rp,] <- rdt
#   PCR_multiab[rp,] <- pcr
# }

