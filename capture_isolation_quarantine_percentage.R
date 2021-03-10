# load a specific condition
library(statnet)
library(dplyr)
library(ergm)
library(reshape2)
library(ggplot2)
library(tidyr)
library(network)
DATE_version <- "20200803"
code_version<-"Multisetting_20200519"

# use color blind friendly scheme
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# set the commonly used color
blue <- cbPalette[6]
red <- cbPalette[7]
green <- cbPalette[4]
orange <- cbPalette[2]
grey <- cbPalette[1]
yellow <- cbPalette[5]
purple <- cbPalette[8]
skyblue <- cbPalette[3]

source(paste0(code_version,"_functions.R"),print.eval =F)

# use Q = 0.5 to indicate close contact quarantined; isolated Q = 1 
Q_update <- function(I, Q, C, O, Z, RDT_n, PCR_n, para, flg_ab, flg_multi_ab, flg_PCR, PCR_need_yesterday){ 
  Q <- Q + 1
  PCR_need_today <-0
  PCR_n_today <-0
  if (sum(!is.na(I))>para$infect_sz) { # take containment intervention only when the infected number passes the threshold 
    # quarantine close contact
    for(i in 1:para$pop_sz){
      if(!is.na(I[i]) & (Z[i] == 2)){ # only when Z == 2 the case is reported at all
        u <- runif(4)
        
        # case detection strategy
        ok_onset_iso <- (I[i]-O[i]) == para$delay_symptom & 
          u[1] < para$onsetiso
        ok_pcr_test <- flg_PCR & 
          (PCR_n_today <= para$pcr_available) & 
          ((I[i]-O[i]) == para$delay_pcr) & 
          u[2] < para$pcr_test_rate 
        ok_ab_test <-  flg_ab & 
          ((I[i]-O[i]) == para$delay_ab) & 
          u[3] < para$ab_test_rate                
        ok_multiab_test <- flg_multi_ab & 
          ((I[i]-O[i]) == 6 | (I[i]-O[i]) == 8 | (I[i]-O[i]) == 10) & 
          u[4] < para$ab_test_rate
        
        # quarantine strategy
        if(ok_onset_iso & is.na(Q[i])){ 
          Q[i] <- 1
          cc <- union(which(C[-i,i] >= 12 - para$tracing_cc_onset), which(C[i,-i] >= 12- para$tracing_cc_onset))
          u <- runif(length(cc))
          # make sure the cc_success_rate for family is always 1
          stopifnot(length(para$family_lable[cc]) == length(cc))
          tracing_rate <- ifelse(para$family_lable[i] == para$family_lable[cc], 1, para$cc_success_symptom)
          #####################################
          # for debug
          # print(length(cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_symptom]))
          #####################################
          Q[cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_symptom]] <- .5
        }
        
        if(ok_pcr_test & is.na(Q[i])){ 
          PCR_need_today <- PCR_need_today + 1
          pcravailable <- ifelse(para$setting==3,
                                 (runif(1)<(para$pcr_available/(PCR_need_yesterday+1)*2)),
                                 T) # in a slum, when PCR is not sufficient, smaller id people have a higher chance of geting the PCR test
          if (pcravailable) {
            PCR_n <- PCR_n + 1
            PCR_n_today <- PCR_n_today +1
            u <- runif(2)
            if(u[1] < para$pcriso & u[2] < para$sensitivity_pcr*(1-para$samplefailure_pcr)){ # if the infected is isolated, quarantine their close contact
              Q[i] <- 1
              cc <- union(which(C[-i,i] >= 12 - para$tracing_cc_pcr), which(C[i,-i] >= 12 - para$tracing_cc_pcr))
              u <- runif(length(cc))
              # make sure the cc_success_rate for family is always 1
              tracing_rate <- ifelse(para$family_lable[i] == para$family_lable[cc], 1, para$cc_success_pcr)
              Q[cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_pcr]] <- .5
              
            }
          }
        }
        
        if(ok_ab_test & is.na(Q[i])){ 
          RDT_n <- RDT_n + 1  
          u <- runif(2)
          seroconvert <- para$ab_rate(I[i],O[i]) # if the test shows positive
          if(u[1] < para$abiso & u[2] < para$sensitivity_ab*seroconvert){ # if the infected is isolated, quarantine their close contact
            Q[i] <- 1
            cc <- union(which(C[-i,i] >= 12 - para$tracing_cc_ab), which(C[i,-i] >= 12 - para$tracing_cc_ab))
            u <- runif(length(cc))
            # make sure the cc_success_rate for family is always 1
            tracing_rate <- ifelse(para$family_lable[i] == para$family_lable[cc], 1, para$cc_success_ab)
            Q[cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_ab]] <- .5
          }
        } else if( ok_multiab_test & is.na(Q[i])){ 
          RDT_n <- RDT_n + 1  
          u <- runif(2)
          seroconvert <- para$ab_rate(I[i],O[i]) # if the test shows positive
          if(u[1] < para$abiso & u[2] < para$sensitivity_ab*seroconvert){ # if the infected is isolated, quarantine their close contact
            Q[i] <- 1
            cc <- union(which(C[-i,i] >= 12 - (I[i]-O[i])), which(C[i,-i] >= 12 - (I[i]-O[i])))
            u <- runif(length(cc))
            tracing_rate <- ifelse(para$family_lable[i] == para$family_lable[cc], 1, para$cc_success_ab)
            Q[cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_ab]] <- 0.5
          }
        }
      }
    }
    # release those after 14 days of quarantine
    for(i in 1:para$pop_sz){
      # we only release those who are do not have onset: if they do have onset,
      # they will be keep quarantined (this is our strategy to calculate medical needs)
      if(!is.na(Q[i]) & Q[i] >= 14){
        if(is.na(I[i])){
          Q[i] <- NA
        }else if(I[i] < O[i]){ # so this guys haven't develop any symptom, and will be free
          Q[i] <- NA
        }
      }
    }
  }
  return(list(Q,RDT_n,PCR_n,PCR_need_today))
}
#######################################
# set working directory
#######################################
setwd(paste0("~/Desktop/Coronavirus/COVID_test_containment/",DATE_version))

#######################################
# strategy 3: Antigen + PCR
#######################################
setting <- 1
v_id <- 3
num_rep <- 10
len_sim <- 100
for(s_id in 1:5){
  load(paste0("set2v", v_id,"s", s_id,".Rdata"))
  para<-SIM[[17]]
  est1 <- SIM[[18]]
  
  # Simulate Scenario 3: symptom & pcr & antigen test
  Rt_ab <- matrix(nrow = num_rep, ncol = len_sim)
  all_ab_new_daily <- matrix(nrow = num_rep, ncol = len_sim)
  cap_ab <- matrix(nrow = num_rep, ncol = len_sim) # store the number of max quarantine
  # only save for this senario, where all three methods are combined
  iso_num <- matrix(nrow = num_rep, ncol = len_sim) # store the number of isoloated via onset
  q_cc_num <- matrix(nrow = num_rep, ncol = len_sim) # store the number of isoloated via close contact tracing
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
    qcc <- rep(NA, len_sim) # store the isolated number 
    ab_new_daily[1] <- para$E_0
    cap[1] <- sum(!is.na(Q_lst[[1]]) & Q_lst[[t]] < 14)
    for(t in 2:len_sim){
      ab_new_daily[t] <- sum(is.na(I_lst[[t-1]]))  - sum(is.na(I_lst[[t]]))
      cap[t] <- sum(!is.na(Q_lst[[t]]) & Q_lst[[t]] < 14)
      iso[t] <- sum(!is.na(Q_lst[[t]]) & (!is.na(I_lst[[t]])) & (Q_lst[[t]] == 1))
      qcc[t] <- sum(!is.na(Q_lst[[t]]) & (!is.na(I_lst[[t]])) & (Q_lst[[t]] == 0.5))
    }
    all_ab_new_daily[rp,] <- ab_new_daily
    cap_ab[rp,] <- cap
    iso_num[rp,] <- iso
    q_cc_num[rp,] <- qcc
    Rt_ab[rp,] <- rt
    RDT_ab[rp,] <- rdt
    PCR_ab[rp,] <- pcr
  }
  CAP_QCC <- list( all_ab_new_daily,cap_ab, Rt_ab, RDT_ab, PCR_ab, para, iso_num, q_cc_num)
  save(CAP_QCC, file = paste0("../Effect_Q_CC_set",setting,"v",v_id,"s",s_id,".Rdata"))
}

# load all files 
num_isolated <- rep(NA, 5)
num_quarantine <- rep(NA, 5)
num_not_q <- rep(NA, 5)
for(s_id in 1:5){
  load(file = paste0("Effect_Q_CC_set",setting,"v",v_id,"s",s_id,".Rdata"))
  all_ab_new_daily <- CAP_QCC[[1]]
  iso_num <- CAP_QCC[[7]]
  q_cc_num <- CAP_QCC[[8]]
  
  num_isolated[s_id] <- mean(rowSums(iso_num, na.rm = T))
  num_quarantine[s_id] <- mean(rowSums(q_cc_num, na.rm = T))
  num_not_q[s_id] <- mean(rowSums(all_ab_new_daily, na.rm = T)) - mean(rowSums(iso_num, na.rm = T)) - mean(rowSums(q_cc_num, na.rm = T))
}
qrate <- c(0.3,0.5,0.7, 0.9, 1)
data_q_percent <- data.frame(qrate, num_isolated, num_quarantine, num_not_q) %>%
  mutate(q_cc_rate = num_quarantine/(num_isolated + num_quarantine + num_not_q)) %>%
  mutate(q_percent = num_isolated/(num_isolated + num_quarantine + num_not_q))
  
plt <- ggplot(data_q_percent,aes(x = qrate)) +
  geom_line(aes(y = q_percent), color = red, linetype = "dashed") + 
  geom_point(aes(y = q_percent, color = "Cases indentified and isolated by testing"), size = 2) + 
  geom_line(aes(y =  q_cc_rate, color = "Cases quarantined by close contact tracing"), linetype = "dashed") +
  geom_point(aes(y = q_cc_rate),color = green, size = 2) + 
  labs(x ="Close contact quarantine consent rate", y = "Percentage of total infected") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,1)) + 
  scale_y_continuous(labels = scales::percent, limits = c(0,0.6)) + 
  scale_color_manual(values = c("Cases indentified and isolated by testing" = red,"Cases quarantined by close contact tracing" = green)) +
  theme(legend.position = "right",panel.background=element_blank())
ggsave(paste0("~/Desktop/Coronavirus/paper_figures/isolation_percentaage_setting",setting, "_antigen_isolation_method.png"),plt, width = 8, height = 4)

#####################
# strategy1: symptom
#####################
setting <- 3
v_id <- 3
setwd(paste0("~/Desktop/Coronavirus/COVID_test_containment/",DATE_version))
num_rep <- 10
len_sim <- 100
for(s_id in 1:5){
  load(paste0("set2v", v_id,"s", s_id,".Rdata"))
  para<-SIM[[17]]
  est1 <- SIM[[18]]
  
  para$qrate_pcr <- para$qrate_symptom
  
  # Simulate Scenario 1: symptom & PCR
  all_Q_new_daily <- matrix(nrow = num_rep, ncol = len_sim)
  cap_Q <- matrix(nrow = num_rep, ncol = len_sim) # store the number of max quarantine
  Rt_Q <- matrix(nrow = num_rep, ncol = len_sim)
  PCR_Q <-matrix(nrow = num_rep, ncol = len_sim)
  
  # only save for this senario, where all three methods are combined
  iso_num <- matrix(nrow = num_rep, ncol = len_sim) # store the number of isoloated via onset
  q_cc_num <- matrix(nrow = num_rep, ncol = len_sim) # store the number of isoloated via close contact tracing

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
    ab_new_daily <- rep(NA, len_sim)
    cap <- rep(NA, len_sim) # store the Q numbers
    iso <- rep(NA, len_sim) # store the isolated number 
    qcc <- rep(NA, len_sim) # store the isolated number 
    ab_new_daily[1] <- para$E_0
    cap[1] <- sum(!is.na(Q_lst[[1]]) & Q_lst[[t]] < 14)
    for(t in 2:len_sim){
      ab_new_daily[t] <- sum(is.na(I_lst[[t-1]]))  - sum(is.na(I_lst[[t]]))
      cap[t] <- sum(!is.na(Q_lst[[t]]) & Q_lst[[t]] < 14)
      iso[t] <- sum(!is.na(Q_lst[[t]]) & (!is.na(I_lst[[t]])) & (Q_lst[[t]] == 1))
      qcc[t] <- sum(!is.na(Q_lst[[t]]) & (!is.na(I_lst[[t]])) & (Q_lst[[t]] == 0.5))
    }
    all_Q_new_daily[rp,] <- ab_new_daily
    cap_Q[rp,] <- cap
    iso_num[rp,] <- iso
    q_cc_num[rp,] <- qcc
    Rt_Q[rp,] <- rt
    PCR_Q[rp,] <- pcr
  }
  CAP_QCC <- list( all_Q_new_daily,cap_Q, Rt_Q,  PCR_Q, para, iso_num, q_cc_num)
  save(CAP_QCC, file = paste0("../Strategy1_symptom_Effect_Q_CC_set",setting,"v",v_id,"s",s_id,".Rdata"))
}

# load all files 
num_isolated <- rep(NA, 5)
num_quarantine <- rep(NA, 5)
num_not_q <- rep(NA, 5)
for(s_id in 1:5){
  load(file = paste0("Strategy1_symptom_Effect_Q_CC_set",setting,"v",v_id,"s",s_id,".Rdata"))
  all_Q_new_daily <- CAP_QCC[[1]]
  iso_num <- CAP_QCC[[6]]
  q_cc_num <- CAP_QCC[[7]]
  
  num_isolated[s_id] <- mean(rowSums(iso_num, na.rm = T))
  num_quarantine[s_id] <- mean(rowSums(q_cc_num, na.rm = T))
  num_not_q[s_id] <- mean(rowSums(all_Q_new_daily, na.rm = T)) - mean(rowSums(iso_num, na.rm = T)) - mean(rowSums(q_cc_num, na.rm = T))
}
qrate <- c(0.1,0.3,0.5, 0.7, 0.9)
data_q_percent <- data.frame(qrate, num_isolated, num_quarantine, num_not_q) %>%
  mutate(q_cc_rate = num_quarantine/(num_isolated + num_quarantine + num_not_q)) %>%
  mutate(q_percent = num_isolated/(num_isolated + num_quarantine + num_not_q))

plt <- ggplot(data_q_percent,aes(x = qrate)) +
  geom_line(aes(y = q_percent), color = red, linetype = "dashed") + 
  geom_point(aes(y = q_percent, color = "Cases indentified and isolated by testing"), size = 2) + 
  geom_line(aes(y =  q_cc_rate, color = "Cases quarantined by close contact tracing"), linetype = "dashed") +
  geom_point(aes(y = q_cc_rate),color = green, size = 2) + 
  labs(x ="Close contact quarantine consent rate", y = "Percentage of total infected") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,1)) + 
  scale_y_continuous(labels = scales::percent, limits = c(0,0.6)) + 
  scale_color_manual(values = c("Cases indentified and isolated by testing" = red,"Cases quarantined by close contact tracing" = green)) +
  theme(legend.position = "right",panel.background=element_blank())
ggsave(paste0("~/Desktop/Coronavirus/paper_figures/isolation_percentaage_setting",setting, "_onset_isolation_method.png"),plt, width = 8, height = 4)

#######################################
# Social distancing + strategy 3: Antigen + PCR
#######################################
setting <- 2
v_id <- 3
setwd(paste0("~/Desktop/Coronavirus/COVID_test_containment/",DATE_version))
num_rep <- 10
len_sim <- 100
for(s_id in 1:5){
  load(paste0("set2v", v_id,"s", s_id,".Rdata"))
  para<-SIM[[17]]
  est1 <- SIM[[18]]
  
  # Simulate Scenario 3: symptom & pcr & antigen test
  Rt_ab <- matrix(nrow = num_rep, ncol = len_sim)
  all_ab_new_daily <- matrix(nrow = num_rep, ncol = len_sim)
  cap_ab <- matrix(nrow = num_rep, ncol = len_sim) # store the number of max quarantine
  # only save for this senario, where all three methods are combined
  iso_num <- matrix(nrow = num_rep, ncol = len_sim) # store the number of isoloated via onset
  q_cc_num <- matrix(nrow = num_rep, ncol = len_sim) # store the number of isoloated via close contact tracing
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
    qcc <- rep(NA, len_sim) # store the isolated number 
    ab_new_daily[1] <- para$E_0
    cap[1] <- sum(!is.na(Q_lst[[1]]) & Q_lst[[t]] < 14)
    for(t in 2:len_sim){
      ab_new_daily[t] <- sum(is.na(I_lst[[t-1]]))  - sum(is.na(I_lst[[t]]))
      cap[t] <- sum(!is.na(Q_lst[[t]]) & Q_lst[[t]] < 14)
      iso[t] <- sum(!is.na(Q_lst[[t]]) & (!is.na(I_lst[[t]])) & (Q_lst[[t]] == 1))
      qcc[t] <- sum(!is.na(Q_lst[[t]]) & (!is.na(I_lst[[t]])) & (Q_lst[[t]] == 0.5))
    }
    all_ab_new_daily[rp,] <- ab_new_daily
    cap_ab[rp,] <- cap
    iso_num[rp,] <- iso
    q_cc_num[rp,] <- qcc
    Rt_ab[rp,] <- rt
    RDT_ab[rp,] <- rdt
    PCR_ab[rp,] <- pcr
  }
  CAP_QCC <- list( all_ab_new_daily,cap_ab, Rt_ab, RDT_ab, PCR_ab,para, iso_num, q_cc_num)
  save(CAP_QCC, file = paste0("../social_dsit_strategy_3_Effect_Q_CC_set",setting,"v",v_id,"s",s_id,".Rdata"))
}

# load all files 
num_isolated <- rep(NA, 5)
num_quarantine <- rep(NA, 5)
num_not_q <- rep(NA, 5)
for(s_id in 1:5){
  load(file = paste0("../social_dsit_strategy_3_Effect_Q_CC_set",setting,"v",v_id,"s",s_id,".Rdata"))
  all_ab_new_daily <- CAP_QCC[[1]]
  iso_num <- CAP_QCC[[7]]
  q_cc_num <- CAP_QCC[[8]]
  
  num_isolated[s_id] <- mean(rowSums(iso_num, na.rm = T))
  num_quarantine[s_id] <- mean(rowSums(q_cc_num, na.rm = T))
  num_not_q[s_id] <- mean(rowSums(all_ab_new_daily, na.rm = T)) - mean(rowSums(iso_num, na.rm = T)) - mean(rowSums(q_cc_num, na.rm = T))
}
qrate <- c(0.3,0.5,0.7, 0.9, 1)
data_q_percent <- data.frame(qrate, num_isolated, num_quarantine, num_not_q) %>%
  mutate(q_cc_rate = num_quarantine/(num_isolated + num_quarantine + num_not_q)) %>%
  mutate(q_percent = num_isolated/(num_isolated + num_quarantine + num_not_q))

plt <- ggplot(data_q_percent,aes(x = qrate)) +
  geom_line(aes(y = q_percent), color = red, linetype = "dashed") + 
  geom_point(aes(y = q_percent, color = "Cases indentified and isolated by testing"), size = 2) + 
  geom_line(aes(y =  q_cc_rate, color = "Cases quarantined by close contact tracing"), linetype = "dashed") +
  geom_point(aes(y = q_cc_rate),color = green, size = 2) + 
  labs(x ="Close contact quarantine consent rate", y = "Percentage of total infected") + 
  scale_x_continuous(labels = scales::percent, limits = c(0,1)) + 
  scale_y_continuous(labels = scales::percent, limits = c(0,0.6)) + 
  scale_color_manual(values = c("Cases indentified and isolated by testing" = red,"Cases quarantined by close contact tracing" = green)) +
  theme(legend.position = "right",panel.background=element_blank())

ggsave(paste0("~/Desktop/Coronavirus/paper_figures/isolation_percentaage_social_distancing_setting",setting, "_antigen_isolation_method.png"),plt, width = 8, height = 4)

####################################################
# perform simulation without close contact tracing
####################################################
source(paste0(code_version,"_functions.R"),print.eval =F)
DATE_version <- "Multisetting_20200519_"
v_id <- 8
s_id <- 1
num_rep <- 200
for(setting in 1:3){
  load(paste0("Dynamics/", DATE_version,"set",setting,"v",v_id,"s",s_id,".Rdata"))
  para<-SIM[[17]]
  est1 <- SIM[[18]]
  # no close contact tracing
  para$qrate_symptom <- 0 
  para$qrate_ab <- 0
  para$qrate_pcr <- 0 
  print(paste("Runing: No close contact tracing","_Setting=",setting,"V=",v_id,"S=",s_id))
  
  source(paste0(code_version,"_base.R"),print.eval =T)
  
  SIM <- list(all_new_daily, all_Q_new_daily, all_ab_new_daily, all_antig_new_daily, 
              cap_Q, cap_ab, cap_antig, Rt, Rt_Q, Rt_ab, Rt_antig, RDT_ab, RDT_antig,
              PCR_Q, PCR_ab, PCR_antig,para, est1, iso_num)
  
  save(SIM, file = paste0("Dynamics/No_contact_tracing_",code_version,"_set",setting,"v",v_id,"s",s_id,".Rdata"))
  
  print(paste("Completed: Setting=",setting,"V=",v_id,"S=",s_id))
}

# compare the results with no close contact tracing
setting <- 3
v_id <- 8
s_id <- 1
len_sim<-100

load(paste0("Dynamics/","No_contact_tracing_", DATE_version, "set",setting,"v",v_id,"s",s_id,".Rdata"))
all_new_daily <- SIM[[1]] # do nothing
all_Q_new_daily <- SIM[[2]] # quarantine on site
# no PCR

all_ab_new_daily <- SIM[[3]] # quarantine + ab testing
all_antig_new_daily <- SIM[[4]] # quarantine + antigen testing
cap_Q <- SIM[[5]]
cap_ab <- SIM[[6]]
cap_antig<-SIM[[7]]
Rt <-SIM[[8]]
Rt_Q <-SIM[[9]]
Rt_ab <-SIM[[10]]
Rt_antig <-SIM[[11]]
RDT_ab <-SIM[[12]]
RDT_antig <-SIM[[13]]
PCR_Q <- SIM[[14]]
PCR_ab <- SIM[[15]]
PCR_antig <-SIM[[16]]
para<-SIM[[17]]
est1 <- SIM[[18]]
iso_num <- SIM[[19]]

err_daily <- apply(all_new_daily, 2, function(x) sd(x)/sqrt(length(x)))
err_Q <- apply(all_Q_new_daily, 2, function(x) sd(x)/sqrt(length(x)))
err_ab <- apply(all_ab_new_daily, 2, function(x) sd(x)/sqrt(length(x)))
err_antig <- apply(all_antig_new_daily, 2, function(x) sd(x)/sqrt(length(x)))

df_plt <- data.frame(day = 1:len_sim, Q_new_daily = colMeans(all_Q_new_daily), new_daily = colMeans(all_new_daily),
                     ab_new_daily = colMeans(all_ab_new_daily),antig_new_daily = colMeans(all_antig_new_daily), 
                     err_daily, err_Q, err_ab,err_antig)
plt <- ggplot(df_plt,aes(x=day)) +
  geom_line(aes( y = Q_new_daily, color = "Strategy 1: symptom + PCR")) +
  geom_line(aes(y = new_daily, color = "Baseline")) +
  geom_line(aes(y = ab_new_daily, color = "Strategy 2: symptom + PCR + antibody RDT")) +
  geom_line(aes(y = antig_new_daily, color = "Strategy 3: antigen RDT + PCR")) +
  geom_ribbon(aes( ymax = Q_new_daily + 2*err_Q, ymin = Q_new_daily - 2*err_Q), fill = blue,alpha = 0.3) +
  geom_ribbon(aes( ymax = new_daily + 2*err_daily, ymin = new_daily - 2*err_daily), fill = red,alpha = 0.3) +
  geom_ribbon(aes( ymax = antig_new_daily + 2*err_antig, ymin = antig_new_daily - 2*err_antig), fill = purple,alpha = 0.3) +
  geom_ribbon(aes( ymax = ab_new_daily + 2*err_ab, ymin = ab_new_daily- 2*err_ab), fill = green,alpha = 0.3) + 
  labs(x= "Days from first importing case", y = "New case per day")+
  scale_color_manual(values = c("Baseline" = red, "Strategy 1: symptom + PCR" = blue, 
                                "Strategy 2: symptom + PCR + antibody RDT" = green, "Strategy 3: antigen RDT + PCR" = purple)) + 
  theme(legend.position = c(0.8,0.8), panel.background=element_blank()) + 
  ylim(c(0,45))

summarise(df_plt, sum(new_daily), sum(Q_new_daily), sum( ab_new_daily), sum(antig_new_daily))/1000

ggsave(paste0("~/Desktop/Coronavirus/paper_figures/No_contact_tracing",setting, "baseline.png"),plt, width = 5, height = 4)

###############################################################
# Get the infection rate for each quarantine consent rate
###############################################################
setting <- 3 # 2
setwd("~/Desktop/Coronavirus/COVID_test_containment/")

v_id <- 3
# first low consent rate 0.1
s_id <- 1
load(paste0("Dynamics/",DATE_version,"set",setting,"v",v_id,"s",s_id,".Rdata"))
all_new_daily <- rowSums(SIM[[1]]) # do nothing
all_Q_new_daily <- rowSums(SIM[[2]])  # quarantine on site
all_ab_new_daily <- rowSums(SIM[[3]])  # quarantine + ab testing
all_antig_new_daily <- rowSums(SIM[[4]])  # quarantine + antigen testing
para<-SIM[[17]]
print(mean(all_Q_new_daily))


# high consent rate 0.9
s_id <- 5
load(paste0("Dynamics/",DATE_version,"set",setting,"v",v_id,"s",s_id,".Rdata"))
social_dist_all_new_daily <- rowSums(SIM[[1]]) # do nothing
social_dist_all_Q_new_daily <- rowSums(SIM[[2]]) # quarantine on site
social_dist_all_ab_new_daily <- rowSums(SIM[[3]])# quarantine + ab testing
social_dist_all_antig_new_daily <- rowSums(SIM[[4]]) # quarantine + antigen testing
para<-SIM[[17]]
print(mean(social_dist_all_Q_new_daily))




