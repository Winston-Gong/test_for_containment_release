t <- try(setwd("C:/Users/hgt/Documents/20200507"),silent=T)
if (!("try-error" %in% class(t))) {
  setwd("C:/Users/hgt/Documents/20200507")
  code_version<-"Multisetting_20200519"
  args <- commandArgs(trailingOnly = TRUE)
  setting <- as.numeric(args[1])
  v_id <- as.numeric(args[2])
  s_id <- as.numeric(args[3])
  ##########################
  # flag for social distancing
  social_distancing_flg <- as.numeric(args[4])
  ##########################
  num_rep <- 200
  senario4_yn <-F # whether no not run senario 4 (multi-ab)
}
t <- try(load(paste0(code_version,"_set",setting,"v",v_id,"s",s_id,".Rdata")),silent=T)
if (!("try-error" %in% class(t))) {
  print(paste("Already exist: Setting=",setting,"V=",v_id,"S=",s_id))
  quit(save = "no", status = 0, runLast = TRUE)
}
###############################################
# simulation for Rural village_bulk simulation
###############################################
library(reshape2)
library(ggplot2)
library(tidyverse)
suppressMessages(library(ergm))

###############################################
# Parameter setting
###############################################
  
len_sim <- 100 # length of simulation

para <- list()
para$setting <- setting # 1=rural 2=affluent 3=slum
# Environment setup
  para$pop_sz <- 1000 # 5000
  para$E_0 <- 2 # number of initial importation
  para$infect_sz <- (-1)*para$pop_sz/1000 # containment intervention only starts after X per 1000 individuals are already infected

R_UK <- 2.7
num_cc_UK <- 13
R0_baseline <- R_UK/num_cc_UK # the R0 is measured in affluent region 
# Transmission parameter & tool availability
if (para$setting==1) {
  para$num_cc <- 7 # set daily close contact to be 7 
  para$family_sz <- 5 # average household size 5
  para$symptom_report <- 0.7 # precentage infected report symptom
  para$R0 <- R0_baseline * para$num_cc # R naught
  para$theta <- 0.1 # quarantine effect: precentage of remaining transmission quarantined individual have
  para$pcr_available <- -1 # daily maximum PCR tests per 1000 population. 
  para$ppe_coef <- 1 # if people wear ppe, then their transmissibility will be different outside their family
  # set the percentage of HH_cc
  para$percent_HH_cc <- .5
}else if (para$setting==2) {
  para$num_cc <- 13 
  para$family_sz <- 3 
  para$symptom_report <- 0.7
  para$R0 <- R0_baseline * para$num_cc
  para$theta <- 0.1 
  para$pcr_available <-1000*para$pop_sz/1000
  para$ppe_coef <- 1 
  # set the percentage of HH_cc
  para$percent_HH_cc <- .3
}else if (para$setting==3) {
  para$num_cc <- 14 
  para$family_sz <- 4 
  para$symptom_report <- 0.6
  para$R0 <- R0_baseline * para$num_cc
  para$theta <- 0.2 
  para$pcr_available <-2*para$pop_sz/1000
  para$ppe_coef <- 1 
  # set the percentage of HH_cc
  para$percent_HH_cc <- .3
}else print ("Parameter setting error")
  
  # subclinical rate
  para$sub_clini_rate <- 0.3
  
  # below is a separate effect for social distancing, meaning when the social distancing is implemented Non_HH_CC_rate < 1, we reduce the number of contacts outside the family
  ################################
  # social distancing kept inposed
  # para$Non_HH_CC_rate <- 0.4 # this will impose a moderate social distancing by cutting 60% of non_HH 
  ################################
  para$Non_HH_CC_rate <- 1 # this term indicate an decreasing rate for contact outside house hold, due to the implementation of social distancing effect. 
  
# Parameters about Infected pts
  para$ab_test_rate <- 0.7 # % accepted ab testing among detected (symptomatic) patients
  para$pcr_test_rate <- 0.8 # % accepted pcr testing among detected (symptomatic) patients

  para$onsetiso <- 0.2 # Isolation compliance rate at onset based on symptom
  para$abiso <-0.9 # Isolation compliance rate among ab testing positive
  para$pcriso <-0.9 # Isolation compliance rate among pcr testing positive
  
  para$delay_symptom <- 1 # days of delay after onset to detect symptomatic patients
  para$delay_ab <- 8 # days of delay after onset to receive ab test and obtain result
  para$delay_pcr <- 5 # days of delay after onset to report pcr test result
  
# Parameters about tracing contects
  para$tracing_cc_onset <- 3 # set how many days we trace close contact back after symptom-based patient detection
  para$tracing_cc_ab <- para$delay_ab # set how many days we trace close contact back after a positive ab_test
  para$tracing_cc_pcr <- para$delay_pcr # set how many days we trace close contact back after a positive ab_test
  
  para$cc_success_symptom <- 0.85 # precentage close contact successfully traced after symptom-based patient detection
  para$cc_success_ab <- 0.75 # precentage close contact successfully traced after positive ab test
  para$cc_success_pcr <- 0.80 # precentage close contact successfully traced after positive pcr test
  
  para$qrate_symptom <- 0.5 # CC quarantine compliance rate based on symptom
  para$qrate_ab <- 0.7 # CC quarantine compliance rate based on positive ab test
  para$qrate_pcr <- 0.7 # CC quarantine compliance rate based on positive pcr test
  
# Parameters about testing tools 
  para$ab_rate <- function(x, t_onset) 1/(1 + exp(7+t_onset-x)) # seroconversion rate from infection day, based on the clinical paper from Yumei Wen
  para$sensitivity_ab <- 0.9 # ab test sensitivity
  para$sensitivity_pcr <- 0.999 # pcr test sensitivity
  para$samplefailure_pcr <- 0.3 # pcr sampling failure

###############################################
# Change parameters for sensitivity analysis
###############################################

if(v_id == 0){          # % reported symptom
    symptom_report <- c(0.3,0.5,0.7,0.9, 1.0)
    para$symptom_report  <- symptom_report[s_id]  
}else if(v_id == 1){    # which day to test antibody = how long we trace close contact for after ab_test
  t_day <- c(5,6,7,8,9,10,11)
  para$tracing_cc_ab <- t_day[s_id]
  para$delay_ab <- t_day[s_id]
}else if(v_id == 2){   # seroconversion rate
  ab_rate <- c( function(x, t_onset) 1/(1 + exp(7+t_onset-x)), # opimistic
                function(x, t_onset) 1/(1 + exp(12-x)),
                function(x, t_onset) 1/(1 + exp((9+t_onset-x)/2))) # conservative
  para$ab_rate <- ab_rate[[s_id]]
  # # plot(para$ab_rate(1:20,10)) # check antibody detection curve
}else if(v_id == 3){  # CC quarantine compliance rate 
  qrate <- c(0.1,0.3,0.5,0.7, 0.9)
  para$qrate_symptom <- qrate[s_id]  
  para$qrate_ab <- ifelse(qrate[s_id]==0.9, 1, qrate[s_id]+0.2)
  para$qrate_pcr <- para$qrate_ab  
}else if(v_id == 4){   # R naught
  r0 <- c(1.5,2.0,2.5,3.0,3.5)
  para$R0 <- r0[s_id]
}else if(v_id == 5){   # set daily close contact 
  num_cc <- c(2,4,7,10,14)
  para$num_cc <- num_cc[s_id]
}else if(v_id == 6){   # number of initial importation
  E_0 <- c(1,2,5,8,15)
  para$E_0 <- E_0[s_id]
}else if(v_id == 7){   # onset sympotom-based isolation sucessuful rate
  onsetiso <- c(0.1,0.2, 0.3,0.4,0.5,0.6, 0.7, 0.8, 0.9)
  para$onsetiso <- onsetiso[s_id]
}else if(v_id == 8){   # social distancing
  Non_HH_CC_rate<- c(1,.8,.6,.4,.2)
  para$Non_HH_CC_rate <- Non_HH_CC_rate[s_id] 
}else if(v_id == 9){  # precentage close contact successfully traced 
  cc_success <- c(0, 0.15,0.35,0.55,0.75,0.85,0.95)
  para$cc_success_symptom <- cc_success[s_id]
  para$cc_success_ab <- cc_success[s_id]-0.1
  para$cc_success_pcr <- cc_success[s_id] -0.05  
}else if(v_id == 10){  # days of delay after onset to report pcr test result
  delay_pcr <- c(0,1,3,5,7,10,14)
  para$delay_pcr <- delay_pcr[s_id]  
  para$tracing_cc_pcr <- delay_pcr[s_id]
}else if(v_id == 11){  # daily maximum PCR tests per 1000 population.
  pcr_available <- c(-1,1,2,10,20,50,1000)
  para$pcr_available <- pcr_available[s_id]  
}else if(v_id == 12){  # simulation size
  pop_sz <- c(100, 1000, 2000)
  para$pop_sz <- pop_sz[s_id]  
  num_rep <- 100
  senario4_yn <-F
}else if(v_id == 13){  # Delayed start of containment intervention when X number is infected
  infect_sz <- c(-1,5,20,50, 100, 250, 500)
  para$infect_sz <- infect_sz[s_id]*para$pop_sz/1000  
}else if(v_id == 14){  # percentage of household contact
  percent_HH_cc <- c(.2,.4,.6,.8)
  para$percent_HH_cc <-  percent_HH_cc[s_id]
}else if(v_id == 15){   # onset sympotom-based isolation sucessuful rate
  onsetiso <- c(0.1,0.2, 0.3,0.4,0.5,0.6, 0.7, 0.8, 0.9)
  para$onsetiso <- onsetiso[s_id]
  para$pcriso <- 0
  para$abiso <- 0
}else if(v_id == 16){   # onset pcr isolation sucessuful rate
  para$onsetiso <- 0
  para$pcr_available <- para$pop_sz
  pcriso <- c(0.1,0.2, 0.3,0.4,0.5,0.6, 0.7, 0.8, 0.9)
  para$pcriso <- pcriso[s_id]
  para$abiso <- 0
}else if(v_id == 17){   # onset rdt isolation sucessuful rate
  para$onsetiso <- 0
  para$pcriso <- 0
  abiso <- c(0.1,0.2, 0.3,0.4,0.5,0.6, 0.7, 0.8, 0.9)
  para$abiso <- abiso[s_id]
}else if(v_id == 18){  # test acceptance rate
  test_rate <- c(0.3,0.5,0.7,0.9)
  para$ab_test_rate <- test_rate[s_id] # % accepted ab testing among detected (symptomatic) patients
  para$pcr_test_rate <- test_rate[s_id] + 0.1 # % accepted pcr testing among detected (symptomatic) patients
}else if(v_id == 19){  # onset to isolation
  delay_pcr <- c(0,1,3,5)
  para$delay_symptom <- 1 + delay_pcr[s_id] # days of delay after onset to detect symptomatic patients
  para$delay_ab <- 8 + delay_pcr[s_id] # days of delay after onset to receive ab test and obtain result
  para$delay_pcr <- 5 + delay_pcr[s_id] # days of delay after onset to report pcr test result
}else if(v_id == 20){   # duration of cc tracing
  tracing_cc <- c(-2,0,2,4)
  para$tracing_cc_onset <- 3 + tracing_cc[s_id] # set how many days we trace close contact back after symptom-based patient detection
  para$tracing_cc_ab <- 8 + tracing_cc[s_id] # set how many days we trace close contact back after a positive ab_test
  para$tracing_cc_pcr <- 5 + tracing_cc[s_id] # set how many days we trace close contact back after a positive ab_test
}else if(v_id == 21){   # qurantine transmission rate
  theta <- c(0.05,0.1,0.15,0.2,0.25)
  para$theta <- theta[s_id]
}
##########################################
# parameter with age
##########################################
  # South africa
  para$age_dist <- c(0.292, 0.193, 0.515)
  # para$age_dist <- c(0.419, 0.195, 0.386) # kenya
  # para$age_dist <- c(0.440, 0.190, 0.370) # nigeria
  # para$age_dist <- c(0.481, 0.203, 0.316) # Uganda

  ##########################################
  para$AGE <- unlist(sapply(1:length(para$age_dist), function(x) rep(x,round(para$age_dist[x] * para$pop_sz))))
  stopifnot(length(para$AGE) == para$pop_sz)
  ############################################################
  # normalized the age matrix
  AGE_matrix <- read.csv("south_africa_age.csv", header = F)
  # AGE_matrix <- read.csv("kenya_age.csv", header = F) # kenya
  # AGE_matrix <- read.csv("Uganda_age.csv", header = F) # Uganda
  # AGE_matrix <- read.csv("Nigeria_age.csv", header = F) # nigeria
  age <- cbind(rowMeans(AGE_matrix[,1:3]),rowMeans(AGE_matrix[,4:5]),rowMeans(AGE_matrix[,6:16]))
  AGE_matrix <- rbind(colMeans(age[1:3,]),colMeans(age[4:5,]),colMeans(age[6:16,]))
  AGE_matrix <- (para$age_dist %*% t(as.matrix(para$age_dist))) * AGE_matrix
  # matrix in the end should be multipl
  AGE_matrix <- (AGE_matrix + t(AGE_matrix))/2
  para$age_mix <- as.matrix(AGE_matrix)[which(upper.tri(AGE_matrix,diag = T))]
  ############################################################
  ############################################################
  # normalized the house-hold age matrix
  HOME_age_matrix <- read.csv("home_south_africa_age.csv", header = F)
  # HOME_age_matrix <- read.csv("home_kenya_age.csv", header = F) # kenya
  # HOME_age_matrix <- read.csv("home_Uganda_age.csv", header = F) # Uganda
  # HOME_age_matrix <- read.csv("home_Nigeria_age.csv", header = F) # Nigeria
  age <- cbind(rowMeans(HOME_age_matrix[,1:3]),rowMeans(HOME_age_matrix[,4:5]),rowMeans(HOME_age_matrix[,6:16]))
  HOME_age_matrix <- rbind(colMeans(age[1:3,]),colMeans(age[4:5,]),colMeans(age[6:16,]))
  HOME_age_matrix <- (para$age_dist %*% t(as.matrix(para$age_dist))) * HOME_age_matrix
  # matrix in the end should be multipl
  HOME_age_matrix <- (HOME_age_matrix + t(HOME_age_matrix))/2
  para$Home_age_mix <- as.matrix(HOME_age_matrix)[which(upper.tri(HOME_age_matrix,diag = T))]
  ############################################################
  # adjust the age matrix to represent the specified household contact rate
  para$age_mix <- para$Home_age_mix + (para$age_mix - para$Home_age_mix) * 
    sum(para$Home_age_mix)/sum(para$age_mix - para$Home_age_mix) *  (1-para$percent_HH_cc)/para$percent_HH_cc
  
#########################################################################
# adjust R0 for 1) young people susceptibility 2) subclinical cases
#########################################################################
  # ajust for the social distancing
  para$num_cc_scdst <- para$num_cc * ((1-para$percent_HH_cc)*para$Non_HH_CC_rate + para$percent_HH_cc) # reduce the number of cc
  para$age_mix_scdst <- para$Home_age_mix + (para$age_mix - para$Home_age_mix) * para$Non_HH_CC_rate
  para$percent_HH_cc_scdst <- para$percent_HH_cc/((1-para$percent_HH_cc)*para$Non_HH_CC_rate + para$percent_HH_cc)
  
  # ajust for R0
  norm1 <- (sum(para$age_mix) - para$age_mix[c(1)]/2- sum(para$age_mix[c(2,4)])/4)/sum(para$age_mix)
  norm2 <- 1 - para$sub_clini_rate * (1-0.2)
  para$R0_adj <- para$R0/(norm1 * norm2)
  ##################################
  # compute Re
  norm_age_mix_scdst <- para$age_mix_scdst/sum(para$age_mix_scdst)
  para$Re <- para$R0_adj/para$num_cc * 
    ((1- norm_age_mix_scdst[c(1)]/2- sum(norm_age_mix_scdst[c(2,4)])/4) * para$num_cc_scdst) * 
    norm2 # approximate Re (haven't taken age-susceptibility into account) 
  ##################################
  
###############################################
# Call master code
###############################################

print(paste("Runing: ",code_version,"_Setting=",setting,"V=",v_id,"S=",s_id))
source(paste0(code_version,"_functions.R"),print.eval =F)
source(paste0(code_version,"_base.R"),print.eval =T)

SIM <- list(all_new_daily, all_Q_new_daily, all_ab_new_daily, all_antig_new_daily, 
            cap_Q, cap_ab, cap_antig, Rt, Rt_Q, Rt_ab, Rt_antig, RDT_ab, RDT_antig,
            PCR_Q, PCR_ab, PCR_antig,para, est1, iso_num)
save(SIM, file = paste0(code_version,"_set",setting,"v",v_id,"s",s_id,".Rdata"))

print(paste("Completed: Setting=",setting,"V=",v_id,"S=",s_id))

