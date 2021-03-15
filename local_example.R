code_version<-"Multisetting_20200519"

############################################
# change following parameters for different setting and social distancing
setting <- 2 # rural simulation; 2 for urban, 3 for slum
country_id <- 1 # we use Uganda demographic for demonstration
social_distancing_flg <- 1 # no social distancing

library(reshape2)
library(ggplot2)
# library(tidyverse)
suppressMessages(library(ergm))

para <- list()
para$setting <- setting # 1=rural 2=affluent 3=slum
# Environment setup
para$pop_sz <- 200 # simulation a smaller community
para$Non_HH_CC_rate <- c(1,.8,.6,.4,.2)[social_distancing_flg] 
community_setting <- c("Rural", "Non-slum urban", "Slum")[setting]
print(paste0("Simulate ",para$pop_sz," individuals in ",community_setting," setting"))
##########################################
# loading country specific variables: 
# age distribution & household
##########################################
if(country_id == 1){
  para$age_dist <- c(0.481, 0.203, 0.316) # Uganda
  para$HH_dist <- c(.11, .22, .27, .4) # Uganda
}else if(country_id == 2){
  para$age_dist <- c(0.292, 0.193, 0.515) # South africa
  para$HH_dist <- c(.27, .35, .23, .15) # South Africa
}else if(country_id == 3){
  para$age_dist <- c(0.419, 0.195, 0.386) # kenya
  para$HH_dist <- c(.19, .28, .3, .23) # kenya
}else if(country_id == 4){
  para$age_dist <- c(0.440, 0.190, 0.370) # nigeria
  para$HH_dist <- c(.16, .26, .26, .32) # nigeria
}

##########################################
# processing demographic information;
# for details: refer to the supplementary methods
##########################################
# para$HH_affluent_dist <- c(.31, .5, .18, .02) # UK

if (para$setting==1) {
  para$num_cc <- 7 # set daily close contact to be 7 
  para$family_sz <- 5 # average household size 5
  # set the percentage of HH_cc
  para$percent_HH_cc <- .5
}else if (para$setting==2) {
  para$num_cc <- 13 
  para$family_sz <- 5
  para$percent_HH_cc <- .23
}else if (para$setting==3) {
  para$num_cc <- 14 
  para$family_sz <- 15 
  para$percent_HH_cc <- .5
  para$HH_dist <- c(.00, .06, .17, .77) # afganistan
}else print ("Parameter setting error")


para$AGE <- unlist(sapply(1:length(para$age_dist), function(x) rep(x,round(para$age_dist[x] * para$pop_sz))))
stopifnot(length(para$AGE) == para$pop_sz)

if(country_id == 1){
  age_mix <- read.csv("Uganda_age_mixture_population_level.csv")
  para$age_mix <- age_mix$total_age
  para$Home_age_mix <- age_mix$home_age
}else{
  if(country_id == 2){
    AGE_matrix <- read.csv("south_africa_age.csv", header = F)
    HOME_age_matrix <- read.csv("home_south_africa_age.csv", header = F)
  }else if(country_id == 3){
    AGE_matrix <- read.csv("kenya_age.csv", header = F) # kenya
    HOME_age_matrix <- read.csv("home_kenya_age.csv", header = F) # kenya
  }
  else if(country_id == 4){
    AGE_matrix <- read.csv("Nigeria_age.csv", header = F) # nigeria
    HOME_age_matrix <- read.csv("home_Nigeria_age.csv", header = F) # Nigeria
  }
  # the age mixing matrix need to collpase into 3*3 since we only have age distribution for 3 groups.
  ##  ##  ##  ## 
  ## this step is important: the contact matrix is not symetrical: row and columns represent different things
  ##  ##  ##  ## 
  # the columns represent contact, we take colSums
  # the rows are the participant, we take the average number of contact for the participant: rowMeans
  age <- cbind(rowMeans(AGE_matrix[,1:3]),rowMeans(AGE_matrix[,4:5]),rowMeans(AGE_matrix[,6:16]))
  AGE_matrix <- rbind(colSums(age[1:3,]),colSums(age[4:5,]),colSums(age[6:16,]))
  # weight each column by the age distribution: we have to change the matrix to reflect the age distribution in the population
  AGE_matrix <- (as.matrix(rep(1,length(para$age_dist))) %*% t(para$age_dist)) * AGE_matrix
  
  # matrix in the end should be the number of contact for each age group pair 
  AGE_matrix <- (AGE_matrix + t(AGE_matrix))/2
  para$age_mix <- as.matrix(AGE_matrix)[which(upper.tri(AGE_matrix,diag = T))]
  
  age <- cbind(rowMeans(HOME_age_matrix[,1:3]),rowMeans(HOME_age_matrix[,4:5]),rowMeans(HOME_age_matrix[,6:16]))
  HOME_age_matrix <- rbind(colSums(age[1:3,]),colSums(age[4:5,]),colSums(age[6:16,]))
  HOME_age_matrix <- (as.matrix(rep(1,length(para$age_dist))) %*% t(para$age_dist)) * HOME_age_matrix
  # matrix in the end should be a triangle one 
  HOME_age_matrix <- (HOME_age_matrix + t(HOME_age_matrix))/2
  para$Home_age_mix <- as.matrix(HOME_age_matrix)[which(upper.tri(HOME_age_matrix,diag = T))]
}

# adjust the age matrix to represent the specified household contact rate
para$age_mix <- para$Home_age_mix + (para$age_mix - para$Home_age_mix) * 
  sum(para$Home_age_mix)/sum(para$age_mix - para$Home_age_mix) *  (1-para$percent_HH_cc)/para$percent_HH_cc

# adjust R0 for 1) young people susceptibility 2) subclinical cases
# adjust for the social distancing
para$num_cc_scdst <- para$num_cc * ((1-para$percent_HH_cc)*para$Non_HH_CC_rate + para$percent_HH_cc) # reduce the number of cc
para$age_mix_scdst <- para$Home_age_mix + (para$age_mix - para$Home_age_mix) * para$Non_HH_CC_rate
para$percent_HH_cc_scdst <- para$percent_HH_cc/((1-para$percent_HH_cc)*para$Non_HH_CC_rate + para$percent_HH_cc)

#####################
# generate synthetic population
source(paste0(code_version,"_functions.R"),print.eval =F)

# generate the contact networks, each one represent a uniquely sampled population (households are formed via sampling the family size distribution)
if(setting ==1 ){
  searched_clustering_number <- 200 
}else if(setting ==2 ){
  searched_clustering_number <- 100
}else if(setting ==3 ){
  searched_clustering_number <- 100
}
print("Construct the synthetic population: this will take a few minutes")
nw_para <- network_generate(para, searched_clustering_number)
searched_clustering_number <- nw_para[[3]] # reduce the searched number if it did not converge
ego.sim100 <- simulate(nw_para[[1]],nsim = 100)
sim.stats <- attr(ego.sim100,"stats")
trgt <- rbind(colMeans(sim.stats), nw_para[[1]]$target.stats)
deviation_target_statistics <- mean(abs(trgt[1,] - trgt[2,])/trgt[2,])
print("Plot one sample of the contact network")
png(filename = paste0("network_example_",community_setting,"_set.png"), width = 1000, height = 1000)
plot(simulate(nw_para[[1]]))
dev.off()

NW_SIM <- list(nw_para, deviation_target_statistics, nw_para[[3]])
save(NW_SIM, file = paste0("network_example_",community_setting,"_set.Rdata"))

###########################################
# simulate a few outbreaks
###########################################
len_sim <- 100 # length of simulation
num_rep <- 20
para$E_0 <- 2 # number of initial importation
para$infect_sz <- (-1)*para$pop_sz/200 # containment intervention only starts after X per 1000 individuals are already infected

R_UK <- 2.7
num_cc_UK <- 13
R0_baseline <- R_UK/num_cc_UK # the R0 is measured in affluent region 
# Transmission parameter & tool availability
if (para$setting==1) {
  para$symptom_report <- 0.7 # precentage infected report symptom
  para$R0 <- R0_baseline * para$num_cc # R naught
  para$theta <- 0.1 # quarantine effect: precentage of remaining transmission quarantined individual have
  para$pcr_available <- -1 # daily maximum PCR tests per 1000 population. 
  para$ppe_coef <- 1 # if people wear ppe, then their transmissibility will be different outside their family
}else if (para$setting==2) {
  para$symptom_report <- 0.7
  para$R0 <- R0_baseline * para$num_cc
  para$theta <- 0.1 
  para$pcr_available <-1000*para$pop_sz/1000
  para$ppe_coef <- 1 
}else if (para$setting==3) {
  para$symptom_report <- 0.6
  para$R0 <- R0_baseline * para$num_cc
  para$theta <- 0.2 
  para$pcr_available <-2*para$pop_sz/1000
  para$ppe_coef <- 1 
}else print ("Parameter setting error")

# subclinical rate
para$sub_clini_rate <- 0.3
para$asym_rate <- 0.2 # transmission rate of asymptomatic patient

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
para$qrate_ab <- 0.9 # CC quarantine compliance rate based on positive ab test
para$qrate_pcr <- 0.9 # CC quarantine compliance rate based on positive pcr test

# Parameters about testing tools 
para$ab_rate <- function(x, t_onset) 1/(1 + exp(7+t_onset-x)) # seroconversion rate from infection day, based on the clinical paper from Yumei Wen
para$sensitivity_ab <- 0.9 # ab test sensitivity
para$sensitivity_pcr <- 0.999 # pcr test sensitivity
para$samplefailure_pcr <- 0.3 # pcr sampling failure
para$sensitivity_antig <- 0.9 # antigen sensitivity is 0.8

#########################################################################
# adjust R0 for 1) young people susceptibility 2) subclinical cases
#########################################################################

# ajust for R0
# norm1 <- (sum(para$age_mix) - para$age_mix[c(1)]/2- sum(para$age_mix[c(2,4)])/4)/sum(para$age_mix)
# using next generation matrix to compute norm1; only kept terms connected with young people susceptibility
Cyy <- para$age_mix[c(1)] # number of comtact for young <--> young
Coy <- sum(para$age_mix[c(2,4)]) # number of comtact for young <--> old
Coo <- sum(para$age_mix[c(3,5,6)]) # number of comtact for old <--> old
Ny <- para$age_dist[1] # number of young people, Cyy/Ny is the average number of young contact for a yong person
No <- sum(para$age_dist[c(2,3)]) # number of young people, Cyy/Ny is the average number of young contact for a yong person
y_sus <- 0.5 # susceptability of young person
NGM <- matrix(c(y_sus * Cyy/Ny,  Coy/Ny, y_sus * Coy/No, Coo/No) , nrow = 2) 
trNGM <- sum(diag(NGM))
detNGM <- det(NGM)
Spectral_radius_half <- trNGM + (trNGM^2 - 4*detNGM )^0.5

y_sus <- 1 # susceptability of young person
NGM <- matrix(c(y_sus * Cyy/Ny,  Coy/Ny, y_sus * Coy/No, Coo/No) , nrow = 2) 
trNGM <- sum(diag(NGM))
detNGM <- det(NGM)
Spectral_radius_1 <- trNGM + (trNGM^2 - 4*detNGM )^0.5

norm1 <- Spectral_radius_half/Spectral_radius_1
# norm2 account for the redution of the low transmission rate of asymptomatic cases
norm2 <- 1 - para$sub_clini_rate * (1 - para$asym_rate)
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

all_new_daily <- matrix(nrow = num_rep, ncol = len_sim)
Rt <- matrix(nrow = num_rep, ncol = len_sim)
for(rp in 1:num_rep){
  est1 <- nw_para[[1]]
  
  print(paste("Running Setting=",community_setting,", Scenario 1, Rep: ",rp))
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
  print(paste("Running Setting=",community_setting,", Scenario 4, Rep: ",rp))
  
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

err_daily <- apply(all_new_daily, 2, function(x) sd(x)/sqrt(length(x)))
err_antig <- apply(all_antig_new_daily, 2, function(x) sd(x)/sqrt(length(x)))

df_plt <- data.frame(day = 1:len_sim, new_daily = colMeans(all_new_daily),
                     antig_new_daily = colMeans(all_antig_new_daily), 
                     err_daily, err_antig)
plt <- ggplot(df_plt,aes(x=day)) +
  geom_line(aes(y = new_daily, color = "Baseline")) +
  geom_line(aes(y = antig_new_daily, color = "Strategy 3: antigen RDT + PCR")) +
  geom_ribbon(aes( ymax = new_daily + 2*err_daily, ymin = new_daily - 2*err_daily), fill = "red",alpha = 0.3) +
  geom_ribbon(aes( ymax = antig_new_daily + 2*err_antig, ymin = antig_new_daily - 2*err_antig), fill = "purple",alpha = 0.3) +
  labs(x= "Days from first importing case", y = "New case per day" ) +#,title = paste0("Physical distancing imposed"))+
  scale_color_manual(values = c("Baseline" = "red",  "Strategy 3: antigen RDT + PCR" = "purple")) + 
  theme(legend.position = "None", panel.background=element_blank()) 


ggsave(paste0("transmission_setting",community_setting, "_baseline.png"),plt, width = 5, height = 4)

