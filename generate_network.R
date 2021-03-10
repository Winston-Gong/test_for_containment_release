
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
  num_rep <- 20
  senario4_yn <-F # whether no not run senario 4 (multi-ab)
}
t <- try(load(paste0("networks_set",setting,"v",v_id,"s",s_id,".Rdata")),silent=T)
if (!("try-error" %in% class(t))) {
  print(paste("Already exist networks_set: Setting=",setting,"V=",v_id,"S=",s_id))
  
  load(paste0("networks_set",setting,"v",v_id,"s",s_id,".Rdata"))
  
  nw <- NW_SIM[[1]] 
  deviation_target_statistics <- NW_SIM[[2]]
  searched_clustering_number <- NW_SIM[[3]]
  
  for(rp in 1:length(nw)){
    NW_SIM <- list(nw[[rp]], deviation_target_statistics[[rp]], searched_clustering_number)
    save(NW_SIM, file = paste0("network_",rp,"_set",setting,"v",v_id,"s",s_id,".Rdata"))
  }
  
  quit(save = "no", status = 0, runLast = TRUE)
}

library(reshape2)
library(ggplot2)
library(tidyverse)
suppressMessages(library(ergm))

###############################################
print(paste("Runing: Setting=",setting,"V=",v_id,"S=",s_id))

para <- list()
para$setting <- setting # 1=rural 2=affluent 3=slum
# Environment setup
para$pop_sz <- 1000 # 5000
para$Non_HH_CC_rate <- 1

##########################################
# parameter with age
##########################################
# South africa
# para$age_dist <- c(0.292, 0.193, 0.515)
# para$age_dist <- c(0.419, 0.195, 0.386) # kenya
# para$age_dist <- c(0.440, 0.190, 0.370) # nigeria
para$age_dist <- c(0.481, 0.203, 0.316) # Uganda

##########################################
# Household size distribution, from UN data
##########################################
para$HH_dist <- c(.11, .22, .27, .4) # Uganda
# para$HH_dist <- c(.19, .28, .3, .23) # kenya
# para$HH_dist <- c(.16, .26, .26, .32) # nigeria
# para$HH_dist <- c(.27, .35, .23, .15) # South Africa
# para$HH_affluent_dist <- c(.31, .5, .18, .02) # UK

if (para$setting==1) {
  para$num_cc <- 7 # set daily close contact to be 7 
  para$family_sz <- 5 # average household size 5
  # set the percentage of HH_cc
  para$percent_HH_cc <- .5
}else if (para$setting==2) {
  para$num_cc <- 13 
  para$family_sz <- 5 
  para$percent_HH_cc <- .3
}else if (para$setting==3) {
  para$num_cc <- 14 
  para$family_sz <- 15 
  para$percent_HH_cc <- .5
  para$HH_dist <- c(.00, .06, .17, .77) # afganistan
}else print ("Parameter setting error")


# Change parameters for sensitivity analysis
###############################################

if(v_id == 5){   # set daily close contact 
  num_cc <- c(2,4,7,10,14)
  para$num_cc <- num_cc[s_id]
}else if(v_id == 8){   # social distancing
  Non_HH_CC_rate<- c(1,.8,.6,.4,.2)
  para$Non_HH_CC_rate <- Non_HH_CC_rate[s_id] 
}else if(v_id == 12){  # simulation size
  pop_sz <- c(100, 1000, 2000)
  para$pop_sz <- pop_sz[s_id]  
}else if(v_id == 14){  # percentage of household contact
  percent_HH_cc <- c(.2,.4,.6,.8)
  para$percent_HH_cc <-  percent_HH_cc[s_id]
}

##########################################
# parameter with age
##########################################
# South africa
# para$age_dist <- c(0.292, 0.193, 0.515)
# para$age_dist <- c(0.419, 0.195, 0.386) # kenya
# para$age_dist <- c(0.440, 0.190, 0.370) # nigeria
para$age_dist <- c(0.481, 0.203, 0.316) # Uganda

##########################################
# Household size distribution, from UN data
##########################################
para$HH_dist <- c(.11, .22, .27, .4) # Uganda
# para$HH_dist <- c(.19, .28, .3, .23) # kenya
# para$HH_dist <- c(.16, .26, .26, .32) # nigeria
# para$HH_dist <- c(.27, .35, .23, .15) # South Africa
# para$HH_affluent_dist <- c(.31, .5, .18, .02) # UK

##########################################
para$AGE <- unlist(sapply(1:length(para$age_dist), function(x) rep(x,round(para$age_dist[x] * para$pop_sz))))
stopifnot(length(para$AGE) == para$pop_sz)
############################################################
# normalized the age matrix
# AGE_matrix <- read.csv("south_africa_age.csv", header = F)
# AGE_matrix <- read.csv("kenya_age.csv", header = F) # kenya
AGE_matrix <- read.csv("Uganda_age.csv", header = F) # Uganda
# AGE_matrix <- read.csv("Nigeria_age.csv", header = F) # nigeria
# the age mixing matrix need to collpase into 3*3 since we only have age distribution for 3 groups. 
age <- cbind(rowMeans(AGE_matrix[,1:3]),rowMeans(AGE_matrix[,4:5]),rowMeans(AGE_matrix[,6:16]))
AGE_matrix <- rbind(colMeans(age[1:3,]),colMeans(age[4:5,]),colMeans(age[6:16,]))
AGE_matrix <- (para$age_dist %*% t(as.matrix(para$age_dist))) * AGE_matrix
# matrix in the end should be the number of contact for each age group pair 
AGE_matrix <- (AGE_matrix + t(AGE_matrix))/2
para$age_mix <- as.matrix(AGE_matrix)[which(upper.tri(AGE_matrix,diag = T))]
############################################################
############################################################
# normalized the house-hold age matrix
# HOME_age_matrix <- read.csv("home_south_africa_age.csv", header = F)
# HOME_age_matrix <- read.csv("home_kenya_age.csv", header = F) # kenya
HOME_age_matrix <- read.csv("home_Uganda_age.csv", header = F) # Uganda
# HOME_age_matrix <- read.csv("home_Nigeria_age.csv", header = F) # Nigeria

age <- cbind(rowMeans(HOME_age_matrix[,1:3]),rowMeans(HOME_age_matrix[,4:5]),rowMeans(HOME_age_matrix[,6:16]))

HOME_age_matrix <- rbind(colMeans(age[1:3,]),colMeans(age[4:5,]),colMeans(age[6:16,]))
HOME_age_matrix <- (para$age_dist %*% t(as.matrix(para$age_dist))) * HOME_age_matrix
# matrix in the end should be a triangle one 
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


source(paste0(code_version,"_functions.R"),print.eval =F)

# generate the contact networks, each one represent a uniquely sampled population (households are formed via sampling the family size distribution)
nw <- list()
deviation_target_statistics <- list()
search_clustering_number_list <- list()
searched_clustering_number <- 500
for(rp in 1:num_rep){
  print(paste0("Vid",v_id,"s_id",s_id, " Repetition for network: ", rp))
  # MCMLE could stuck, regenerate the network if that is the case, and reduce the geographical clustering 
  nw_para <- network_generate(para, searched_clustering_number)
  searched_clustering_number <- nw_para[[3]] # reduce the searched number if it did not converge
  ###############################################
  # check the target.stats
  ego.sim100 <- simulate(nw_para[[1]],nsim = 100)
  sim.stats <- attr(ego.sim100,"stats")
  trgt <- rbind(colMeans(sim.stats), nw_para[[1]]$target.stats)
  deviation_target_statistics[[rp]] <- mean(abs(trgt[1,] - trgt[2,])/trgt[2,])
  nw[[rp]] <- nw_para
  search_clustering_number_list[[rp]] <- searched_clustering_number
}

# save all the files in separate files
for(rp in 1:length(nw)){
  NW_SIM <- list(nw[[rp]], deviation_target_statistics[[rp]], search_clustering_number_list[[rp]])
  save(NW_SIM, file = paste0("network_",rp,"_set",setting,"v",v_id,"s",s_id,".Rdata"))
}

print(paste("Completed: Setting=",setting,"V=",v_id,"S=",s_id))



