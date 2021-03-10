###############################################
# Preparation functions
###############################################

# incubation time: shape = 7.92, scale = 0.69
incub <- function(x) floor(pmax(pmin(rgamma(x, shape = 7.92, scale = 0.69), 14),1)) # maximum incubation time is 14 days

# here we use the shape of weibull from the science paper.
# The onset should be the mode of weibull, so we compute the scale from the onset .
# weibull has the property of keeping the shape, and mode change with scale
# physical property of weibull also make sense: single event happen with an rate of t^k

g_transmissibility <- function(t_onset) para$R0_adj/para$num_cc * dweibull(1:(len_sim+1), shape = 2.826, scale = t_onset/0.8568)
# plot(g_transmissibility(20))
# weibull distribution: shape: 2.8 scale 5.7
# mean 5.665*(1-1/2.826)^(1/2.826)
# normalized transmissibility per contact para$R0_adj/para$num_cc

###############################################
# load an existing the contact network
###############################################
load_sample_nw <- function(setting,v_id,s_id, para){
  NW_SIM <- NA
  while (is.na(NW_SIM)){
    try({
      idx <- sample(1:20, 1)
      if(v_id == 5 | v_id == 8 | v_id ==  12 | v_id ==  14){
        load(paste0("Networks/network_",idx,"_set",setting,"v",v_id,"s",s_id,".Rdata"))
      }else{
        social_dst_idx <- 1+(1-para$Non_HH_CC_rate)/0.2
        load(paste0("Networks/network_",idx,"_set",setting,"v8s", social_dst_idx,".Rdata"))
      }
    })
  }
  nw <- NW_SIM[[1]]
  para_nw <- nw[[2]]
  para$family_lable <- para_nw$family_lable
  para$clustering_x <- para_nw$clustering_x
  para$clustering_y <- para_nw$clustering_y
  para$clustering_effect <- para_nw$clustering_effect
  
  # load all the network parameter, include AGE 
  para$pop_sz <- para_nw$pop_sz
  para$AGE <- para_nw$AGE
  para$age_dist <- para_nw$age_dist
  para$family_sz <- para_nw$family_sz
  para$HH_dist <- para_nw$HH_dist
  para$num_cc_scdst <- para_nw$num_cc_scds
  para$percent_HH_cc_scdst <- para_nw$percent_HH_cc_scdst 
  para$age_mix_scdst <- para_nw$age_mix_scdst 
  # non-social-distanced para
  para$age_mix <- para_nw$age_mix 
  para$Home_age_mix <- para_nw$Home_age_mix
  para$percent_HH_cc <- para_nw$percent_HH_cc
  para$num_cc <- para_nw$num_cc
  
  return(list(nw[[1]], para))
}
###############################################
# #  generate the contact network
###############################################
network_generate <- function(para, searched_clustering_number){
  nw <- network.initialize(n = para$pop_sz , directed = FALSE)
  # specify house hold
  adult_idx <- which(para$AGE == 3)
  family_core <- sample(adult_idx, para$pop_sz/para$family_sz, replace =F)# a family should have at least one adult
  # assign other adults randomly 
  non_core_adult_idx <- adult_idx[! adult_idx %in% family_core]
  # family_adults <- split(non_core_adult_idx, sample(para$pop_sz/para$family_sz, length(non_core_adult_idx), replace = TRUE) )
  
  non_core_idx <- sample(c(which(para$AGE <= 2), non_core_adult_idx))
  # comupte the split point for family size
  family_sz_label <- cumsum(para$HH_dist)
  
  # assigning family labels
  family_lable <- rep(NA, para$pop_sz)
  family_lable[family_core] <- 1:(para$pop_sz/para$family_sz)
  
  # using a rolling rotation to creat family based on the household size distribution
  for(sz_gp in 1:2){
    # for each group 2-3/4-5, I assign half of family group 2-3 to be 2 and the other half to be 3
    idx1 <- floor(family_sz_label[sz_gp] * para$pop_sz/para$family_sz):(para$pop_sz/para$family_sz)
    family_lable[non_core_idx[1:(idx1[length(idx1)] - idx1[1] + 1)]] <-  idx1
    non_core_idx <- non_core_idx[(1+length(idx1)):length(non_core_idx)]
    # pick half of thsese families to assiagn another member: the result will be in the family of size (2-3, that's the data we have unfortunately), 
    # half of then are havign 2 members and the other half have 3 members
    idx2 <- floor((family_sz_label[sz_gp] + family_sz_label[sz_gp+1])/2 * para$pop_sz/para$family_sz):(para$pop_sz/para$family_sz)
    family_lable[non_core_idx[1:(idx2[length(idx2)] - idx2[1] + 1)]] <-  idx2
    non_core_idx <- non_core_idx[(1+length(idx2)):length(non_core_idx)]
  }
  # assign the remaining people randomly to those family with more than 5 people
  family_lable[non_core_idx] <- sample(floor(family_sz_label[3] * para$pop_sz/para$family_sz + 1):(para$pop_sz/para$family_sz),length(non_core_idx), replace = T)
  
  # create social clusters, for school and work space, increase clustering coefficient
  clustering_x <- runif(para$pop_sz/para$family_sz)# rnorm(para$pop_sz/para$family_sz, mean = 0, sd = 1)
  # clustering_x <-  rnorm(para$pop_sz/para$family_sz)
  clustering_x <- clustering_x[family_lable] # make sure a family share the same clustering 
  
  clustering_y <- runif(para$pop_sz/para$family_sz)# rnorm(para$pop_sz/para$family_sz, mean = 0, sd = 1)
  #clustering_y <-  rnorm(para$pop_sz/para$family_sz)
  clustering_y <- clustering_y[family_lable] # make sure a family share the same clustering 
  
  # this family and age label will be add to the population
  nw <- set.vertex.attribute(nw, "family", family_lable)
  nw <- set.vertex.attribute(nw, "age", para$AGE)
  # nw <- set.vertex.attribute(nw, "clustering", c(clustering_x,clustering_y))
  nw <- set.vertex.attribute(nw, "clustering_x", clustering_x)
  nw <- set.vertex.attribute(nw, "clustering_y", clustering_y)
  
  est1 <- NA
  grid_id <- 1
  # perform grid search to get the local clustering coefficient
  while(is.na(est1)){
    clustering_effect <- (para$pop_sz/2*para$num_cc_scdst)*(1 - para$percent_HH_cc_scdst)/searched_clustering_number #  14 (25) 23 (20) 
    target.stats <- c(para$pop_sz/2 * para$num_cc_scdst, para$pop_sz/2 * para$num_cc_scdst * para$percent_HH_cc_scdst,
                      para$age_mix_scdst/sum(para$age_mix_scdst) * para$pop_sz/2 * para$num_cc_scdst, clustering_effect, clustering_effect)
    try(suppressMessages( est1 <- ergm(nw ~ edges  + nodematch("family") + mm("age") + absdiff("clustering_x", pow=2) + absdiff("clustering_y", pow=2),
                                       target.stats = target.stats,control = control.ergm(MCMLE.maxit = 400, SAN.maxit = 200)) ))
    searched_clustering_number <- max(6, searched_clustering_number* 0.9) # reduce the searched number if it did not converge; 6 is when there is no geographical clustering effect
    print(paste0("search ", grid_id, " for geographical clustering"))
    grid_id <- grid_id + 1
  }

  para$family_lable <- family_lable
  para$clustering_x <- clustering_x
  para$clustering_y <- clustering_y
  para$clustering_effect <- clustering_effect
  ###################################################################
  # generate a number of networks for analysis (should be an input)
  ###################################################################
  return(list(est1,para,searched_clustering_number/0.9))
}

###############################################
# Transmission update functions
###############################################
# update C
# the function below provide a way to update daily contact 
C_update <- function(C, Q, est_nw, para){
  # one day pass; C is 12 when just infected, then C decrease by 1 everyday pass
  C <- pmax(C-1, 0)
  # then update close contact of this day. we assume the close contact to happen between the people around
  # cc store the distance from the close contact to index case for each close contact of a particular day
  edge_lst <- as.edgelist(simulate(est_nw))
  prob_cc <- ifelse(!is.na(Q[edge_lst[,1]]), para$theta, 1) * ifelse(!is.na(Q[edge_lst[,2]]),para$theta, 1)
  true_cc <- (runif(dim(edge_lst)[1]) < prob_cc)
  C[cbind(edge_lst[true_cc,1], edge_lst[true_cc,2])] <- 12
  return(C)
}
# update the Infected individual & keep an record of the incubation time
I_O_update <- function(I, Q, C, O, Z, trace_inf_n, para){
  I <- I + 1
  for(i in 1:para$pop_sz){
    # i is the infection source, idx is the infected
    # cc_idx find the infected individual and he has contact someone that day 
      # (our model is almost certain to have contact, unless he is quarantined)
      # C[-i,i] "-i" is to avoid count the index case himself
    if( !is.na(I[i]) & (length(which(C[-i,i] == 12)) | length(which(C[i,-i] == 12)) ) ){
      cc_idx <- union(which(C[-i,i] == 12), which(C[i,-i] == 12))
      u <- runif(length(cc_idx))
      # get the individual transmissibility, it need to be normalized by the number of individuals
      transmissibility <- g_transmissibility(O[i]) * ifelse(Z[i] == 1, para$asym_rate, 1)
      # pick the index of individual 1) who is not yet infected & 2) transmited
      # the expected transimisstion rate need to be normalized by daily contact 2*num_cc
      # reduce the sesceptability for young people (AGE == 1) by half
      susceptibility <- ifelse(para$AGE[cc_idx] == 1, .5, 1)
      u1 <- runif(length(cc_idx))
      idx <- cc_idx[(u < transmissibility[I[i]]) & is.na(I[cc_idx]) & u1 < susceptibility]
      I[idx] <- 0 # consider the incubation time
      u1 <- runif(length(idx))
      # if Z=2, the infected case is detectable, if Z=1, subclinical, if Z=F, not reported and infected
      Z[idx] <- sapply(u1, function(x) if(x<para$symptom_report*(1-para$sub_clini_rate)) 2 else if(x<(1-para$sub_clini_rate)) F else 1) 
      O[idx] <- incub(length(idx)) # incubation duration
      # trace the Rt
      trace_inf_n[i] <- trace_inf_n[i] + length(idx)
    }
  }
  return(list(I, O, Z,trace_inf_n))
}
########################################################
# Quarantine functions
########################################################
# this function will quarantine after ab tests 
# flg_ab <- FALSE
# flg_multi_ab <- TRUE
# flg_PCR <- TRUE
########################################
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
          Q[cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_symptom]] <- 1
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
              Q[cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_pcr]] <- 1
              
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
            Q[cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_ab]] <- 1
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
            Q[cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_ab]] <- 1
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

