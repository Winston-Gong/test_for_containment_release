# this file contains code for performing genetic analysis on the topic loadings 


loadings %>% 
  left_join(data, by = c("eid" = "V1")) %>% 
  write.table("environment_covariates.txt", sep="\t", col.names = FALSE, row.names = FALSE)