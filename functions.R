##############################################
#                 functions                  #
#  Firkowski et al. - Ecological Monographs  #
##############################################

# loop id function
ijloop <- function(idloop=loop_id, Comnum=numCom){
  ((idloop*Comnum)-(Comnum-1)):(idloop*Comnum)
}

# beta diversity function
# reference: Ricotta 2017
betaeve <- function(dataset){
  library(vegan)
  n_col <- ncol(dataset)
  n_row <- nrow(dataset)
  total <- apply(dataset, 1, sum)
  species_weights <- total/sum(dataset)
  rel_abu <- sweep(dataset, 1, total, "/")
  rel_abu[is.na(rel_abu)] <- 0
  h_shannon <- vegan::diversity(rel_abu, index = "shannon", MARGIN = 1)
  pielou <- h_shannon/log(n_col)
  pielou_com <- 1-sum(pielou/n_row)
  weighted_pielou <- 1-sum(pielou*species_weights)
  output <- list("Average Beta"=pielou_com,"Species weights"=species_weights, "Weighted Beta" = weighted_pielou)
  return(output)
}

# synchrony function
# reference: function "synch_onerep" modified from codyn package
synchrony <- function(df, time.var, species.var, abundance.var, metric = "Gross"){
  metric = match.arg(metric, choices = c("Gross")) # removed "Loreau" from choices
  # remove any species that were never present
  df <- subset(df, abundance.var > 0)
  # fill in 0
  spplist <- unique(df[species.var])
  # added dimension of species lists #
  sppsize <- dim(spplist)[1]
  # change ends here #
  yearlist <- unique(df[time.var])
  fulllist <- expand.grid(species.var = spplist[[species.var]],
                          time.var = yearlist[[time.var]])
  # recapture original names
  names(fulllist) = c(species.var, time.var)
  df2 <- merge(df[c(species.var, time.var, abundance.var)], fulllist, all.y = T)
  df2[is.na(df2)] <- 0
  # removed if statement for metric="Loreau" #
  corout <- as.data.frame(cbind(species.var = as.character(),
                                "sppcor" =  as.numeric()))
  # check to see if there are species which do not vary within a subplot
  nonvary <- apply(transpose_community(df, time.var, species.var, abundance.var), 2, sd)
  # added NA to 0 #
  #nonvary[is.na(nonvary)] <- 0
  # change ends here #
  if(any(nonvary == 0 | is.na(nonvary))) {
    warning("One or more species has non-varying abundance within a subplot and has been omitted")
    # remove non-varying species from "spplist" for this subplot
    spplist <- data.frame(species =
                            spplist[[1]][is.na(match(as.character(unlist(spplist)), names(nonvary)[which(nonvary == 0)]))]
    )
  }
  # added #
  # return NA
  if(dim(spplist)[1]==0){
    corout <- as.data.frame(cbind(1:sppsize, NA))
  } else {
  # change ends here #  
    for (i in 1:nrow(spplist)){
      myspp <- as.character(spplist[[1]][i])
      focalspp <- df2[which(df2[species.var] == myspp),]
      com.focalspp <- transpose_community(focalspp, time.var, species.var, abundance.var)
      otherspp <- df2[which(df2[species.var] != myspp),]
      com.otherspp <- transpose_community(otherspp, time.var, species.var, abundance.var)
      agg.otherspp <- rowSums(com.otherspp)
      sppcor <- stats::cor(agg.otherspp, com.focalspp)
      subout <- as.data.frame(cbind(myspp, sppcor))
      names(subout) = c(species.var, "sppcor")
      subout$sppcor <- as.numeric(as.character(subout$sppcor))
      corout <- rbind(corout, subout)
    }
  }
  # average correlation for the community
  synchrony <- mean(corout$sppcor)
  return(synchrony)
}
# transpose community
transpose_community <- function(df, time.var, species.var, abundance.var){
  df <- as.data.frame(df)
  # remove unused levels if species is a factor
  df[species.var] <- if(is.factor(df[[species.var]]) == TRUE){factor(df[[species.var]])} else {df[species.var]}
  # sort by time and species
  df <- df[order(df[[time.var]], df[[species.var]]),]
  # cast as a species x time data frame; insert NA for 0
  comdat <- tapply(df[[abundance.var]], list(df[[time.var]], as.vector(df[[species.var]])), sum)
  comdat[is.na(comdat)] <- 0
  comdat <- as.data.frame(comdat)
  # results
  return(comdat)
}

# variability function
# reference: Wang et al. 2019
var.partition <- function(metacomm_tsdata){
  ts_metacom <- apply(metacomm_tsdata,2,sum)
  ts_patch <-  apply(metacomm_tsdata,c(2,3),sum)
  ts_species <- apply(metacomm_tsdata,c(1,2),sum)
  sd_metacom <- sd(ts_metacom)
  sd_patch_k <- apply(ts_patch,2,sd)
  sd_species_i <- apply(ts_species,1,sd)
  sd_species_patch_ik <- apply(metacomm_tsdata,c(1,3),sd)
  mean_metacom <- mean(ts_metacom)
  CV_S_L <- sum(sd_species_patch_ik)/mean_metacom
  CV_C_L <- sum(sd_patch_k)/mean_metacom
  CV_S_R <- sum(sd_species_i)/mean_metacom
  CV_C_R <- sd_metacom/mean_metacom
  phi_S_L2R <- CV_S_R/CV_S_L
  phi_C_L2R <- CV_C_R/CV_C_L
  phi_S2C_L <- CV_C_L/CV_S_L
  phi_S2C_R <- CV_C_R/CV_S_R
  partition_3level <- c(CV_S_L=CV_S_L, CV_C_L=CV_C_L, CV_S_R=CV_S_R, CV_C_R=CV_C_R,
                        phi_S_L2R=phi_S_L2R, phi_C_L2R=phi_C_L2R, phi_S2C_L=phi_S2C_L,
                        phi_S2C_R=phi_S2C_R)
  return(partition_3level)
}

# weighted means
weighted_mean <- function(x, w, ..., na.rm = FALSE){
  if(na.rm){
    df_omit <- na.omit(data.frame(x, w))
    return(weighted.mean(df_omit$x, df_omit$w, ...))
  } 
  weighted.mean(x, w, ...)
}

# Gross' synchrony & variability partitioning function 
syn_CV_partitioning <- function(X_saved){
  # data frame for all results
  synchrony_data <- data.frame(replicates=as.numeric(), environment=as.numeric(), 
                               patch=as.numeric(), species=as.character(),
                               scale=as.character(), trophic=as.character(),
                               synchrony=as.numeric(), variability=as.numeric())
  ### local species synchrony ###
  # note:
  # local species needs to be weighted
  # parameters
  environment <- ePeriod
  patches <- numCom
  # data
  # all species
  all_data <- X_saved
  # number of species per trophic level
  num_pred <- nSp-2 # change
  num_prey <- nSp-4 # change
  num_nprey <- nSp-3 # change
  # data frame
  localspp_out <- data.frame(replicates=as.numeric(), environment=as.numeric(), 
                             patch=as.numeric(), species=as.character(),
                             scale=as.character(), trophic=as.character(),
                             synchrony=as.numeric(), variability=as.numeric())
  # loop
  for(p in 1:patches){
    # define patch
    pp <- as.numeric(p)
    # subset data
    d1 <- subset(all_data, patch==pp) 
    data <- subset(d1, select=c(time, replicate, (6:14)))
    # define community
    community_all <- data.frame(time=rep(NA,(dim(data)[1]*nSp)), species=NA, abundances=NA, replicates=NA)
    community_pred <- data.frame(time=rep(NA,(dim(data)[1]*num_pred)), species=NA, abundances=NA, replicates=NA)
    community_prey <- data.frame(time=rep(NA,(dim(data)[1]*num_prey)), species=NA, abundances=NA, replicates=NA)
    community_nprey <- data.frame(time=rep(NA,(dim(data)[1]*num_nprey)), species=NA, abundances=NA, replicates=NA)
    #
    community_all$time <- community_pred$time <- community_prey$time <- community_nprey$time <- data$time
    community_all$replicates <- community_pred$replicates <- community_prey$replicates <- community_nprey$replicates <- data$replicate                   
    #
    community_all$species <- rep(c("V1","V2","V3","V4","V5","V6","V7","V8","V9"), each=dim(data)[1])
    community_all$abundances <- c(data$V1,data$V2,data$V3,data$V4,data$V5,data$V6,data$V7,data$V8,data$V9)
    #
    #community_pred$species <- rep(c("V8","V9"), each=dim(data)[1])
    #community_pred$abundances <- c(data$V8,data$V9)
    community_pred$species <- rep(c("V1","V2","V3","V4","V5","V6","V7"), each=dim(data)[1]) # change
    community_pred$abundances <- c(data$V1,data$V2,data$V3,data$V4,data$V5,data$V6,data$V7) # change
    #
    #community_prey$species <- rep(c("V4","V5","V6","V7"), each=dim(data)[1])
    #community_prey$abundances <- c(data$V4,data$V5,data$V6,data$V7)
    community_prey$species <- rep(c("V1","V2","V3","V8","V9"), each=dim(data)[1]) # change
    community_prey$abundances <- c(data$V1,data$V2,data$V3,data$V8,data$V9) # change
    #
    #community_nprey$species <- rep(c("V1","V2","V3"), each=dim(data)[1])
    #community_nprey$abundances <- c(data$V1,data$V2,data$V3)
    community_nprey$species <- rep(c("V4","V5","V6","V7","V8","V9"), each=dim(data)[1]) # change
    community_nprey$abundances <- c(data$V4,data$V5,data$V6,data$V7,data$V8,data$V9) # change
    # calculate synchrony
    # for all species
    Synchrony_all <- synchrony(community_all, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    # for predators
    Synchrony_pred <- synchrony(community_pred, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    # for prey
    Synchrony_prey <- synchrony(community_prey, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    # for non-prey
    Synchrony_nprey <- synchrony(community_nprey, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    # outputs
    out_all <- data.frame(replicates=r, environment=ePeriod, 
                          patch=pp, species=NA, scale="localspp", trophic="all",
                          synchrony=Synchrony_all, variability=NA)
    out_pred <- data.frame(replicates=r, environment=ePeriod, 
                           patch=pp, species=NA, scale="localspp", trophic="pred",
                           synchrony=Synchrony_pred, variability=NA)
    out_prey <- data.frame(replicates=r, environment=ePeriod, 
                           patch=pp, species=NA, scale="localspp", trophic="prey",
                           synchrony=Synchrony_prey, variability=NA)
    out_nprey <- data.frame(replicates=r, environment=ePeriod, 
                            patch=pp, species=NA, scale="localspp", trophic="nprey",
                            synchrony=Synchrony_nprey, variability=NA)
    # outputs combined
    localspp_out <- rbind(localspp_out, out_all, out_pred, out_prey, out_nprey)
  }
  # calculate weights for synchrony
  localspp_all_weights_S <- as.vector(unlist(tapply(localspp_out[localspp_out$trophic=="all",7], localspp_out[localspp_out$trophic=="all",3], function(x) (abs(x)/sum(abs(x))))))
  localspp_pred_weights_S <- as.vector(unlist(tapply(localspp_out[localspp_out$trophic=="pred",7], localspp_out[localspp_out$trophic=="pred",3], function(x) (abs(x)/sum(abs(x))))))
  localspp_prey_weights_S <- as.vector(unlist(tapply(localspp_out[localspp_out$trophic=="prey",7], localspp_out[localspp_out$trophic=="prey",3], function(x) (abs(x)/sum(abs(x))))))
  localspp_nprey_weights_S <- as.vector(unlist(tapply(localspp_out[localspp_out$trophic=="nprey",7], localspp_out[localspp_out$trophic=="nprey",3], function(x) (abs(x)/sum(abs(x))))))
  # weight data for synchrony
  localspp_all_weighted_S <- weighted_mean(localspp_out[localspp_out$trophic=="all",7], localspp_all_weights_S, na.rm=T)
  localspp_pred_weighted_S <- weighted_mean(localspp_out[localspp_out$trophic=="pred",7], localspp_pred_weights_S, na.rm=T)
  localspp_prey_weighted_S <- weighted_mean(localspp_out[localspp_out$trophic=="prey",7], localspp_prey_weights_S, na.rm=T)
  localspp_nprey_weighted_S <- weighted_mean(localspp_out[localspp_out$trophic=="nprey",7], localspp_nprey_weights_S, na.rm=T)
  # final output for local species results
  localspp_data <- data.frame(replicates=r, environment=ePeriod, 
                              patch=NA, species=NA, scale="localspp", trophic=c("all","pred","prey","nprey"),
                              synchrony=c(localspp_all_weighted_S, localspp_pred_weighted_S, localspp_prey_weighted_S, localspp_nprey_weighted_S),
                              variability=NA)
  # combine to return object
  synchrony_data <- rbind(synchrony_data, localspp_data)
  
  
  ### regional species synchrony ###
  # create metacommunity abundance data
  reg_list <- split(X_saved[,6:14], list(X_saved$time, X_saved$fluctuation, X_saved$replicate))
  reg_sum <- as.data.frame(do.call(rbind,(lapply(names(reg_list), function(x) apply(reg_list[[x]], 2, sum)))))
  treats <- data.frame(str_split_fixed(names(reg_list), "[.]", 3))
  reg_abund <- cbind(treats, reg_sum)
  colnames(reg_abund)[colnames(reg_abund)=="X1"] <- "time"
  colnames(reg_abund)[colnames(reg_abund)=="X2"] <- "environment"
  colnames(reg_abund)[colnames(reg_abund)=="X3"] <- "r"
  # data
  data_reg <- reg_abund
  # "loop"
  # define community
  community_all <- data.frame(time=rep(NA,dim(data_reg)[1]*nSp), species=NA, abundances=NA, replicates=NA)
  community_pred <- data.frame(time=rep(NA,dim(data_reg)[1]*num_pred), species=NA, abundances=NA, replicates=NA)
  community_prey <- data.frame(time=rep(NA,dim(data_reg)[1]*num_prey), species=NA, abundances=NA, replicates=NA)
  community_nprey <- data.frame(time=rep(NA,dim(data_reg)[1]*num_nprey), species=NA, abundances=NA, replicates=NA)
  #
  community_all$time <- community_pred$time <- community_prey$time <- community_nprey$time <- as.numeric(data_reg$time)
  community_all$replicates <- community_pred$replicates <- community_prey$replicates <- community_nprey$replicates <- data_reg$r                
  #
  community_all$species    <- rep(c("V1","V2","V3","V4","V5","V6","V7","V8","V9"), each=dim(data_reg)[1])
  community_all$abundances <- c(data_reg$V1,data_reg$V2,data_reg$V3,data_reg$V4,data_reg$V5,data_reg$V6,data_reg$V7,data_reg$V8,data_reg$V9)
  #
  #community_pred$species    <- rep(c("V8","V9"), each=dim(data_reg)[1])
  #community_pred$abundances <- c(data_reg$V8,data_reg$V9)
  community_pred$species <- rep(c("V1","V2","V3","V4","V5","V6","V7"), each=dim(data_reg)[1]) # change
  community_pred$abundances <- c(data_reg$V1,data_reg$V2,data_reg$V3,data_reg$V4,data_reg$V5,data_reg$V6,data_reg$V7) # change
  #
  #community_prey$species    <- rep(c("V4","V5","V6","V7"), each=dim(data_reg)[1])
  #community_prey$abundances <- c(data_reg$V4,data_reg$V5,data_reg$V6,data_reg$V7)
  community_prey$species <- rep(c("V1","V2","V3","V8","V9"), each=dim(data_reg)[1]) # change
  community_prey$abundances <- c(data_reg$V1,data_reg$V2,data_reg$V3,data_reg$V8,data_reg$V9) # change
  #
  #community_nprey$species    <- rep(c("V1","V2","V3"), each=dim(data_reg)[1])
  #community_nprey$abundances <- c(data_reg$V1,data_reg$V2,data_reg$V3)
  community_nprey$species <- rep(c("V4","V5","V6","V7","V8","V9"), each=dim(data_reg)[1]) # change
  community_nprey$abundances <- c(data_reg$V4,data_reg$V5,data_reg$V6,data_reg$V7,data_reg$V8,data_reg$V9) # change
  # calculate synchrony
  # for all species
  Synchrony_all <- synchrony(community_all, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  # for predators
  Synchrony_pred <- synchrony(community_pred, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  # for prey
  Synchrony_prey <- synchrony(community_prey, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  # for non-prey
  Synchrony_nprey <- synchrony(community_nprey, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  # outputs
  out_all <- data.frame(replicates=r, environment=ePeriod, 
                        patch=NA, species=NA, scale="regionalspp", trophic="all",
                        synchrony=Synchrony_all, variability=NA)
  out_pred <- data.frame(replicates=r, environment=ePeriod, 
                         patch=NA, species=NA, scale="regionalspp", trophic="pred",
                         synchrony=Synchrony_pred, variability=NA)
  out_prey <- data.frame(replicates=r, environment=ePeriod, 
                         patch=NA, species=NA, scale="regionalspp", trophic="prey",
                         synchrony=Synchrony_prey, variability=NA)
  out_nprey <- data.frame(replicates=r, environment=ePeriod, 
                          patch=NA, species=NA, scale="regionalspp", trophic="nprey",
                          synchrony=Synchrony_nprey, variability=NA)
  # combine to return object
  synchrony_data <- rbind(synchrony_data, out_all, out_pred, out_prey, out_nprey)
  
  
  ### species spatial synchrony ###
  # note:
  # species spatial needs to be weighted
  # setting parameters
  species <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9")
  # dataframe
  sppspatial_out <- data.frame(replicates=as.numeric(), environment=as.numeric(), 
                               patch=as.numeric(), species=as.character(),
                               scale=as.character(), trophic=as.character(),
                               synchrony=as.numeric(), variability=as.numeric())
  # loop
  for(spp in 6:14){
    # community dataframe
    community <- data.frame(time=rep(NA,(dim(data)[1]*5)), species=NA, abundances=NA, replicates=NA)
    # subset data
    data <- subset(all_data, select=c(time, patch, replicate, spp))
    # define community
    community$time       <- data$time
    community$species    <- data$patch
    community$abundances <- data[,4]
    community$replicates <- data$replicate                   
    # calculate synchrony
    Synchrony_spp <- synchrony(community, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    # output
    out_spp <- data.frame(replicates=r, environment=ePeriod, 
                          patch=NA, species=spp, scale="sppspatial", trophic=NA,
                          synchrony=Synchrony_spp, variability=NA)
    # combine outputs
    sppspatial_out <- rbind(sppspatial_out, out_spp)
  }
  # calculate weights for synchrony
  sppspatial_all_weights_S <- as.vector(unlist(tapply(sppspatial_out[,7], sppspatial_out[,4], function(x) (abs(x)/sum(abs(x))))))
  #sppspatial_pred_weights_S <- as.vector(unlist(tapply(sppspatial_out[c(8,9),7], sppspatial_out[c(8,9),4], function(x) (abs(x)/sum(abs(x))))))
  #sppspatial_prey_weights_S <- as.vector(unlist(tapply(sppspatial_out[c(4:7),7], sppspatial_out[c(4:7),4], function(x) (abs(x)/sum(abs(x))))))
  #sppspatial_nprey_weights_S <- as.vector(unlist(tapply(sppspatial_out[c(1:3),7], sppspatial_out[c(1:3),4], function(x) (abs(x)/sum(abs(x))))))
  sppspatial_pred_weights_S <- as.vector(unlist(tapply(sppspatial_out[c(1:7),7], sppspatial_out[c(1:7),4], function(x) (abs(x)/sum(abs(x)))))) # change
  sppspatial_prey_weights_S <- as.vector(unlist(tapply(sppspatial_out[c(1:3,8,9),7], sppspatial_out[c(1:3,8,9),4], function(x) (abs(x)/sum(abs(x)))))) # change
  sppspatial_nprey_weights_S <- as.vector(unlist(tapply(sppspatial_out[c(4:9),7], sppspatial_out[c(4:9),4], function(x) (abs(x)/sum(abs(x)))))) # change
  # weight data for synchrony
  sppspatial_all_weighted_S <- weighted_mean(sppspatial_out[,7], sppspatial_all_weights_S, na.rm=T)
  #sppspatial_pred_weighted_S <- weighted_mean(sppspatial_out[c(8,9),7], sppspatial_pred_weights_S, na.rm=T)
  #sppspatial_prey_weighted_S <- weighted_mean(sppspatial_out[c(4:7),7], sppspatial_prey_weights_S, na.rm=T)
  #sppspatial_nprey_weighted_S <- weighted_mean(sppspatial_out[c(1:3),7], sppspatial_nprey_weights_S, na.rm=T)
  sppspatial_pred_weighted_S <- weighted_mean(sppspatial_out[c(1:7),7], sppspatial_pred_weights_S, na.rm=T) # change
  sppspatial_prey_weighted_S <- weighted_mean(sppspatial_out[c(1:3,8,9),7], sppspatial_prey_weights_S, na.rm=T) # change
  sppspatial_nprey_weighted_S <- weighted_mean(sppspatial_out[c(4:9),7], sppspatial_nprey_weights_S, na.rm=T) # change
  # final output for local species results
  sppspatial_data <- data.frame(replicates=r, environment=ePeriod, 
                                patch=NA, species=NA, scale="sppspatial", trophic=c("all","pred","prey","nprey"),
                                synchrony=c(sppspatial_all_weighted_S, sppspatial_pred_weighted_S, sppspatial_prey_weighted_S, sppspatial_nprey_weighted_S),
                                variability=NA)
  # combine to return object
  synchrony_data <- rbind(synchrony_data, sppspatial_data)
  
  
  ### community spatial synchrony ###
  # data
  all_data$comm <- apply(X_saved[,6:14], 1, sum)
  #all_data$comm_pred <- apply(X_saved[,13:14], 1, sum)
  #all_data$comm_prey <- apply(X_saved[,9:12], 1, sum)
  #all_data$comm_nprey <- apply(X_saved[,6:8], 1, sum)
  all_data$comm_pred <- apply(X_saved[,6:12], 1, sum) # change
  all_data$comm_prey <- apply(X_saved[,c(6:8,13,14)], 1, sum) # change
  all_data$comm_nprey <- apply(X_saved[,9:14], 1, sum) # change
  # community data frame
  community_all <- data.frame(time=rep(NA,(dim(all_data)[1])), species=NA, abundances=NA, replicates=NA)
  community_pred <- data.frame(time=rep(NA,(dim(all_data)[1])), species=NA, abundances=NA, replicates=NA)
  community_prey <- data.frame(time=rep(NA,(dim(all_data)[1])), species=NA, abundances=NA, replicates=NA)
  community_nprey <- data.frame(time=rep(NA,(dim(all_data)[1])), species=NA, abundances=NA, replicates=NA)
  # subset data
  data_all <- subset(all_data, select=c(time, patch, replicate, comm))
  data_pred <- subset(all_data, select=c(time, patch, replicate, comm_pred))
  data_prey <- subset(all_data, select=c(time, patch, replicate, comm_prey))
  data_nprey <- subset(all_data, select=c(time, patch, replicate, comm_nprey))
  # define community
  community_all$time <- community_pred$time <- community_prey$time <- community_nprey$time <- data_all$time
  community_all$species <- community_pred$species <- community_prey$species <- community_nprey$species <- data_all$patch
  community_all$replicates <- community_pred$replicates <- community_prey$replicates <- community_nprey$replicates <- data_all$replicate
  #
  community_all$abundances <- data_all[,4]
  community_pred$abundances <- data_pred[,4]
  community_prey$abundances <- data_prey[,4]
  community_nprey$abundances <- data_nprey[,4]
  # calculate synchrony
  # for all species
  Synchrony_all <- synchrony(community_all, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  # for predator
  Synchrony_pred <- synchrony(community_pred, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  # for prey
  Synchrony_prey <- synchrony(community_prey, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  # for non-prey
  Synchrony_nprey <- synchrony(community_nprey, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  # outputs
  out_all <- data.frame(replicates=r, environment=ePeriod, 
                        patch=NA, species=NA, scale="commspatial", trophic="all",
                        synchrony=Synchrony_all, variability=NA)
  out_pred <- data.frame(replicates=r, environment=ePeriod, 
                         patch=NA, species=NA, scale="commspatial", trophic="pred",
                         synchrony=Synchrony_pred, variability=NA)
  out_prey <- data.frame(replicates=r, environment=ePeriod, 
                         patch=NA, species=NA, scale="commspatial", trophic="prey",
                         synchrony=Synchrony_prey, variability=NA)
  out_nprey <- data.frame(replicates=r, environment=ePeriod, 
                          patch=NA, species=NA, scale="commspatial", trophic="nprey",
                          synchrony=Synchrony_nprey, variability=NA)
  # combine to return object
  synchrony_data <- rbind(synchrony_data, out_all, out_pred, out_prey, out_nprey)
  
  ###
  
  # variability
  
  # transform abundance into biomass
  X_biom <- X_saved
  #for(i in 1:dim(X_biom)[1]){
  #  X_biom[i,6:14] <- X_biom[i,6:14]*mass
  #}
  
  # create empty array of dimensions
  # row=spp, column=time(10), array=patch(5)
  empty_array_all <- array(NA, dim=c(nSp, length(unique(X_biom$time)), 5))
  empty_array_pred <- array(NA, dim=c(num_pred, length(unique(X_biom$time)), 5))
  empty_array_prey <- array(NA, dim=c(num_prey, length(unique(X_biom$time)), 5))
  empty_array_nprey <- array(NA, dim=c(num_nprey, length(unique(X_biom$time)), 5))
  
  # populate array with community data
  for(pp in 1:5){
    comm_all <- t(subset(X_biom, X_biom$patch==pp)[,6:14])
    #comm_pred <- t(subset(X_biom, X_biom$patch==pp)[,13:14])
    #comm_prey <- t(subset(X_biom, X_biom$patch==pp)[,9:12])
    #comm_nprey <- t(subset(X_biom, X_biom$patch==pp)[,6:8])
    comm_pred <- t(subset(X_biom, X_biom$patch==pp)[,6:12]) # change
    comm_prey <- t(subset(X_biom, X_biom$patch==pp)[,c(6:8,13,14)]) # change
    comm_nprey <- t(subset(X_biom, X_biom$patch==pp)[,9:14]) # change
    empty_array_all[,,pp] <- comm_all
    empty_array_pred[,,pp] <- comm_pred
    empty_array_prey[,,pp] <- comm_prey
    empty_array_nprey[,,pp] <- comm_nprey
  }
  
  # calculate synchrony and CV partitioning (Wang et al. 2019)
  part_calc_all <- as.data.frame(var.partition(empty_array_all))
  part_calc_pred <- as.data.frame(var.partition(empty_array_pred))
  part_calc_prey <- as.data.frame(var.partition(empty_array_prey))
  part_calc_nprey <- as.data.frame(var.partition(empty_array_nprey))
  
  # organize output
  out_CV_all <- data.frame(replicates=r, environment=ePeriod, 
                           patch=NA, species=NA, 
                           scale=c("Population","Community","Metapopulation","Metacommunity"), 
                           trophic="all", synchrony=NA, 
                           variability=c(part_calc_all[1,1], part_calc_all[2,1],
                                         part_calc_all[3,1], part_calc_all[4,1]))
  out_CV_pred <- data.frame(replicates=r, environment=ePeriod, 
                            patch=NA, species=NA, 
                            scale=c("Population","Community","Metapopulation","Metacommunity"), 
                            trophic="pred", synchrony=NA, 
                            variability=c(part_calc_pred[1,1], part_calc_pred[2,1],
                                          part_calc_pred[3,1], part_calc_pred[4,1]))
  out_CV_prey <- data.frame(replicates=r, environment=ePeriod, 
                            patch=NA, species=NA, 
                            scale=c("Population","Community","Metapopulation","Metacommunity"), 
                            trophic="prey", synchrony=NA, 
                            variability=c(part_calc_prey[1,1], part_calc_prey[2,1],
                                          part_calc_prey[3,1], part_calc_prey[4,1]))
  out_CV_nprey <- data.frame(replicates=r, environment=ePeriod, 
                             patch=NA, species=NA, 
                             scale=c("Population","Community","Metapopulation","Metacommunity"), 
                             trophic="nprey", synchrony=NA, 
                             variability=c(part_calc_nprey[1,1], part_calc_nprey[2,1],
                                           part_calc_nprey[3,1], part_calc_nprey[4,1]))
  
  # combine all results 
  synchrony_data <- rbind(synchrony_data, out_CV_all, out_CV_pred, out_CV_prey, out_CV_nprey)
  
  # return
  return(synchrony_data)
}

# Gross' synchrony & variability partitioning function for <1 % abundance
syn_CV_partitioning_abund <- function(X_saved){
  # data frame for all results
  synchrony_cv_data <- data.frame(replicates=as.numeric(), environment=as.numeric(), 
                                  patch=as.numeric(), species=as.character(),
                                  scale=as.character(), trophic=as.character(),
                                  synchrony=as.numeric(), CV=as.numeric())
  ### local species synchrony ###
  # note:
  # local species needs to be weighted
  # parameters
  environment <- ePeriod
  patches <- 5
  # data
  # all species
  all_data <- X_saved <- D2
  # data frame
  localspp_out <- data.frame(replicates=as.numeric(), environment=as.numeric(), 
                             patch=as.numeric(), species=as.character(),
                             scale=as.character(), trophic=as.character(),
                             synchrony=as.numeric(), variability=as.numeric())
  # loop
  for(p in 1:patches){
    # define patch
    pp <- as.numeric(p)
    # subset data
    d1 <- subset(all_data, patch==pp) 
    data <- subset(d1, select=c(time, replicate, (6:10)))
    # define community
    community_all <- data.frame(time=rep(NA,(dim(data)[1]*5)), species=NA, abundances=NA, replicates=NA)
    community_all$time <- data$time
    community_all$replicates <- data$replicate                   
    community_all$species <- rep(c("V4","V5","V7","V8","V9"), each=dim(data)[1])
    community_all$abundances <- c(data$V4,data$V5,data$V7,data$V8,data$V9)
    # calculate synchrony
    Synchrony_all <- synchrony(community_all, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    # outputs
    out_all <- data.frame(replicates=r, environment=ePeriod, 
                          patch=pp, species=NA, scale="localspp", trophic="all",
                          synchrony=Synchrony_all, variability=NA)
    # outputs combined
    localspp_out <- rbind(localspp_out, out_all)
  }
  # calculate weights for synchrony
  localspp_all_weights_S <- as.vector(unlist(tapply(localspp_out[localspp_out$trophic=="all",7], localspp_out[localspp_out$trophic=="all",3], function(x) (abs(x)/sum(abs(x))))))
  # weight data for synchrony
  localspp_all_weighted_S <- weighted_mean(localspp_out[localspp_out$trophic=="all",7], localspp_all_weights_S, na.rm=T)
  # final output for local species results
  localspp_data <- data.frame(replicates=r, environment=ePeriod, 
                              patch=NA, species=NA, scale="localspp", trophic="all",
                              synchrony=localspp_all_weighted_S,
                              variability=NA)
  # combine to return object
  synchrony_cv_data <- rbind(synchrony_cv_data, localspp_data)
  
  
  ### regional species synchrony ###
  # create metacommunity abundance data
  reg_list <- split(X_saved[,6:10], list(X_saved$time, X_saved$fluctuation, X_saved$replicate))
  reg_sum <- as.data.frame(do.call(rbind,(lapply(names(reg_list), function(x) apply(reg_list[[x]], 2, sum)))))
  treats <- data.frame(str_split_fixed(names(reg_list), "[.]", 3))
  reg_abund <- cbind(treats, reg_sum)
  colnames(reg_abund)[colnames(reg_abund)=="X1"] <- "time"
  colnames(reg_abund)[colnames(reg_abund)=="X2"] <- "environment"
  colnames(reg_abund)[colnames(reg_abund)=="X3"] <- "r"
  # data
  data_reg <- reg_abund
  # "loop"
  # define community
  community_all <- data.frame(time=rep(NA,dim(data_reg)[1]*5), species=NA, abundances=NA, replicates=NA)
  community_all$time <- as.numeric(data_reg$time)
  community_all$replicates <- data_reg$r
  community_all$species    <- rep(c("V4","V5","V7","V8","V9"), each=dim(data_reg)[1])
  community_all$abundances <- c(data_reg$V4,data_reg$V5,data_reg$V7,data_reg$V8,data_reg$V9)
  # calculate synchrony
  # for all species
  Synchrony_all <- synchrony(community_all, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  # outputs
  out_all <- data.frame(replicates=r, environment=ePeriod, 
                        patch=NA, species=NA, scale="regionalspp", trophic="all",
                        synchrony=Synchrony_all, variability=NA)
  # combine to return object
  synchrony_cv_data <- rbind(synchrony_cv_data, out_all)
  
  
  ### species spatial synchrony ###
  # note:
  # species spatial needs to be weighted
  # setting parameters
  species <- c("V4","V5","V7","V8","V9")
  # dataframe
  sppspatial_out <- data.frame(replicates=as.numeric(), environment=as.numeric(), 
                               patch=as.numeric(), species=as.character(),
                               scale=as.character(), trophic=as.character(),
                               synchrony=as.numeric(), variability=as.numeric())
  # loop
  for(spp in 6:10){
    # community dataframe
    community <- data.frame(time=rep(NA,(dim(data)[1]*5)), species=NA, abundances=NA, replicates=NA)
    # subset data
    data <- subset(all_data, select=c(time, patch, replicate, spp))
    # define community
    community$time       <- data$time
    community$species    <- data$patch
    community$abundances <- data[,4]
    community$replicates <- data$replicate                   
    # calculate synchrony
    Synchrony_spp <- synchrony(community, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    # output
    out_spp <- data.frame(replicates=r, environment=ePeriod, 
                          patch=NA, species=spp, scale="sppspatial", trophic=NA,
                          synchrony=Synchrony_spp, variability=NA)
    # combine outputs
    sppspatial_out <- rbind(sppspatial_out, out_spp)
  }
  # calculate weights for synchrony
  sppspatial_all_weights_S <- as.vector(unlist(tapply(sppspatial_out[,7], sppspatial_out[,4], function(x) (abs(x)/sum(abs(x))))))
  # weight data for synchrony
  sppspatial_all_weighted_S <- weighted_mean(sppspatial_out[,7], sppspatial_all_weights_S, na.rm=T)
  # final output for local species results
  sppspatial_data <- data.frame(replicates=r, environment=ePeriod, 
                                patch=NA, species=NA, scale="sppspatial", trophic="all",
                                synchrony=sppspatial_all_weighted_S,
                                variability=NA)
  # combine to return object
  synchrony_cv_data <- rbind(synchrony_cv_data, sppspatial_data)
  
  
  ### community spatial synchrony ###
  # data
  all_data$comm <- apply(X_saved[,6:10], 1, sum)
  # community data frame
  community_all <- data.frame(time=rep(NA,(dim(all_data)[1])), species=NA, abundances=NA, replicates=NA)
  # subset data
  data_all <- subset(all_data, select=c(time, patch, replicate, comm))
  # define community
  community_all$time <- data_all$time
  community_all$species <- data_all$patch
  community_all$replicates <- data_all$replicate
  community_all$abundances <- data_all[,4]
  # calculate synchrony
  # for all species
  Synchrony_all <- synchrony(community_all, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  # outputs
  out_all <- data.frame(replicates=r, environment=ePeriod, 
                        patch=NA, species=NA, scale="commspatial", trophic="all",
                        synchrony=Synchrony_all, variability=NA)
  # combine to return object
  synchrony_cv_data <- rbind(synchrony_cv_data, out_all)
  
  ###
  
  # variability
  
  # transform abundance into biomass
  X_biom <- X_saved
  #for(i in 1:dim(X_biom)[1]){
  #  X_biom[i,6:14] <- X_biom[i,6:14]*mass
  #}
  
  # create empty array of dimensions
  # row=spp, column=time(10), array=patch(5)
  empty_array_all <- array(NA, dim=c(5, length(unique(X_biom$time)), 5))
  
  # populate array with community data
  for(pp in 1:5){
    comm_all <- t(subset(X_biom, X_biom$patch==pp)[,6:10])
    empty_array_all[,,pp] <- comm_all
  }
  
  # calculate synchrony and CV partitioning (Wang et al. 2019)
  part_calc_all <- as.data.frame(var.partition(empty_array_all))
  
  # organize output
  out_CV_all <- data.frame(replicates=r, environment=ePeriod, 
                           patch=NA, species=NA, 
                           scale=c("Population","Community","Metapopulation","Metacommunity"), 
                           trophic="all", synchrony=NA, 
                           variability=c(part_calc_all[1,1], part_calc_all[2,1],
                                         part_calc_all[3,1], part_calc_all[4,1]))
  
  # combine all results 
  synchrony_cv_data <- rbind(synchrony_cv_data, out_CV_all)
  
  # return
  return(synchrony_cv_data)
}



