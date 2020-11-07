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

# gross' synchrony & variability (CV) partitioning function
syn_CV_partitioning <- function(X_saved){
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
  patches <- numCom
  # data
  # all species
  all_data <- X_saved
  # number of species per trophic level
  num_pred <- 2
  num_prey <- 4
  num_nprey <- 3
  # data frame
  localspp_out <- data.frame(replicates=as.numeric(), environment=as.numeric(), 
                             patch=as.numeric(), species=as.character(),
                             scale=as.character(), trophic=as.character(),
                             synchrony=as.numeric(), CV=as.numeric())
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
    community_pred$species <- rep(c("V8","V9"), each=dim(data)[1])
    community_pred$abundances <- c(data$V8,data$V9)
    #
    community_prey$species <- rep(c("V4","V5","V6","V7"), each=dim(data)[1])
    community_prey$abundances <- c(data$V4,data$V5,data$V6,data$V7)
    #
    community_nprey$species <- rep(c("V1","V2","V3"), each=dim(data)[1])
    community_nprey$abundances <- c(data$V1,data$V2,data$V3)
    # calculate synchrony & CV
    # for all species
    Synchrony_all <- synchrony(community_all, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    CV_all <- sd(community_all$abundances, na.rm=T)/mean(community_all$abundances, na.rm=T)
    # for predators
    Synchrony_pred <- synchrony(community_pred, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    CV_pred <- sd(community_pred$abundances, na.rm=T)/mean(community_pred$abundances, na.rm=T)
    # for prey
    Synchrony_prey <- synchrony(community_prey, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    CV_prey <- sd(community_prey$abundances, na.rm=T)/mean(community_prey$abundances, na.rm=T)
    # for non-prey
    Synchrony_nprey <- synchrony(community_nprey, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    CV_nprey <- sd(community_nprey$abundances, na.rm=T)/mean(community_nprey$abundances, na.rm=T)
    # outputs
    out_all <- data.frame(replicates=r, environment=ePeriod, 
                          patch=pp, species=NA, scale="localspp", trophic="all",
                          synchrony=Synchrony_all, CV=CV_all)
    out_pred <- data.frame(replicates=r, environment=ePeriod, 
                           patch=pp, species=NA, scale="localspp", trophic="pred",
                           synchrony=Synchrony_pred, CV=CV_pred)
    out_prey <- data.frame(replicates=r, environment=ePeriod, 
                           patch=pp, species=NA, scale="localspp", trophic="prey",
                           synchrony=Synchrony_prey, CV=CV_prey)
    out_nprey <- data.frame(replicates=r, environment=ePeriod, 
                            patch=pp, species=NA, scale="localspp", trophic="nprey",
                            synchrony=Synchrony_nprey, CV=CV_nprey)
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
  # calculate weights for CV
  localspp_all_weights_CV <- as.vector(unlist(tapply(localspp_out[localspp_out$trophic=="all",8], localspp_out[localspp_out$trophic=="all",3], function(x) (abs(x)/sum(abs(x))))))
  localspp_pred_weights_CV <- as.vector(unlist(tapply(localspp_out[localspp_out$trophic=="pred",8], localspp_out[localspp_out$trophic=="pred",3], function(x) (abs(x)/sum(abs(x))))))
  localspp_prey_weights_CV <- as.vector(unlist(tapply(localspp_out[localspp_out$trophic=="prey",8], localspp_out[localspp_out$trophic=="prey",3], function(x) (abs(x)/sum(abs(x))))))
  localspp_nprey_weights_CV <- as.vector(unlist(tapply(localspp_out[localspp_out$trophic=="nprey",8], localspp_out[localspp_out$trophic=="nprey",3], function(x) (abs(x)/sum(abs(x))))))
  # weight data for CV
  localspp_all_weighted_CV <- weighted_mean(localspp_out[localspp_out$trophic=="all",8], localspp_all_weights_CV, na.rm=T)
  localspp_pred_weighted_CV <- weighted_mean(localspp_out[localspp_out$trophic=="pred",8], localspp_pred_weights_CV, na.rm=T)
  localspp_prey_weighted_CV <- weighted_mean(localspp_out[localspp_out$trophic=="prey",8], localspp_prey_weights_CV, na.rm=T)
  localspp_nprey_weighted_CV <- weighted_mean(localspp_out[localspp_out$trophic=="nprey",8], localspp_nprey_weights_CV, na.rm=T)
  # final output for local species results
  localspp_data <- data.frame(replicates=r, environment=ePeriod, 
                              patch=NA, species=NA, scale="localspp", trophic=c("all","pred","prey","nprey"),
                              synchrony=c(localspp_all_weighted_S, localspp_pred_weighted_S, localspp_prey_weighted_S, localspp_nprey_weighted_S), 
                              CV=c(localspp_all_weighted_CV, localspp_pred_weighted_CV, localspp_prey_weighted_CV, localspp_nprey_weighted_CV))
  # combine to return object
  synchrony_cv_data <- rbind(synchrony_cv_data, localspp_data)
  
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
  community_all$time <- community_pred$time <- community_prey$time <- community_nprey$time <- data_reg$time
  community_all$replicates <- community_pred$replicates <- community_prey$replicates <- community_nprey$replicates <- data_reg$r                
  #
  community_all$species    <- rep(c("V1","V2","V3","V4","V5","V6","V7","V8","V9"), each=dim(data_reg)[1])
  community_all$abundances <- c(data_reg$V1,data_reg$V2,data_reg$V3,data_reg$V4,data_reg$V5,data_reg$V6,data_reg$V7,data_reg$V8,data_reg$V9)
  #
  community_pred$species    <- rep(c("V8","V9"), each=dim(data_reg)[1])
  community_pred$abundances <- c(data_reg$V8,data_reg$V9)
  #
  community_prey$species    <- rep(c("V4","V5","V6","V7"), each=dim(data_reg)[1])
  community_prey$abundances <- c(data_reg$V4,data_reg$V5,data_reg$V6,data_reg$V7)
  #
  community_nprey$species    <- rep(c("V1","V2","V3"), each=dim(data_reg)[1])
  community_nprey$abundances <- c(data_reg$V1,data_reg$V2,data_reg$V3)
  # calculate synchrony & CV
  # for all species
  Synchrony_all <- synchrony(community_all, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  CV_all <- sd(community_all$abundances, na.rm=T)/mean(community_all$abundances, na.rm=T)
  # for predators
  Synchrony_pred <- synchrony(community_pred, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  CV_pred <- sd(community_pred$abundances, na.rm=T)/mean(community_pred$abundances, na.rm=T)
  # for prey
  Synchrony_prey <- synchrony(community_prey, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  CV_prey <- sd(community_prey$abundances, na.rm=T)/mean(community_prey$abundances, na.rm=T)
  # for non-prey
  Synchrony_nprey <- synchrony(community_nprey, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  CV_nprey <- sd(community_nprey$abundances, na.rm=T)/mean(community_nprey$abundances, na.rm=T)
  # outputs
  out_all <- data.frame(replicates=r, environment=ePeriod, 
                        patch=NA, species=NA, scale="regionalspp", trophic="all",
                        synchrony=Synchrony_all, CV=CV_all)
  out_pred <- data.frame(replicates=r, environment=ePeriod, 
                         patch=NA, species=NA, scale="regionalspp", trophic="pred",
                         synchrony=Synchrony_pred, CV=CV_pred)
  out_prey <- data.frame(replicates=r, environment=ePeriod, 
                         patch=NA, species=NA, scale="regionalspp", trophic="prey",
                         synchrony=Synchrony_prey, CV=CV_prey)
  out_nprey <- data.frame(replicates=r, environment=ePeriod, 
                          patch=NA, species=NA, scale="regionalspp", trophic="nprey",
                          synchrony=Synchrony_nprey, CV=CV_nprey)
  # combine to return object
  synchrony_cv_data <- rbind(synchrony_cv_data, out_all, out_pred, out_prey, out_nprey)
  
  ### species spatial synchrony ###
  # note:
  # species spatial needs to be weighted
  # setting parameters
  species <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9")
  # dataframe
  sppspatial_out <- data.frame(replicates=as.numeric(), environment=as.numeric(), 
                               patch=as.numeric(), species=as.character(),
                               scale=as.character(), trophic=as.character(),
                               synchrony=as.numeric(), CV=as.numeric())
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
    # calculate synchrony & CV
    Synchrony_spp <- synchrony(community, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    CV_spp <- sd(community$abundances, na.rm=T)/mean(community$abundances, na.rm=T)
    # output
    out_spp <- data.frame(replicates=r, environment=ePeriod, 
                          patch=NA, species=spp, scale="sppspatial", trophic=NA,
                          synchrony=Synchrony_spp, CV=CV_spp)
    # combine outputs
    sppspatial_out <- rbind(sppspatial_out, out_spp)
  }
  # calculate weights for synchrony
  sppspatial_all_weights_S <- as.vector(unlist(tapply(sppspatial_out[,7], sppspatial_out[,4], function(x) (abs(x)/sum(abs(x))))))
  sppspatial_pred_weights_S <- as.vector(unlist(tapply(sppspatial_out[c(8,9),7], sppspatial_out[c(8,9),4], function(x) (abs(x)/sum(abs(x))))))
  sppspatial_prey_weights_S <- as.vector(unlist(tapply(sppspatial_out[c(4:7),7], sppspatial_out[c(4:7),4], function(x) (abs(x)/sum(abs(x))))))
  sppspatial_nprey_weights_S <- as.vector(unlist(tapply(sppspatial_out[c(1:3),7], sppspatial_out[c(1:3),4], function(x) (abs(x)/sum(abs(x))))))
  # weight data for synchrony
  sppspatial_all_weighted_S <- weighted_mean(sppspatial_out[,7], sppspatial_all_weights_S, na.rm=T)
  sppspatial_pred_weighted_S <- weighted_mean(sppspatial_out[c(8,9),7], sppspatial_pred_weights_S, na.rm=T)
  sppspatial_prey_weighted_S <- weighted_mean(sppspatial_out[c(4:7),7], sppspatial_prey_weights_S, na.rm=T)
  sppspatial_nprey_weighted_S <- weighted_mean(sppspatial_out[c(1:3),7], sppspatial_nprey_weights_S, na.rm=T)
  # calculate weights for CV
  sppspatial_all_weights_CV <- as.vector(unlist(tapply(sppspatial_out[,8], sppspatial_out[,4], function(x) (abs(x)/sum(abs(x))))))
  sppspatial_pred_weights_CV <- as.vector(unlist(tapply(sppspatial_out[c(8,9),8], sppspatial_out[c(8,9),4], function(x) (abs(x)/sum(abs(x))))))
  sppspatial_prey_weights_CV <- as.vector(unlist(tapply(sppspatial_out[c(4:7),8], sppspatial_out[c(4:7),4], function(x) (abs(x)/sum(abs(x))))))
  sppspatial_nprey_weights_CV <- as.vector(unlist(tapply(sppspatial_out[c(1:3),8], sppspatial_out[c(1:3),4], function(x) (abs(x)/sum(abs(x))))))
  # weight data for CV
  sppspatial_all_weighted_CV <- weighted_mean(sppspatial_out[,8], sppspatial_all_weights_CV, na.rm=T)
  sppspatial_pred_weighted_CV <- weighted_mean(sppspatial_out[c(8,9),8], sppspatial_pred_weights_CV, na.rm=T)
  sppspatial_prey_weighted_CV <- weighted_mean(sppspatial_out[c(4:7),8], sppspatial_prey_weights_CV, na.rm=T)
  sppspatial_nprey_weighted_CV <- weighted_mean(sppspatial_out[c(1:3),8], sppspatial_nprey_weights_CV, na.rm=T)
  # final output for local species results
  sppspatial_data <- data.frame(replicates=r, environment=ePeriod, 
                                patch=NA, species=NA, scale="sppspatial", trophic=c("all","pred","prey","nprey"),
                                synchrony=c(sppspatial_all_weighted_S, sppspatial_pred_weighted_S, sppspatial_prey_weighted_S, sppspatial_nprey_weighted_S), 
                                CV=c(sppspatial_all_weighted_CV, sppspatial_pred_weighted_CV, sppspatial_prey_weighted_CV, sppspatial_nprey_weighted_CV))
  # combine to return object
  synchrony_cv_data <- rbind(synchrony_cv_data, sppspatial_data)
  
  ### community spatial synchrony ###
  # data
  all_data$comm <- apply(X_saved[,6:14], 1, sum)
  all_data$comm_pred <- apply(X_saved[,13:14], 1, sum)
  all_data$comm_prey <- apply(X_saved[,9:12], 1, sum)
  all_data$comm_nprey <- apply(X_saved[,6:8], 1, sum)
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
  CV_all <- NA
  # for predator
  Synchrony_pred <- synchrony(community_pred, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  CV_pred <- NA
  # for prey
  Synchrony_prey <- synchrony(community_prey, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  CV_prey <- NA
  # for non-prey
  Synchrony_nprey <- synchrony(community_nprey, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  CV_nprey <- NA
  # outputs
  out_all <- data.frame(replicates=r, environment=ePeriod, 
                        patch=NA, species=NA, scale="commspatial", trophic="all",
                        synchrony=Synchrony_all, CV=CV_all)
  out_pred <- data.frame(replicates=r, environment=ePeriod, 
                         patch=NA, species=NA, scale="commspatial", trophic="pred",
                         synchrony=Synchrony_pred, CV=CV_pred)
  out_prey <- data.frame(replicates=r, environment=ePeriod, 
                         patch=NA, species=NA, scale="commspatial", trophic="prey",
                         synchrony=Synchrony_prey, CV=CV_prey)
  out_nprey <- data.frame(replicates=r, environment=ePeriod, 
                          patch=NA, species=NA, scale="commspatial", trophic="nprey",
                          synchrony=Synchrony_nprey, CV=CV_nprey)
  # combine to return object
  synchrony_cv_data <- rbind(synchrony_cv_data, out_all, out_pred, out_prey, out_nprey)
  
  ### population CV ###
  # note: population CV needs to be weighted
  # data frame
  popCV_out <- data.frame(replicates=as.numeric(), environment=as.numeric(), 
                          patch=as.numeric(), species=as.character(),
                          scale=as.character(), trophic=as.character(),
                          synchrony=as.numeric(), CV=as.numeric())
  # loop
  for(spp in 6:14){
    for(p in 1:5){
      # define patch
      pp <- as.numeric(p)  
      # community data frame
      community <- data.frame(time=rep(NA,(dim(data)[1])), species=NA, abundances=NA, replicates=NA)
      # subset data
      d1 <- subset(all_data, patch==pp)
      data <- subset(d1, select=c(time, patch, replicate, spp))
      # define community
      community$time <- data$time
      community$species <- paste(data$patch, spp, sep="")
      community$abundances <- data[,4]
      community$replicates <- data$replicate                   
      # calculate CV
      Synchrony_spp_p <- NA
      CV_spp_p <- sd(community$abundances, na.rm=T)/mean(community$abundances, na.rm=T)
      # output
      out_spp_p <- data.frame(replicates=r, environment=ePeriod, 
                              patch=pp, species=spp, scale="popCV", trophic=NA,
                              synchrony=Synchrony_spp_p, CV=CV_spp_p)
      # combine outputs
      popCV_out <- rbind(popCV_out, out_spp_p)
    }
  }
  # calculate weights for synchrony
  popCV_all_weights_S <- as.vector(unlist(tapply(popCV_out[,7], popCV_out[,4], function(x) (abs(x)/sum(abs(x))))))
  popCV_pred_weights_S <- as.vector(unlist(tapply(popCV_out[popCV_out$species==13 | popCV_out$species==14,7], popCV_out[popCV_out$species==13 | popCV_out$species==14,4], function(x) (abs(x)/sum(abs(x))))))
  popCV_prey_weights_S <- as.vector(unlist(tapply(popCV_out[popCV_out$species==9 | popCV_out$species==10 | popCV_out$species==11 | popCV_out$species==12,7], 
                                                  popCV_out[popCV_out$species==9 | popCV_out$species==10 | popCV_out$species==11 | popCV_out$species==12,4], function(x) (abs(x)/sum(abs(x))))))
  popCV_nprey_weights_S <- as.vector(unlist(tapply(popCV_out[popCV_out$species==6 | popCV_out$species==7 | popCV_out$species==8,7], 
                                                   popCV_out[popCV_out$species==6 | popCV_out$species==7 | popCV_out$species==8,4], function(x) (abs(x)/sum(abs(x))))))
  # weight data for synchrony
  popCV_all_weighted_S <- weighted_mean(popCV_out[,7], popCV_all_weights_S, na.rm=T)
  popCV_pred_weighted_S <- weighted_mean(popCV_out[popCV_out$species==13 | popCV_out$species==14,7], popCV_pred_weights_S, na.rm=T)
  popCV_prey_weighted_S <- weighted_mean(popCV_out[popCV_out$species==9 | popCV_out$species==10 | popCV_out$species==11 | popCV_out$species==12,7], popCV_prey_weights_S, na.rm=T)
  popCV_nprey_weighted_S <- weighted_mean(popCV_out[popCV_out$species==6 | popCV_out$species==7 | popCV_out$species==8,7], popCV_nprey_weights_S, na.rm=T)
  # calculate weights for CV
  popCV_all_weights_CV <- as.vector(unlist(tapply(popCV_out[,8], popCV_out[,4], function(x) (abs(x)/sum(abs(x))))))
  popCV_pred_weights_CV <- as.vector(unlist(tapply(popCV_out[popCV_out$species==13 | popCV_out$species==14,8], popCV_out[popCV_out$species==13 | popCV_out$species==14,4], function(x) (abs(x)/sum(abs(x))))))
  popCV_prey_weights_CV <- as.vector(unlist(tapply(popCV_out[popCV_out$species==9 | popCV_out$species==10 | popCV_out$species==11 | popCV_out$species==12,8], 
                                                   popCV_out[popCV_out$species==9 | popCV_out$species==10 | popCV_out$species==11 | popCV_out$species==12,4], function(x) (abs(x)/sum(abs(x))))))
  popCV_nprey_weights_CV <- as.vector(unlist(tapply(popCV_out[popCV_out$species==6 | popCV_out$species==7 | popCV_out$species==8,8], 
                                                    popCV_out[popCV_out$species==6 | popCV_out$species==7 | popCV_out$species==8,4], function(x) (abs(x)/sum(abs(x))))))
  # weight data for CV
  popCV_all_weighted_CV <- weighted_mean(popCV_out[,8], popCV_all_weights_CV, na.rm=T)
  popCV_pred_weighted_CV <- weighted_mean(popCV_out[popCV_out$species==13 | popCV_out$species==14,8], popCV_pred_weights_CV, na.rm=T)
  popCV_prey_weighted_CV <- weighted_mean(popCV_out[popCV_out$species==9 | popCV_out$species==10 | popCV_out$species==11 | popCV_out$species==12,8], popCV_prey_weights_CV, na.rm=T)
  popCV_nprey_weighted_CV <- weighted_mean(popCV_out[popCV_out$species==6 | popCV_out$species==7 | popCV_out$species==8,8], popCV_nprey_weights_CV, na.rm=T)
  # output
  popCV_data <- data.frame(replicates=r, environment=ePeriod, 
                           patch=NA, species=NA, scale="popCV", trophic=c("all","pred","prey","nprey"),
                           synchrony=c(popCV_all_weighted_S, popCV_pred_weighted_S, popCV_prey_weighted_S, popCV_nprey_weighted_S), 
                           CV=c(popCV_all_weighted_CV, popCV_pred_weighted_CV, popCV_prey_weighted_CV, popCV_nprey_weighted_CV))
  # combine to return object
  synchrony_cv_data <- rbind(synchrony_cv_data, popCV_data)
  # return
  return(synchrony_cv_data)
}

