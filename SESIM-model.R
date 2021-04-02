##############################################
#                 SESIM model                #
#  Firkowski et al. - Ecological Monographs  #
##############################################

library(ggplot2)
library(reshape2)
library(vegan)
library(codyn)
library(stringr)
library(RColorBrewer)
library(ggbeeswarm)

# functions
source("./functions.R")

# model scenario
# choose "species" for multi-trophic food-web 
# or "competitive" for purely competitive interactions 
scenario <- "species" # "species" or "competitive"

# input parameters
# ignore warning
source("./input.R")

# simulation parameters
source("./simulation.R")

# output parameters
source("./output.R")

# model #
# experimental treatment count
treatment_id = 1
# loop count
loop_id = 1
# for loop of experimental design
for(ePeriod in eP){
  for(r in 1:reps){
    
    # sampling count          
    sampling = 1
    # set seed          
    set.seed(xseed[r])
              
    # initial abundances
    high <- base <- rep(100, nSp)
    X <- cbind(high, base, high, high, base)
    Xd <- t(X)
              
    # growth rate matrix
    C <- cbind(r.high, r.base, r.high, r.high, r.base)
              
    # species environmental optima
    Env_Opt <- matrix(runif(nSp, 0, 1), nSp, nEnv)
    
    # environmental states
    if(ePeriod==0) envt.v <- c(1, 1, 1, 1, 1)
    if(ePeriod!=0) envt.v <- c(1, 0, 1, 0, 1)
    env1 <- c(1, 0, 1, 0, 1)
    env2 <- c(0, 1, 0, 1, 0)
              
    # environmental fluctuation regime
    if(ePeriod==0) fluc.time <- 0
    if(ePeriod!=0) fluc.time <- seq(0, Tmax, by=ePeriod)
       
    # interaction strengths matrices
    # positive interactions
    for(j in 1:ncol(BB_pos)){
      for(i in 1:nrow(BB_pos)){
        set.seed(xseed[r])
        #pos[i,j] <- runif(1, min=min(0, BB_pos[i,j]), max=max(0, BB_pos[i,j]))
        pos[i,j] <- runif(1, min=min(BB_pos[i,j]/2, BB_pos[i,j]), max=max(BB_pos[i,j]/2, BB_pos[i,j]))
      }
    }
    pos[pos > 0] <- 0
    # negative interactions
    for(j in 1:ncol(BB_neg)){
      for(i in 1:nrow(BB_neg)){
        set.seed(xseed[r])
        #neg[i,j] <- runif(1, min=min(0, BB_neg[i,j]), max=max(0, BB_neg[i,j]))
        neg[i,j] <- runif(1, min=min(BB_neg[i,j]/2, BB_neg[i,j]), max=max(BB_neg[i,j]/2, BB_neg[i,j]))
       }
    }
    
    # for loop of species dynamics
    for(l in 1:Tmax){
      
      # environment
      Env <- matrix(rep(envt.v, each=nSp), nSp, numCom)
                
      # environmental match
      enviro <- t(apply(1-abs(Env_Opt[,1] - Env), 1, function(x) x/enveff))
      
      # negative interactions
      neg_interactions <- neg %*% (enviro*X)    
      # positive interactions
      pos_interactions <- (pos %*% X) * enviro  
      # total interactions
      interactions <- pos_interactions + neg_interactions
                
      # growth
      growth <- C * enviro
                
      # migration
      migrants <- apply(X, 2, function(x) x*(a*disp))
                
      # immigration
      immigrants <- matrix(NA, nSp, numCom)
      for(i in 1:nSp){
        immigrants[i,] <- dispersal_m%*%X[i,]*(a*disp[i])
      }
      
      # generalized lotka-volterra equation
      Xt <- X * exp(growth + interactions) + immigrants - migrants
                
      # extinctions
      Xt[Xt < a] <- 0
      X <- Xt
      Xd <- t(X)
                
      # sample
      if(any(sampleV==l) && l<=Tmax){
        X_save[,,sampling] <- Xd
        sampling <- sampling+1
      }
                
      # environmental fluctuation
      if(ePeriod!=0 && any(l==fluc.time)){
        if(envt.v[1]==1){
          envt.v <- env2
        }else{
          envt.v <- env1
        }
      }
                
      } # end of for loop for species dynamics
              
    # output data frames
    # experimental design
    experiment$environment[treatment_id] <- ePeriod
    # results
    results$treatment[ijloop()] <- treatment_id
    results$environment[ijloop()] <- ePeriod
    results$replicate[ijloop()] <- r
    # diversity
    results$alpha[ijloop()] <- vegan::diversity(log(Xd+10), index="shannon")
    results$gamma[ijloop()] <- vegan::diversity(apply(log(Xd+10), 2, sum), index="shannon")
    results$beta[ijloop()] <- betaeve(log(Xd+10))$"Weighted Beta"
    # richness
    results$comrich[ijloop()] <- specnumber(Xd)
    results$metarich[ijloop()] <- specnumber(apply(Xd, 2, sum))
    # evenness
    results$comeve[ijloop()] <- vegan::diversity(log(Xd+10), index="shannon")/log(specnumber(Xd)+10)
    results$metaeve[ijloop()] <- vegan::diversity(apply(log(Xd+10), 2, sum), index="shannon")/log(specnumber(apply(Xd, 2, sum))+10)
    # productivity
    productivity <- t(log(Xd+10)) * mass
    results$comprod[ijloop()] <- apply(productivity, 2, sum)
    results$metaprod[ijloop()] <- sum(productivity)
    # temporal sampling
    treatments <- data.frame(time = rep(c(sampleV), each=dim(X_save)[1]),
                             fluctuation = ePeriod,
                             replicate = r,
                             patch = c(1:5),
                             resource = c("high","base","high","high","base"))
    # abundances
    X_saved <- cbind(treatments, as.data.frame(apply(X_save[,,], 2, cbind)))
    # sampled abundances per replicate and environmental fluctuation period
    #assign(paste("rep_sampled", ePeriod, r, sep="_"), X_saved)
    # synchrony & variability results
    sync_cv_results <- rbind(sync_cv_results, syn_CV_partitioning(X_saved))
            
    # loop count
    loop_id <- loop_id+1
              
    } # end of for loop for replicates
  
  # count treatment          
  treatment_id  <- treatment_id+1
  
} # end of for loop for experimental design

# end
