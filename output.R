##############################################
#              output data frames            #
#  Firkowski et al. - Ecological Monographs  #
##############################################

# experimental design
experiment <- data.frame(treatment=rep(1:(length(eP))), environment=NA)
# sampled abundances
X_save <- array(data=NA, dim=c(numCom, nSp, length(sampleV)))
# diversity results
results <- data.frame(treatment=rep(NA, dim(experiment)[1]*numCom*reps),
                      environment=NA, replicate=NA, patch=rep(1:numCom),
                      alpha=NA, gamma=NA, 
                      comrich=NA, metarich=NA, 
                      comeve=NA, metaeve=NA, 
                      comprod=NA, metaprod=NA)
# synchrony & variability results
sync_cv_results <- data.frame(replicates=as.numeric(), environment=as.numeric(), 
                              patch=as.numeric(), species=as.character(),
                              scale=as.character(), trophic=as.character(),
                              synchrony=as.numeric(), CV=as.numeric())
