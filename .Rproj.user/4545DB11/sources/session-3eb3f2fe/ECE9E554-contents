library(huge)
source("scripts/JGNsc_simulation_functions.R")

num_reps=1 # how many "repetitions" of simulations to run
true_graph_m=1 # controls the number of edges in the ground truth network
rdseed <- 200 # controls the seed
nsample = 300 # number of cells
sim_list <- c(1) # What simulation settings to run? options: 1-5 for 2-condition settings, 1-2 for 3-condition settings
a_3=9
b_3=2.5
# a_3 and b_3 control the number of zeros in the resulting expression matrix. I used 9, 2.5 for 20% overall 0s in the expression; 6.1, 3.5 for 30%; 4.4, 4.8 for 40% (can adjust)
# Note that higher a_3, lower b_3 lead to fewer zeros

## 2 Conditions
str_list <- list("Identical S, Identical W", "Diff S, Identical W", "Diff S, Identical W", "Diff S, Identical W", "Change Weights")
diffblk_list <- list(NULL, list(7:8,7:8), list(5:8,5:8), list(3:8,3:8), NULL)
nivec.list.diff <- list(nivec= rep(25,8),nivec2 = rep(25,8))

## 3 Conditions
# NOTE: change weights only changes the second group of 4 modules (So if 8 modules, 5-8 have weights changed; if 9 modules, 6-9 have weights changed)
#str_list <- list(
#  c("Change Weights","3C: Change Weights/Diff S, Identical W"),
#  c("Change Weights", "3C: Change Weights/Diff S, Identical W") #"Diff S, Identical W")
#)
#diffblk_list <- list(
#  list(diffblk1=NULL,diffblk2=c(6:9)), # for C1-C2, change weights for modules 6-9; for C1-C3, change weights for modules 2-5, change modules 6-9 (module 1 is housekeeping)
#  list(diffblk1=NULL,diffblk2=c(2:9)) # for C1-C2, change weights for modules 6-9; for C1-C3, change modules 2-9 (module 1 is housekeeping)
 #)
# nivec.list.diff <- list(nivec= rep(25,9),nivec2 = rep(25,9), nivec3=rep(25,9))

# The following code simulates raw count data from 2 or 3 conditions, depending on the settings above



for (i in sim_list){
  str = str_list[[i]]
  diffblk = diffblk_list[[i]]
  tnet <- lapply(1:num_reps,c)
  simu_data <- lapply(1:num_reps,c)
  scrmats.true <- lapply(1:num_reps,c)

  for (rep in 1:num_reps){
    cormt <- list()
    pcsmt <- list()
    pij <- list()
    theta <- list()
    set.seed(rdseed*rep) # sets a different seed for each repetition
    tempout <-  generateSigmaList(
      nivec.list = nivec.list.diff, structure = str, diffblk = diffblk, blk_runif_threshold=edge_threshold, true_net="power_law", pos_def="diagplus1", add_m=true_graph_m)
    sigma.list.1 <- tempout$sigma.list

     # Create a list of count matrices
    countlist.1 <- getCountList(sigma.list = sigma.list.1, nvec = rep(nsample,2), a3=a_3, b3=b_3) # for higher, 47, a3=0.95 and b3=1.6, default is 3 1 i think
    theta <- lapply(countlist.1, function(c) c$thetaj)
    cormt <- lapply(countlist.1, function(c) c$sigma)
    pcsmt <- lapply(countlist.1, function(c) c$precision)
    adjt <- tempout$adj.list
    pij <- lapply(countlist.1, function(c) c$pij)

    for (index in 1:length(theta)){
      pcsmt[[index]][abs(pcsmt[[index]])<1e-8]=0
    }

    if (length(theta)==2){
    # true screening matrix is 1 if the absolute difference between the true precision matrices is < the threshold, 0 otherwise
      scrmats.true[[rep]] <- abs(pcsmt[[1]]-pcsmt[[2]])<0.01
    } else if (length(theta)==3){
        scrmats.true[[rep]] <- lapply(1:3,c)
         # true screening matrix is 1 if the absolute difference between the true precision matrices is < the threshold, 0 otherwise
        scrmats.true[[rep]][[1]] <- abs(pcsmt[[1]]-pcsmt[[2]])<0.01
        scrmats.true[[rep]][[2]] <- abs(pcsmt[[1]]-pcsmt[[3]])<0.01
        scrmats.true[[rep]][[3]] <- abs(pcsmt[[2]]-pcsmt[[3]])<0.01
        }

      simu_data[[rep]] <- lapply(1:length(countlist.1), function(i) huge.npn(countlist.1[[i]]$count))
      simu_data[[rep]] <- lapply(simu_data[[rep]], function(c) apply(c, 2, function(g) g - min(g))) # make min value 0 ## ?
      observed.list <- list(t(countlist.1[[1]]$count), t(countlist.1[[2]]$count))

    # find and remove columns that have one unique value, 0
    const_cols <- c()
    for (condit in 1:length(simu_data[[rep]])){
      if (sum(simu_data[[rep]][[condit]] > 20) > 0){
        print("warning: abnormally large values in simulated data, condition:")
        print(condit)
        print("max value:")
        print(max(simu_data[[rep]][[condit]]))
        q()
      }

      zero_cols <- sum(apply(simu_data[[rep]][[condit]], 2, function(c) length(unique(c))==1))
      if (zero_cols > 0){
        print(paste0("condition:", condit))
        print(paste0("number of constant columns: ", zero_cols))
        const_cols <- c(const_cols, which(apply(simu_data[[rep]][[condit]], 2, function(c) length(unique(c))==1)))
      }
    }

    if (length(const_cols) > 0){
      print("REMOVING CONSTANT COLUMNS IN ALL CONDITIONS")
      const_cols <- unique(const_cols)
      print("indices of const_cols:"); print(const_cols)
      simu_data[[rep]] <- lapply(simu_data[[rep]], function(m) m[,-(const_cols)])
      pcsmt <- lapply(pcsmt, function(m) m[-(const_cols),-(const_cols)])
      cormt <- lapply(cormt, function(m) m[-(const_cols),-(const_cols)])
      adjt <- lapply(adjt, function(m) m[-(const_cols),-(const_cols)])
      scrmats.true[[rep]] <- scrmats.true[[rep]][-(const_cols),-(const_cols)] # why did we have to change this???
      countlist.drop <- lapply(countlist.1, function(c) c$count[,-(const_cols)])
    } else{
      countlist.drop <- lapply(countlist.1, function(c) c$count)
    }

    tnet[[rep]]$pcsm <- pcsmt
    tnet[[rep]]$corm <- cormt
    tnet[[rep]]$matA <- adjt
    tnet[[rep]]$pij <- pij
    tnet[[rep]]$theta <- theta
  }


  save(file=paste("data/simu_data_", i, ".RData",sep=""),
       a_3, b_3, simu_data, tnet, scrmats.true, nivec.list.diff, str, diffblk, rdseed)
}
