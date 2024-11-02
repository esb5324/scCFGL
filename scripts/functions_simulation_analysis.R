# This file has functions for calculating TPR, precision.

curve_base_func <- function(truth, est){
  TP <- 0; FP <- 0; TN <- 0; FN <- 0
  for (ix in length(truth)){
    truth[[ix]] <- truth[[ix]][upper.tri(truth[[ix]])]
    est[[ix]] <- est[[ix]][upper.tri(est[[ix]])]
    TP <- TP + sum((truth[[ix]] > 0 & est[[ix]] > 0) | (truth[[ix]] < 0 & est[[ix]] < 0))
    TN <- TN + sum(truth[[ix]] == 0 & est[[ix]] == 0)
    FP <- FP + sum(truth[[ix]] == 0 & est[[ix]] != 0)
    FN <- FN + sum((truth[[ix]] < 0 & est[[ix]] > 0) | (truth[[ix]] > 0 & est[[ix]] < 0) | (truth[[ix]] != 0 & est[[ix]] == 0))
  }
  tpr <- TP/(TP + FN)
  fpr <- FP/(FP+TN)
  prec <- TP/(TP + FP)
  return(list(TPR=tpr,FPR=fpr,precision=prec))
}

nets_res <- function(s){
  load(paste0("data/networks_", s, ".RData"))
  num_reps <- length(simu_data)
  scLink_Sigma <- lapply(1:num_reps,c)
  TPR <- lapply(1:4,c)
  FPR <- lapply(1:4,c)
  SSE <- lapply(1:4,c)
  precision <- lapply(1:4,c)
  edges <- lapply(1:4,c)
  num_cond <- length(simu_data[[1]])
  CFGL_list <- list(CFGL_out, CFGL_r_out, CFGL_rp_out)

  if(length(simu_data[[1]])==2) {n_vec <- c(nrow(simu_data[[1]][[1]]),nrow(simu_data[[1]][[2]]))}
  else  if(length(simu_data[[1]])==3) {n_vec <- c(nrow(simu_data[[1]][[1]]),nrow(simu_data[[1]][[2]]),nrow(simu_data[[1]][[3]]))}

  if (length(scLink_out[[1]])>1){
    # extract/rearrange scLink precision matrices
    for (rep in 1:num_reps){
      scLink_Sigma[[rep]] <- lapply(1:length(scLink_lda), c)
      for (sci in 1:length(scLink_lda)){
        scLink_Sigma[[rep]][[sci]] <- lapply(1:length(scLink_out[[rep]]), c) #repi
        for (scc in 1:num_cond){
          temp <- scLink_out[[rep]][[scc]]$summary[[sci]]$Sigma #repi
          scLink_Sigma[[rep]][[sci]][[scc]] <- temp
        }
      }
    }

    TPR[[1]] <- lapply(1:length(scLink_lda), function(n) 0)
    FPR[[1]] <- lapply(1:length(scLink_lda), function(n) 0)
    precision[[1]] <- lapply(1:length(scLink_lda), function(n) 0)
    SSE[[1]] <- lapply(1:length(scLink_lda), function(n) 0)
    edges[[1]] <- lapply(1:length(scLink_lda),c)

    for (scl in 1:length(scLink_lda)){
      if (num_cond==2){edges[[1]][[scl]] <- data.frame(C1=0,C2=0,C12=0,C1_only=0,C2_only=0)}
      else if (num_cond==3){edges[[1]][[scl]] <- data.frame(C1=0,C2=0,C3=0,C12_only=0,C13_only=0,C23_only=0,C123=0,C1_only=0,C2_only=0,C3_only=0)}

      for (rep in 1:num_reps){
        s_pcsm <- scLink_Sigma[[rep]][[scl]]
        pcsm_t_list <- tnet[[rep]]$pcsm

        pcsm_ <- scLink_Sigma[[rep]][[scl]]
        pcsm_t <- pcsm_t_list

        out <- curve_base_func(truth=pcsm_t, est=pcsm_)
        out$precision[is.na(out$precision)] <- 0
        TPR[[1]][[scl]] <- TPR[[1]][[scl]] + out$TPR
        FPR[[1]][[scl]] <- FPR[[1]][[scl]] + out$FPR
        precision[[1]][[scl]] <- precision[[1]][[scl]] + out$precision
      }
      TPR[[1]][[scl]] <- TPR[[1]][[scl]]/num_reps
      FPR[[1]][[scl]] <- FPR[[1]][[scl]]/num_reps
      precision[[1]][[scl]] <- precision[[1]][[scl]]/num_reps
      edges[[1]][[scl]] <- edges[[1]][[scl]]/num_reps
    }
  }

  # CFGL-based methods
  for (a in 1:3){
    if (length(CFGL_list[[a]][[1]])>1){
      # here we index so that each CFGL variation is indexed by lambda combinations
      TPR[[a+1]] <- lapply(1:length(lda[[2]][[a]]), c)
      FPR[[a+1]] <- lapply(1:length(lda[[2]][[a]]), c)
      precision[[a+1]] <- lapply(1:length(lda[[2]][[a]]), c)
      SSE[[a+1]] <- lapply(1:length(lda[[2]][[a]]), c)
      edges[[a+1]] <- lapply(1:length(lda[[2]][[a]]), c)
      count <- 0
      for (j in 1:length(lda[[2]][[a]])){
        TPR[[a+1]][[j]] <- lapply(1:length(lda[[1]][[a]]),function(n) 0)
        FPR[[a+1]][[j]] <- lapply(1:length(lda[[1]][[a]]),function(n) 0)
        precision[[a+1]][[j]] <- lapply(1:length(lda[[1]][[a]]), function(n) 0)
        SSE[[a+1]][[j]] <- lapply(1:length(lda[[1]][[a]]), function(n) 0)
        edges[[a+1]][[j]] <- lapply(1:length(lda[[1]][[a]]), function(n) 0)

        for (i in 1:length(lda[[1]][[a]])){
          count <- count + 1
          if (num_cond==2){edges[[a+1]][[j]][[i]] <- data.frame(C1_only=0,C2_only=0,C12=0,C1=0,C2=0)}
          else if (num_cond==3){edges[[a+1]][[j]][[i]] <- data.frame(C1_only = 0, C2_only = 0, C3_only = 0, C12_only = 0, C13_only = 0, C23_only = 0, C123 = 0, C1 = 0, C2 = 0, C3=0)}
          for (rep in 1:num_reps){
            pcsm_t_list <- tnet[[rep]]$pcsm
            C_pcsm <- CFGL_list[[a]][[rep]][[i]][[j]]$theta

            pcsm_t <- pcsm_t_list # true precision matrix

            out <- curve_base_func(truth=pcsm_t, est=C_pcsm)
            TPR[[a+1]][[j]][[i]] <- TPR[[a+1]][[j]][[i]] + out$TPR
            FPR[[a+1]][[j]][[i]] <- FPR[[a+1]][[j]][[i]] + out$FPR
            out$precision[is.na(out$precision)] <- 0
            precision[[a+1]][[j]][[i]] <- precision[[a+1]][[j]][[i]] + out$precision
          }
          TPR[[a+1]][[j]][[i]] <- TPR[[a+1]][[j]][[i]]/num_reps
          FPR[[a+1]][[j]][[i]] <- FPR[[a+1]][[j]][[i]]/num_reps
          precision[[a+1]][[j]][[i]] <- precision[[a+1]][[j]][[i]]/num_reps
          edges[[a+1]][[j]][[i]] <- edges[[a+1]][[j]][[i]]/num_reps
        }
      }
    }
  }

  p <- dim(pcsm_t_list[[1]])[1] # number of genes

  temp <- c()
  for (z in 1:length(lda[[2]][[1]])){

    # data frames to make ROC/PR curves: original, averaged TPR/precision
    if (length(CFGL_list[[1]][[1]])>1){
      Ctemp <-
        data.frame(
          TPR = unlist(TPR[[2]][[z]]), precision = unlist(precision[[2]][[z]]),
          FPR = unlist(FPR[[2]][[z]]),
          lda2 = lda[[2]][[1]][z],
          SSE = unlist(SSE[[2]][[z]]), edges=unlist(lapply(edges[[2]][[z]], function(a) sum(a$C1,a$C2))), edges1 = unlist(lapply(edges[[2]][[z]], function(a) a$C1)), edges2 = unlist(lapply(edges[[2]][[z]], function(a) a$C2)), lda1 = unlist(lda[[1]][[1]]), method = "CFGL")
    }

    if (length(CFGL_list[[2]][[1]])>1){
      Crtemp <-
        data.frame(
          TPR = unlist(TPR[[3]][[z]]), precision = unlist(precision[[3]][[z]]),
          FPR = unlist(FPR[[3]][[z]]),
          lda2 = lda[[2]][[2]][z],
          SSE = unlist(SSE[[3]][[z]]), edges=unlist(lapply(edges[[3]][[z]], function(a) sum(a$C1,a$C2))), edges1 = unlist(lapply(edges[[3]][[z]], function(a) a$C1)), edges2 = unlist(lapply(edges[[3]][[z]], function(a) a$C2)), lda1 = unlist(lda[[1]][[2]]), method = "CFGL-r")
    }

    if (length(CFGL_list[[3]][[1]])>1){
      Crptemp <-
        data.frame(
          TPR = unlist(TPR[[4]][[z]]), precision = unlist(precision[[4]][[z]]),
          FPR = unlist(FPR[[4]][[z]]),
          lda2 = lda[[2]][[3]][z],
          SSE = unlist(SSE[[4]][[z]]), edges=unlist(lapply(edges[[4]][[z]], function(a) sum(a$C1,a$C2))), edges1 = unlist(lapply(edges[[4]][[z]], function(a) a$C1)), edges2 = unlist(lapply(edges[[4]][[z]], function(a) a$C2)), lda1 = unlist(lda[[1]][[3]]), method = "CFGL-rp")
    }


    if (length(CFGL_list[[1]][[1]])>1){ temp <- rbind(temp, Ctemp)}
    if (length(CFGL_list[[2]][[1]])>1){ temp <- rbind(temp, Crtemp)}
    if (length(CFGL_list[[3]][[1]])>1){ temp <- rbind(temp, Crptemp)}
  }

  if (length(scLink_out[[1]])>1){
    stemp <-
      data.frame(
        TPR = unlist(TPR[[1]]), precision = unlist(precision[[1]]),
        FPR = unlist(FPR[[1]]),
        lda2 = 0,
        SSE = unlist(SSE[[1]]),
        edges=unlist(lapply(edges[[1]], function(a) sum(a$C1,a$C2))), edges1 = unlist(lapply(edges[[1]], function(a) a$C1)), edges2 = unlist(lapply(edges[[1]], function(a) a$C2)),
        lda1 = unlist(scLink_lda), method = "scLink")
  }

  if (length(scLink_out[[1]])>1){ temp <- rbind(temp, stemp)}

  ROC_data <- temp

  ##---------------------------------------------------------------------------
   save(list=c("ROC_data","TPR","FPR","precision","edges","lda"),file=paste0("data/simu_res_", s,".RData"))
}

