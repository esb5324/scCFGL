## NOTE: to run these graphs, the lda2 values of all methods must be equal

l<- c(1) # the simulation settings you want to run graphs with. Note that the files must be labeled as "simu_res_1.RData" for S1, for example
library(pracma)
library(ggplot2)
source("scripts/functions_simulation_analysis.R")

auc <- data.frame()
for (s in l){
  nets_res(s)
  load(paste0("data/simu_res_", s, ".RData"))
  if (lda[[2]][[1]]!=lda[[2]][[2]] | lda[[2]][[1]] != lda[[2]][[3]] | lda[[2]][[2]] != lda[[2]][[3]]){
    print("For these simulation graphs, the lambda2 values of all methods (CFGL, CFGL-r, CFGL-rp) must be equal")
    q()
  }
  cbpalette <- c("#999999", "#56B4E9", "#009E73", "#CC79A7")

  for (i in 1:length(lda[[2]][[1]])) { # for each lambda2
    ROC_data_sub <- ROC_data[(ROC_data$lda2==lda[[2]][[1]][i] | ROC_data$method=="scLink"),]

    ROC_data_sc <- ROC_data[ROC_data$method == "scLink",]
    temp <- data.frame(setting=s, lda2=lda[[2]][[1]][i], method = "scLink", AUPRC=trapz(ROC_data_sc$TPR, ROC_data_sc$precision), AUROC=trapz(ROC_data_sc$FPR, ROC_data_sc$TPR))
    auc <- rbind(auc, temp); rm(temp)

    ROC_data_C <- ROC_data_sub[ROC_data_sub$method=="CFGL",]
    temp <- data.frame(setting=s, lda2=lda[[2]][[1]][i], method = "CFGL", AUPRC=trapz(ROC_data_C$TPR, ROC_data_C$precision), AUROC=trapz(ROC_data_C$FPR, ROC_data_C$TPR))
    auc <- rbind(auc, temp); rm(temp)

    ROC_data_Cr <- ROC_data_sub[ROC_data_sub$method=="CFGL-r",]
    temp <- data.frame(setting=s, lda2=lda[[2]][[1]][i], method = "CFGL-r", AUPRC=trapz(ROC_data_Cr$TPR, ROC_data_Cr$precision), AUROC=trapz(ROC_data_Cr$FPR, ROC_data_Cr$TPR))
    auc <- rbind(auc, temp); rm(temp)

    ROC_data_Crp <- ROC_data_sub[ROC_data_sub$method=="CFGL-rp",]
    temp <- data.frame(setting=s, lda2=lda[[2]][[1]][i], method = "CFGL-rp", AUPRC=trapz(ROC_data_Crp$TPR, ROC_data_Crp$precision), AUROC=trapz(ROC_data_Crp$FPR, ROC_data_Crp$TPR))
    auc <- rbind(auc, temp); rm(temp)

    # each graph is for a particular lambda2, and contains points for different lambda1's
    assign(paste0("S", s, "_PR_", lda[[2]][[1]][i]), ggplot(ROC_data_sub,aes(TPR, precision, color=method))
           + geom_line(linewidth=2,linetype=1)
           + coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + labs(x = NULL, y = NULL)
           + scale_color_manual(values=cbpalette)
           + theme(legend.position = "none")
    )

    assign(paste0("S", s, "_ROC_", lda[[2]][[1]][i]), ggplot(ROC_data_sub, aes(FPR, TPR, color = method))
           +geom_line(linewidth=2,linetype=1) + coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + labs(x = NULL, y = NULL)
           + scale_color_manual(values=cbpalette)
           + theme(legend.position = "none")
    )
  }
  rm(ROC_data)
}

save(list=c("S1_PR_0.01", "S1_ROC_0.01", # changed based on simulation settings run and lambda2's used
            "auc"),
     file=paste0("data/graphs.RData"))
