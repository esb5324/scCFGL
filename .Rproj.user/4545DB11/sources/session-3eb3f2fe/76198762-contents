#' single-cell Condition-adaptive fused graphical lasso.
#' The function jointly construct gene co-expression network for multiple class using single-cell Condition-adaptive Fused Graphical Lasso. Pairwise screening matrices are required to adjust between-condition lasso penalty.
#' @param Y A list of expression data which are n*p matrices. all matrices should have the same p. Values should be 0 or above.
#' @import matrixcalc
#' @import parallel
#' @import glmnet
#' @param lambda1 The tuning parameter for the graphical lasso penalty.
#' @param lambda2 The tuning parameter for the between condition group lasso penalty.
#' @param btc.screening A list of screening matrices (p*p) for between condition penalty. Can be obtained using the function \code{get_scr_mat}.
#' @param rho Step size parameter for ADMM algorithm. Large values decrease the step size.
#' @param rho.increment Adjustment for rho. In each ADMM iteration, rho will be updated as rho=rho*rho.increment.
#' @param maxiter The maximum number of ADMM iterations.
#' @param tol The criterion for ADMM convergence.
#' @param truncate All value in the estimated inverse convenience below this number will be set to 0.
#' @param loglik.trace Store trace of the likelihood of estimation in each iteration.
#' @param nthre An integer specifying a threshold on the number of complete observations.
#' @param dthre A number specifying the threshold on dropout probabilities.
#' @param ncores: set to 1 if using Windows, of >1 can allow for parallelization when calculating the robust correlation
#' @param change_penalty: if TRUE, alters the first penalty term to include the multiplicative term
#' @param robust_S: if TRUE, uses the robust S instead of the typical estimated covariance
#' @details To run CFGL, set change_penalty=FALSE and robust_S=FALSE. To run CFGL-r, set change_penalty=FALSE and robust_S=TRUE. To run CFGL-rp, set change_penalty=TRUE and robust_S=TRUE. Do not set change_penalty=TRUE and robust_S=FALSE. Note that this function has adaptations from scLink and CFGL.
#' @return \code{sccfgl} produces a list that contains estimated precision matrices and other components.
#' \itemize{
#'  \item{$theta} {The estimation of inverse matrices}
#'  \item{$cor} {The correlation used; Pearson correlation for CFGL and robust correlation otherwise}
#'  \item{$iters} {The number of ADMM iterations}
#'  \item{$loglik.trace} {Trace of log-likelihood}
#' }
#' @export
sccfgl <- function(Y, lambda1, lambda2, btc.screening=NULL,
                    rho=1, rho.increment=1, maxiter=2000, tol=1e-4, truncate=1e-05,
                   loglik.trace=FALSE, nthre=10,dthre=.5, ncores=4, change_penalty=FALSE, robust_S=TRUE){
  K <- length(Y) # number of conditions
  p <- c() # number of genes
  n <- c() # number of cells
  S <- list() # est cov matrix list
  corr <- list()
  weights_p <- list() # these are the weights for the adjusted penalty term, if we adjust it (CFGL-rp)
  for (k in 1:K) {
    p[k] <- dim(Y[[k]])[2]
    n[k] <- dim(Y[[k]])[1]
    x.q = apply(Y[[k]], 2, sd)
    cor_pearson = cor(Y[[k]])
    cor = est_pearson_scimpute.c(Y[[k]], ncores = ncores, cor.p = cor_pearson,
                                           thre_no = nthre, dthre = dthre)
    weights_p[[k]] = 1 - abs(cor) # for adjusted penalty term if we adjust it
    if  (!robust_S){ # if CFGL
      corr[[k]] <- cor_pearson
      S[[k]] <- cov(Y[[k]])* ((n[k] - 1)/n[k])
    } else {
      corr[[k]] <- cor
      cov = diag(x.q) %*% cor %*% diag(x.q)
      S[[k]] <-  (easy.psd(cov,method="perturb"))
    }
  }
  if (!(K==2||K==3)) stop("K must be 2 or 3 in this version")
  if (p[1]!=p[2]) stop("p must be the same for each k")
  p <- p[1]
  weights <- rep(1, K)
  lam1.m <- get_lam_mat(lambda1, p, penalize.diagonal = FALSE)
  lam2.m0 <- get_lam_mat(lambda2, p, penalize.diagonal = FALSE)
  if ((K==2)&&(!is.null(btc.screening))) lam2.m = lam2.m0 * btc.screening else lam2.m = lam2.m0
  if (K==3){
    if (!is.null(btc.screening)){
      lam2.m <- list()
      lam2.m[[1]] <- btc.screening[[1]] * lam2.m0 #1-2
      lam2.m[[2]] <- btc.screening[[2]] * lam2.m0 #1-3
      lam2.m[[3]] <- btc.screening[[3]] * lam2.m0 #2-3
    }
    else{
      lam2.m <- list()
      lam2.m[[1]] <- lam2.m0 #1-2
      lam2.m[[2]] <- lam2.m0 #1-3
      lam2.m[[3]] <- lam2.m0 #2-3
    }
  }
  if (change_penalty){
    lam1.weights <- weights_p
  }
  else{
    lam1.weights <- list()
    for (k in 1:K){
      lam1.weights[[k]] <- matrix(1,nrow=nrow(lam1.m),ncol=nrow(lam1.m))
      #
    }
  }
  # admm start
  out.admm = admm.iter(S = S,n = n,lam1.m = lam1.m,lam2.m = lam2.m, weights = weights, rho = rho,rho.increment = rho.increment,maxiter = maxiter,tol = tol,loglik.trace = loglik.trace,lam1.weights)
  diff_theta_z = 0
  for (k in 1:K) {
    diff_theta_z = diff_theta_z + sum(abs(out.admm$theta[[k]] - out.admm$Z[[k]]))
  }
  # round down
  theta = list()
  for (k in 1:K) {
    rounddown = abs(out.admm$Z[[k]]) < truncate
    diag(rounddown) = FALSE
    theta[[k]] = out.admm$Z[[k]] * (1 - rounddown)
  }
  out.scCFGL= list(theta = theta,iters = out.admm$iter)
  if (loglik.trace) {
    out.V$loglik.trace = out.admm$loglik.trace
  }
  out.scCFGL$cor=corr
  return(out.scCFGL)
}

## -----------------------------------------------------------------------------

# This function runs the methods on the simulated data. See scripts folder, file estimate_networks.R for an example of how to run it
run_nets <- function(scLink_lda=NA, CFGL_lda1=NA, CFGLr_lda1=NA, CFGLrp_lda1=NA, CFGL_lda2, CFGLr_lda2, CFGLrp_lda2,
                     n_cores, SM, data_mat, txt=""){
    if(length(scLink_lda)!=length(data_mat[[1]])){
    print("The length of the lambda1 list for scLink (scLink_Lda) must equal the number of conditions")
    q()
  }
  if (is.na(SM)){
  if (length(data_mat[[1]])==3){
    scr_mat <- lapply(1:length(data_mat),c)
    for (i in 1:length(data_mat)){
      temp1_2<-get_scr_mat(expr1=data_mat[[i]][[1]],expr2=data_mat[[i]][[2]])$scr.mat
      temp1_3<-get_scr_mat(expr1=data_mat[[i]][[1]],expr2=data_mat[[i]][[3]])$scr.mat
      temp2_3<-get_scr_mat(expr1=data_mat[[i]][[2]],expr2=data_mat[[i]][[3]])$scr.mat
      scr_mat[[i]] <- list(temp1_2, temp1_3, temp2_3)
    }
  } else if (length(data_mat[[1]])==2){
    ## 2 Conditions
    scr_mat <- mclapply(data_mat, function(x) get_scr_mat(expr1=x[[1]], expr2 = x[[2]])$scr.mat,mc.cores=20)
  }
  }
  if (!is.na(scLink_lda[[1]][1])){
    scLink_out <- mclapply(data_mat, function(x) mclapply(x, function(y) sclink_net(y, ncores = n_cores, lda = scLink_lda, nthre=10, dthre=0.5),mc.cores=n_cores,mc.preschedule=FALSE),mc.cores=n_cores,mc.preschedule=FALSE)
    for (a in 1:length(scLink_out)){
          for (b in 1:length(scLink_out[[a]])){
            for (c in 1:length(scLink_out[[a]][[b]]$summary)){
              temp <- scLink_out[[a]][[b]]$summary[[c]]$adj
              temp[lower.tri(temp)]<-t(temp)[lower.tri(temp)]
              temp[abs(temp)<1e-5] <- 0
            }
          }
        }
  } else{
    scLink_out <- lapply(data_mat, function(a) NA)
  }
  if (!is.na(CFGL_lda1[1])){
     CFGL_out <- lapply(1:length(data_mat), function(k) { # NOTE: EDITED ncores to n_cores for mc.cores!!!!!
       mclapply(CFGL_lda1, function(x) lapply(CFGL_lda2, function(y) sccfgl(data_mat[[k]], btc.screening = scr_mat[[k]], lambda1=x, lambda2=y, nthre=10, dthre = 0.5, change_penalty = FALSE, robust_S=FALSE, ncores=n_cores)), mc.cores = n_cores, mc.preschedule = FALSE)
     })
  } else{
    CFGL_out <- lapply(data_mat, function(a) NA)
  }
  if (!is.na(CFGLr_lda1[1])){
    CFGL_r_out <- lapply(1:length(data_mat), function(k) { # NOTE: EDITED ncores to n_cores for mc.cores!!!!!
      mclapply(CFGLr_lda1, function(x) mclapply(CFGLr_lda2, function(y) sccfgl(data_mat[[k]], btc.screening = scr_mat[[k]], lambda1=x, lambda2=y, nthre=10, dthre = 0.5, change_penalty = FALSE, robust_S=TRUE, ncores=n_cores),mc.cores=n_cores,mc.preschedule=FALSE), mc.cores = n_cores, mc.preschedule = FALSE)
    })
  } else{
    CFGL_r_out <- lapply(data_mat, function(a) NA)
  }
  if (!is.na(CFGLrp_lda1[1])){
    CFGL_rp_out <- lapply(1:length(data_mat), function(k) {
      mclapply(CFGLrp_lda1, function(x) mclapply(CFGLrp_lda2, function(y) sccfgl(data_mat[[k]], btc.screening = scr_mat[[k]], lambda1=x, lambda2=y, nthre=10, dthre = 0.5, change_penalty = TRUE, robust_S=TRUE, ncores=n_cores),mc.cores=n_cores,mc.preschedule=FALSE), mc.cores = n_cores, mc.preschedule = FALSE)
    })
  } else{
    CFGL_rp_out <- lapply(data_mat, function(a) NA)
  }
  lda <- list()
  lda[[1]] <- list(CFGL=CFGL_lda1, CFGLr=CFGLr_lda1, CFGLrp=CFGLrp_lda1)
  lda[[2]] <- list(CFGL=CFGL_lda2, CFGLr=CFGLr_lda2, CFGLrp=CFGLrp_lda2)
  ## -----------------------------------------------------------------------------
  save(list=c(
"scLink_out", "CFGL_out","CFGL_r_out","CFGL_rp_out", "scLink_lda", "lda"),file=paste0("data/networks_", txt, ".RData"))
}
