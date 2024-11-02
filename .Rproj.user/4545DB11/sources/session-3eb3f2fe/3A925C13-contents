# functions from CFGL which may be altered, see: https://github.com/Yafei611/CFGL for original functions and information

get_lam_mat <- function (lambda, p, penalize.diagonal) {
  if (is.matrix(lambda)) {
    if (sum(lambda != t(lambda)) > 0) stop("error: penalty matrix is not symmetric")
    if (sum(abs(dim(lambda) - p)) != 0) stop("error: penalty matrix has wrong dimension")
  }
  if (length(lambda) == 1) lambda = matrix(lambda, p, p)
  if (!penalize.diagonal) diag(lambda) = 0
  return(lambda)
} # convert lambda to matrices

soft_thresholding <- function (a, lam.m) {
  out <- sign(a) * pmax(0, abs(a) - lam.m)
  return(out) # used to solve for updated Z in functions below
}

cal_z_L1_t2 <- function (A, rho, lam1.m, lam2.m, lam1.weights){
  S1 <- abs(A[[1]] - A[[2]]) <= 2 * lam2.m/rho
  Z1_1 <- (A[[1]] + A[[2]])/2
  Z2_1 <- Z1_1
  S2 <- (A[[1]] > A[[2]] + 2 * lam2.m/rho) # error
  Z1_2 <- A[[1]] - lam2.m/rho
  Z2_2 <- A[[2]] + lam2.m/rho
  S3 <- (A[[2]] > A[[1]] + 2 * lam2.m/rho)
  Z1_3 <- A[[1]] + lam2.m/rho
  Z2_3 <- A[[2]] - lam2.m/rho
  Z1 <- soft_thresholding(a = S1 * Z1_1 + S2 * Z1_2 + S3 * Z1_3, lam.m = hadamard.prod(lam1.m/rho,lam1.weights[[1]])) # edited: multiplying lam1 penalty by weights if we change penalty
  Z2 <- soft_thresholding(a = S1 * Z2_1 + S2 * Z2_2 + S3 * Z2_3, lam.m = hadamard.prod(lam1.m/rho,lam1.weights[[2]])) # edited
  return(list(Z1, Z2))
}
# function above and below update Z
cal_z_L1_t3 <- function(A, rho, lam1.m, lam2.m, lam1.weights){
  p <- dim(A[[1]])[1]
  vl <- dim(A[[1]])[1]^2
  # turn matrix list to vector
  A_vec <- cbind(as.vector(A[[1]]),as.vector(A[[2]]),as.vector(A[[3]])) #1,2,3
  P_vec <- cbind(as.vector(lam2.m[[3]]),as.vector(lam2.m[[2]]),as.vector(lam2.m[[1]])) / rho  # 23,13,12
  condi <- rep(0,vl)
  Z_vec <- matrix(0,nc=3,nr=vl)
  # cond 1 z1=z2=z3
  fcd.temp <- find_condi_1(av = A_vec,pv = P_vec )
  condi[fcd.temp$id] <- 1
  Z_vec[fcd.temp$id,] <- fcd.temp$zv
  # cond 2 z1>z2>z3
  index1 <- matrix(c(1,2,3,1,3,2,
                     2,1,3,2,3,1,
                     3,1,2,3,2,1),nc=3,byrow = T)
  for (indexi in 1:6){
    inputi <- which(condi==0)
    fcd.temp <- find_condi_2(av = A_vec[inputi,],pv = P_vec[inputi,], ind = index1[indexi,] )
    if (sum(condi[inputi][fcd.temp$id]>0)>0) stop("error")  # for test
    condi[inputi][fcd.temp$id] <- 2
    Z_vec[inputi[fcd.temp$id],] <- fcd.temp$zv
  }
  # cond 3 z1=z2>z3
  index2 <- matrix(c(1,2,3,2,3,1,3,1,2),nc=3,byrow = T)
  for (indexi in 1:3){
    inputi <- which(condi==0)
    fcd.temp <- find_condi_3(av = A_vec[inputi,],pv = P_vec[inputi,], ind = index2[indexi,] )
    if (sum(condi[inputi][fcd.temp$id]>0)>0) stop("error")  # for test
    condi[inputi][fcd.temp$id] <- 3
    Z_vec[inputi[fcd.temp$id],] <- fcd.temp$zv
  }
  # cond 4 z1=z2<z3
  for (indexi in 1:3){
    inputi <- which(condi==0)
    fcd.temp <- find_condi_4(av = A_vec[inputi,],pv = P_vec[inputi,], ind = index2[indexi,] )
    if (sum(condi[inputi][fcd.temp$id]>0)>0) stop("error")  # for test
    condi[inputi][fcd.temp$id] <- 4
    Z_vec[inputi[fcd.temp$id],] <- fcd.temp$zv
  }
  if (sum(condi==0)){ print("A");print(str(A_vec)); print(range(A_vec)); print("P"); print(str(P_vec)); print(range(P_vec)); stop("error: condi==0") }
  Z <- list(matrix(Z_vec[,1],nc=p,byrow = F),
            matrix(Z_vec[,2],nc=p,byrow = F),
            matrix(Z_vec[,3],nc=p,byrow = F))
  Z[[1]] <- soft_thresholding(Z[[1]], lam.m = hadamard.prod(lam1.m/rho,lam1.weights[[1]])) # edited
  Z[[2]] <- soft_thresholding(Z[[2]], lam.m = hadamard.prod(lam1.m/rho,lam1.weights[[2]])) # edited
  Z[[3]] <- soft_thresholding(Z[[3]], lam.m = hadamard.prod(lam1.m/rho,lam1.weights[[3]])) # edited
  return(Z)
}

find_condi_1 <- function(av,pv){
  av <- matrix(av,ncol=3)
  pv <- matrix(pv,ncol=3)
  z <- (av[,1]+av[,2]+av[,3])/3
  c_temp <- which( (z<=av[,1]+pv[,2]+pv[,3])&(z>=av[,1]-pv[,2]-pv[,3])&
                     (z<=av[,2]+pv[,1]+pv[,3])&(z>=av[,2]-pv[,1]-pv[,3])&
                     (z<=av[,3]+pv[,1]+pv[,2])&(z>=av[,3]-pv[,1]-pv[,2]))
  zv <- cbind(z[c_temp],z[c_temp],z[c_temp])
  return(list(zv=zv,id=c_temp))
}

find_condi_2 <- function(av,pv,ind){
  av <- matrix(av,ncol=3)
  pv <- matrix(pv,ncol=3)
  c_temp <- which(((av[,ind[1]]-pv[,ind[2]]-pv[,ind[3]])>(av[,ind[2]]+pv[,ind[3]]-pv[,ind[1]])) &
                    ((av[,ind[2]]+pv[,ind[3]]-pv[,ind[1]])>(av[,ind[3]]+pv[,ind[1]]+pv[,ind[2]]))==T)
  zv <- cbind(av[c_temp,ind[1]]-pv[c_temp,ind[2]]-pv[c_temp,ind[3]],
              av[c_temp,ind[2]]+pv[c_temp,ind[3]]-pv[c_temp,ind[1]],
              av[c_temp,ind[3]]+pv[c_temp,ind[1]]+pv[c_temp,ind[2]])
  zv[,c(ind[1],ind[2],ind[3])] <- zv
  return(list(zv=zv,id=c_temp))
}

find_condi_3 <- function(av,pv,ind){
  av <- matrix(av,ncol=3)
  pv <- matrix(pv,ncol=3)
  c_temp <- which( (abs(av[,ind[1]]-av[,ind[2]]-pv[,ind[2]]+pv[,ind[1]])<=2*pv[,ind[3]]) &
                     ((av[,ind[1]]+av[,ind[2]]-pv[,ind[2]]-pv[,ind[1]])/2>=av[,ind[3]]+pv[,ind[1]]+pv[,ind[2]]) )
  temp <- (av[c_temp,ind[1]]+av[c_temp,ind[2]]-pv[c_temp,ind[2]]-pv[c_temp,ind[1]])/2
  zv <- cbind(temp,
              temp,
              av[c_temp,ind[3]]+pv[c_temp,ind[1]]+pv[c_temp,ind[2]])
  zv[,c(ind[1],ind[2],ind[3])] <- zv
  return(list(zv=zv,id=c_temp))
}

find_condi_4 <- function(av,pv,ind){
  av <- matrix(av,ncol=3)
  pv <- matrix(pv,ncol=3)
  c_temp <- which( (abs(av[,ind[1]]-av[,ind[2]]+pv[,ind[2]]-pv[,ind[1]])<=2*pv[,ind[3]]) &
                     ((av[,ind[1]]+av[,ind[2]]+pv[,ind[2]]+pv[,ind[1]])/2<=av[,ind[3]]-pv[,ind[1]]-pv[,ind[2]]) )
  temp <- (av[c_temp,ind[1]]+av[c_temp,ind[2]]+pv[c_temp,ind[2]]+pv[c_temp,ind[1]])/2
  zv <- cbind(temp,
              temp,
              av[c_temp,ind[3]]-pv[c_temp,ind[1]]-pv[c_temp,ind[2]])
  zv[,c(ind[1],ind[2],ind[3])] <- zv
  return(list(zv=zv,id=c_temp))
}

admm.iter <- function(S,n,lam1.m,lam2.m,weights=NULL,rho=1,rho.increment=1,maxiter=500,tol=1e-05,loglik.trace=FALSE,lam1.weights){
  K <- length(S)
  p <- dim(S[[1]])[1]
  theta <- list()
  Z <- list()
  U <- list()
  for (k in 1:K) { # initialization
    theta[[k]] <- diag(1/diag(S[[k]]))
    Z[[k]] <- matrix(0, p, p)
    U[[k]] <- matrix(0, p, p)
  }
  iter <- 0
  diff_value <- 1
  loglik.tr <- rep(0,maxiter)
  DiffVal.tr <- rep(0,maxiter)
  while (iter < maxiter && diff_value > tol) {
    theta.prev <- theta
    iter <- iter + 1
    if (loglik.trace) loglik.tr[iter] <- cal_loglik(theta = Z,S = S,n = n,lam1.m = lam1.m,lam2.m = lam2.m)
    # step a update theta
    for (k in 1:K) {
      if (!is.positive.definite(theta[[k]])){print("theta not pd! quitting");q()}
    }
    # step b update Z
    A <- list()
    for (k in 1:K) {
      A[[k]] <- theta[[k]] + U[[k]]
    }
    if (K == 2) Z <- cal_z_L1_t2(A, rho, lam1.m, lam2.m,lam1.weights) # edited
    if (K == 3) Z <- cal_z_L1_t3(A, rho, lam1.m, lam2.m,lam1.weights) # edited
    # step c update U
    for (k in 1:K) {
      U[[k]] <- U[[k]] + (theta[[k]] - Z[[k]])
    }
    # difference between est
    diff_value <- 0
    for (k in 1:K) {
      diff_value <- diff_value + sum(abs(theta[[k]] - theta.prev[[k]]))/sum(abs(theta.prev[[k]]))
    }
    if (loglik.trace) DiffVal.tr[iter] <- diff_value
    rho <- rho * rho.increment
  }
  out.admm = list(theta = theta, Z = Z, iters = iter)
  if (loglik.trace) {
    out.admm$loglik.trace = loglik.tr[1:iter]
    out.admm$DiffVal.trace = DiffVal.tr[1:iter]
  }
  return(out.admm)
}

# selecting tuning parameter s automatically
s_selection <- function(expr1,expr2,ss=seq(0.1,2,0.1),verbose=F){
  p <- dim(expr1)[2]
  temp1 <- pnorm(sqrt(log(p)))
  d <- matrix(0,nr=length(ss),nc=10)
  for (i in 1:length(ss)){
    W <- get_diff_W(expr1,expr2,s = ss[i], tri = T)
    for (j in 1:10){
      temp2 <- (1-temp1)*j/10
      nomi <- sum(abs(W)>qnorm(1-temp2))
      deno <- temp2*p*(p-1)
      d[i,j] <- (nomi/deno-1)^2
    }
    if(verbose) print(paste("ss =",ss[i],"/",tail(ss,1),"   d =",sum(d[i,])))
  }
  s.slected <- ss[which.min(rowSums(d))]
  if(verbose) print(paste("selected s =",s.slected))
  return(s.slected)
}

# estimate the different between precision matrixes
# note: fixed a couple of erroneous lines from CFGL
get_diff_W <- function(expr1,expr2,s=2,tri=F){
  n1 <- dim(expr1)[1]
  n2 <- dim(expr2)[1]
  p <- dim(expr1)[2]
  covm1 <- cov(expr1)
  covm2 <- cov(expr2)
  sigma1 <- diag(covm1)
  sigma2 <- diag(covm2)
  expr1.t <- t(expr1)
  expr2.t <- t(expr2)
  b1 <- matrix(0,nr=p-1,nc=p)
  b2 <- matrix(0,nr=p-1,nc=p)
  c1 <- matrix(0,nr=n1,nc=p)
  c2 <- matrix(0,nr=n2,nc=p)
  T1 <- matrix(0,nr=p,nc=p)
  T2 <- matrix(0,nr=p,nc=p)
  W <- matrix(0,nr=p,nc=p)
  for (i in 1:p){
    y1 <- expr1[,i]
    x1 <- expr1[,-i]
    lam1 <- s*sqrt(sigma1[i]*log(p)/n1)
    temp1 <- glmnet(x1,y1,family = "gaussian",alpha = 1, lambda = lam1)
    b1[,i] <- as.matrix(coefficients(temp1)[-1])
    y2 <- expr2[,i]
    x2 <- expr2[,-i]
    lam2 <- s*sqrt(sigma2[i]*log(p)/n2)
    temp2 <- glmnet(x2,y2,family = "gaussian",alpha = 1, lambda = lam2)
    b2[,i] <- as.matrix(coefficients(temp2)[-1])
    c1[,i] <- y1-mean(y1) - (x1-rep(1,nrow(x1))%*%t(colMeans(x1))) %*%b1[,i]
    c2[,i] <- y2-mean(y2) - (x2-rep(1,nrow(x2))%*%t(colMeans(x2))) %*%b2[,i]
  }
  r1 <- t(c1)%*%c1/n1
  r2 <- t(c2)%*%c2/n2
  s1 <- colMeans(c1^2)
  s2 <- colMeans(c2^2)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      T1[i,j] <- (r1[i,j]+s1[i]*b1[i,j]+s1[j]*b1[j-1,i])/(r1[i,i]*r1[j,j])
      T2[i,j] <- (r2[i,j]+s2[i]*b2[i,j]+s2[j]*b2[j-1,i])/(r2[i,i]*r2[j,j])
      W[i,j]  <- (T1[i,j]-T2[i,j])/
        sqrt(  (1+b1[i,j]^2*r1[i,i]/r1[j,j]) / (r1[i,i]*r1[j,j]*n1) + (1+b2[i,j]^2*r2[i,i]/r2[j,j]) / (r2[i,i]*r2[j,j]*n2) )
    }
  }
  if (!tri) W <- W+t(W)
  return(W)
}

# multiple testing procedure by controlling FPR
get_W_theshold <- function(W,alpha){
  p <- dim(W)[1]
  t.upper <- 2*sqrt(log(p))
  t0 <- abs(W[upper.tri(W)])
  t1 <- t0[which(t0<t.upper)]
  t2 <- sort(t1,decreasing = T)
  temp <- (p^2-p)/2
  thes <- NULL
  use.t.upper <- F
  x=NULL
  for(i in 1:length(t2)){
    x[i] <- 2*(1-pnorm(t2[i]))*temp/i
    if (x[i]>=alpha) {
      if (i>1) {thes<-t2[i-1]}
      if (i==1) {thes<-t.upper; use.t.upper<-T}
      break
    }
  }
  W_thes <- abs(W)>=thes
  return(list(W_thes=W_thes, thes=thes, use.t.upper=use.t.upper))
}


#' Determine screening matrix by testing differential entries between 2 precision matrices.
#' The function estimates differences between 2 precision matrices. Then a multiple testing procedure with false discovery rate control will be applied to determine different entries between 2 matrices. The testing result will be turned into a binary matrix which can be used as screening matrix.
#' @param expr1 A n*p matrix or data frame of normalized gene expression data. The rows correspond samples (n), the columns correspond genes (p).
#' @param expr2 The second gene expression data, should be in the same format, size as expr1.
#' @param s The tuning parameter for matrices differences estimation, leave it as NULL to automatically select.
#' @param s.seq The candidates for s selection.
#' @param alpha Pre-specified level of false discovery rate. A relatively loose criterion is suggested for determines screening matrix.
#' @param verbose Set verbose to TURE to show details of s selection.
#' @export
get_scr_mat <- function(expr1,expr2,s=NULL,s.seq=seq(0.2,2,0.2),alpha=0.4,verbose=F){
  if (is.null(s)){
    s.selected <- s_selection(expr1,expr2,ss = s.seq,verbose = verbose)
    if (verbose) print("S was not assigned, will be auto-selected")
  }
  W <- get_diff_W(expr1 = expr1,expr2 = expr2,s = s.selected)
  temp <- get_W_theshold(W,alpha = alpha)
  W.diff.m <- !temp$W_thes
  return(list(scr.mat=W.diff.m,s=s.selected,t.threshold=temp$thes,W.test=W))
}
