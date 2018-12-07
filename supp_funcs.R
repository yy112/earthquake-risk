# supp_funcs.R
# supplementary functions to use
#######################################3
# prelec and tversky
prelec <- function(p, psi, psi2 = 1){
  w = ifelse(p!=0, exp(- psi2 * (-log(p)) ^ psi), 0)
  return(w)
}

tversky <- function(p, psi = 1){
  w <- ifelse(p!=0, p^psi / (p^psi + (1-p)^psi)^(1/psi), 0)
  return(w)
}

#########################################
# derivatives
Z_prelec <- function(p, psi, log = FALSE){
  if(log == FALSE){
    z <- ifelse(p!=0, exp(-  (-log(p)) ^ psi) * (- ((-log(p)) ^ psi)) * log(-log(p)),
                0)
  } else{
    z <- ifelse(p!=0, (- ((-log(p)) ^ psi)) * log(-log(p)),
                0)
  }
  return(z)
}
Z_tversky <- function(p, psi, log = FALSE){
  if(log == FALSE){
    z <- ifelse(p!=0,  p^psi / (p^psi + (1-p)^psi)^(1/psi) * 
                  (log(p) + 1/psi^2 * (log(p^psi + (1-p)^psi)) -
                     1/psi * (p^psi * log(p) + (1-p)^psi * log(1-p) )/(p^psi + (1-p)^psi)),
                0)
  } else{
    z <- ifelse(p!=0,  log(p) + 1/psi^2 * (log(p^psi + (1-p)^psi)) -
                  1/psi * (p^psi * log(p) + (1-p)^psi * log(1-p) )/(p^psi + (1-p)^psi),
                0)
  }
  return(z)
}
#############################################################################3
# X_invOmega_Y.R
# This function calculates the functions in the form 
# X'*Inv(Omega)*y. Number of periods T, number of areas N

X_invOmega_Y <- function(X, Y, N, Tn, p, 
                         O1inv, O2inv, O3inv, O4inv){
  # X and Y are of dimensions (NTp)*k, k can be different
  # Deltas are of dimension p*p
  # number of cols of X and Y
  k1 <- ifelse(length(ncol(X))==0, 1, ncol(X))
  k2 <- ifelse(length(ncol(Y))==0, 1, ncol(Y))
  X <- as.matrix(X, ncol=k1)
  Y <- as.matrix(Y, ncol=k2)
  # iFail <- 1
  # if(nrow(X) == N*Tn*p & nrow(Y) == N*Tn*p & 
  #    nrow(O1inv) == p & ncol(O1inv) == p & 
  #    nrow(O2inv) == p & ncol(O2inv) == p & 
  #    nrow(O3inv) == p & ncol(O3inv) == p & 
  #    nrow(O4inv) == p & ncol(O4inv) == p ){
  #   iFail <- 0
  # } else {print("Dimension mismatch")}
  # 
  
  term1 <- matrix(0, nrow = k1, ncol = k2)
  if(!all((O1inv - O2inv - O3inv + O4inv)==0)){
    sum_it_X <- matrix(0, nrow = p, ncol = k1)
    sum_it_Y <- matrix(0, nrow = p, ncol = k2)
    for(it in 1:(N*Tn)){
      sum_it_X <- sum_it_X + X[((it-1)*p+1):(it*p), ]
      sum_it_Y <- sum_it_Y +  Y[((it-1)*p+1):(it*p), ]
    }
    term1 <- t(sum_it_X) %*% (O1inv - O2inv - O3inv + O4inv) %*%
      as.matrix(sum_it_Y)
  }
  
  
  
  term2 <- matrix(0, nrow = k1, ncol = k2)
  if(!all((O2inv - O4inv)==0)){
    for(i in 1:N){
      sumtXDsumtY <- matrix(0, nrow = k1, ncol = k2)
      sumtX <- matrix(0, nrow = p, ncol = k1)
      sumtY <- matrix(0, nrow = p, ncol = k2)
      for(t in 1:Tn){
        sumtX <- sumtX + 
          (X[(((t-1)*N+i-1)*p+1):(((t-1)*N+i)*p), ]) 
        sumtY <- sumtY + 
          (Y[(((t-1)*N+i-1)*p+1):(((t-1)*N+i)*p), ]) 
      }
      sumtXDsumtY <- t(sumtX) %*% (O2inv - O4inv) %*% 
        as.matrix(sumtY)
      term2 <- term2 + sumtXDsumtY
    }
  }
  
  
  
  term3 <- matrix(0, nrow = k1, ncol = k2)
  if(!all((O3inv - O4inv)==0)){
    for(t in 1:Tn){
      sumiXDsumiY <- matrix(0, nrow = k1, ncol = k2)
      sumiX <- matrix(0, nrow = p, ncol = k1)
      sumiY <- matrix(0, nrow = p, ncol = k2)
      for(i in 1:N){
        sumiX <- sumiX + 
          X[(((t-1)*N+i-1)*p+1):(((t-1)*N+i)*p), ]
        sumiY <- sumiY + 
          Y[(((t-1)*N+i-1)*p+1):(((t-1)*N+i)*p), ]
      }
      sumiXDsumiY <- t(sumiX) %*%
        (O3inv - O4inv) %*% 
        as.matrix(sumiY)
      term3 <- term3 + sumiXDsumiY
    }
  }
  
  
  term4 <- matrix(0, nrow = k1, ncol = k2)
  if(!all(O4inv==0) & length(O4inv)!=1){
    for(it in 1:(N*Tn)){
      term4 <- term4 + t(X[((it-1)*p+1):(it*p), ]) %*% 
        O4inv %*% as.matrix(Y[((it-1)*p+1):(it*p), ])
    }
  }else if(!all(O4inv==0) & length(O4inv)==1){
    for(it in 1:(N*Tn)){
      term4 <- term4 + as.numeric(O4inv) * 
        as.matrix(X[((it-1)*p+1):(it*p), ]) %*% t(as.matrix(Y[((it-1)*p+1):(it*p), ]))
    }
  }
  
  
  XOY <- (1/N)*(1/Tn)*term1 + (1/Tn)*term2 + 
    (1/N)*term3 + term4
  return(XOY)
}




####################################################################################
gradient <- function(par, include, N, Tn, HX, Hy, p = 3){
  # This function is the gradient function of the full model 
  # (with 3*p*(p+1)/2 parameters).
  npar <- length(par)
  full.par <- rep(0, 3*p*(p+1)/2)
  full.par[which(include!=0)] <- par
  full.par <- matrix(full.par, ncol=3, byrow = F)
  tmp <- matrix(0, p, p)
  eye <- diag(1, p)
  if(npar != sum(include) | length(include)!=(3*p*(p+1)/2)){
    print("number of parameters not match")
  } else if(p==3){
    L.zeta <- tmp
    L.zeta[lower.tri(tmp, diag = T)] <- full.par[, 1]
    diag(L.zeta) <- exp(diag(L.zeta)) * include[c(1, 4, 6)]
    L.eta <- tmp
    L.eta[lower.tri(tmp, diag = T)] <- full.par[, 2]
    diag(L.eta) <- exp(diag(L.eta)) * include[c(7, 10, 12)]
    L.eps <- tmp
    L.eps[lower.tri(tmp, diag = T)] <- full.par[, 3]
    diag(L.eps) <- exp(diag(L.eps)) * include[c(13, 16, 18)]
    # generate sigma_zeta, sigma_eta, sigma_eps
    sigmazeta <-  tcrossprod(L.zeta)
    sigmaeta <- tcrossprod(L.eta)
    sigmaeps <- tcrossprod(L.eps)
    # generate Omega matrices
    Omega1 <- sigmaeps + Tn * sigmazeta + N * sigmaeta
    Omega2 <- sigmaeps + Tn * sigmazeta
    Omega3 <- sigmaeps +  N * sigmaeta
    Omega4 <- sigmaeps 
    O1inv <- solve(Omega1)
    O2inv <- solve(Omega2)
    O3inv <- solve(Omega3)
    O4inv <- solve(Omega4)
    XOX <- X_invOmega_Y(HX, HX, N, Tn, p, 
                        O1inv, O2inv, O3inv, O4inv)
    XOy <- X_invOmega_Y(HX, Hy, N, Tn, p, 
                        O1inv, O2inv, O3inv, O4inv)
    beta <- solve(XOX, XOy)
    vit <- Hy - HX %*% beta
    vOX <- X_invOmega_Y(vit, HX, N, Tn, p, 
                        O1inv, O2inv, O3inv, O4inv)
    gr <- rep(0, 3*p*(p+1)/2)
    for(m in 1:(3*p*(p+1)/2)){
      if(include[m]!=0){
        n.mat <- ceiling(m/(p*(p+1)/2))
        j <- which((m - (p*(p+1)/2)*(n.mat-1)) <= 
                     cumsum(seq(p,1,-1)))[1]
        i <- 3 + m - (p*(p+1)/2)*(n.mat-1) - 
          cumsum(seq(p,1,-1))[j]
        ei <- matrix(0, nrow = p, ncol = 1) 
        ei[i, 1] <- 1
        ej <- matrix(0, nrow = p, ncol = 1) 
        ej[j, 1] <- 1
        # term1: tr(inv(Omega) dOmega)
        if(n.mat==1){
          eeL <-  ei %*% t(ej) %*%  t(L.zeta) + 
            L.zeta %*% ej %*% t(ei)
          term1 <- (Tn*sum(diag(O1inv %*% eeL+(N-1)*O2inv %*% eeL)))* 
            ifelse(i==j, L.zeta[i, j], 1)
          dE1 <- Tn * O1inv %*% eeL %*% O1inv * 
            ifelse(i==j, L.zeta[i, j], 1)
          dE2 <- Tn * O2inv %*% eeL %*% O2inv * 
            ifelse(i==j, L.zeta[i, j], 1)
          dE3 <- 0
          dE4 <- 0
        } else if(n.mat==2) {
          eeL <-  ei %*% t(ej) %*%  t(L.eta) + 
            L.eta %*% ej %*% t(ei)
          term1 <- (N*sum(diag(O1inv %*% eeL+(Tn-1)*O3inv %*% eeL)))* 
            ifelse(i==j, L.eta[i, j], 1)
          dE1 <- N * O1inv %*% eeL %*% O1inv * 
            ifelse(i==j, L.eta[i, j], 1)
          dE2 <- 0
          dE3 <- N * O3inv %*% eeL %*% O3inv * 
            ifelse(i==j, L.eta[i, j], 1)
          dE4 <- 0
        } else if(n.mat==3){
          eeL <-  ei %*% t(ej) %*%  t(L.eps) + 
            L.eps %*% ej %*% t(ei)
          term1 <- (sum(diag(O1inv %*% eeL + (N - 1)*O2inv %*% eeL +
                               (Tn - 1)*O3inv %*% eeL + 
                               (N - 1)*(Tn -1 )*O4inv %*% eeL))) *
            ifelse(i==j, L.eps[i, j], 1)
          dE1 <- O1inv %*% eeL %*% O1inv * 
            ifelse(i==j, L.eps[i, j], 1)
          dE2 <- O2inv %*% eeL %*% O2inv * 
            ifelse(i==j, L.eps[i, j], 1)
          dE3 <- O3inv %*% eeL %*% O3inv * 
            ifelse(i==j, L.eps[i, j], 1)
          dE4 <- O4inv %*% eeL %*% O4inv * 
            ifelse(i==j, L.eps[i, j], 1)
        }
        term2 <- vOX %*% solve(XOX) %*% 
          X_invOmega_Y(HX, vit, N, Tn, p, dE1, dE2, dE3, dE4) 
        term3 <- X_invOmega_Y(vit, vit, N, Tn, p, 
                              dE1, dE2, dE3, dE4) 
        gr[m] <-  1/2 * term1 + term2 - 1/2 * term3 
        
      } 
    }
  }else if(p==1){
    L.zeta <- exp(full.par[, 1]) * include[1]
    L.eta <- exp(full.par[, 2]) * include[2]
    L.eps <- exp(full.par[, 3]) * include[3]
    # generate sigma_zeta, sigma_eta, sigma_eps
    sigmazeta <-  tcrossprod(L.zeta)
    sigmaeta <- tcrossprod(L.eta)
    sigmaeps <- tcrossprod(L.eps)
    # generate Omega matrices
    Omega1 <- sigmaeps + Tn * sigmazeta + N * sigmaeta
    Omega2 <- sigmaeps + Tn * sigmazeta
    Omega3 <- sigmaeps +  N * sigmaeta
    Omega4 <- sigmaeps 
    O1inv <- solve(Omega1)
    O2inv <- solve(Omega2)
    O3inv <- solve(Omega3)
    O4inv <- solve(Omega4)
    XOX <- X_invOmega_Y(HX, HX, N, Tn, p, 
                        O1inv, O2inv, O3inv, O4inv)
    XOy <- X_invOmega_Y(HX, Hy, N, Tn, p, 
                        O1inv, O2inv, O3inv, O4inv)
    beta <- solve(XOX, XOy)
    vit <- Hy - HX %*% beta
    vOX <- X_invOmega_Y(vit, HX, N, Tn, p, 
                        O1inv, O2inv, O3inv, O4inv)
    gr <- rep(0, 3*p*(p+1)/2)
    for(m in 1:(3*p*(p+1)/2)){
      if(include[m]!=0){
        n.mat <- ceiling(m/(p*(p+1)/2))
        j <- which((m - (p*(p+1)/2)*(n.mat-1)) <= 
                     cumsum(seq(p,1,-1)))[1]
        i <- 3 + m - (p*(p+1)/2)*(n.mat-1) - 
          cumsum(seq(p,1,-1))[j]
        ei  <- 1
        ej  <- 1
        # term1: tr(inv(Omega) dOmega)
        if(n.mat==1){
          eeL <-  ei %*% t(ej) %*%  t(L.zeta) + 
            L.zeta %*% ej %*% t(ei)
          term1 <- (Tn*sum(diag(O1inv %*% eeL+(N-1)*O2inv %*% eeL)))* 
            ifelse(i==j, L.zeta[i, j], 1)
          dE1 <- Tn * O1inv %*% eeL %*% O1inv * 
            ifelse(i==j, L.zeta[i, j], 1)
          dE2 <- Tn * O2inv %*% eeL %*% O2inv * 
            ifelse(i==j, L.zeta[i, j], 1)
          dE3 <- 0
          dE4 <- 0
        } else if(n.mat==2) {
          eeL <-  ei %*% t(ej) %*%  t(L.eta) + 
            L.eta %*% ej %*% t(ei)
          term1 <- (N*sum(diag(O1inv %*% eeL+(Tn-1)*O3inv %*% eeL)))* 
            ifelse(i==j, L.eta[i, j], 1)
          dE1 <- N * O1inv %*% eeL %*% O1inv * 
            ifelse(i==j, L.eta[i, j], 1)
          dE2 <- 0
          dE3 <- N * O3inv %*% eeL %*% O3inv * 
            ifelse(i==j, L.eta[i, j], 1)
          dE4 <- 0
        } else if(n.mat==3){
          eeL <-  ei %*% t(ej) %*%  t(L.eps) + 
            L.eps %*% ej %*% t(ei)
          term1 <- (sum(diag(O1inv %*% eeL + (N - 1)*O2inv %*% eeL +
                               (Tn - 1)*O3inv %*% eeL + 
                               (N - 1)*(Tn -1 )*O4inv %*% eeL))) *
            ifelse(i==j, L.eps[i, j], 1)
          dE1 <- O1inv %*% eeL %*% O1inv * 
            ifelse(i==j, L.eps[i, j], 1)
          dE2 <- O2inv %*% eeL %*% O2inv * 
            ifelse(i==j, L.eps[i, j], 1)
          dE3 <- O3inv %*% eeL %*% O3inv * 
            ifelse(i==j, L.eps[i, j], 1)
          dE4 <- O4inv %*% eeL %*% O4inv * 
            ifelse(i==j, L.eps[i, j], 1)
        }
        term2 <- vOX %*% solve(XOX) %*% 
          X_invOmega_Y(HX, vit, N, Tn, p, dE1, dE2, dE3, dE4) 
        term3 <- X_invOmega_Y(vit, vit, N, Tn, p, 
                              dE1, dE2, dE3, dE4) 
        gr[m] <-  1/2 * term1 + term2 - 1/2 * term3 
        
      } 
    }
  }
  gr <- gr[which(include!=0)]
  return(gr)
} 


#####################################################################
negloglik <- function(par, include, N, Tn, HX, Hy, p = 3){
  # include is the dummy for whether to include non-zero 
  # theta parameters, in order of L.zeta, L.eta, L.eps, by col  
  npar <- length(par)
  full.par <- rep(0, 3*p*(p+1)/2)
  full.par[which(include!=0)] <- par
  full.par <- matrix(full.par, ncol=3, byrow = F)
  tmp <- matrix(0, p, p)
  if(npar != sum(include) | length(include)!=(3*p*(p+1)/2)){
    print("number of parameters not match")
  } else if(p==3){
    L.zeta <- tmp
    L.zeta[lower.tri(tmp, diag = T)] <- full.par[, 1]
    diag(L.zeta) <- exp(diag(L.zeta)) * include[c(1, 4, 6)]
    L.eta <- tmp
    L.eta[lower.tri(tmp, diag = T)] <- full.par[, 2]
    diag(L.eta) <- exp(diag(L.eta)) * include[c(7, 10, 12)]
    L.eps <- tmp
    L.eps[lower.tri(tmp, diag = T)] <- full.par[, 3]
    diag(L.eps) <- exp(diag(L.eps)) * include[c(13, 16, 18)]
    # generate sigma_zeta, sigma_eta, sigma_eps
    sigmazeta <-  tcrossprod(L.zeta)
    sigmaeta <- tcrossprod(L.eta)
    sigmaeps <- tcrossprod(L.eps)
    # generate Omega matrices
    Omega1 <- sigmaeps + Tn * sigmazeta + N * sigmaeta
    Omega2 <- sigmaeps + Tn * sigmazeta
    Omega3 <- sigmaeps +  N * sigmaeta
    Omega4 <- sigmaeps 
    O1inv <- solve(Omega1)
    O2inv <- solve(Omega2)
    O3inv <- solve(Omega3)
    O4inv <- solve(Omega4)
    XOX <- X_invOmega_Y(HX, HX, N, Tn, p, 
                        O1inv, O2inv, O3inv, O4inv)
    XOy <- X_invOmega_Y(HX, Hy, N, Tn, p, 
                        O1inv, O2inv, O3inv, O4inv)
    #print(solve(XOX))
    beta <- solve(XOX, XOy)
    # vectors of residuals, each column is for a different "type"
    # this is a (N*Tn) * p matrix
    # inside each col, {t1, i1}, {t1, i2}, {t1, i3}, ...{t2, i1}, {t2, i2}, {t2, i3}..
    vit <- Hy - HX %*% beta
    negll <- 1/2 * (log(det(Omega1)) + (N-1) * log(det(Omega2)) + 
                      (Tn-1) * log(det(Omega3)) + (N-1) * (Tn-1) * log(det(Omega4)) + 
                      X_invOmega_Y(vit, vit, N, Tn, p, 
                                   O1inv, O2inv, O3inv, O4inv))
    
    return(negll)
  } else if(p==1){
    L.zeta <- exp(full.par[, 1]) * include[1]
    L.eta <- exp(full.par[, 2]) * include[2]
    L.eps <- exp(full.par[, 3]) * include[3]
    
    # generate sigma_zeta, sigma_eta, sigma_eps
    sigmazeta <-  tcrossprod(L.zeta)
    sigmaeta <- tcrossprod(L.eta)
    sigmaeps <- tcrossprod(L.eps)
    # generate Omega matrices
    Omega1 <- sigmaeps + Tn * sigmazeta + N * sigmaeta
    Omega2 <- sigmaeps + Tn * sigmazeta
    Omega3 <- sigmaeps +  N * sigmaeta
    Omega4 <- sigmaeps 
    O1inv <- solve(Omega1)
    O2inv <- solve(Omega2)
    O3inv <- solve(Omega3)
    O4inv <- solve(Omega4)
    XOX <- X_invOmega_Y(HX, HX, N, Tn, p, 
                        O1inv, O2inv, O3inv, O4inv)
    XOy <- X_invOmega_Y(HX, Hy, N, Tn, p, 
                        O1inv, O2inv, O3inv, O4inv)
    #print(solve(XOX))
    beta <- solve(XOX, XOy)
    # vectors of residuals, each column is for a different "type"
    # this is a (N*Tn) * p matrix
    # inside each col, {t1, i1}, {t1, i2}, {t1, i3}, ...{t2, i1}, {t2, i2}, {t2, i3}..
    vit <- Hy - HX %*% beta
    negll <- 1/2 * (log(det(Omega1)) + (N-1) * log(det(Omega2)) + 
                      (Tn-1) * log(det(Omega3)) + (N-1) * (Tn-1) * log(det(Omega4)) + 
                      X_invOmega_Y(vit, vit, N, Tn, p, 
                                   O1inv, O2inv, O3inv, O4inv))
    
    return(negll)
    
  }
}

####################################################################
# vectorize.R
# This function reads the data frames containing y, X,
# and information for i(district), t(time period) and p(type), 
# and take the average over each i-t-p combination.
# The output is a vector y, matrix X, and H (number of observations for each cell).
# All data should be numericals except for the 
# columns containing i, t and p.


vectorize <- function(data, colName.i, colName.t, colName.p, colName.y, colName.X){
  # retrieve the unique names of district, time and type
  distr <- sort(unique(data[, colName.i]))
  time <- sort(unique(data[, colName.t]))
  type <- unique(data[, colName.p])
  # generate a tabulated dataframe 
  tab <- expand.grid(type, distr, time)[, c(3, 2, 1)]
  colnames(tab) <- c(colName.t, colName.i, colName.p)
  # factorise the i/t/p columns so that they are ordered. 
  # this is important for generating the averages.
  data[, colName.i] <- factor(x = data[, colName.i], levels = distr)
  data[, colName.t] <- factor(x = data[, colName.t], levels = time)
  data[, colName.p] <- factor(x = data[, colName.p], levels = type)
  # group data by i-t-p and take the means of each variable
  # result <- tab
  # result$H <- 0
  # result[, colName.y] <- 0
  # result[, colName.X] <- 0
  # for(i in 1:nrow(tab)){
  #   rows <- data[(data[, colName.t]==tab$t[i])&
  #     (data[, colName.i]==tab$Area.Ward.City[i])&
  #     (data[, colName.p]==tab$Type[i]), , drop=FALSE]
  #   result$H[i] <- nrow(rows)
  #   means <- sapply(Filter(is.numeric, rows), mean)
  #   result[i, colName.y] <- means[colName.y]
  #   result[i, colName.X] <- means[colName.X]
  #   if(i %% 10000 ==0){
  #     print(i)
  #   }
  
  sumstat <- data %>%
    dplyr:::group_by_(colName.t, colName.i, colName.p) %>%
    dplyr:::summarise_at(c(colName.X, colName.y), mean) 
  sumH <- data %>%
    dplyr:::group_by_(colName.t, colName.i, colName.p) %>%
    dplyr:::summarise(n = n()) 
  # merge the name table "tab" with averages in table "sumstat"
  # and counts in table "sumH"
  result <- merge(tab, sumstat, all=TRUE)
  result <- merge(result, sumH, all=TRUE)
  # for(i in 1:nrow(sumcode)){
  #   rows <- which(result$StationCity==sumcode$StationCity[i])
  #   result$City.Town.Ward.Village.code[rows] <- sumcode$City.Town.Ward.Village.code[i] 
  # }
  # for(i in 1:nrow(result)){
  #   if(!is.na(result[i, "n"])){
  #     xrows <- which(data[, colName.i]==result[i, colName.i] & 
  #                      data[, colName.t]==result[i, colName.t])
  #     result[i, col.code] <- unique(data[xrows, col.code])[1]
  #   }
  # }
  y <- result[, colName.y]
  H <- result[, "n"]
  X <- result[, c(colName.t, colName.i, colName.p, colName.X)]
  # replace the NA rows with zero
  H[is.na(H)] <- 0
  y[is.na(y)] <- 0
  X[is.na(X)] <- 0
  # 
  return(list("district" = distr, "time" = time, "type" = type,
              "H" = H, "y" = y, "X" = X))
  
}

###################################################
# ec_reg.R
# error components regression 
ec_reg <- function(data_X, data_y, data_H, 
                   colName.i = "i", colName.t = "t",  colName.p = "Type",
                   district, time, type, 
                   var,
                   par.include = rep(1, 18),
                   par.init = rep(0.5, 18)){
  # The main optimization procedure for the regression model with multivariate error components. 
  N <- length(district)
  Tn <- length(time)
  X.num <- data_X[, var]
  k <- length(var)
  p <- length(type)
  # supplementary matrices
  JT <- matrix(1, ncol=Tn, nrow=Tn)/Tn
  JN <- matrix(1, ncol=N, nrow=N)/N
  Xnum <- X.num
  for(i in 1:ncol(Xnum)){
    Xnum[,i]<-as.numeric(Xnum[,i])
  }
  Hy <- as.vector(sqrt(data_H)) * data_y
  Hy[data_H==0] <- 0
  HX <- as.vector(sqrt(data_H)) * Xnum
  HX[data_H==0, ] <- 0
  HX <- as.matrix(HX)
  # do optimization
  results <- data.frame(var = character(k),
                        par.wgr = numeric(k),
                        count.wgr  = numeric(k),
                        gradient = numeric(k),
                        val = numeric(k),
                        coef = numeric(k),
                        t.stat = numeric(k))

  init <- par.init[par.include==1]
  result.wgr <- optim(par = init, fn=negloglik, gr=gradient,
                      include = par.include, N = N, Tn = Tn, p = p,HX = HX, Hy = Hy,
                      method="L-BFGS-B", hessian = T,
                      control = list(trace = 2, maxit = 500))
  results$var <- var
  results$par.wgr[1:sum(par.include)] <- result.wgr$par
  results$count.wgr[1] <- result.wgr$counts[1]
  results$val[1] <- result.wgr$value[1]
  par <- result.wgr$par
  npar <- length(par)
  full.par <- rep(0, 3*p*(p+1)/2)
  full.par[which(par.include!=0)] <- par
  full.par <- matrix(full.par, ncol=3, byrow = F)
  tmp <- matrix(0, p, p)
  L.zeta <- tmp
  L.zeta[lower.tri(tmp, diag = T)] <- full.par[, 1]
  diag(L.zeta) <- exp(diag(L.zeta)) * par.include[c(1, 4, 6)]
  L.eta <- tmp
  L.eta[lower.tri(tmp, diag = T)] <- full.par[, 2]
  diag(L.eta) <- exp(diag(L.eta)) * par.include[c(7, 10, 12)]
  L.eps <- tmp
  L.eps[lower.tri(tmp, diag = T)] <- full.par[, 3]
  diag(L.eps) <- exp(diag(L.eps)) * par.include[c(13, 16, 18)]
  # generate sigma_zeta, sigma_eta, sigma_eps
  sigmazeta <-  tcrossprod(L.zeta)
  sigmaeta <- tcrossprod(L.eta)
  sigmaeps <- tcrossprod(L.eps)
  Omega1 <- sigmaeps + Tn * sigmazeta + N * sigmaeta
  Omega2 <- sigmaeps + Tn * sigmazeta
  Omega3 <- sigmaeps +  N * sigmaeta
  Omega4 <- sigmaeps
  O1inv <- solve(Omega1)
  O2inv <- solve(Omega2)
  O3inv <- solve(Omega3)
  O4inv <- solve(Omega4)
  XOX <- X_invOmega_Y(HX, HX, N, Tn, p,
                      O1inv, O2inv, O3inv, O4inv)
  XOy <- X_invOmega_Y(HX, Hy, N, Tn, p,
                      O1inv, O2inv, O3inv, O4inv)
  V11inv <- solve(XOX)
  beta <- solve(XOX, XOy)
  results$coef <- as.vector(beta)
  var_beta <- V11inv
  t_beta <- beta/sqrt(diag(var_beta))
  results$t.stat <- as.vector(t_beta)
  return(results)
}

###################################################################
# reg_psi.R
# When one of the regressors depend on a parameter (psi), 
# Adjust the corresponding variance of beta and calculate the variance of psi.
# Here we assume Xpsi(psi) is a function of Xpsi_obj with parameter psi. 
# The functional form can be a Prelec function (method="p") 
# or a Tversky-Kahneman function (method="t"). 
# functions prelec, tversky, and their derivative Z_prelec, Z_tversky are needed.
reg_psi <- function(data_X, data_y, data_H, 
                    colName.i = "i", colName.t = "t",  colName.p = "Type",
                    psi = 1, method = "p",
                    district, time, type, 
                    var, 
                    par.include = rep(1, 18),
                    par.init = rep(0.5, 18)){
  if(! "Xpsi"%in%var){
    print("Xpsi not in list of regressors!")
    return(0)
  } else{
    data_X$Xpsi <- 0
    if(method == "p"){
      data_X$Xpsi <- ifelse(data_X$Xpsi_obj==0, 0, prelec(data_X$Xpsi_obj, psi = psi))
      z_psi <-  Z_prelec(p = (data_X$Xpsi_obj), psi = psi, log = FALSE)
      Hz <- as.vector(sqrt(data_H)) * z_psi
    } else if(method == "t"){
      data_X$Xpsi <- ifelse(data_X$Xpsi_obj==0, 0, tversky(data_X$Xpsi_obj, psi = psi))
      z_psi <-  Z_tversky(p = (data_X$Xpsi_obj), psi = psi, log = FALSE)
      Hz <- as.vector(sqrt(data_H)) * z_psi
    }
    include_2error <- c(rep(1, 6), rep(0,6), rep(1,6))
    # include_3error <- rep(1, 18)
    initpar <- c(-1.5, 0.1, -0.1, -2.1, -0.01, -1.1,
                 -2.5, 0.02, -0.1, -3.5,  0.1, -3.5,
                 -0.9, 0.008, 0.005, -0.9, 0.01, -0.9)
    results <- ec_reg(data_X,  data_y, data_H, 
                      colName.i = "Area.Ward.City", colName.t = "t", 
                      colName.p = "Type", 
                      var = var,
                      district = district, time = time, type = type,
                      par.include = include_2error, 
                      par.init = initpar)
    # given the regression results, reconstruct the variances 
    N <- length(district)
    Tn <- length(time)
    p <- length(type)
    Hy <- as.vector(sqrt(data_H)) * data_y
    Hy[data_H==0] <- 0
    HX <- as.vector(sqrt(data_H)) * data_X[, var]
    HX[data_H==0, ] <- 0
    HX <- as.matrix(HX)
    par.include <- include_2error
    npar <- sum(par.include!=0)
    par <- results$par.wgr[1:npar]
    full.par <- rep(0, 3*p*(p+1)/2)
    full.par[which(par.include!=0)] <- par
    full.par <- matrix(full.par, ncol=3, byrow = F)
    tmp <- matrix(0, p, p)
    L.zeta <- tmp
    L.zeta[lower.tri(tmp, diag = T)] <- full.par[, 1]
    diag(L.zeta) <- exp(diag(L.zeta)) * par.include[c(1, 4, 6)]
    L.eta <- tmp
    L.eta[lower.tri(tmp, diag = T)] <- full.par[, 2]
    diag(L.eta) <- exp(diag(L.eta)) * par.include[c(7, 10, 12)]
    L.eps <- tmp
    L.eps[lower.tri(tmp, diag = T)] <- full.par[, 3]
    diag(L.eps) <- exp(diag(L.eps)) * par.include[c(13, 16, 18)]
    # generate sigma_zeta, sigma_eta, sigma_eps
    sigmazeta <-  tcrossprod(L.zeta)
    sigmaeta <- tcrossprod(L.eta)
    sigmaeps <- tcrossprod(L.eps)
    Omega1 <- sigmaeps + Tn * sigmazeta + N * sigmaeta
    Omega2 <- sigmaeps + Tn * sigmazeta
    Omega3 <- sigmaeps +  N * sigmaeta
    Omega4 <- sigmaeps
    O1inv <- solve(Omega1)
    O2inv <- solve(Omega2)
    O3inv <- solve(Omega3)
    O4inv <- solve(Omega4)
    XOX <- X_invOmega_Y(HX, HX, N, Tn, p,
                        O1inv, O2inv, O3inv, O4inv)
    XOy <- X_invOmega_Y(HX, Hy, N, Tn, p,
                        O1inv, O2inv, O3inv, O4inv)
    V11inv <- solve(XOX)
    V12 <- X_invOmega_Y(HX, Hz, N, Tn, p,
                        O1inv, O2inv, O3inv, O4inv)
    V21 <- X_invOmega_Y(Hz, HX, N, Tn, p,
                        O1inv, O2inv, O3inv, O4inv)
    V22 <- X_invOmega_Y(Hz, Hz, N, Tn, p,
                        O1inv, O2inv, O3inv, O4inv)
    beta <- solve(XOX, XOy)
    beta_k <- results$coef[which(results$var=="Xpsi")]
    var_beta <- V11inv + V11inv %*% V12 %*% V21 %*% V11inv /
      as.numeric((V22 - V21 %*% V11inv %*% V12))
    t_beta <- beta/sqrt(diag(var_beta))
    results$t.stat <- as.vector(t_beta)
    var_psi <- 1/ ( beta_k^2 * (V22 - V21 %*% V11inv %*% V12) )
      tstat_psi <- (psi-1) / sqrt(var_psi)
      results_psi <- rbind(results, c(-results$val[1], 0,0,0,0, psi, tstat_psi))
  }
  return(results_psi)
}


#################################################################################
# opt_psi.R
# given a list of psi's,  this function finds the psi that maximizes likelihood
opt_psi <- function(data_X, data_y, data_H, 
                    list_psi = 1, method = "p",
                    district, time, type,
                    var,
                    par.include = rep(1, 18),
                    par.init = rep(0.5, 18),
                    plot = FALSE){
  list_logL <- numeric(length(list_psi))
  for(i in 1:length(list_psi)){
    res <- reg_psi(data_X, data_y, data_H,
            psi = list_psi[i], method = "p",
            district = district, time = time, type = type,
            var = Xnames_base,
            par.include,
            par.init)
    list_logL[i] <- -res$val[1]
  }
  if(plot==TRUE){
    plot(list_psi, list_logL)
  }
  return(list_psi[which.max(list_logL)])
}
