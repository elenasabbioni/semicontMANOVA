###############################################################################
#### functions needed to run CHEN et al. method on a simulated scenario
#### with semicontinuous high-dimensional data
###############################################################################


###### theta_lambda ######
## function that returns the values of \hat{Theta}_1 and \hat{Theta}_2 described in (7) and (8) of Chen et al.
#' @param lambda , scalar, regularization parameter
#' @param S, variance and covariance matrix
#' @param n , vector with the number of observations for each group
# It returns a vector with two components,\hat{Theta}_1 and \hat{Theta}_2
theta_lambda <- function(lambda, S, n){
    n <- sum(n)
    p <- dim(S)[1] 
    m <- 1/p*sum(diag((solve(S + lambda*diag(1, p)))))
    m1 <- sum(diag((solve((S + lambda*diag(1, p))%*%(S + lambda*diag(1, p))))))/p
    
    theta1_lambda <- (1 - lambda*m)/(1 - p/n*(1 - lambda*m))

    theta2_lambda <- (1 - lambda*m)/(1 - p/n + p/n*lambda*m)^3 - lambda*(m - lambda*m1)/(1 - p/n + p/n*lambda*m)^4

    return(c(theta1_lambda, theta2_lambda))
}
########################


###### Dpart_lambda ######
## function that returns difference between the asymptotic (1 − \alpha) quantile of Dependent Null and Independent Null 
#' @param lambda , scalar, regularization parameter
#' @param S , variance and covariance matrix
#' @param theta , vector with \hat{Theta}_1 and \hat{Theta}_2
#' @param p , number of components
#' @param n , vector with the number of observations for each group
#' @param alpha , target significance level of the test
# It returns a scalar with D_alpha(\lambda)
Dpart_lambda <- function(lambda, S = NA, theta = NA, p, n, alpha = 0.05){
    if(is.na(theta)){
        theta <- theta_lambda(lambda, S, n)
    }

    eps <- qnorm(1 - alpha)
    D <- (theta[1] + sqrt(2)*eps/sqrt(p)*sqrt(theta[2])) 
    return(D)
}
########################




###### estimationCHEN ######
## function that returns the estimated parameters of the model 
## and run the regularized MANOVA test
#' @param x , n x p matrix with the data
#' @param n , vector with the number of observations for each group
#' @param Lambda , vector of regularization parameters in which we search for the optimal lambda
#' @param B , number of permutations used in the permutation test
#' @param delta , parameter used to check to condition on lambda. According to the authors it should be a small positive value
#' @param alpha , target significance level of the test
# It outputs the estimated p-value of the permutation test

estimationCHEN <- function(x, n, Lambda = seq(0.01, 10, 0.05), B = 100, delta = 0.5, alpha = 0.05){
    
    p <- dim(x)[2] # numer of components
    x[which(x == 0)] <- NA 

    # the approach is for two groups'testing (describe din the Supplementary material of Chen)
    dati1 <- x[1:n[1],]
    dati2 <- x[(n[1] + 1):(sum(n)),]
    
    # sample means
    mean1 <- apply(dati1, FUN = mean, MARGIN = 2, rm.na = TRUE)
    mean2 <- apply(dati2, FUN = mean, MARGIN = 2, rm.na = TRUE)
    # deal with missing data 
    mean1[which(is.na(mean1))] <- 0
    mean2[which(is.na(mean2))] <- 0

    # sample variance and covariance matrix
    tmp1 <- matrix(0, nrow = p, ncol = p)
    for(i in 1:n[1]){
        tt <- (dati1[i, ] - mean1)%*%t(dati1[i, ] - mean1)
        if(is.null(dim(tt))){
            tt <- matrix(0, nrow = p, ncol = p)
        }
        tt[which(is.na(tt))] <- 0
        tmp1 <- tmp1 + tt
    }

    tmp2 <- matrix(0, nrow = p, ncol = p)
    for(i in 1:n[2]){
        tt <- (dati2[i, ] - mean2)%*%t(dati2[i, ] - mean2)
        if(is.null(dim(tt))){
            tt <- matrix(0, nrow = p, ncol = p)
        }
        tt[which(is.na(tt))] <- 0
        tmp2 <- tmp2 + tt
    }
    S <- 1/(sum(n) - 2)*(tmp1 + tmp2)
    
    # decomposition of S
    r <- eigen(S)
    U <- r$vectors
    D <- diag(r$values)

    # difference between the asymptotic (1-alpha) quantile of DN and IN for the different values of Lambda
    D_orig <- c()
    for(lambda in Lambda){
        D_orig <- c(D_orig, Dpart_lambda(lambda, S = S, theta = NA, p, n, alpha = alpha))
    }

    # resampling
    resample <- array(0, dim = c(dim(x), B))
    lambda_B <- rep(0, B)
    m1 <- matrix(0, nrow = p, ncol = B)
    m2 <- matrix(0, nrow = p, ncol = B)
    S_r <- array(0, dim = c(p, p, B))

    for(i in 1:B){
        resample[,,i] <- x[sample(sum(n)),]
    
        r1 <- resample[1:n[1], , i]
        r2 <- resample[(n[1] + 1):(sum(n)), , i]
        m1[, i] <- apply(r1, FUN = mean, MARGIN = 2, rm.na = TRUE)
        m2[, i] <- apply(r2, FUN = mean, MARGIN = 2, rm.na = TRUE)

        m1[which(is.na(m1[,i])), i] <- 0
        m2[which(is.na(m2[,i])), i] <- 0

        tmp1 <- matrix(0, nrow = p, ncol = p)
        for(j in 1:n[1]){
            tt <- (r1[j, ] - m1[, i])%*%t(r1[j, ] - m1[, i])
            if(is.null(dim(tt))){
                tt <- matrix(0, nrow = p, ncol = p)
            }
            tt[which(is.na(tt))] <- 0
            tmp1 <- tmp1 + tt
        }

        tmp2 <- matrix(0, nrow = p, ncol = p)
        for(j in 1:n[2]){
            tt <- (r2[j, ] - m2[, i])%*%t(r2[j, ] - m2[, i])
            if(is.null(dim(tt))){
                tt <- matrix(0, nrow = p, ncol = p)
            }
            tt[which(is.na(tt))] <- 0
            tmp2 <- tmp2 + tt
        }
        
        S_r[, , i] <- 1/(sum(n) - 2)*(tmp1 + tmp2)

        cond <- TRUE
        k <- 1
        while(cond){
            D_res <- Dpart_lambda(Lambda[k], S = S_r[,,i], theta = NA, p, n, alpha = alpha)

            if(!is.na(D_orig[k]) & !is.na(D_res) & (abs(D_orig[k] - D_res) < delta)){
                cond <- FALSE
                lambda_B[i] <- Lambda[k]
            }else{
                k <- k + 1
            }
        }
    }

    lambda_sel <- median(lambda_B)

    # test statistic on the original sample
    RHT <- n[1]*n[2]/sum(n)*(mean1 - mean2)%*%U%*%solve(D + lambda_sel * diag(1, p))%*%t(U)%*%t(matrix(mean1 - mean2, nrow = 1))

    RHT_B <- rep(0, B)
    for(i in 1:B){
        r <- eigen(S_r[,,i])
        U <- r$vectors
        D <- diag(r$values)

        # test statistics on the resampled data
        RHT_B[i] <- n[1]*n[2]/sum(n)*(m1[,i] - m2[,i])%*%U%*%solve(D + lambda_sel * diag(1, p))%*%t(U)%*%t(matrix(m1[,i] - m2[,i], nrow = 1))
    }

    # p-value
    pval <- 1/B*length(which(RHT_B > as.numeric(RHT)))

    return(pval)
}
########################


######## simulaCHEN ########
## auxiliary function that do all the Monte Carlo repetitions of the Chen et al. approach
## on simulated data
#' @param dorep , number of Monte Carlo repetitions
#' @param n , vector with the number of observations for each group
#' @param p , number of components
#' @param Lambda ,  vector of regularization parameters in which we search for the optimal lambda
#' @param pmiss , probability \pi_{j1} of a missing data in the first group 
#' @param rho , covariance 
#' @param B , number of permutations used in the permutation test
#' @param delta , parameter used to check to condition on lambda. According to the authors it should be a small positive value
#' @param alpha , target significance level of the test
#' @param file.name , name used to save the files
#' @param mc.cores , number of cores used for the parallelization of the Monte Carlo repetitions
#' @param tySpecific , name of the simulation we are performing, contains information about n, p, pmiss, rho
#' @param ty , define if we are under the null hypothesis or the alternative hypothesis
#' @param c1c2 , values of c_1 and c_2
#' 
simulaCHEN <- function(dorep, n, p, Lambda = seq(0.01, 10, 0.05), pmiss = NA, rho = 0, B = 1000, delta = 0.5, alpha = 0.05, file.name = "", mc.cores = 1, tySpecific, ty, c1c2){

  res <- parallel::mclapply(X=dorep,
    FUN=function(i)  simula.internalCHEN(i, n = n, p = p, Lambda = Lambda, pmiss = pmiss, rho = rho, B = B, file.name = file.name, temp = TRUE, tySpecific = tySpecific, c1c2 = c1c2),  mc.cores = mc.cores)
  
  len <- 2 # length output of estimation 
  res <- matrix(0, nrow=length(dorep), ncol=len)
  k <- 0
  K <- length(n) # number of groups

  # merge together all the output of the Monte Carlo repetitions
  for (irep in dorep){
    print(irep)
    file.irep <- paste("./temp", K, "-", tySpecific, "/", irep, "-", file.name, sep="")
    if (file.exists(file.irep)) {
      k <- k + 1
      tmp <- data.matrix(read.table(file=file.irep, sep=";", header=FALSE, blank.lines.skip = FALSE))
      if(length(tmp) <= len){
        res[k,1:length(tmp)] <- tmp
      }else{
        res[k,] <- c(k, rep(NA, len - 1)) 
      }
    }
  }


  write.table(file=paste("./Results_CHEN/K", K, "/", ty, "/", file.name, sep = ""), x=res[1:k,], append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE,sep=";")
  system(paste("rm ./temp",K, "-",  tySpecific, "/*-",file.name,sep=""))
}
########################


#### simula.internalCHEN ####
## function that perform one single iteration of the Monte Carlo repetition: 
## it simulates the data according to c_1 and c_2 and then performs the 
## test described by Chen et al. through the "estimationCHEN" function
#' @param i , index of the Monte Carlo repetition
#' @param n , vector with the number of observations for each group
#' @param p , number of components
#' @param Lambda ,  vector of regularization parameters in which we search for the optimal lambda
#' @param pmiss , probability \pi_{j1} of a missing data in the first group 
#' @param rho , covariance 
#' @param B , number of permutations used in the permutation test
#' @param delta , parameter used to check to condition on lambda. According to the authors it should be a small positive value
#' @param alpha , target significance level of the test
#' @param file.name , name used to save the files
#' @param temp , put the file with the output of the test in the tmeporary folder
#' @param tySpecific , name of the simulation we are performing, contains information about n, p, pmiss, rho
#' @param c1c2 , values of c_1 and c_2
## it produces a file as output with the result of the test: 
## the file contains:
## - the index of the iteration (used to set the seed)
## - p-value of the permutation test,

simula.internalCHEN <- function(i, n, p, Lambda = seq(0.01, 10, 0.05), pmiss = NA, rho = 0, B = 1000, delta = 0.5, alpha = 0.05, file.name = "", temp = TRUE, tySpecific, c1c2) {

  set.seed(1235+i)

  # number of groups
  K <- length(n)

  # simulate the data
  if(c1c2[1] == 0 & c1c2[2] == 0){ # H0
    mu <- NULL
  }else if(c1c2[1] > 0 & c1c2[2] == 0){ # H1_1
    mu <- matrix(0, nrow=K, ncol=p)
    mu[1:K, ] <- mu[1:K, ] + matrix(c1c2[1]*seq(0, K-1)/(K-1), nrow =K, ncol = p)
  }else if(c1c2[1] == 0 & c1c2[2] > 0){ # H1_2
    mu <- NULL
    pmiss <- rep(pmiss, K) + matrix(c1c2[2]*seq(0, K-1)/(K-1), nrow =K, ncol = p)
  }
  dati <- scMANOVAsimulation(n = n, p = p, pmiss = pmiss, rho = rho, mu = mu) 

  # estimate the parameters and run the test
  res <- estimationCHEN(x = dati, n = n, Lambda = Lambda, B = B, delta = delta, alpha = alpha)
  res <- c(i, res)

  # write the result of the test in a file
  if (temp) {
    file.name <- paste("./temp", K, "-", tySpecific, "/",i, "-", file.name, sep="")
  }
  write.table(file = file.name, x=matrix(res, nrow=1), append=FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=";")
  return(i)
}
########################
