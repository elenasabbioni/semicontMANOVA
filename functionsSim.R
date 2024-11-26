###############################################################################
#### functions needed to run a regularized MANOVA test on a simulated scenario
#### with semicontinuous high-dimensional data
###############################################################################


###### ESTIMATION ######
## function that returns the estimated parameters of the model 
## and run the regularized MANOVA test
#' @param x , n x p matrix with the data
#' @param n , vector with the number of observations for each group
#' @param lambda ,  minimum and maximum values between we search for the optimal lambda
#' @param lambda0 , minimum and maximum values between we search for the optimal lambda0
#' @param ncpus , number of cpus used for the parallelization of the permutation test
#' @param B , number of permutations used in the permutation test
#' @param penalty , penalty used in the Information criteria
#' @param ident , use or not the identity function in the regularization
# It outputs the estimated parameters of the test and the p-value of the permutation test

estimation <- function(x, n, lambda = c(0, 100), lambda0 = c(0, 100), ncpus = 1L, B = 1000, penalty = function(n,p) log(n) + 0.5*log(p), ident = TRUE) {
  # estimate the parameters and run the test
  res <- try(scMANOVA(x = x, n = n, lambda = lambda, lambda0 = lambda0, tol = 1e-8, p.value.perm = TRUE, fixed.lambda = FALSE, parallel = 'multicore', ncpus = ncpus, B = B, penalty = penalty, ident = ident))
  if (class(res)[1] =='try-error') {
    print('Try-error')
    rep(NA, 12)
  } else { 
    sigmaRidge <- res$sigmaRidge
    sigma0Ridge <- res$sigma0Ridge
    logLik <- res$logLik
    logLik0 <- res$logLik0
    lambda <- res$lambda
    lambda0 <- res$lambda0
    df <- res$df
    df0 <- res$df0
    aic <- res$aic
    aic0 <- res$aic0
    statistic <- res$statistic
    p.valuePERM <- res$p.value
    
    c(logLik, logLik0, lambda, lambda0, df, df0, aic, aic0, statistic, p.valuePERM, dim(sigmaRidge)[1], dim(sigma0Ridge)[1])
  }
}
########################


######## SIMULA ########
## auxiliary function that do all the Monte Carlo repetitions of the MANOVA test
## on simulated data
#' @param dorep , number of Monte Carlo repetitions
#' @param n , vector with the number of observations for each group
#' @param p , number of components
#' @param lambda ,  minimum and maximum values between we search for the optimal lambda
#' @param lambda0 , minimum and maximum values between we search for the optimal lambda0
#' @param pmiss , probability \pi_{j1} of a missing data in the first group 
#' @param rho , covariance 
#' @param ncpus , number of cpus used for the parallelization of the permutation test
#' @param B , number of permutations used in the permutation test
#' @param penalty , penalty used in the Information criteria
#' @param ident , use or not the identity function in the regularization
#' @param file.name , name used to save the files
#' @param mc.cores , number of cores used for the parallelization of the Monte Carlo repetitions
#' @param tySpecific , name of the simulation we are performing, contains information about n, p, pmiss, rho
#' @param ty , define if we are under the null hypothesis or the alternative hypothesis
#' @param c1c2 , values of c_1 and c_2

simula <- function(dorep, n, p, lambda = c(0, 100), lambda0 = c(0, 100), pmiss = NA, rho = 0, ncpus = 1L, B = 1000, penalty = function(n,p) log(n) + 0.5*log(p), ident = TRUE, file.name = "", mc.cores = 1, tySpecific, ty, c1c2) {
  # do all the Monte Carlo repetition of the test
  res <- parallel::mclapply(X=dorep,
                            FUN=function(i)  simula.internal(i, n = n, p = p, lambda = lambda, lambda0 = lambda0, pmiss = pmiss, rho = rho, ncpus = ncpus, B = B, penalty = penalty, ident = ident, file.name = file.name, temp = TRUE, tySpecific = tySpecific, c1c2 = c1c2),  mc.cores = mc.cores)
  
  len <- 12 + 1 + 5 # length output: estimation's output + i + time
  res <- matrix(0, nrow=length(dorep), ncol=len)
  k <- 0
  K <- length(n) # number of groups
  # merge together all the output of the Monte Carlo repetitions
  for (irep in dorep){
    file.irep <- paste("./temp",  K, "-", tySpecific, "/", irep, "-", file.name, sep="")
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
  write.table(file=paste("./Results_semicontinuousMANOVA/K", K, "/", ty, "/", file.name, sep = ""), x=res[1:k,], append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE,sep=";")
  system(paste("rm ./temp", K, "-",  tySpecific, "/*-",file.name,sep=""))  
}
#########################


#### SIMULA.INTERNAL ####
## function that perform one single iteration of the Monte Carlo repetition: 
## it simulates the data according to c_1 and c_2 and then performs the 
## regularized MANOVA test through the "estimation" function
#' @param i , index of the Monte Carlo repetition
#' @param n , vector with the number of observations for each group
#' @param p , number of components
#' @param lambda ,  minimum and maximum values between we search for the optimal lambda
#' @param lambda0 , minimum and maximum values between we search for the optimal lambda0
#' @param pmiss , probability \pi_{j1} of a missing data in the first group 
#' @param rho , covariance 
#' @param ncpus , number of cpus used for the parallelization of the permutation test
#' @param B , number of permutations used in the permutation test
#' @param penalty , penalty used in the Information criteria
#' @param ident , use or not the identity function in the regularization
#' @param file.name , name used to save the files
#' @param temp , put the file with the output of the test in the tmeporary folder
#' @param tySpecific , name of the simulation we are performing, contains information about n, p, pmiss, rho
#' @param c1c2 , values of c_1 and c_2
## it produces a file as output with the result of the test: 
## the file contains:
## - the index of the iteration (used to set the seed)
## - the total user and system CPU time of the current process and the real elapsed time since the process was started (5 entrances) 
## - results of the test: loglikelihood under no hypothesis, loglikelihood under the null hypothesis, 
##   selected value of lambda, selected value of lambda0, model complexity measure under no hypothesis, 
##   model complexity measure under H0, Information criteria under no hypothesis, 
##   Information criteria under the null hypothesis, test statistic, p-value of the permutation test,
##   number of component of the final estimator under no hypothesis, 
##   number of component of the final estimator under H0 (equal to the number of components under no hypothesis)
simula.internal <- function(i, n, p, lambda = c(0, 100), lambda0 = c(0, 100), pmiss = NA, rho = 0, ncpus = 1, B = 1000, penalty = function(n,p) log(n) + 0.5*log(p), ident = TRUE, file.name = "", temp = TRUE, tySpecific, c1c2) {
  # set the seed
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
  dati <- scMANOVAsimulation(n = n, p = p, pmiss = pmiss, rho = rho, mu= mu) 
  
  # estimate the parameters and run the test
  time <- system.time({
    res <- estimation(x = dati, n = n, lambda = lambda, lambda0 = lambda0, ncpus = ncpus, B = B, penalty = penalty, ident = ident)
  })
  res <- c(i, time, res)
  
  # write the result of the test in a file
  if (temp) {
    file.name <- paste("./temp", K, "-", tySpecific, "/",i, "-", file.name, sep="")
  }
  write.table(file=file.name, x=matrix(res, nrow=1), append=FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=";")
  return(i)
}
#########################
