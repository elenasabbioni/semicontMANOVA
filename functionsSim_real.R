
######## SIMULA_REAL ########
## aulxiliary function that do all the Monte Carlo repetitions of the MANOVA test
## on data simulated starting from sample parameters derived from real data
#' @param dorep , number of Monte Carlo repetitions
#' @param n , vector with the number of observations for each group
#' @param lambda ,  minimum and maximum values between we search for the optimal lambda
#' @param lambda0 , minimum and maximum values between we search for the optimal lambda0
#' @param rho , covariance 
#' @param ncpus , number of cpus used for the parallelization of the permutation test
#' @param B , number of permutations used in the permutation test
#' @param penalty , penalty used in the Information criteria
#' @param ident , use or not the identity function in the regularization
#' @param file.name , name used to save the files
#' @param mc.cores , number of cores used for the parallelization of the Monte Carlo repetitions
#' @param ty , define if we are under the null hypothesis or the alternative hypothesis
#' @param mu0 , mean parameter under H0
#' @param sigma0, variance and covariance matrix under H0
#' @param pmiss0, probability of a missing value under H0
#' @param mu1, mean parameter under H1
#' @param sigma1, variance and covariance matrix under H1
#' @param pmiss1, probability of a missing value under H1

simula_real <- function(dorep, n, p, lambda = c(0, 100), lambda0 = c(0, 100), rho, ncpus = 50L, B, penalty, ident, file.name, mc.cores, ty, mu0 = NA, sigma0 = NA, pmiss0 = NA, mu1 = NA, sigma1 = NA, pmiss1 = NA) {
    # do all the Monte Carlo repetition of the test
    res <- parallel::mclapply(X=dorep,
        FUN=function(i)  simula.internal_real(i, n = n, p = p, lambda = lambda, lambda0 = lambda0, rho = rho, B = B, ncpus = ncpus, penalty = penalty, ident = ident, file.name = file.name, temp = TRUE, ty = ty, mu0 = mu0, sigma0 = sigma0, pmiss0 = pmiss0, mu1 = mu1, sigma1 = sigma1, pmiss1 = pmiss1),  mc.cores = mc.cores)
    
    len <- 12 + 1 + 5 # length output: estimation's output + i + time
    res <- matrix(0, nrow=length(dorep), ncol=len)
    k <- 0
    # merge together all the output of the Monte Carlo repetitions
    for (irep in dorep){
        file.irep <- paste("./temp", ty, "/", irep, "-", file.name, sep="")
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
    write.table(file=paste("./Results/", ty, "/", file.name, sep = ""), x=res[1:k,], append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE,sep=";")
}
#########################




#### SIMULA.INTERNAL_REAL ####
## function that perform one single iteration of the Monte Carlo repetition: 
## it simulates the data according to c_1 and c_2 and then performs the 
## regularized MANOVA test through the "estimation" function
#' @param i , index of the Monte Carlo repetition
#' @param n , vector with the number of observations for each group
#' @param p , number of components
#' @param lambda ,  minimum and maximum values between we search for the optimal lambda
#' @param lambda0 , minimum and maximum values between we search for the optimal lambda0
#' @param rho , covariance 
#' @param ncpus , number of cpus used for the parallelization of the permutation test
#' @param B , number of permutations used in the permutation test
#' @param penalty , penalty used in the Information criteria
#' @param ident , use or not the identity function in the regularization
#' @param file.name , name used to save the files
#' @param temp , put the file with the output of the test in the tmeporary folder
#' @param ty , define if we are under the null hypothesis or the alternative hypothesis
#' @param mu0 , mean parameter under H0
#' @param sigma0, variance and covariance matrix under H0
#' @param pmiss0, probability of a missing value under H0
#' @param mu1, mean parameter under H1
#' @param sigma1, variance and covariance matrix under H1
#' @param pmiss1, probability of a missing value under H1
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
simula.internal_real <- function(i, n, p, lambda = c(0, 100), lambda0 = c(0, 100), rho, ncpus, B, penalty, ident, file.name, temp = TRUE, ty, mu0 = NA, sigma0 = NA, pmiss0 = NA, mu1 = NA, sigma1 = NA, pmiss1 = NA) {
    # set the seed
    set.seed(1235+i)

    # number of groups
    K <- length(n)

    # simulate the data
    if(grepl("H0", ty)){
    mu <- mu0
    sigma <- sigma0
    pmiss <- rep(pmiss0, K)
    }else if(grepl("H1", ty)){
    mu <- mu1
    sigma <- sigma1
    pmiss <- pmiss1
    }

    dati <- scMANOVAsimulation(n = n, p = p, pmiss = pmiss, sigma = sigma, mu= mu) 

    # estimate the parameters and run the test
    time <- system.time({
    res <- estimation(x = dati, n=n, lambda = lambda, lambda0 = lambda0, ncpus = ncpus, B = B, penalty = penalty, ident = ident)
    })
    res <- c(i, time, res)

    # write the result of the test in a file
    if (temp) {
    file.name <- paste("./temp", ty, "/",i, "-", file.name, sep="")
    }
    write.table(file=file.name, x=matrix(res, nrow=1), append=FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=";")
    return(i)
}
#########################