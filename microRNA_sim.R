###############################################################################
#### Resularized MANOVA test for semicontinuous high-dimensional data 
#### SIMULATED DATASETS similar to the real microRNA dataset
###############################################################################

method <- 'semicontinuousMANOVA' # Estimation procedure

## Packages
library("semicontMANOVA")
library("parallel")

## Working directory
setwd("./")

## Functions needed for the simulations
source("./functionsSim.R")
source("./functionsSim_real.R")


## Sample quantities obtained from the real dataset of microRNA differential expression: 
## we have sample mean, sample mean under H0, sample covariance matrix, sample covariance matrix under H0, probability of a missing value, probability of a missing values under H0
load("microRNA_sim.RData")


## Input missing values in the sample means and sample covariance matrices, as described in Section 3 of the manuscript
for(k in 1:length(n)){
  mu1[k, which(is.na(mu1[k,]))] <- mu0[which(is.na(mu1[k,]))]
}

sigma0[which(is.na(diag(sigma0))), which(is.na(diag(sigma0)))] <- 1
sigma0[which(is.na(sigma0))] <- 0
diag(sigma0) <- diag(sigma0) + 0.1
sigma1[which(is.na(diag(sigma1))), which(is.na(diag(sigma1)))] <- 1
sigma1[which(is.na(sigma1))] <- 0
diag(sigma1) <- diag(sigma1) + 0.1



## Number of permutations in the permutation test
B <- 1000

## Number of Monte Carlo simulations
nrep <- 1000

## Number of cores available for the parallelization of the different simulations
mc.cores <- 10L

## Number of cpus available for the parallelization of the permutation test
ncpus <- 50L

## Penalty to be used in the Information Criteria
penalty <- function(n,p) log(n) + 0.5*log(p)

## Use or not the identity function in the regularization of \Sigma and \Sigma_0
ident <- TRUE

## Minimum and Maximum lambda in which we look for the optimal lambda 
lambda <- c(0, 100)
## Minimum and Maximum lambda_0 in which we look for the optimal lambda_0
lambda0 <- c(0, 100)


## Create the Results folder
dir.create(paste("./Results", sep = ""), showWarnings = FALSE)


## Explored values of c_1 and c_2
typeC <- list()
typeC[["H0_microRNA"]] <- list("mu" = mu0, "sigma" = sigma0, "pmiss" = pmiss0)
typeC[["H1_microRNA"]] <- list("mu" = mu1,  "sigma" = sigma1,  "pmiss" = pmiss1)



## Run all the simulations
for(ty in names(typeC)){
  dir.create(paste("./Results/",ty, sep = ""), showWarnings = FALSE)
  ## Create a folder to store temporary output files
  dir.create(paste("./temp",ty, sep = ""), showWarnings = FALSE)

  cat(ty, " --- n: ", n, ", p:", p, "\n", sep = " ")

  file.name <- paste(method,"-", ty, "-", 2*n[1], "-", p,'.txt', sep='')
  ## Check if some output for this simulation setting 
  ## In case some result is already present, create the corresponding 
  ## temporary files 
  newdoRep <- 1:nrep

    
  if(file.exists(paste("./Results/", ty, "/", file.name, sep = ""))){
    
    old <- read.table(paste("./Results/",ty,"/", file.name, sep = ""), sep = ";", header = FALSE)
    newdoRep <- setdiff(1:nrep, old[,1])
    for(i in 1:dim(old)[1]){
      file.name.i <- paste("./temp", ty, "/", i, "-", file.name, sep="")
      write.table(file=file.name.i, x=matrix(old[i,], nrow=1), append=FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=";")
    }
  }

  
  dorep <- 1:nrep 
  len <- 12 + 1 + 5 # length of the output of "estimation" function
  res <- matrix(0, nrow=length(dorep), ncol=len)
  k <- 0

  ## If some temporary files are already present, store their result in res
  for (irep in dorep) {
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

  dorep <- newdoRep
  ndorep <- length(dorep)
  print(paste("ndorep ",ndorep))

  
  if (ndorep==1) {
    simula.internal_real(dorep, n, p, lambda = lambda, lambda0 = lambda0, rho, ncpus = 50L, B,penalty = penalty, ident = ident, file.name, temp=FALSE, ty = ty,  mu0 = mu0, sigma0 = sigma0, pmiss0 = pmiss0, mu1 = mu1, sigma1 = sigma1, pmiss1 = pmiss1)
  } else if (ndorep > 1) {
    simula_real(dorep, n, p, lambda = lambda, lambda0 = lambda0, rho, ncpus = 50L, B, penalty = penalty, ident = ident, file.name, mc.cores, ty = ty, mu0 = mu0, sigma0 = sigma0, pmiss0 = pmiss0, mu1 = mu1, sigma1 = sigma1, pmiss1 = pmiss1)
  }
}

system(paste("rm -r ./temp*", sep =""))



sessionInfo()



###############################################################################
#### create files with the summaries of the results of the test 
###############################################################################
## Working directory
setwd("./Results/")

## nominal size
alpha <- 0.05 

# Create files with the summary of the simulation results
for(ty in names(typeC)){
  k <- 1

  res <- matrix(NA, nrow = 1, ncol = 8)
  colnames(res) <- c( "n", "p", "rejPERM", "n_simEnded", "meanDim", "meanDim0", "meanLambda", "meanLambda0")
  
  tySpecific <- paste(ty, 2*n[1], p,sep = "-")
  cat(tySpecific, " --- n: ", n, ", p:", p, "\n", sep = " ")
  file.name <- paste(method,"-", tySpecific, '.txt', sep='')
  
  # read the file with the results of the specific scenario of n, p, pmiss, rho
  tmp <- data.matrix(read.table(file=paste(ty,"/", file.name, sep = ""), sep=";", header=FALSE, blank.lines.skip = FALSE))
  
  if(dim(tmp)[2] == 1){
    tmp <- t(tmp)
  }
  
  colnames(tmp) <- c("index", "time1", "time2", "time3", "time4", "time5", "logLik", "logLik0", "lambda", "lambda0", "df", "df0", "aic", "aic0", "statistic", "p.valuePERM", "dimSigma", "dimSigma0")
  rownames(tmp) <- tmp[,1]
  
  colPvaluePERM <- which(colnames(tmp) == "p.valuePERM")
  colmeanDim <- which(colnames(tmp) == "dimSigma")
  colmeanDim0 <- which(colnames(tmp) == "dimSigma0")
  colLambda <- which(colnames(tmp) == "lambda")
  colLambda0 <- which(colnames(tmp) == "lambda0")
  
  # Consider only the simulations that end (in some cases p^* = 0 and the simulation is discared)
  simEnded <- which(!is.na(tmp[, colPvaluePERM]))
  
  # Consider how many simulations has a p.value of the permutation test smaller than alpha
  rejPERM <- round(sum(tmp[simEnded, colPvaluePERM] < alpha)/length(simEnded), 3)
  
  meanDim <- round(mean(tmp[simEnded, colmeanDim]), 2)
  meanDim0 <- round(mean(tmp[simEnded, colmeanDim0]),2)
  meanLambda <- round(mean(tmp[simEnded, colLambda]), 2)
  meanLambda0 <- round(mean(tmp[simEnded, colLambda0]), 2)
  
  
  res[k, ] <- c( sum(n), p, rejPERM, length(simEnded), meanDim, meanDim0, meanLambda, meanLambda0)
  
  write.table(file=paste(ty, "_res", ".csv", sep = ""), x=res, append=FALSE, row.names=FALSE, col.names=TRUE, quote=FALSE,sep=";")
}





