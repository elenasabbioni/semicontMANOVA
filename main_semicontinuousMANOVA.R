###############################################################################
#### Regularized MANOVA test for semicontinuous high-dimensional data 
#### SIMULATION STUDY
###############################################################################

method <- 'semicontinuousMANOVA' # Estimation procedure

## Packages
library("semicontMANOVA")
library("parallel")

## Working directory
currWD <- getwd()

## Functions needed for the simulations
source("./functionsSim.R")

## Explored values of c_1 and c_2
typeC <- list()
typeC[["H0"]] <- c(0, 0)
typeC[["H1_1-1"]] <- c(1, 0)
typeC[["H1_1-5"]] <- c(5, 0)
typeC[["H1_2-0.15"]] <- c(0, 0.15)
typeC[["H1_2-0.3"]]  <- c(0, 0.3)


## Simulated numbers of components
allp <- c(50, 100, 150, 200)

## Simulated probabilities of a missing data in the first group \pi_j1
allpmiss <- c(0.2, 0.5, 0.8)

## Simulated covariances between two components
allrho <- c(0, 0.4)

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
dir.create(paste0("./Results_", method), showWarnings = FALSE)




# different number of groups
for(K in c(2, 4)){
  dir.create(paste0("./Results_", method, "/K", K), showWarnings = FALSE)

  # Simulated sample sizes 
  n_sizes <- c(5, 10)*K # total sample size in the two groups 
  alln <- list()
  for(n in n_sizes){
    alln[[as.character(n)]] <- rep(n/K, K) # sample size in each group
  } 


  ## Run all the simulations for the different combination of n, p, pmiss, rho
  for(ty in names(typeC)){
    dir.create(paste("./Results_", method, "/K", K, "/", ty, sep = ""), showWarnings = FALSE)
    
    for (n in alln){
      for (p in allp){
        for (pmiss in allpmiss){
          if( (typeC[[ty]][2] + pmiss) < 1){ # check if we have an admissible probability \pi_j2
            for (rho in allrho){
              tySpecific <- paste(ty, K*n[1], p, pmiss, rho, sep = "-")
              
              ## Create a folder to store temporary output files
              dir.create(paste("./temp", K, "-",tySpecific, sep = ""), showWarnings = FALSE)
              
              cat(tySpecific, " --- n: ", n, ", p:", p, ", pmiss: ", pmiss, ", rho: ", rho, "\n", sep = " ")
              file.name <- paste(method,"-", ty, "-", K*n[1], "-", p, "-", pmiss, '-', rho, '.txt', sep='')
              
              ## Check if some output for this simulation setting 
              ## In case some result is already present, create the corresponding 
              ## temporary files 
              newdoRep <- 1:nrep
              if(file.exists(paste("./Results_", method, "/K", K, "/", ty, "/", file.name, sep = ""))){
                old <- read.table(paste("./Results_", method, "/K", K, "/", ty, "/", file.name, sep = ""), sep = ";", header = FALSE)
                newdoRep <- setdiff(1:nrep, old[,1])
                for(i in 1:dim(old)[1]){
                  file.name.i <- paste("./temp", K, "-", tySpecific, "/", i, "-", file.name, sep="")
                  write.table(file=file.name.i, x=matrix(old[i,], nrow=1), append=FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=";")
                }
              }
              
              dorep <- 1:nrep 
              len <- 12 + 1 + 5 # length of the output of "estimation" function
              res <- matrix(0, nrow=length(dorep), ncol=len)
              k <- 0
              
              ## If some temporary files are already present, store their result in res
              for (irep in dorep) {
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
    
              
              dorep <- newdoRep
              ndorep <- length(dorep)
              print(paste("ndorep ",ndorep))
              if (ndorep==1) {
                simula.internal(dorep, n, p, lambda = lambda, lambda0 = lambda0, pmiss = pmiss, rho = rho, ncpus = ncpus, B = B, penalty = penalty, ident = ident, file.name = file.name, temp = FALSE, tySpecific = tySpecific, c1c2 = typeC[[ty]])
              } else if (ndorep > 1) {
                simula(dorep, n, p, lambda = lambda, lambda0 = lambda0, pmiss = pmiss, rho = rho, ncpus = ncpus, B = B, penalty = penalty, ident = ident, file.name = file.name, mc.cores = mc.cores, tySpecific = tySpecific, ty = ty, c1c2 = typeC[[ty]])
              }
            }
          }
        }
      }
    }
    system(paste("rm -r ./temp", K, "-*", sep =""))
  }
}


sessionInfo()
