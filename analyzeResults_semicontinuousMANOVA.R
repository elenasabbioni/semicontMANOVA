#####################################################################################
#### create files with the summaries of the results of the regularized MANOVA test 
#### for the simulation study
####################################################################################
method <- 'semicontinuousMANOVA' # Estimation procedure

## Working directory
currWD <- getwd()

# Explored values of c_1 and c_2
typeC <- list()
typeC[["H0"]] <- c(0, 0)
typeC[["H1_1-1"]] <- c(1, 0)
typeC[["H1_1-5"]] <- c(5, 0)
typeC[["H1_2-0.15"]] <- c(0, 0.15)
typeC[["H1_2-0.3"]]  <- c(0, 0.3)


# Simulated numbers of components
allp <- c(50, 100, 150, 200)

# Simulated probabilities of a missing data in the first group \pi_j1
allpmiss <- c(0.2, 0.5, 0.8)

# Simulated covariances between two components
allrho <- c(0, 0.4)

# nominal size
alpha <- 0.05 


# different number of groups
for(K in c(2, 4)){
  setwd(paste0(currWD, "/Results_", method, "/K", K))

  # Simulated sample sizes 
  n_sizes <- c(5, 10)*K # total sample size in the two groups 
  alln <- list()
  for(n in n_sizes){
    alln[[as.character(n)]] <- rep(n/K, K) # sample size in each group
  } 

  # Create files with the summary of the simulation results
  for(ty in names(typeC)){
    res <- matrix(NA, nrow = length(alln)*length(allp)*length(allpmiss)*length(allrho), ncol = 10)
    colnames(res) <-  c( "n", "p", "pmiss", "rho", "rejPERM", "n_simEnded", "meanDim", "meanDim0", "meanLambda", "meanLambda0")

    k <- 0
    for (n in alln){
      for (p in allp){
        for (pmiss in allpmiss){
          if( (typeC[[ty]][2] + pmiss) < 1){ # check if we have an admissible probability \pi_j2
            for (rho in allrho){
              k <- k + 1
              tySpecific <- paste(ty, K*n[1], p, pmiss, rho, sep = "-")
              cat(tySpecific, " --- n: ", n, ", p:", p, ", pmiss: ", pmiss, ", rho: ", rho, "\n", sep = " ")
              file.name <- paste(method,"-", tySpecific, '.txt', sep='')
              
              # read the file with the results of the specific scenario of n, p, pmiss, rho
              tmp <- data.matrix(read.table(file=paste(ty,"/", file.name, sep = ""), sep=";", header=FALSE, blank.lines.skip = FALSE))
              
              if(dim(tmp)[2] == 1){
                tmp <- t(tmp)
              }
              
              colnames(tmp) <-  c("index", "time1", "time2", "time3", "time4", "time5", "logLik", "logLik0", "lambda", "lambda0", "df", "df0", "aic", "aic0", "statistic", "p.valuePERM", "dimSigma", "dimSigma0")
              rownames(tmp) <- tmp[,1]
              
              colPvaluePERM <- which(colnames(tmp) == "p.valuePERM")          
              colmeanDim <- which(colnames(tmp) == "dimSigma")
              colmeanDim0 <- which(colnames(tmp) == "dimSigma0")
              colLambda <- which(colnames(tmp) == "lambda")
              colLambda0 <- which(colnames(tmp) == "lambda0")

              # Consider only the simulations that end (in some cases p^* = 0 and the simulation is discared)
              simEnded <- which(!is.na(tmp[, colPvaluePERM]))

              meanDim <- round(mean(tmp[simEnded, colmeanDim]), 2)
              meanDim0 <- round(mean(tmp[simEnded, colmeanDim0]),2)
              meanLambda <- round(mean(tmp[simEnded, colLambda]), 2)
              meanLambda0 <- round(mean(tmp[simEnded, colLambda0]), 2)  
              
              
              # Consider how many simulations has a p.value of the permutation test smaller than alpha
              rejPERM <- round(sum(tmp[simEnded, colPvaluePERM] < alpha)/length(simEnded), 3)
              

              res[k, ] <- c( sum(n), p, pmiss, rho, rejPERM, length(simEnded), meanDim, meanDim0, meanLambda, meanLambda0)
            }
          }
        }
      }
    }
    
    res <- res[which(!is.na(res[,1])), ]
    
    write.table(file=paste(ty, "_res", ".csv", sep = ""), x=res, append=FALSE, row.names=FALSE, col.names=TRUE, quote=FALSE,sep=";")
  }
}

setwd(currWD)
