#####################################################################################
#### create files with the summaries of the results of the test of CHEN et al.
#### for the simulation study
####################################################################################
K <- 2 # number of groups

method <- 'CHEN' # Estimation procedure



## Working directory
currWD <- getwd()
setwd(paste0("./Results_", method, "/K", K))


# Explored values of c_1 and c_2
typeC <- list()
typeC[["H0"]] <- c(0, 0)
typeC[["H1_1-1"]] <- c(1, 0)
typeC[["H1_1-5"]] <- c(5, 0)
typeC[["H1_2-0.15"]] <- c(0, 0.15)
typeC[["H1_2-0.3"]]  <- c(0, 0.3)

# Simulated sample sizes 
n_sizes <- c(5, 10)*K # total sample size in the two groups 
alln <- list()
for(n in n_sizes){
  alln[[as.character(n)]] <- rep(n/K, K) # sample size in each group
} 

# Simulated numbers of components
allp <- c(50, 100, 150, 200)

# Simulated probabilities of a missing data in the first group \pi_j1
allpmiss <- c(0.2, 0.5, 0.8)

# Simulated covariances between two components
allrho <- c(0, 0.4)

# nominal size
alpha <- 0.05 

# Create files with the summary of the simulation results
for(ty in names(typeC)){
  res <- matrix(NA, nrow = length(alln)*length(allp)*length(allpmiss)*length(allrho), ncol = 6)
  colnames(res) <- c( "n", "p", "pmiss", "rho", "rejPERM", "n_simEnded")
  
  k <- 0
  for (n in alln){
    for (p in allp){
      for (pmiss in allpmiss){
        if( (typeC[[ty]][2] + pmiss) < 1){ # check if we have an admissible probability \pi_j2
          for (rho in allrho){
            k <- k + 1
            tySpecific <- paste(ty, 2*n[1], p, pmiss, rho, sep = "-")
            cat(tySpecific, " --- n: ", n, ", p:", p, ", pmiss: ", pmiss, ", rho: ", rho, "\n", sep = " ")
            file.name <- paste(method,"-", tySpecific, '.txt', sep='')
            
            # read the file with the results of the specific scenario of n, p, pmiss, rho
            tmp <- data.matrix(read.table(file=paste("./", ty,"/", file.name, sep = ""), sep=";", header=FALSE, blank.lines.skip = FALSE))
            
            if(dim(tmp)[2] == 1){
              tmp <- t(tmp)
            }
            
            colnames(tmp) <- c("index", "p.valuePERM")
            rownames(tmp) <- tmp[,1]
            
            colPvaluePERM <- which(colnames(tmp) == "p.valuePERM")
            
            # Consider only the simulations that end 
            simEnded <- which(!is.na(tmp[, colPvaluePERM]))
            
            # Consider how many simulations has a p.value of the permutation test smaller than alpha
            rejPERM <- round(sum(tmp[simEnded, colPvaluePERM] < alpha)/length(simEnded), 3)
            
            res[k, ] <- c( sum(n), p, pmiss, rho, rejPERM, length(simEnded))
          }
        }
      }
    }
  }
  
  res <- res[which(!is.na(res[,1])), ]
  
  write.table(file=paste(ty, "_res", ".csv", sep = ""), x=res, append=FALSE, row.names=FALSE, col.names=TRUE, quote=FALSE,sep=";")
}


setwd(currWD)
