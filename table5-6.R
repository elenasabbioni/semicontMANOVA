###############################################################################
#### produce Table 5 e 6
###############################################################################

method <- 'semicontinuousMANOVA' # Estimation procedure

## Packages
library(xtable)

## Create files with the summaries of the simulation study
source("analyzeResults_semicontinuousMANOVA.R")

# Current working directory
currWD <- getwd()


# Explored values of c_1 and c_2
typeC <- list()
typeC[["H0"]] <- c(0, 0)
typeC[["H1_1_1"]] <- c(1, 0)
typeC[["H1_1_5"]] <- c(5, 0)
typeC[["H1_2_0.15"]] <- c(0, 0.15)
typeC[["H1_2_0.3"]]  <- c(0, 0.3)


## Read the files with the summaries
for(K in c(2, 4)){
  ## Working directory
    setwd(paste0(currWD, "./Results_", method, "/K", K))

    H0 <- data.frame(read.table(file=paste("H0", "_res", ".csv", sep = ""),  sep=";", header = TRUE))
    H1_1_1 <- read.table(file=paste("H1_1-1", "_res", ".csv", sep = ""),  sep=";", header = TRUE)
    H1_1_5 <- read.table(file=paste("H1_1-5", "_res", ".csv", sep = ""),  sep=";", header = TRUE)
    H1_2_0.15 <- read.table(file=paste("H1_2-0.15", "_res", ".csv", sep = ""),  sep=";", header = TRUE)
    H1_2_0.3 <- read.table(file=paste("H1_2-0.3", "_res", ".csv", sep = ""),  sep=";", header = TRUE)

    df <- data.frame(typeH = rep(NA,28), c1 = NA, c2 = NA,
                    pmiss10.50  = NA, rho10.50 = NA, dim10.50 = NA, lambda10.50 = NA, lambda010.50 = NA, 
                    pmiss10.100 = NA, rho10.100 = NA, dim10.100 = NA, lambda10.100 = NA, lambda010.100 = NA,
                    pmiss10.150 = NA, rho10.150 = NA, dim10.150 = NA, lambda10.150 = NA, lambda010.150 = NA,
                    pmiss10.200 = NA, rho10.200 = NA, dim10.200 = NA, lambda10.200 = NA, lambda010.200 = NA,
                    pmiss20.50  = NA, rho20.50 = NA, dim20.50 = NA, lambda20.50 = NA, lambda020.50 = NA, 
                    pmiss20.100 = NA, rho20.100 = NA, dim20.100 = NA, lambda20.100 = NA, lambda020.100 = NA,
                    pmiss20.150 = NA, rho20.150 = NA, dim20.150 = NA, lambda20.150 = NA, lambda020.150 = NA,
                    pmiss20.200 = NA, rho20.200 = NA, dim20.200 = NA, lambda20.200 = NA, lambda020.200 = NA, rownames = FALSE)

    k1050 <- 1
    k10100 <- 1
    k10150 <- 1
    k10200 <- 1
    k2050 <- 1
    k20100 <- 1
    k20150 <- 1
    k20200 <- 1
    for(ty in names(typeC)){
      h <- get(ty)
      h <- h[, -c(which(colnames(h) %in% c("n_simEnded", "meanDim0")))]
      
      for(i in 1:dim(h)[1]){
        if(h$n[i] == 5*K & h$p[i] == 50){
          df$typeH[k1050] <- ty
          
          df$pmiss10.50[k1050] <- h$pmiss[i]
          df$rho10.50[k1050] <- h$rho[i]
          df$dim10.50[k1050] <- h$meanDim[i]
          df$lambda10.50[k1050] <- h$meanLambda[i]
          df$lambda010.50[k1050] <- h$meanLambda0[i]
          
          k1050 <- k1050 + 1
        }else if(h$n[i] == 5*K & h$p[i] == 100){
          df$typeH[k10100] <- ty
          
          df$pmiss10.100[k10100] <- h$pmiss[i]
          df$rho10.100[k10100] <- h$rho[i]
          df$dim10.100[k10100] <- h$meanDim[i]
          df$lambda10.100[k10100] <- h$meanLambda[i]
          df$lambda010.100[k10100] <- h$meanLambda0[i]
          
          k10100 <- k10100 + 1
        }else if(h$n[i] == 5*K & h$p[i] == 150){
          df$typeH[k10150] <- ty
          
          df$pmiss10.150[k10150] <- h$pmiss[i]
          df$rho10.150[k10150] <- h$rho[i]
          df$dim10.150[k10150] <- h$meanDim[i]
          df$lambda10.150[k10150] <- h$meanLambda[i]
          df$lambda010.150[k10150] <- h$meanLambda0[i]
          
          k10150 <- k10150 + 1
        }else if(h$n[i] == 5*K & h$p[i] == 200){
          df$typeH[k10200] <- ty
          
          df$pmiss10.200[k10200] <- h$pmiss[i]
          df$rho10.200[k10200] <- h$rho[i]
          df$dim10.200[k10200] <- h$meanDim[i]
          df$lambda10.200[k10200] <- h$meanLambda[i]
          df$lambda010.200[k10200] <- h$meanLambda0[i]
          
          k10200 <- k10200 + 1
        }else if(h$n[i] == 10*K & h$p[i] == 50){
          df$typeH[k2050] <- ty
          
          df$pmiss20.50[k2050] <- h$pmiss[i]
          df$rho20.50[k2050] <- h$rho[i]
          df$dim20.50[k2050] <- h$meanDim[i]
          df$lambda20.50[k2050] <- h$meanLambda[i]
          df$lambda020.50[k2050] <- h$meanLambda0[i]
          
          k2050 <- k2050 + 1
        }else if(h$n[i] == 10*K & h$p[i] == 100){
          df$typeH[k20100] <- ty
          
          df$pmiss20.100[k20100] <- h$pmiss[i]
          df$rho20.100[k20100] <- h$rho[i]
          df$dim20.100[k20100] <- h$meanDim[i]
          df$lambda20.100[k20100] <- h$meanLambda[i]
          df$lambda020.100[k20100] <- h$meanLambda0[i]
          
          k20100 <- k20100 + 1
        }else if(h$n[i] == 10*K & h$p[i] == 150){
          df$typeH[k20150] <- ty
          
          df$pmiss20.150[k20150] <- h$pmiss[i]
          df$rho20.150[k20150] <- h$rho[i]
          df$dim20.150[k20150] <- h$meanDim[i]
          df$lambda20.150[k20150] <- h$meanLambda[i]
          df$lambda020.150[k20150] <- h$meanLambda0[i]
          
          k20150 <- k20150 + 1
        }else if(h$n[i] == 10*K & h$p[i] == 200){
          df$typeH[k20200] <- ty
          
          df$pmiss20.200[k20200] <- h$pmiss[i]
          df$rho20.200[k20200] <- h$rho[i]
          df$dim20.200[k20200] <- h$meanDim[i]
          df$lambda20.200[k20200] <- h$meanLambda[i]
          df$lambda020.200[k20200] <- h$meanLambda0[i]
          
          k20200 <- k20200 + 1
        }
      }
      
      df$c1[which(df$typeH == ty)] <- round(typeC[[ty]][1])
      df$c2[which(df$typeH == ty)] <- typeC[[ty]][2]
    }




    df <- df[, -1]


    df$lambda10.50 <- paste(df$lambda10.50, df$lambda010.50, sep = "-")
    df$lambda10.100 <- paste(df$lambda10.100, df$lambda010.100, sep = "-")
    df$lambda10.150 <- paste(df$lambda10.150, df$lambda010.150, sep = "-")
    df$lambda10.200 <- paste(df$lambda10.200, df$lambda010.200, sep = "-")
    df$lambda20.50 <- paste(df$lambda20.50, df$lambda020.50, sep = "-")
    df$lambda20.100 <- paste(df$lambda20.100, df$lambda020.100, sep = "-")
    df$lambda20.150 <- paste(df$lambda20.150, df$lambda020.150, sep = "-")
    df$lambda20.200 <- paste(df$lambda20.200, df$lambda020.200, sep = "-")

    df$lambda10.50 <- paste(df$dim10.50, df$lambda10.50, sep = "; ")
    df$lambda10.100 <- paste(df$dim10.100, df$lambda10.100, sep = "; ")
    df$lambda10.150 <- paste(df$dim10.150, df$lambda10.150, sep = "; ")
    df$lambda10.200 <- paste(df$dim10.200, df$lambda10.200, sep = "; ")
    df$lambda20.50 <- paste(df$dim20.50, df$lambda20.50, sep = "; ")
    df$lambda20.100 <- paste(df$dim20.100, df$lambda20.100, sep = "; ")
    df$lambda20.150 <- paste(df$dim20.150, df$lambda20.150, sep = "; ")
    df$lambda20.200 <- paste(df$dim20.200, df$lambda20.200, sep = "; ")

    df$lambda10.50 <- paste("(", df$lambda10.50, ")", sep = "")
    df$lambda10.100 <- paste("(", df$lambda10.100, ")", sep = "")
    df$lambda10.150 <- paste("(", df$lambda10.150, ")", sep = "")
    df$lambda10.200 <- paste("(", df$lambda10.200, ")", sep = "")
    df$lambda20.50 <- paste("(", df$lambda20.50, ")", sep = "")
    df$lambda20.100 <- paste("(", df$lambda20.100, ")", sep = "")
    df$lambda20.150 <- paste("(", df$lambda20.150, ")", sep = "")
    df$lambda20.200 <- paste("(", df$lambda20.200, ")", sep = "")


    if(K == 2){
      table5 <- df[, c("c1", "c2", "pmiss10.50", "rho10.50", "lambda10.50", "lambda10.100", "lambda10.150", "lambda10.200", "lambda20.50", "lambda20.100", "lambda20.150",  "lambda20.200")]
      table5 <- rbind(c(rep("", 4), rep("n = 10", 4), rep("n = 20", 4)), table5)
      colnames(table5) <- c("$c_1$", "$c_2$",  "$\\pi_{j1}$", "$\\rho$",rep(c("$p = 50$", "$p = 100$", "$p = 150$", "$p = 200$"), 2))
    }else{
      table6 <- df[, c("c1", "c2", "pmiss10.50", "rho10.50", "lambda10.50", "lambda10.100", "lambda10.150", "lambda10.200", "lambda20.50", "lambda20.100", "lambda20.150",  "lambda20.200")]
      table6 <- rbind(c(rep("", 4), rep("n = 10", 4), rep("n = 20", 4)), table6)
      colnames(table6) <- c("$c_1$", "$c_2$",  "$\\pi_{j1}$", "$\\rho$",rep(c("$p = 50$", "$p = 100$", "$p = 150$", "$p = 200$"), 2))
    }
}


# Table 5
print(xtable(table5, digits=c(0, 0, 2, 1, 1, NA, NA, NA, NA, NA, NA, NA, NA)), include.rownames = FALSE)

# Table 6
print(xtable(table6, digits=c(0, 0, 2, 1, 1, NA, NA, NA, NA, NA, NA, NA, NA)), include.rownames = FALSE)
