###############################################################################
#### produce Table4
###############################################################################
K <- 4 # number of groups


method <- 'semicontinuousMANOVA' # Estimation procedure

## Packages
library(xtable)

## Working directory
setwd(paste0("./Results_", method, "/K", K))

## Create files with the summaries of the simulation study
source("analyzeResults.R")


# Explored values of c_1 and c_2
typeC <- list()
typeC[["H0"]] <- c(0, 0)
typeC[["H1_1_1"]] <- c(1, 0)
typeC[["H1_1_5"]] <- c(5, 0)
typeC[["H1_2_0.15"]] <- c(0, 0.15)
typeC[["H1_2_0.3"]]  <- c(0, 0.3)


## Read the files with the summaries
H0 <- data.frame(read.table(file=paste("H0", "_res", ".csv", sep = ""),  sep=";", header = TRUE))
H1_1_1 <- read.table(file=paste("H1_1-1", "_res", ".csv", sep = ""),  sep=";", header = TRUE)
H1_1_5 <- read.table(file=paste("H1_1-5", "_res", ".csv", sep = ""),  sep=";", header = TRUE)
H1_2_0.15 <- read.table(file=paste("H1_2-0.15", "_res", ".csv", sep = ""),  sep=";", header = TRUE)
H1_2_0.3 <- read.table(file=paste("H1_2-0.3", "_res", ".csv", sep = ""),  sep=";", header = TRUE)



df <- data.frame(typeH = rep(NA,28), c1 = NA, c2 = NA,
pmiss10.50  = NA, rho10.50 = NA, propRej10.50 = NA, 
pmiss10.100 = NA, rho10.100 = NA, propRej10.100 = NA,
pmiss10.150 = NA, rho10.150 = NA, propRej10.150 = NA,
pmiss10.200 = NA, rho10.200 = NA, propRej10.200 = NA,
pmiss20.50  = NA, rho20.50 = NA, propRej20.50 = NA,
pmiss20.100 = NA, rho20.100 = NA, propRej20.100 = NA,
pmiss20.150 = NA, rho20.150 = NA, propRej20.150 = NA,
pmiss20.200 = NA, rho20.200 = NA, propRej20.200 = NA, rownames = FALSE)
 
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
            df$propRej10.50[k1050] <- h$rejPERM[i]
            

            k1050 <- k1050 + 1
        }else if(h$n[i] == 5*K & h$p[i] == 100){
            df$typeH[k10100] <- ty

            df$pmiss10.100[k10100] <- h$pmiss[i]
            df$rho10.100[k10100] <- h$rho[i]
            df$propRej10.100[k10100] <- h$rejPERM[i]
           

            k10100 <- k10100 + 1
        }else if(h$n[i] == 5*K & h$p[i] == 150){
            df$typeH[k10150] <- ty

            df$pmiss10.150[k10150] <- h$pmiss[i]
            df$rho10.150[k10150] <- h$rho[i]
            df$propRej10.150[k10150] <- h$rejPERM[i]
            

            k10150 <- k10150 + 1
        }else if(h$n[i] == 5*K & h$p[i] == 200){
            df$typeH[k10200] <- ty

            df$pmiss10.200[k10200] <- h$pmiss[i]
            df$rho10.200[k10200] <- h$rho[i]
            df$propRej10.200[k10200] <- h$rejPERM[i]
           

            k10200 <- k10200 + 1
        }else if(h$n[i] == 10*K & h$p[i] == 50){
            df$typeH[k2050] <- ty

            df$pmiss20.50[k2050] <- h$pmiss[i]
            df$rho20.50[k2050] <- h$rho[i]
            df$propRej20.50[k2050] <- h$rejPERM[i]
           

            k2050 <- k2050 + 1
        }else if(h$n[i] == 10*K & h$p[i] == 100){
            df$typeH[k20100] <- ty

            df$pmiss20.100[k20100] <- h$pmiss[i]
            df$rho20.100[k20100] <- h$rho[i]
            df$propRej20.100[k20100] <- h$rejPERM[i]
          

            k20100 <- k20100 + 1
        }else if(h$n[i] == 10*K & h$p[i] == 150){
            df$typeH[k20150] <- ty

            df$pmiss20.150[k20150] <- h$pmiss[i]
            df$rho20.150[k20150] <- h$rho[i]
            df$propRej20.150[k20150] <- h$rejPERM[i]
           

            k20150 <- k20150 + 1
        }else if(h$n[i] == 10*K & h$p[i] == 200){
            df$typeH[k20200] <- ty

            df$pmiss20.200[k20200] <- h$pmiss[i]
            df$rho20.200[k20200] <- h$rho[i]
            df$propRej20.200[k20200] <- h$rejPERM[i]
           

            k20200 <- k20200 + 1
        }
    }

    df$c1[which(df$typeH == ty)] <- round(typeC[[ty]][1])
    df$c2[which(df$typeH == ty)] <- typeC[[ty]][2]
}



table <- df[, c(which(colnames(df) %in% c("c1", "c2", "pmiss10.50", "rho10.50")), which(grepl("propRej", colnames(df))))]
table[, c(5:12)] <- round(table[, c(5:12)], 3)
table <- rbind(c(rep("", 4), rep("$n_k = 5$", 4), rep("$n_k = 10$", 4)), table)
colnames(table) <- c("$c_1$", "$c_2$", "$\\pi_{j1}$", "$\\rho$", rep(c("$p = 50$", "$p = 100$", "$p = 150$", "$p = 200$"), 2) )


# Table 4 
print(xtable(table), include.rownames = FALSE)


