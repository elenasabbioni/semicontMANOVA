###############################################################################
#### produce Table1, 2 and 3
###############################################################################
K <- 2 # number of groups

## Packages
library(xtable)

###########
## Create files with the summaries of the simulation study of semicontinuous_MANOVA
source("analyzeResults_semicontinuousMANOVA.R")

## Create files with the summaries of the simulation study of CHEN
source("analyzeResults_CHEN.R")

# Explored values of c_1 and c_2
typeC <- list()
typeC[["H0"]] <- c(0, 0)
typeC[["H1_1_1"]] <- c(1, 0)
typeC[["H1_1_5"]] <- c(5, 0)
typeC[["H1_2_0.15"]] <- c(0, 0.15)
typeC[["H1_2_0.3"]]  <- c(0, 0.3)

## Read the files with the summaries
pathSemMANOVA <- "./Results_semicontinuousMANOVA/K2/"
H0_semMANOVA <- data.frame(read.table(file=paste(pathSemMANOVA, "H0", "_res", ".csv", sep = ""),  sep=";", header = TRUE))
H1_1_1_semMANOVA <- read.table(file=paste(pathSemMANOVA, "H1_1-1", "_res", ".csv", sep = ""),  sep=";", header = TRUE)
H1_1_5_semMANOVA <- read.table(file=paste(pathSemMANOVA, "H1_1-5", "_res", ".csv", sep = ""),  sep=";", header = TRUE)
H1_2_0.15_semMANOVA <- read.table(file=paste(pathSemMANOVA, "H1_2-0.15", "_res", ".csv", sep = ""),  sep=";", header = TRUE)
H1_2_0.3_semMANOVA <- read.table(file=paste(pathSemMANOVA, "H1_2-0.3", "_res", ".csv", sep = ""),  sep=";", header = TRUE)


pathCHEN <- "./Results_CHEN/K2/"
H0_CHEN <- data.frame(read.table(file=paste(pathCHEN, "H0", "_res", ".csv", sep = ""),  sep=";", header = TRUE))
H1_1_1_CHEN <- read.table(file=paste(pathCHEN, "H1_1-1", "_res", ".csv", sep = ""),  sep=";", header = TRUE)
H1_1_5_CHEN <- read.table(file=paste(pathCHEN, "H1_1-5", "_res", ".csv", sep = ""),  sep=";", header = TRUE)
H1_2_0.15_CHEN <- read.table(file=paste(pathCHEN, "H1_2-0.15", "_res", ".csv", sep = ""),  sep=";", header = TRUE)
H1_2_0.3_CHEN <- read.table(file=paste(pathCHEN, "H1_2-0.3", "_res", ".csv", sep = ""),  sep=";", header = TRUE)




df <- data.frame(typeH = rep(NA,28*2), c1 = NA, c2 = NA,
pmiss10.50  = NA, rho10.50 = NA, propRej10.50 = NA, 
pmiss10.100 = NA, rho10.100 = NA, propRej10.100 = NA,
pmiss10.150 = NA, rho10.150 = NA, propRej10.150 = NA, 
pmiss10.200 = NA, rho10.200 = NA, propRej10.200 = NA, 
pmiss20.50  = NA, rho20.50 = NA, propRej20.50 = NA, 
pmiss20.100 = NA, rho20.100 = NA, propRej20.100 = NA, 
pmiss20.150 = NA, rho20.150 = NA, propRej20.150 = NA, 
pmiss20.200 = NA, rho20.200 = NA, propRej20.200 = NA, 
type = rep(c("scMAN", "CHEN"), 28),
rownames = FALSE)
 



k1050 <- 1
k10100 <- 1
k10150 <- 1
k10200 <- 1
k2050 <- 1
k20100 <- 1
k20150 <- 1
k20200 <- 1
for(ty in names(typeC)){
    h.semMANOVA <- get(paste0(ty, "_semMANOVA"))
    h.semMANOVA <- h.semMANOVA[, -c(which(colnames(h.semMANOVA) %in% c("n_simEnded", "meanDim0")))]

    h.CHEN <- get(paste0(ty, "_CHEN"))
    h <- h.semMANOVA
    h$rejPERM.CHEN <- h.CHEN$rejPERM

    for(i in 1:dim(h)[1]){
        if(h$n[i] == 5*K & h$p[i] == 50){
            df$typeH[k1050] <- ty
            df$pmiss10.50[k1050] <- h$pmiss[i]
            df$rho10.50[k1050] <- h$rho[i]
            df$propRej10.50[k1050] <- h$rejPERM[i]

            df[k1050 + 1, c("typeH", "pmiss10.50", "rho10.50")] <- df[k1050, c("typeH", "pmiss10.50", "rho10.50")]
            df$propRej10.50[k1050 + 1] <- h$rejPERM.CHEN[i]

            k1050 <- k1050 + 2
        }else if(h$n[i] == 5*K & h$p[i] == 100){
            df$typeH[k10100] <- ty

            df$pmiss10.100[k10100] <- h$pmiss[i]
            df$rho10.100[k10100] <- h$rho[i]
            df$propRej10.100[k10100] <- h$rejPERM[i]

            df[k10100 + 1, c("typeH", "pmiss10.100", "rho10.100")] <- df[k10100, c("typeH", "pmiss10.100", "rho10.100")]

            df$propRej10.100[k10100 + 1] <- h$rejPERM.CHEN[i]           

            k10100 <- k10100 + 2
        }else if(h$n[i] == 5*K & h$p[i] == 150){
            df$typeH[k10150] <- ty

            df$pmiss10.150[k10150] <- h$pmiss[i]
            df$rho10.150[k10150] <- h$rho[i]
            df$propRej10.150[k10150] <- h$rejPERM[i]

            df[k10150 + 1, c("typeH", "pmiss10.150", "rho10.150")] <- df[k10150, c("typeH", "pmiss10.150", "rho10.150")]

            df$propRej10.150[k10150 + 1] <- h$rejPERM.CHEN[i]           
            
            k10150 <- k10150 + 2
        }else if(h$n[i] == 5*K & h$p[i] == 200){
            df$typeH[k10200] <- ty

            df$pmiss10.200[k10200] <- h$pmiss[i]
            df$rho10.200[k10200] <- h$rho[i]
            df$propRej10.200[k10200] <- h$rejPERM[i]

            df[k10200 + 1, c("typeH", "pmiss10.200", "rho10.200")] <- df[k10200, c("typeH", "pmiss10.200", "rho10.200")]

            df$propRej10.200[k10200 + 1] <- h$rejPERM.CHEN[i]           
           

            k10200 <- k10200 + 2
        }else if(h$n[i] == 10*K & h$p[i] == 50){
            df$typeH[k2050] <- ty

            df$pmiss20.50[k2050] <- h$pmiss[i]
            df$rho20.50[k2050] <- h$rho[i]
            df$propRej20.50[k2050] <- h$rejPERM[i]

            df[k2050 + 1, c("typeH", "pmiss20.50", "rho20.50")] <- df[k2050, c("typeH", "pmiss20.50", "rho20.50")]

            df$propRej20.50[k2050 + 1] <- h$rejPERM.CHEN[i]           

           
            k2050 <- k2050 + 2
        }else if(h$n[i] == 10*K & h$p[i] == 100){
            df$typeH[k20100] <- ty

            df$pmiss20.100[k20100] <- h$pmiss[i]
            df$rho20.100[k20100] <- h$rho[i]
            df$propRej20.100[k20100] <- h$rejPERM[i]

            df[k20100 + 1, c("typeH", "pmiss20.100", "rho20.100")] <- df[k20100, c("typeH", "pmiss20.100", "rho20.100")]

            df$propRej20.100[k20100 + 1] <- h$rejPERM.CHEN[i]           
          
            k20100 <- k20100 + 2
        }else if(h$n[i] == 10*K & h$p[i] == 150){
            df$typeH[k20150] <- ty

            df$pmiss20.150[k20150] <- h$pmiss[i]
            df$rho20.150[k20150] <- h$rho[i]
            df$propRej20.150[k20150] <- h$rejPERM[i]

            df[k20150 + 1, c("typeH", "pmiss20.150", "rho20.150")] <- df[k20150, c("typeH", "pmiss20.150", "rho20.150")]

            df$propRej20.150[k20150 + 1] <- h$rejPERM.CHEN[i]           

            k20150 <- k20150 + 2
        }else if(h$n[i] == 10*K & h$p[i] == 200){
            df$typeH[k20200] <- ty

            df$pmiss20.200[k20200] <- h$pmiss[i]
            df$rho20.200[k20200] <- h$rho[i]
            df$propRej20.200[k20200] <- h$rejPERM[i]

            df[k20200 + 1, c("typeH", "pmiss20.200", "rho20.200")] <- df[k20200, c("typeH", "pmiss20.200", "rho20.200")]

            df$propRej20.200[k20200 + 1] <- h$rejPERM.CHEN[i]           
           

            k20200 <- k20200 + 2
        }
    }

    df$c1[which(df$typeH == ty)] <- round(typeC[[ty]][1])
    df$c2[which(df$typeH == ty)] <- typeC[[ty]][2]
}



table <- df[, c(which(colnames(df) %in% c("c1", "c2", "pmiss10.50", "rho10.50", "type")), which(grepl("propRej", colnames(df))))]
table[, c(6:13)] <- round(table[, c(6:13)], 3)
table <- rbind(c(rep("", 5), rep("$n_k = 5$", 4), rep("$n_k = 10$", 4)), table)
colnames(table) <- c("$c_1$", "$c_2$", "$\\pi_{j1}$", "$\\rho$", "test", rep(c("$p = 50$", "$p = 100$", "$p = 150$", "$p = 200$"), 2))

table[which(table$test == "CHEN"), 1:4] <- rep("", 4)


table1 <- table[1:13, ]
table2 <- table[c(1, 14:37), ]
table3 <- table[c(1, 38:57), ]

# Table 1
print(xtable(table1), include.rownames = FALSE)

# Table 2
print(xtable(table2), include.rownames = FALSE)

# Table 3
print(xtable(table3), include.rownames = FALSE)


