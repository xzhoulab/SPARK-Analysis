##------------------------------------------------
## Visualization of the pattern summarized from SpatialDE result
##------------------------------------------------

rm(list = ls())
source("./funcs/funcs.R")
load("./output/MOB_Pattern_SpatialDE.rds")

tmplist <- datlist
for (i in 1:5) {
    tmplist[[i]][, 2] <- relative_func(datlist[[i]][, 2])
}

patterns = c("I", "II", "III")
## three major pattern were used for simulation
df <- setNames(cbind.data.frame(tmplist[[1]][, 1], do.call(cbind, sapply(tmplist[c(5, 
    1, 4)], "[", 2))), c("xy", paste0("Pattern ", patterns)))
pp <- lapply(1:3, function(x) {
    pattern_plot2(df, x, xy = F, main = T, titlesize = 1.5)
})


grid.arrange(grobs = pp, ncol = 3)


##------------------------------------------------
## Data Generation
##------------------------------------------------
## Generate the countdata based on the info from the realdata and
## patterns summarized from the SpatialDE result

## Fold Change
##-----------------
rm(list = ls())

load(paste0("./output/Rep11_MOB_spark.rds"))
info <- spark@location
info$total_counts <- spark@lib_size

beta <- sapply(1:length(spark@res_vc), function(x) {
    spark@res_vc[[x]]$coefficients
})
nb <- sapply(1:4, function(x) {
    log(x * exp(median(beta)))
})
tau2 <- sapply(1:length(spark@res_vc), function(x) {
    spark@res_vc[[x]]$theta[2]
})

load("./output/MOB_Pattern_SpatialDE.rds")

itau = 35

for (ifc in 2:4) {
    newN <- info$total_counts
    for (ipt in 2:4) {
        pattern <- datlist[[paste0("pattern", ipt)]]
        grp <- as.numeric(pattern[, 2] > mean(pattern[, 2])) + 1
        uu <- c(nb[1], nb[ifc], nb[1])
        numSignal <- 1000
        numGene <- 10000
        for (ipow in 1:10) {
            set.seed(ipow)
            beta0 <- uu[grp]
            lambda0 <- sapply(1:numSignal, function(x) {
                exp(beta0 + rnorm(length(beta0), 0, itau/100))
            })
            newCt0 <- lapply(1:numSignal, function(x) {
                rpois(length(lambda0[, x]), newN * lambda0[, x])
            })
            
            beta1 <- rep(uu[3], nrow(pattern))
            lambda1 <- sapply(1:(numGene - numSignal), function(x) {
                exp(beta1 + rnorm(length(beta1), 0, itau/100))
            })
            newCt1 <- lapply(1:(numGene - numSignal), function(x) {
                rpois(length(lambda1[, x]), newN * lambda1[, x])
            })
            
            countdata <- data.frame(rbind(do.call(rbind, newCt0), do.call(rbind, 
                newCt1)))
            rownames(countdata) <- paste0("gene", 1:nrow(countdata))
            colnames(countdata) <- pattern[, 1]
            
            write.csv(t(countdata), file = paste0("./data/sim_MOB_pattern", 
                ipt, "_fc", ifc, "_tau", itau, "_count_power", ipow, ".csv"), 
                row.names = T)
        }
    }
}

write.csv(info, file = "./data/Rep11_MOB_info_spark.csv", row.names = T)


## Noise Change
##-----------------

rm(list = ls())
load(paste0("./output/Rep11_MOB_spark.rds"))
info <- spark@location
info$total_counts <- spark@lib_size

beta <- sapply(1:length(spark@res_vc), function(x) {
    spark@res_vc[[x]]$coefficients
})
nb <- sapply(1:4, function(x) {
    log(x * exp(median(beta)))
})
tau2 <- sapply(1:length(spark@res_vc), function(x) {
    spark@res_vc[[x]]$theta[2]
})

load("./output/MOB_Pattern_SpatialDE.rds")

for (itau in c(20, 60)) {
    for (ifc in 3:3) {
        newN <- info$total_counts
        for (ipt in 2:4) {
            pattern <- datlist[[paste0("pattern", ipt)]]
            grp <- as.numeric(pattern[, 2] > mean(pattern[, 2])) + 1
            uu <- c(nb[1], nb[ifc], nb[1])
            numSignal <- 1000
            numGene <- 10000
            
            for (ipow in 1:10) {
                set.seed(ipow)
                beta0 <- uu[grp]
                lambda0 <- sapply(1:numSignal, function(x) {
                  exp(beta0 + rnorm(length(beta0), 0, itau/100))
                })
                newCt0 <- lapply(1:numSignal, function(x) {
                  rpois(length(lambda0[, x]), newN * lambda0[, x])
                })
                
                beta1 <- rep(uu[3], nrow(pattern))
                lambda1 <- sapply(1:(numGene - numSignal), function(x) {
                  exp(beta1 + rnorm(length(beta1), 0, itau/100))
                })
                newCt1 <- lapply(1:(numGene - numSignal), function(x) {
                  rpois(length(lambda1[, x]), newN * lambda1[, x])
                })
                
                countdata <- data.frame(rbind(do.call(rbind, newCt0), do.call(rbind, 
                  newCt1)))
                rownames(countdata) <- paste0("gene", 1:nrow(countdata))
                colnames(countdata) <- pattern[, 1]
                
                write.csv(t(countdata), file = paste0("./data/sim_MOB_pattern", 
                  ipt, "_fc", ifc, "_tau", itau, "_count_power", ipow, 
                  ".csv"), row.names = T)
            }
        }
    }
}



## Analyze with SPARK
##-----------------
rm(list = ls())

library(SPARK)
itau = 35
ifc = 3
ipt = 2
ipow = 1

info <- read.csv("./processed_data/Rep11_MOB_info_spark.csv", row.names = 1)

countdata <- t(read.csv(paste0("./processed_data/sim_MOB_pattern", ipt, "_fc", ifc, 
    "_tau", itau, "_count_power", ipow, ".csv"), row.names = 1))
spark <- CreateSPARKObject(counts = countdata, location = info[, 1:2], 
    percentage = 0.1, min_total_counts = 10)

spark@lib_size <- info$total_counts

t1 <- proc.time()
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
    num_core = 1, verbose = T, fit.maxiter = 500)
spark <- spark.test(spark, check_positive = T, verbose = T)
time_comp <- proc.time() - t1

