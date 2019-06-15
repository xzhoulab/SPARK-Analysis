##-------------------------------------------------------------
##  Mouse OB Data Analysis  
##-------------------------------------------------------------

rm(list = ls())

# load the R package
library(SPARK)

# read the raw counts
counts <- read.table("./raw_data/Rep11_MOB_count_matrix-1.tsv", check.names = F)
rn <- rownames(counts)
info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), "[", 1)), 
                         y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

spark <- CreateSPARKObject(counts = t(counts), location = info[, 1:2], 
                           percentage = 0.1, min_total_counts = 10)

spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
                  num_core = 10, verbose = T, fit.maxiter = 500)
spark <- spark.test(spark, check_positive = T, verbose = T)

save(spark, file = "./output/Rep11_MOB_spark.rds")

##-------------------------------------------------------------
## Spatial Distribution of Representative Genes 
##-------------------------------------------------------------
rm(list = ls())
source("./funcs/funcs.R")

counts <- read.table("./raw_data/Rep11_MOB_count_matrix-1.tsv", check.names = F)
rn <- rownames(counts)
info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), "[", 1)), 
                         y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

gene_plot <- c("Reln", "Cldn5", "Gng4", "Doc2g", "Kctd12", "Penk")
vst_ct <- var_stabilize(t(counts)) # R function in funcs.R
sig_vst_ct <- vst_ct[gene_plot, ]
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)


pltdat <- cbind.data.frame(info[,1:2],rel_vst_ct)
genetitle <- c(expression("Reln "*(7.6 %*% 10^{-09})),
               expression("Cldn5 "*(8.4 %*% 10^{-13})),
               expression("Gng4 "*(5.6 %*% 10^{-17})),
               expression("Doc2g "*(1.6 %*% 10^{-10})),
               expression("Kctd12 "*(5.6 %*% 10^{-17})),
               expression("Penk "*(5.6 %*% 10^{-17})))

pp <- lapply(1:(ncol(pltdat)-2),
			function(x){pattern_plot2(pltdat,x,main=T,titlesize=1.5,title=genetitle[x])})
grid.arrange(grobs=pp, ncol=3)




##-------------------------------------------------------------
## Summarize Pattern using Hclust
##-------------------------------------------------------------
rm(list = ls())
source("./funcs/funcs.R")

LMReg <- function(ct, T) {
  return(lm(ct ~ T)$residuals)
}

load("./output/Rep11_MOB_spark.rds")

counts <- spark@counts
info <- spark@location
vst_count <- var_stabilize(counts) # R function in funcs.R
sig_vst_count <- vst_count[which(spark@res_mtest$adjusted_pvalue < 0.05), 
                           ]
sig_vst_res <- t(apply(sig_vst_count, 1, LMReg, T = log(spark@lib_size)))


library(amap)
hc <- hcluster(sig_vst_res, method = "euc", link = "ward", nbproc = 1, 
               doubleprecision = TRUE)
numC <- 3
memb <- cutree(hc, k = numC)

cent <- NULL
for (k in 1:numC) {
  cent <- cbind(cent, colMeans(sig_vst_res[memb == k, , drop = FALSE]))
}

position_cord <- info[, 1:2]
rownames(position_cord) <- rownames(cent)

rel_cent <- t(apply(cent, 1, relative_func))
pd <- setNames(cbind.data.frame(position_cord, rel_cent), c("x", "y", paste0("Pattern", c("III", "II", "I"))))
# pd <-
# setNames(cbind.data.frame(position_cord,rel_cent),c('x','y',paste0('Pattern',c('I','II','III','IV','V'))))
MBP <- lapply(1:numC, function(x) {
  pattern_plot2(pd, x, xy = T, main = T, titlesize = 1.5)
})

grid.arrange(grobs = MBP[3:1], ncol = numC)


