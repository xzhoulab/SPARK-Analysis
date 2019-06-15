##-------------------------------------------------------------
## Breast Cancer Data Analysis
##-------------------------------------------------------------

rm(list = ls())

library(SPARK)

counts <- read.table("./raw_data/Layer2_BC_count_matrix-1.tsv", check.names = F)
rn <- rownames(counts)
info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), 
    "[", 1)), y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

spark <- CreateSPARKObject(counts = t(counts), location = info[, 1:2], 
    percentage = 0.1, min_total_counts = 10)

spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
    num_core = 10, verbose = T, fit.maxiter = 500)
spark <- spark.test(spark, check_positive = T, verbose = T)

save(spark, file = "./output/Layer2_BC_spark.rds")


##-------------------------------------------------------------
## Ten Times Permutation
##-------------------------------------------------------------

rm(list = ls())

library(SPARK)
counts <- read.table("./raw_data/Layer2_BC_count_matrix-1.tsv", check.names = F)
rn <- rownames(counts)
info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), 
    "[", 1)), y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

for (iper in 1:10) {
    set.seed(iper)
    ran_index <- sample(1:nrow(info))
    perm_info <- info[ran_index, ]
    rownames(perm_info) <- rownames(info)
    
    spark <- CreateSPARKObject(counts = t(counts), location = perm_info[, 
        1:2], percentage = 0.1, min_total_counts = 10)
    
    spark@lib_size <- apply(spark@counts, 2, sum)
    spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
        num_core = 10, verbose = T)
    spark <- spark.test(spark, check_positive = T, verbose = T)
    save(spark, file = paste0("./output/Layer_BC_perm", iper, "_spark.rds"))
    rm(spark)
}


##-------------------------------------------------------------
## Spatial Distribution of Representative Genes
##-------------------------------------------------------------
rm(list = ls())
source("./funcs/funcs.R")

counts <- read.table("./raw_data/Layer2_BC_count_matrix-1.tsv", check.names = F)
rn <- rownames(counts)
info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), 
    "[", 1)), y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

gene_plot <- c("HLA-B", "EEF1A1", "ERBB2", "MMP14", "CD44")

vst_ct <- var_stabilize(t(counts))
sig_vst_ct <- vst_ct[gene_plot, ]
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)


pltdat <- cbind.data.frame(info[, 1:2], rel_vst_ct)
genetitle <- c(expression("HLA-B " * (5.6 %*% 10^{-17})), 
				expression("EEF1A1 " * (3 %*% 10^{-7})), 
				expression("ERBB2 " * (2.9 %*% 10^{-6})), 
				expression("MMP14 " * (2.5 %*% 10^{-5})), 
				expression("CD44 " * (1.9 %*% 10^{-4})))

pp <- lapply(1:(ncol(pltdat) - 2), function(x) {
    pattern_plot2(pltdat, x, main = T, titlesize = 1.5, title = genetitle[x])
})
grid.arrange(grobs = pp, ncol = 3)


