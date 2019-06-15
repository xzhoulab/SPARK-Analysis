##-------------------------------------------------------------
## Extract counts from the raw data
##-------------------------------------------------------------

rm(list = ls())

library(openxlsx)
dat <- read.xlsx("./raw_data/mmc6.xlsx", sheet = 1, colNames = F, rowNames = T)
field <- read.xlsx("./raw_data/mmc6.xlsx", sheet = 2, colNames = F)
allinfo <- read.xlsx("./raw_data/mmc6.xlsx", sheet = 3, colNames = F)

count43 <- dat[, which(field[1, ] == 43)]
info43 <- setNames(allinfo[which(field[1, ] == 43), ], nm = c("x", "y", 
    "z"))
xidx <- which(info43$x > 203 & info43$x < 822)
yidx <- which(info43$y > 203 & info43$y < 822)
info <- info43[intersect(xidx, yidx), ]
counts <- count43[, intersect(xidx, yidx)]
rownames(info) <- colnames(counts) <- paste0("C", 1:nrow(info))

write.csv(info, file = "./processed_data/seqFISH_field43_info.csv", row.names = T)
write.csv(t(counts), file = "./processed_data/seqFISH_field43_countdata.csv", 
    row.names = T)


##-------------------------------------------------------------
## SeqFISH Data Analysis
##-------------------------------------------------------------
rm(list = ls())

library(SPARK)

counts <- t(read.csv("./data/seqFISH_field43_countdata.csv", row.names = 1, 
    check.names = F))
info <- read.csv("./data/seqFISH_field43_info.csv", row.names = 1)
spark <- CreateSPARKObject(counts = counts, location = info[, 1:2], 
    percentage = 0.1, min_total_counts = 10)

spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
    num_core = 10, verbose = T, fit.maxiter = 500)
spark <- spark.test(spark, check_positive = T, verbose = T)

save(spark, file = "./output/seqFISH_field43_spark.rds")

##-------------------------------------------------------------
## One Hundred Times Permutation
##-------------------------------------------------------------
rm(list = ls())

library(SPARK)

counts <- read.csv("./processed_data/seqFISH_field43_countdata.csv", row.names = 1, 
    check.names = F)
info <- read.csv("./processed_data/seqFISH_field43_info.csv", row.names = 1)

for (iper in 1:100) {
    set.seed(iper)
    ran_index <- sample(1:nrow(info))
    perm_info <- info[ran_index, ]
    rownames(perm_info) <- rownames(info)
    
    spark <- CreateSPARKObject(counts = t(counts), location = perm_info[, 
        1:2], percentage = 0, min_total_counts = 10)
    
    spark@lib_size <- apply(spark@counts, 2, sum)
    spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
        num_core = 10, verbose = T)
    spark <- spark.test(spark, check_positive = T, verbose = T)
    save(spark, file = paste0("./output/MERFISH_Animal18_Bregma0.11_perm", 
        iper, "_spark.rds"))
    rm(spark)
}


##-------------------------------------------------------------
## Spatial Distribution of Representative Genes
##-------------------------------------------------------------

rm(list = ls())
source("./funcs/funcs.R")

load("./output/seqFISH_field43_spark.rds")
counts <- spark@counts
info <- spark@location

newrn <- sapply(strsplit(rownames(counts), split = "'"), "[", 2)
rownames(counts) <- newrn

gene_plot <- c("mog", "myl14", "lyve", "ndnf", "foxj1", "sst", "ctss", "xdh", 
				"gad1", "camk2", "mfge8")
acatp <- spark@res_mtest$combined_pvalue
names(acatp) <- newrn

glsm_res <- t(sapply(1:nrow(spark@counts), function(x) {
    spark@res_vc[[x]]$residuals
}))

rownames(glsm_res) <- rownames(counts) <- sapply(strsplit(rownames(spark@counts), split = "'"), "[[", 2)

sig_glsm_res <- glsm_res[gene_plot, ]
rel_glsm_res <- apply(sig_glsm_res, 1, relative_func)
pltdat <- cbind.data.frame(info[, 1:2], rel_glsm_res)


colnames(pltdat) <- c("x", "y", paste0(rownames(sig_glsm_res), " (", format(acatp[gene_plot], 
    digits = 1), ")"))
pp <- lapply(1:(ncol(pltdat) - 2), function(x) {
    pattern_plot3(pltdat, x, main = T, titlesize = 1.5, pointsize = 3, 
        min.pand = 0.9, max.pand = 1.05)
})
grid.arrange(grobs = pp, ncol = 4)


