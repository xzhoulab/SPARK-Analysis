rm(list = ls())
library(data.table)
library(trendsceek)

counts <- t(read.csv("./processed_data/MERFISH_Animal18_Bregma0.11_countdata.csv", 
    row.names = 1))
info <- read.csv("./processed_data/MERFISH_Animal18_Bregma0.11_info.csv", row.names = 1)

## Filter genes on being expressed
min.ncells.expr = 3
min.expr = 5
counts_filt = genefilter_exprmat(counts, min.expr, min.ncells.expr)
dim(counts_filt)

xy <- info[, 1:2]
pp <- pos2pp(xy)

## Set marks as the logged normalized gene expression
log.fcn = log10
pp = set_marks(pp, counts_filt, log.fcn = log.fcn)

## Run trendsceek
nrand = 10000
ncores = 1

set.seed(1)
trendstat_list = trendsceek_test(pp, nrand, ncores)
save(trendstat_list, time_comp, file = paste0("./output/MERFISH_Animal18_Bregma0.11_sid1_tsc.rds"))

