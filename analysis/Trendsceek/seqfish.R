
rm(list = ls())
## process with the seqFISH
winsorize <- function(x, n.win = 1) {
    n.vals = length(x)
    win.sorted.ind = order(x)
    x[win.sorted.ind[(n.vals - n.win + 1):n.vals]] = x[win.sorted.ind[n.vals - 
        n.win]]
    return(x)
}

library(trendsceek)
library(data.table)

counts <- t(read.csv("./processed_data/seqFISH_field43_countdata.csv", row.names = 1))
info <- read.csv("./processed_data/seqFISH_field43_info.csv", row.names = 1)
win.counts <- t(apply(counts, 1, winsorize, 4))

## Filter genes on being expressed
min.ncells.expr = 3
min.expr = 5
counts_filt = genefilter_exprmat(win.counts, min.expr, min.ncells.expr)
dim(counts_filt)

xy <- info[, 1:2]
pp <- pos2pp(xy)

## Set marks as the logged normalized gene expression
log.fcn = log10
pp = set_marks(pp, counts_filt, log.fcn = log.fcn)

## Run trendsceek
nrand = 10000
ncores = 1

for (isid in 1:10) {
    set.seed(isid)
    trendstat_list = trendsceek_test(pp, nrand, ncores)
    ## one sided winsorize
    save(trendstat_list, file = paste0("./output/seqFISH_field43_all_log10_sid",isid, "_osw_tsc.rds"))
    rm(trendstat_list)
}


