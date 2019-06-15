rm(list = ls())

library(data.table)
library(trendsceek)

counts <- t(read.table("./raw_data/Rep11_MOB_count_matrix-1.tsv", check.names = F))
tot <- apply(counts, 2, sum)

## Filter genes on being expressed
min.ncells.expr = 3
min.expr = 5
counts_filt = genefilter_exprmat(counts, min.expr, min.ncells.expr)
dim(counts_filt)

## Calculate gene variability stats
quantile.cutoff = 0.9  ## filter out the most lowly expressed genes from the fitting
method = "rlm"
vargenes_stats = calc_varstats(counts_filt, counts_filt, quant.cutoff = quantile.cutoff, 
    method = method)

## Select subset with the top-ranked most variable genes
n.topvar = 500
topvar.genes = rownames(vargenes_stats[["real.stats"]])[1:n.topvar]

## Subset normalized counts on the most variable genes
counts_sub = counts_filt[topvar.genes, ]
dim(counts_sub)

## Create point-pattern using positions as spatial distribution and
## expression levels as mark distribution
cn <- colnames(counts_sub)
xy <- matrix(as.numeric(do.call(rbind, strsplit(cn, split = "x"))), ncol = 2)
pp <- pos2pp(xy)


## Set marks as the logged normalized gene expression
log.fcn = log10
pp = set_marks(pp, counts_sub, log.fcn = log.fcn)

topvar.genes = rownames(vargenes_stats[["real.stats"]])[1:n.topvar]
pp2plot = pp_select(pp, topvar.genes)

## Run trendsceek
nrand = 10000
ncores = 1

for (isid in 1:10) {
    set.seed(isid)
    trendstat_list = trendsceek_test(pp2plot, nrand, ncores)
    save(trendstat_list, time_comp, file = paste0("./output/Rep11_MOB_rlm_top500_log10_sid",isid,".rds"))
}

