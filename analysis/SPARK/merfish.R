##-------------------------------------------------------------
## Extract counts from the raw data
##-------------------------------------------------------------

rm(list = ls())
library(openxlsx)
genelist <- read.xlsx("./raw_data/aau5324_Moffitt_Table-S6.xlsx", startRow = 2)

library(data.table)
newcount <- fread("./raw_data/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv")

single_animal <- newcount[newcount$Animal_ID == 18, ]
animal18_11 <- as.data.frame(single_animal[single_animal$Bregma == 0.11, 
    ])

select_gene <- c(paste0("Blank_", 1:5), genelist$Gene.name)
geneidx <- which(colnames(animal18_11) %in% select_gene)

counts <- t(animal18_11[-which(animal18_11$Cell_class %in% c("Ambiguous")), 
    geneidx]) * 1000
info <- animal18_11[-which(animal18_11$Cell_class %in% c("Ambiguous")), 
    c("Centroid_X", "Centroid_Y", "Cell_class")]
colnames(info) <- c("x", "y", "Cell_class")

write.csv(info, file = paste0("./processed_data/MERFISH_Animal18_Bregma0.11_info.csv"), 
    row.names = T)
write.csv(t(counts), file = paste0("./processed_data/MERFISH_Animal18_Bregma0.11_countdata.csv"), 
    row.names = T)

##-------------------------------------------------------------
## MERFISH Data Analysis
##-------------------------------------------------------------

rm(list = ls())
library(SPARK)
counts <- t(read.csv("./processed_data/MERFISH_Animal18_Bregma0.11_countdata.csv", 
    row.names = 1))
info <- read.csv("./processed_data/MERFISH_Animal18_Bregma0.11_info.csv", row.names = 1)
spark <- CreateSPARKObject(counts = counts, location = info[, 1:2], 
    percentage = 0, min_total_counts = 10)

spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
    num_core = 10, verbose = T, fit.maxiter = 500)
spark <- spark.test(spark, check_positive = T, verbose = T)


save(spark, file = paste0("./output/MERFISH_Animal18_Bregma0.11_spark.rds"))



##-------------------------------------------------------------
## One Hundred Times Permutation
##-------------------------------------------------------------

rm(list = ls())

library(SPARK)

counts <- read.csv("./processed_data/MERFISH_Animal18_Bregma0.11_countdata.csv", 
    row.names = 1, check.names = F)
info <- read.csv("./processed_data/MERFISH_Animal18_Bregma0.11_info.csv", row.names = 1)

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
## Cell Class Plot
##-------------------------------------------------------------

## all in one
##-----------------------
rm(list = ls())

info <- read.csv("./processed_data/MERFISH_Animal18_Bregma0.11_info.csv", row.names = 1)
tmp1 <- sapply(strsplit(as.character(info$Cell_class), split = " "), "[",1)
tmp2 <- sapply(strsplit(as.character(info$Cell_class), split = " "), "[",2)

idx <- which(tmp1 == "OD")
cell_label <- tmp1
cell_label[idx] <- paste0(tmp2[idx], " ", tmp1[idx])
info$cell_label <- cell_label

pd <- cbind.data.frame(x = info$x, y = info$y, cell_label = info$cell_label)

## set the name Pericytes to Mural, consistent with the original paper
cellname1 <- c("Inhibitory", "Excitatory", "Mature OD", "Immature OD", 
    "Astrocyte", "Microglia", "Ependymal", "Endothelial", "Pericytes")
cellname2 <- c("Inhibitory", "Excitatory", "Mature OD", "Immature OD", 
    "Astrocyte", "Microglia", "Ependymal", "Endothelial", "Mural")

tmp3 <- LETTERS[1:9]
names(tmp3) <- cellname1

pd$label_idx <- tmp3[as.character(pd$cell_label)]

min.pand = 0.99
max.pand = 1.01
pointsize = 1
titlesize = 1

library(ggplot2)
library(RColorBrewer)
pal <- c("#BC3C29FF", "#5CB85CFF", "#E18727FF", "#1B191999", "#7876B1FF", 
    "#6F99ADFF", "#0072B5FF", "#EE4C9799", "#9632B8FF")

class_dist <- ggplot(pd, aes(x = x, y = y, color = label_idx)) + geom_point(size = pointsize) + 
    scale_colour_manual(name = "", values = pal, labels = cellname2) + 
    scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0, 
    1)) + expand_limits(x = c(min(pd$x) * min.pand, max(pd$x) * max.pand), 
    y = c(min(pd$y) * min.pand, max(pd$y) * max.pand)) + labs(title = "Bregma 0.11 (mm)", 
    x = NULL, y = NULL) + theme_bw() + theme(legend.position = "none", 
    plot.title = element_text(hjust = 0.5, size = rel(titlesize)))



## individual class plots
##-----------------------
rm(list = ls())

library(ggplot2)
library(RColorBrewer)

cell_class_plot <- function(pltdat, iclass, xy = T, main = F, titlesize = 2, 
    pointsize = 3, min.pand = 0.99, max.pand = 1.01, title = NULL) {
    pd <- pltdat
    pal <- c("#BC3C29FF", "#5CB85CFF", "#E18727FF", "#1B191999", "#7876B1FF", 
        "#6F99ADFF", "#0072B5FF", "#EE4C9799", "#9632B8FF")
    
    cellname1 <- c("Inhibitory", "Excitatory", "Mature OD", "Immature OD", 
        "Astrocyte", "Microglia", "Ependymal", "Endothelial", "Pericytes")
    gpt <- ggplot(pd, aes(x = x, y = y, color = (cell_label == cellname1[iclass]))) + 
        geom_point(size = pointsize) + scale_colour_manual(name = "", values = c("gray90", 
        pal[iclass])) + scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0, 
        1)) + expand_limits(x = c(min(pd$x) * min.pand, max(pd$x) * max.pand), 
        y = c(min(pd$y) * min.pand, max(pd$y) * max.pand)) + theme_bw()
    
    
    if (main) {
        if (is.null(title)) {
            title = colnames(pd)[igene + 2]
        }
        out = gpt + labs(title = title, x = NULL, y = NULL) + theme(legend.position = "none", 
            plot.title = element_text(hjust = 0.5, size = rel(titlesize)))
    } else {
        out = gpt + labs(title = NULL, x = NULL, y = NULL) + theme(legend.position = "none")
    }
    return(out)
}


info <- read.csv("./processed_data/MERFISH_Animal18_Bregma0.11_info.csv", row.names = 1)
tmp1 <- sapply(strsplit(as.character(info$Cell_class), split = " "), "[", 1)
tmp2 <- sapply(strsplit(as.character(info$Cell_class), split = " "), "[",2)

idx <- which(tmp1 == "OD")
cell_label <- tmp1
cell_label[idx] <- paste0(tmp2[idx], " ", tmp1[idx])
info$cell_label <- cell_label
pd <- cbind.data.frame(x = info$x, y = info$y, cell_label = info$cell_label)

cellname1 <- c("Inhibitory", "Excitatory", "Mature OD", "Immature OD", 
    "Astrocyte", "Microglia", "Ependymal", "Endothelial", "Pericytes")
cellname2 <- c("Inhibitory", "Excitatory", "Mature OD", "Immature OD", 
    "Astrocyte", "Microglia", "Ependymal", "Endothelial", "Mural")
tmp3 <- LETTERS[1:9]
names(tmp3) <- cellname1
pd$label_idx <- tmp3[as.character(pd$cell_label)]

min.pand = 0.99
max.pand = 1.01
pointsize = 1
titlesize = 1

pp <- lapply(1:9, function(x) {
    cell_class_plot(pltdat = pd, iclass = x, main = T, title = cellname2[x], 
        pointsize = 1, titlesize = 1.5)
})
grid.arrange(grobs = pp, ncol = 3)




##-------------------------------------------------------------
## Spatial Distribution of Representative Genes
##-------------------------------------------------------------

rm(list = ls())

library(SPARK)

source("./funcs/funcs.R")

counts <- read.csv("./processed_data/MERFISH_Animal18_Bregma0.11_countdata.csv", 
    row.names = 1)
info <- read.csv("./processed_data/MERFISH_Animal18_Bregma0.11_info.csv", row.names = 1)

spark <- CreateSPARKObject(counts = t(counts), location = info[, 1:2], 
    percentage = 0.1, min_total_counts = 10)

rm(counts, info)

info <- spark@location
info$total_counts <- apply(spark@counts, 2, sum)
counts <- spark@counts


gene_plot <- c("Gad1", "Mbp", "Cd24a", "Myh11")
sig_ct <- counts[gene_plot, ]
sig_vst_ct <- var_stabilize(sig_ct)
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)
pltdat <- cbind.data.frame(info[, 1:2], rel_vst_ct)

# acatp[mainGene]
genetitle <- c(expression("Gad1 " * (5.6 %*% 10^{-17})), 
				expression("Mbp " * (5.6 %*% 10^{-17})), 
				expression("Cd24a " * (5.6 %*% 10^{-17})), 
				expression("Myh11 " * (5.6 %*% 10^{-17})))
pp <- lapply(1:(ncol(pltdat) - 2), function(x) {
    pattern_plot3(pltdat, x, main = T, titlesize = 1.5, pointsize = 1, 
        title = genetitle[x])
})

grid.arrange(grobs = pp, ncol = 2)



