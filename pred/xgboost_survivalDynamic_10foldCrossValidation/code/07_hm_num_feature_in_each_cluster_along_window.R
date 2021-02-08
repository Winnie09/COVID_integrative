rm(list=ls())
setwd('/Users/wenpinhou/Dropbox/covid/')
rdir <- 'pred/xgboost_survivalDynamic/result/'
pdir <- 'pred/xgboost_survivalDynamic/plot/'
source('/Users/wenpinhou/Dropbox/trajectory_variability/function/mySTIP2.R')
library(pheatmap)
d.pt <- readRDS('/Users/wenpinhou/Dropbox/covid/pred/data/data/sample_info_pb_opt.rds')
rownames(d.pt) <- d.pt$Patient
str(d.pt)
d.pt <- d.pt[d.pt[,'Clinical.Outcome'] %in% c('Recovered','Deceased'),]

stu <- d.pt$Study
sel <- names(which(sapply(unique(stu),function(s) {
  tab <- table(d.pt[stu==s,'Clinical.Outcome'])
  length(tab) == 2 & min(tab) > 1
})))
d.pt <- d.pt[stu %in% sel,]
pt <- d.pt$Pseudotime
qt <- quantile(pt,seq(0,1,length.out=21))
type.splie <- sapply(1:(length(qt)-5),function(i) {
  table(d.pt[pt > qt[i] & pt < qt[i+5],'Clinical.Outcome'])
})
pd <- do.call(cbind, type.splie)
colnames(pd) <- paste0('window', 1:ncol(pd))
write.csv(pd, paste0(pdir, 'window_sample_size.csv'))

sp.splie <- pat <- sapply(1:(length(qt)-5),function(i) {
  d.pt[pt > qt[i] & pt < qt[i+5],'Patient']
},simplify = F)

## pseudobulk data
d <- readRDS('pred/data/data/pbmcnorm_combat.rds')
# > sapply(d, ncol)
#  c1  c2  c3  c4  c5  c6  c7  c8  c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 
# 429 429 429 429 429 418 425 419 262 333 314  12 419 414 415 391 366 254  20 
# c20 c21 c22 c23 c24 c25 c26 c27 c28 c29 c30 c31 
#  79   6   3  14 428 421 352 262 158 406   5 359 
id <- which(sapply(d, ncol) > 100)
d <- d[id]

## --------------------------
## read in candidate features
## --------------------------
# meta <- readRDS('pred/data/data/meta.rds')
cellProp.mat <- readRDS('pred/data/prop/prop.rds')
str(cellProp.mat)
cyto <- readRDS('pred/data/data/cytokine_receptor.rds')
str(cyto)
cytokine.v <- unique(as.character(cyto[,1]))
receptor.v <- unique(as.vector(apply(cyto[,2:ncol(cyto)], 2, as.character)))
str(cytokine.v)
str(receptor.v)
tf.v <- readRDS('pred/data/TF/TFonly.rds')
str(tf.v)

## ----------------
## define functions: Deceased = 0, Recovered = 1
## ----------------
af <- list.files(rdir, pattern = 'impt_')
impt <- t(sapply(af, function(f){
  impt <- readRDS(paste0(rdir, f))
  a <- cbind(window = gsub('_.*', '', sub('impt_', '', f)), type = strsplit(f,'_')[[1]][[3]], fea = impt[,1])
}))
impt <- do.call(rbind, impt)

for (type in unique(impt[,2])){
  print(type)
  impt.tmp <- impt[impt[,2] == type, ]
  impt.tmp <- cbind(impt.tmp, cluster = sub('cluster', '', sub(';.*', '', impt.tmp[,3])))
  num <- sapply(unique(impt.tmp[,1]), function(i){
    tab <- table(impt.tmp[impt.tmp[,1] == i, 4])
    v <- rep(0, 31)
    names(v) <- paste0('c', 1:31)
    v[paste0('c',names(tab))] <- tab
    v
  })
  pdf(paste0(pdir, 'hm_num_', type, '_along_pseudotime.pdf'), width = 2.8, height = 5)
  pheatmap(num, cluster_rows = FALSE, cluster_cols = FALSE)
  dev.off()
}
  
