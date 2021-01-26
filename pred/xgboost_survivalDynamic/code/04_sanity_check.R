rm(list=ls())
setwd('/Users/wenpinhou/Dropbox/covid/')
rdir <- 'pred/xgboost_survivalDynamic/result/'
pdir <- 'pred/xgboost_survivalDynamic/plot/'
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
str(id)
d <- d[id]
str(d)
length(d)

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

library(ggpubr)

plist1 <- plist2 <- list()
for (window in 11:16){
  print(window)
  sp.select <- sp.splie[[window]]
  
  pbmc.s <- intersect(colnames(d[[1]]), sp.select)
  
  pd <- data.frame(sample = pbmc.s, sex = d.pt[pbmc.s, 'Sex'], age = d.pt[pbmc.s, 'Age'], outcome = d.pt[pbmc.s, 'Clinical.Outcome'],stringsAsFactors = F )
  
  library(ggplot2)
  plist1[[paste0('window',window)]] <- ggplot(data = pd,aes(x = outcome, y = age, fill = outcome)) + geom_boxplot() + theme_classic() + ggtitle(paste0('window',window)) +
    theme(legend.position = 'bottom') +
    stat_compare_means(method = "t.test")
    
  plist2[[paste0('window',window)]] <- ggplot(data = pd) + geom_bar(aes(x = outcome, fill = sex),position='dodge') + theme_classic() + ggtitle(paste0('window',window)) +
    theme(legend.position = 'bottom')
}
  
library(gridExtra)
pdf(paste0(pdir, 'age.pdf'), width = 8, height = 8)
grid.arrange(grobs = plist1, nrow = 2)
dev.off()

pdf(paste0(pdir, 'sex.pdf'), width = 8, height = 8)
grid.arrange(grobs = plist2, nrow = 2)
dev.off()

  
  