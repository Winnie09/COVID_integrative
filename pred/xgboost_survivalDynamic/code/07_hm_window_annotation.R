rm(list=ls())
setwd('/Users/wenpinhou/Dropbox/covid/')
rdir <- 'pred/xgboost_survivalDynamic/result/'
pdir <- 'pred/xgboost_survivalDynamic/plot/'
d.pt <- readRDS('/Users/wenpinhou/Dropbox/covid/pred/data/data/sample_info_pb_opt.rds')
rownames(d.pt) <- d.pt$Patient
str(d.pt)
d.pt <- d.pt[d.pt[,'Clinical.Outcome'] %in% c('Recovered','Deceased'),]
d.pt <- d.pt[d.pt[,'Phenotype'] %in% c('Mi', 'Mod', 'Se'),]

stu <- d.pt$Study
sel <- names(which(sapply(unique(stu),function(s) {
  tab <- table(d.pt[stu==s,'Clinical.Outcome'])
  length(tab) == 2 & min(tab) > 1
})))
d.pt <- d.pt[stu %in% sel,]
pt <- d.pt$Pseudotime
qt <- quantile(pt,seq(0,1,length.out=21))

deceased.id = which(d.pt$Clinical.Outcome == 'Deceased')
windowcut <- pt[deceased.id[1:length(deceased.id)%%5 == 0]]
windowcut[length(windowcut)] <- max(pt)
windowcut <- c(deceased.id[1]-1, windowcut)

sp.splie <- sapply(1:(length(windowcut)-2),function(i) {
 d.pt[pt > windowcut[i] & pt <= windowcut[i+2],'Patient']
  
})
names(sp.splie) <- paste0('window', 1:length(sp.splie))
cut <- d.pt$Pseudotime[which(d.pt$Clinical.Outcome=='Deceased')[3]]
sp.splie[[1]] <- setdiff(sp.splie[[1]],d.pt$Patient[d.pt$Pseudotime < cut])

d.pt <- readRDS('/Users/wenpinhou/Dropbox/covid/pred/data/data/sample_info_pb_opt.rds')
rownames(d.pt) <- d.pt$Patient
str(d.pt)

mat <- t(sapply(1:length(sp.splie), function(i){
  (d.pt[,'Patient'] %in% sp.splie[[i]]) + 0
}))
dimnames(mat) <- list(paste0('window', 1:nrow(mat)), d.pt[,'Patient'])
colann <- data.frame(pseudotime = d.pt[,'Pseudotime'],
                     severity = sub('-.*', '',d.pt[,'Subject.type']),
                     survival = d.pt[,'Clinical.Outcome'],
                     age = d.pt[,'Age'],
                     sex = d.pt[,'Sex'])
rownames(colann) <- colnames(mat)
library(pheatmap)
pdf(paste0(pdir, 'hm_window_annotation.pdf'), width = 7, height = 6)
pheatmap(mat, cluster_cols = F, cluster_rows = F,
         border_color = NA,
         annotation_col = colann, 
         show_colnames = F)
dev.off()

