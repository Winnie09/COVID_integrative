rm(list=ls())
setwd('/Users/wenpinhou/Dropbox/covid/')
rdir <- 'pred/xgboost_survivalDynamic/result/'
pdir <- 'pred/xgboost_survivalDynamic/plot/'
af <- list.files(rdir, pattern = 'AUC_window')
auc <- sapply(af, function(f){
  auc <- readRDS(paste0(rdir, f))
  auc <- colMeans(auc)
})
colnames(auc) <- gsub('.rds', '', sub('AUC_', '', colnames(auc)))
saveRDS(auc, paste0(rdir, 'AUC_summary.rds'))
write.csv(auc, paste0(pdir, 'AUC_summary.csv'))


pd <- reshape2::melt(auc)
str(pd)
colnames(pd) <- c('featureType', 'window', 'AUC')
pd[,1] <- as.character(pd[,1])
pd[,2] <- as.character(pd[,2])
library(pheatmap)
breaksList <- seq(0.5, max(auc), by = 0.01)
pdf(paste0(pdir, 'AUC_feature_by_window.pdf'), width = 6, height = 5)
pheatmap(auc, cluster_cols = FALSE, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)
dev.off()

pdf(paste0(pdir, 'AUC_feature_by_window_rank.pdf'), width = 6, height = 5)
pheatmap(apply(auc, 2, rank), cluster_cols = FALSE)
dev.off()

