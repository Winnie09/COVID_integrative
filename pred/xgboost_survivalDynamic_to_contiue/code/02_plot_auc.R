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

str(auc)
window.n <- as.numeric(sub('window', '', colnames(auc)))
auc <- auc[, order(window.n)]
library(pheatmap)
pdf(paste0(pdir, 'AUC_feature_by_window.pdf'), width = 6, height = 5)
pheatmap(auc, cluster_cols = FALSE)
dev.off()


pdf(paste0(pdir, 'AUC_feature_by_window_rank.pdf'), width = 6, height = 5)
pheatmap(apply(auc, 2, rank), cluster_cols = FALSE)
dev.off()

