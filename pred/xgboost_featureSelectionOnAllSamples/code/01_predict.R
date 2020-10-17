require(xgboost)  
require(Matrix)
require(data.table)
if (!require('vcd')) install.packages('vcd')
library(sva)

d <- readRDS('/dcs01/gcode/zji4/covid/analysis/pbmc/diffgene/pb/res/pbmcnorm.rds')
meta <- readRDS('/dcl02/hongkai/data/covid/data/200916/meta.rds')
cytokine <- readLines('/dcl02/hongkai/data/whou/covid/pred/xgboost/data/cytokine.txt')
cytokine <- sapply(cytokine, function(i) sub('     ', '', i), USE.NAMES = F)
cytokine <- sapply(cytokine, function(i) sub('   ', '', i), USE.NAMES = F)
cytokine <- sapply(cytokine, function(i) sub('  ', '', i), USE.NAMES = F)
cytokine <- sapply(cytokine, function(i) sub(' ', '', i), USE.NAMES = F)


# within each cluster, combat, leave one study out, xgboost predict 
# aim: see how much each cluster constribute to severity levels
for (i in 1:length(d)){
  d[[i]] <- ComBat(d[[i]],sub('.*-','',colnames(d[[i]])))
}
saveRDS(d,'/dcl02/hongkai/data/whou/covid/pred/data/data/pbmcnorm_combat.rds')

trans <- function(k) {
    a <- rep(0,length(k))
    a[k=='Mi'] <- 1
    a[k=='Se'] <- 2
    a
}

revtrans <- function(k){
  a <- rep('HD', length(k))
  a[k == 1] <- 'Mi'
  a[k == 2] <- 'Se'
  a
}

getcontmat <- function(pred, true){ # rows are true, col are pred
  contmat <- matrix(0, nrow = length(unique(true)), ncol = length(unique(pred)))
  rownames(contmat) <- unique(true)
  colnames(contmat) <- unique(pred)
  for (predi in unique(pred)){
      for (truei in unique(true)){
        tab <- table(pred[true == truei])
        contmat[truei, names(tab)] <- tab
    }
  }
  rownames(contmat) <- paste0('true_', rownames(contmat))
  colnames(contmat) <- paste0('pred_', colnames(contmat))
  return(contmat)
}
# pred
# true
# getcontmat(pred, true)
    

remove.s <- meta[which(meta[,'type'] == 'IAV' | meta[,'type'] == 'Rec'), 'Library Sample Code']

fea <- sapply(cytokine, function(i){
  cluid <- as.numeric(sub('.*cluster', '', i))
  g <- sub('.cluster.*', '', i)
  mat <- d[[cluid]]
  mat <- mat[, !colnames(mat) %in% remove.s]
  mat[g, ]
})
feamat <- matrix(0, nrow=length(d), ncol = ncol(d[[1]]))
rownames(feamat) <- names(fea)
colnames(feamat) <- names(fea[[1]])

for (i in 1:length(fea)){
  feamat[i, names(fea[[i]])] <- fea[[i]]
}


tmp <- d[[1]]
tmp <- tmp[, !colnames(tmp) %in% remove.s]
as <- meta[match(colnames(tmp), meta[,'Library Sample Code']), 'Paper code']



contmatlist <- list()
for (j in 1:length(unique(as))) {
  leaveid <- which(as == unique(as)[j]) 
  trainid <- setdiff(1:ncol(feamat), leaveid)
  y <- meta[match(colnames(feamat)[trainid], meta[, 'Library Sample Code']), 'type']

  model <- xgboost(data = t(feamat[, trainid]), label = trans(y), nrounds = 10,objective='multi:softprob',num_class=3)
  pred <- predict(model,t(feamat[,leaveid]))
  pred <- matrix(pred,nrow=3)
  id <- apply(pred,2,which.max) - 1
  predy <- revtrans(id)
  truey <- meta[match(colnames(feamat)[leaveid], meta[, 'Library Sample Code']), 'type']
  contmatlist[[j]] <- getcontmat(pred = predy, true = truey)
  num.correct[j] <- mean(id==trans(truey))
}
  
names(contmatlist) <- names(num.correct) <- unique(as)  

saveRDS(num.correct, '/dcl02/hongkai/data/whou/covid/pred/xgboost/res/num.correct.pred.rds')

saveRDS(contmatlist, '/dcl02/hongkai/data/whou/covid/pred/xgboost/res/contegency_table.rds')

library(pheatmap)
for (i in 1:length(contmatlist)){
  if (nrow(contmatlist[[i]]) > 1){
    fn <- sub('\\..*', '', names(contmatlist)[i])
    pdf(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost/plot/contegency_', fn, '.pdf'), width =4, height = 4)
    print(pheatmap(contmatlist[[i]], scale = 'none', cluster_rows = F, cluster_cols = F, main = fn))
    dev.off()
  }
}
  






