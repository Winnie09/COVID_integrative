require(xgboost)  
require(Matrix)
require(data.table)
if (!require('vcd')) install.packages('vcd')
library(sva)

d <- readRDS('/dcs01/gcode/zji4/covid/analysis/pbmc/diffgene/pb/res/pbmcnorm.rds')
meta <- readRDS('/dcl02/hongkai/data/covid/data/200916/meta.rds')

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

remove.s <- meta[which(meta[,'type'] == 'IAV' | meta[,'type'] == 'Rec'), 'Library Sample Code']

num.correct.list <- list()  
for (i in 1:length(d)){
  print(paste0('cluster ', i))
  mat <- d[[i]]
  mat <- mat[, !colnames(mat) %in% remove.s]
  as <- meta[match(colnames(mat), meta[,'Library Sample Code']), 'Paper code']
  res <- NULL
  
  num.correct <- sapply(1:length(unique(as)), function(j){
      leaveid <- which(as == unique(as)[j]) 
      trainid <- setdiff(1:ncol(mat), leaveid)
      y <- meta[match(colnames(mat)[trainid], meta[, 'Library Sample Code']), 'type']
    
      model <- xgboost(data = t(mat[, trainid]), label = trans(y), nrounds = 10,objective='multi:softprob',num_class=3)
      pred <- predict(model,t(mat[,leaveid]))
      pred <- matrix(pred,nrow=3)
      id <- apply(pred,2,which.max) - 1
      truey <- meta[match(colnames(mat)[leaveid], meta[, 'Library Sample Code']), 'type']
      mean(id==trans(truey))
    })
    names(num.correct) <- unique(as)  
    num.correct.list[[i]] <- num.correct
}

names(num.correct.list) <- paste0('cluster', 1:length(num.correct.list)) 
saveRDS(num.correct.list, '/dcl02/hongkai/data/whou/covid/pred/xgboost/res/num.correct.pred.rds')





