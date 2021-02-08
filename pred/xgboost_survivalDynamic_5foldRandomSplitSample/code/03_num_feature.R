rm(list=ls())
setwd('/Users/wenpinhou/Dropbox/covid/')
rdir <- 'pred/xgboost_survivalDynamic/result/'
pdir <- 'pred/xgboost_survivalDynamic/plot/'
auc <- readRDS(paste0(rdir, 'AUC_summary.rds'))

af <- list.files(rdir, pattern = 'impt_')

impt <- t(sapply(af, function(f){
  impt <- readRDS(paste0(rdir, f))
  c(window = gsub('_.*', '', sub('impt_', '', f)), type = strsplit(f,'_')[[1]][[3]], num = nrow(impt))
}))
str(impt)
head(impt)
impt <- data.frame(impt)
impt <- tapply(as.numeric(impt$num),list(paste0(impt$window,'_',impt$type)),sum)
impt <- data.frame(window = sub('_.*', '', names(impt)), featureType = sub('.*_', '', names(impt)),
                   num = impt)


window = colnames(auc)[1]
numfea <- sapply(colnames(auc), function(window){
  fealen <- c(impt[impt[,1] == window, 3], 2, 1)
  names(fealen) <- c(impt[impt[,1] == window, 2], 'base', 'cellProp')
  sapply(rownames(auc), function(i){
   tmp <- strsplit(i, '_') [[1]]
   tmp[tmp == 'rec'] <- 'receptor'
   tmp[tmp == 'cyt'] <- 'cytokine'
   if ('all' %in% tmp){
    return(sum(fealen)) 
   } else if ('gene' %in% tmp){
    tmp <- tmp[!tmp %in% c('cytokine', 'receptor', 'tf')]   
    return(sum(fealen[tmp]) )
   } else {
    return(sum(fealen[tmp]) )
   } 
  })
})
numfea <- numfea[order(rowSums(numfea)), ]
pdf(paste0(pdir, 'number_of_clusterFeature_by_window.pdf'), width = 6, height = 5)
pheatmap(numfea, cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

write.csv(numfea, paste0(pdir, 'number_of_clusterFeature_by_window.csv'))

######### unique features
impt.unique <- lapply(af, function(f){
  impt <- readRDS(paste0(rdir, f))
  a <- cbind(window = gsub('_.*', '', sub('impt_', '', f)), type = strsplit(f,'_')[[1]][[3]], impt)
})
impt <- do.call(rbind, impt.unique)
impt <- cbind(impt, windowType = paste0(impt[,1],';', impt[,2]))

a <- t(sapply(unique(impt$windowType), function(i){
  tmp <- impt[impt$windowType == i, ]
  c(window = unique(tmp[,1]), type = unique(tmp[,2]), num = length(unique(tmp[,3])))
}))  

str(a)
window = colnames(auc)[1]
j = rownames(auc)[1]

numfea <- sapply(colnames(auc), function(window){
  fealen <- c(as.numeric(a[a[,1] == window, 3]), 2, 1)
  names(fealen) <- c(a[a[,1] == window, 2], 'base', 'cellProp')
  sapply(rownames(auc), function(i){
   tmp <- strsplit(i, '_') [[1]]
   tmp[tmp == 'rec'] <- 'receptor'
   tmp[tmp == 'cyt'] <- 'cytokine'
   if ('all' %in% tmp){
    return(sum(fealen)) 
   } else {
    return(sum(fealen[tmp]) )
   } 
  })
})
numfea <- numfea[order(rowSums(numfea)), ]
pdf(paste0(pdir, 'number_of_uniqueFeature_by_window.pdf'), width = 6, height = 5)
pheatmap(numfea, cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

write.csv(numfea, paste0(pdir, 'number_of_uniqueFeature_by_window.csv'))
