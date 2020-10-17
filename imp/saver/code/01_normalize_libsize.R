library(Matrix)
library(parallel)
setwd('/dcl02/hongkai/data/covid/data/200916/')
rdir <- '/dcl02/hongkai/data/whou/covid/imp/data/norm/'

pbmc <- readRDS('pbmc.rds')
ap <- sub(':.*', '', colnames(pbmc))
libsize <- colSums(pbmc)
libsize <- libsize/median(libsize)

nn <- sapply(unique(ap), function(p){
  print(p)
  tmp <- pbmc[, ap == p, drop = FALSE]
  tmp <- sweep(tmp, 2, libsize[ap == p], '/')
  saveRDS(tmp, paste0(rdir, p, '.rds'))
  return(0)
})




