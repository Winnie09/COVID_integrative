setwd('/dcl02/hongkai/data/whou/covid/imp/saver/res/sample')
af = list.files(getwd())
am <- sapply(af, function(f){
  tmp <- readRDS(f)$estimate
  tmp <- log2(tmp + 1)
})
saveRDS(am, '/dcl02/hongkai/data/whou/covid/imp/saver/res/all/log2norm.rds')
