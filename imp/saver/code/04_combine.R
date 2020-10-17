a <- readRDS('/dcl02/hongkai/data/whou/covid/imp/saver/res/all/log2norm.rds')
b = do.call(cbind, a)
saveRDS(b, '/dcl02/hongkai/data/whou/covid/imp/saver/res/all/log2norm_matrix.rds')
