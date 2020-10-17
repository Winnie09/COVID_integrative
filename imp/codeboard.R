# normalize
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

# ---------
# SAVER impute 
# ---------
library(SAVER)
library(Matrix)
library(parallel)
packageVersion('SAVER')
setwd('/dcl02/hongkai/data/whou/covid/imp/data/norm/')
rdir <- '/dcl02/hongkai/data/whou/covid/imp/res/'
f <- as.character(commandArgs(trailingOnly = T)[[1]])
print(f)
d <- readRDS(f)
d.saver <- saver(d, ncores = detectCores(), size.factor = 1)
saveRDS(d.saver, paste0(rdir, f))

#
ml R
Rscript 02_imputation.R $1

#
for i in `ls /dcl02/hongkai/data/whou/covid/imp/data/norm/`
do
qsub run02.sh $i
done

# rerun uncompleted samples
getf = sub('.rds', '', list.files('/dcl02/hongkai/data/whou/covid/imp/res/'))
af <- list.files('/dcl02/hongkai/data/whou/covid/imp/data/sample/')
runf <- paste0(setdiff(af, getf), '.rds')

fileConn<-file("/dcl02/hongkai/data/whou/covid/imp/code/runsample.txt")
writeLines(runf, fileConn)
close(fileConn)

# ------------
# saverx impute
# ------------
##
library(SAVERX)
library(Matrix)
packageVersion('SAVERX')
setwd('/dcl02/hongkai/data/whou/covid/imp/data/norm/')
rdir <- '/dcl02/hongkai/data/whou/covid/imp/saverx/res/'

d.saverx <- saverx(commandArgs(trailingOnly = T)[[1]], data.species = 'Human')


##
af = list.files('/dcl02/hongkai/data/whou/covid/imp/data/norm/')
for (f in af){
  print(f)
  rdir = paste0('/dcl02/hongkai/data/whou/covid/imp/data/', sub('.rds', '', f), '/')
  dir.create(rdir, recursive = T, showWarnings = F)
  system(paste0('cp /dcl02/hongkai/data/whou/covid/imp/data/norm/', f, ' ', rdir, 'norm.rds'))
}

# ------------
# MAGIC impute
# ------------
ddir <- '/dcl02/hongkai/data/whou/covid/imp/data/'
af = list.files(ddir)

dlist = list()
for (f in af){
  dlist[[f]] <- readRDS(paste0(ddir, f, '/norm.rds'))
}
d <- do.call(cbind, dlist)
dir.create('/dcl02/hongkai/data/whou/covid/imp/data/all/', recursive = T)
saveRDS(d, '/dcl02/hongkai/data/whou/covid/imp/data/all/norm.rds')

# run MAGIC
library(Rmagic)
library(Matrix)
ddir <- '/dcl02/hongkai/data/whou/covid/imp/data/sample/'
f <- as.character(commandArgs(trailingOnly = T)[1])
d <- readRDS(paste0(ddir, f, '/norm.rds'))
d <- log2(d + 1)
d.magic <- magic(d)
saveRDS(d.magic, paste0('/dcl02/hongkai/data/whou/covid/imp/magic/sample/', f, '.rds'))

use_python("/jhpce/shared/jhpce/core/python/2.7.9/bin/python")








