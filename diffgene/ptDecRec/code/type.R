library(parallel)
library(limma)
library(splines)
library(ggplot2)
library(mgcv)
library(reshape2)
meta <- readRDS('/home-4/zji4@jhu.edu/work-zfs/covid/data/current/meta.rds')
d <- readRDS('/home-4/zji4@jhu.edu/scratch/covid/info/pbmcnorm_combat.rds')
pt <- readRDS('/home-4/zji4@jhu.edu/scratch/covid/info/sample_info_pb_opt.rds')
pt <- pt[!pt$Patient %in% pt[which(pt$Clinical.Outcome=='Deceased' & pt$Pseudotime < 200),'Patient'],]
meta <- meta[meta[,'Library Sample Code'] %in% pt$Patient & meta[,'Clinical Outcome(COVID-19)'] %in% c('Recovered','Deceased') & meta$type %in% c('Mi','Mod','Se'),]
pb <- readRDS('/home-4/zji4@jhu.edu/work-zfs/covid/data/current/pbmc/pb_norm.rds')

stu <- sub('.*-','',meta[,'Library Sample Code'])
sel <- names(which(sapply(unique(stu),function(s) {
  tab <- table(meta[stu==s,'Clinical Outcome(COVID-19)'])
  length(tab) == 2 & min(tab) > 1
})))
meta <- meta[stu %in% sel,]

af <- list.files('/home-4/zji4@jhu.edu/scratch/covid/analysis/pbmc/diffgene/withincluster/deceasept/fullres')
type <- mclapply(af,function(f) {
  k <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/covid/analysis/pbmc/diffgene/withincluster/deceasept/fullres/',f))
  i <- sub('.rds','',f)
  m <- d[[i]]
  m <- m[,intersect(colnames(m),meta[,'Library Sample Code'])]
  selg <- rownames(k)[k[,'fdr'] < 0.05]
  
  m <- m[selg,]
  group <- meta[match(colnames(m),meta[,'Library Sample Code']),'Clinical Outcome(COVID-19)']
  names(group) <- colnames(m)
  p <- pt[match(colnames(m),pt$Patient),'Pseudotime']
  names(p) <- colnames(m)
  p <- sort(p)
  r1 <- range(p[names(which(group=='Deceased'))])
  r2 <- range(p[names(which(group=='Recovered'))])
  fit <- sapply(rownames(m),function(g) {
    set.seed(12345)
    m <- mgcv::gam(e~s(p,k=3,by=group),data=data.frame(e=m[g,],p=p[colnames(m)],group=group))
    pred1 <- predict(m,data.frame(p=seq(min(r1),max(r1),length.out = 1000),group='Deceased'))
    pred2 <- predict(m,data.frame(p=seq(min(r1),max(r1),length.out = 1000),group='Recovered'))
    #pred3 <- predict(m,data.frame(p=min(r2):max(r2),group='Recovered'))
    names(pred1) <- paste0('Deceased:',seq(min(r1),max(r1),length.out = 1000))
    names(pred2) <- paste0('Recovered:',seq(min(r1),max(r1),length.out = 1000))
    diff <- pred1-pred2
    names(diff) <- paste0('Difference:',sub('.*:','',names(diff)))
    #names(pred3) <- paste0('3:',min(r2):max(r2))
    #c(pred1,pred2,pred3)
    c(diff,pred1,pred2)
  })
  
  tmpgroup <- sub(':.*','',rownames(fit))
  for (u in unique(tmpgroup))
    for (j in 1:ncol(fit))
      fit[tmpgroup==u,j] <- scale(fit[tmpgroup==u,j])
  
  
  diffit <- t(fit[tmpgroup=='Difference',])
  zp <- apply(diffit[,1:(ncol(diffit)-1)]*diffit[,2:ncol(diffit)] < 0,1,which)
  zpnum <- sapply(zp,length)
  
  type <- rep(NA,length(zpnum))
  type[zpnum==1 & diffit[,1] > 0] <- 'D'
  type[zpnum==1 & diffit[,1] < 0] <- 'I'
  type[zpnum==2 & diffit[,1] > 0] <- 'DI'
  type[zpnum==2 & diffit[,1] < 0] <- 'ID'
  p1 <- sapply(zp,function(i) i[1])
  p2 <- sapply(zp,function(i) i[2])
  tab <- data.frame(type,p1,p2,stringsAsFactors = F)
  rownames(tab) <- names(zp)
  tab
},mc.cores=40)
names(type) <- sub('.rds','',af)
saveRDS(type,file='/home-4/zji4@jhu.edu/scratch/covid/analysis/pbmc/diffgene/withincluster/deceasept/fulltype/type.rds')

