library(parallel)
library(limma)
library(splines)
library(ggplot2)
library(mgcv)
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

cid <- as.numeric(commandArgs(trailingOnly = T))
i <- names(d)[cid]

print(i)
m <- d[[i]]
m <- m[,intersect(colnames(m),meta[,'Library Sample Code'])]
group <- meta[match(colnames(m),meta[,'Library Sample Code']),'Clinical Outcome(COVID-19)']

tpb <- pb[[i]][,colnames(m)]
m <- m[names(which(rowMeans(tpb > 0) > 0.1)),]

if (length(group) > 0 && min(table(group)) >= 10) {
  p <- pt[match(colnames(m),pt$Patient),'Pseudotime']

  pval <- mclapply(rownames(m),function(g) {
    set.seed(12345)
    reald <- mgcv::gam(e~s(p,k=3,by=group),data=data.frame(e=m[g,],p=p,group=group))$deviance
    permud <- sapply(1:1000,function(i) mgcv::gam(e~s(p,k=3,by=group),data=data.frame(e=m[g,],p=p,group=sample(group)))$deviance)
    den <- density(permud,bw='SJ')$bw
    c(mean(pnorm(reald, permud, sd=den)),(mean(permud)-reald)/sd(permud))
  },mc.cores=detectCores())
  names(pval) <- rownames(m)
  pval <- do.call(rbind,pval)
  pval <- data.frame(pval)
  colnames(pval) <- c('pvalue','zscore')
  pval$fdr <- p.adjust(pval$pvalue,method='fdr')
  saveRDS(pval,file=paste0('/home-4/zji4@jhu.edu/scratch/covid/analysis/pbmc/diffgene/withincluster/deceasept/res/',i,'.rds'))
}



