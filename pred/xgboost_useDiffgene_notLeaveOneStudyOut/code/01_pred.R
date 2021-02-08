rm(list=ls())
library(here)
here()
setwd(here())
## conda_R/4.0.x
# within each cluster, combat, leave one study out
# use the training samples to calculate correlation between feasutures and severity
# xgboost predict 
# aim: see how much each cluster constribute to severity levels
library(data.table)
require(xgboost)  
require(Matrix)
require(data.table)
# if (!require('vcd')) install.packages('vcd')
library(sva)  ## not available (for R version 4.0.2)
library(pROC)
library(ggplot2)
library(RColorBrewer)
set.seed(12345)
## -------------------------
## combat to rm batch effect
## -------------------------
# d <- readRDS('/dcl02/hongkai/data/covid/data/current/pbmc/pb_norm.rds')
# for (i in 1:length(d)){
#   d[[i]] <- ComBat(d[[i]],sub('.*-','',colnames(d[[i]])))
# }
# saveRDS(d,'/dcl02/hongkai/data/whou/covid/pred/data/data/pbmcnorm_combat.rds')
d <- readRDS('pred/data/data/pbmcnorm_combat.rds')
# > sapply(d, ncol)
#  c1  c2  c3  c4  c5  c6  c7  c8  c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 
# 429 429 429 429 429 418 425 419 262 333 314  12 419 414 415 391 366 254  20 
# c20 c21 c22 c23 c24 c25 c26 c27 c28 c29 c30 c31 
#  79   6   3  14 428 421 352 262 158 406   5 359 
id <- which(sapply(d, ncol) > 100)
str(id)
d <- d[id]
str(d)
length(d)
## --------------------------
## read in candidate features
## --------------------------
meta <- readRDS('pred/data/data/meta.rds')
table(meta$type)
rownames(meta) <- meta$`Library Sample Code`
cellProp.mat <- readRDS('pred/data/prop/prop.rds')
str(cellProp.mat)

af <- list.files('pred/data/diffgene/reprolist/')
f <- af[1]
allgene <- list()
gene <- lapply(c('MiltHD.csv', 'MigtHD.csv'), function(f){
  tb <- read.csv(paste0('pred/data/diffgene/reprolist/', f))
  str(tb)
  paste0(gsub('c','cluster',tb$Cluster), ';', tb$Gene)
})
gene <- unlist(gene)
gene <- unique(gene)
str(gene)  
allgene[['HD_Mi']] <- gene
gene <- lapply(c('SeltMi.csv', 'SegtMi.csv'), function(f){
  tb <- read.csv(paste0('pred/data/diffgene/reprolist/', f))
  str(tb)
  paste0(gsub('c','cluster',tb$Cluster), ';', tb$Gene)
})
gene <- unlist(gene)
gene <- unique(gene)
str(gene)  
allgene[['Mi_Se']] <- gene
gene <- lapply(c('SeltHD.csv', 'SegtHD.csv'), function(f){
  tb <- read.csv(paste0('pred/data/diffgene/reprolist/', f))
  str(tb)
  paste0(gsub('c','cluster',tb$Cluster), ';', tb$Gene)
})
gene <- unlist(gene)
gene <- unique(gene)
str(gene)  
allgene[['HD_Se']] <- gene
str(allgene)

cyto <- readRDS('pred/data/data/cytokine_receptor.rds')
cytokine.v <- unique(as.character(cyto[,1]))
receptor.v <- unique(as.vector(apply(cyto[,2:ncol(cyto)], 2, as.character)))
tf.v <- readRDS('pred/data/TF/TFonly.rds')



## ----------------
## define functions
## ----------------
trans <- function(k) {
  a <- rep(0,length(k))
  a[k=='Mi'] <- 1
  a[k=='Mod'] <- 2
  a[k=='Se'] <- 3
  a[k=='Rec'] <- 4
  a[k=='Flu'] <- 5
  a
}
revtrans <- function(k){
  a <- rep('HD', length(k))
  a[k==1] <- 'Mi'
  a[k==2] <- 'Mod'
  a[k==3] <- 'Se'
  a[k==4] <- 'Rec'
  a[k==5] <- 'Flu'
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

## ---------------
## select features
## ---------------
unitype <- c('HD','Mi','Se')
typecom <- expand.grid(1:length(unitype), 1:length(unitype))
typecom <- typecom[typecom[,1] < typecom[,2], ]
typecom[,1] <- unitype[typecom[,1]]
typecom[,2] <- unitype[typecom[,2]]
# typecom <- typecom[-1,]
# for (arg in list(c('HD','Mi'), c('HD', 'Se'), c('Se', 'Mi'))){
for (arg in lapply(1:nrow(typecom), function(i) typecom[i, ])){
  print(arg)
  pbmc.s <- intersect(colnames(d[[1]]),meta[, 'Library Sample Code'][meta$type %in% arg])
  str(pbmc.s)
  as <- sub('.*-', '', pbmc.s)
  table(as)
  v <- NULL
  for (s in unique(as)){
    # v[s] <- length(table(meta[match(meta[, 'Library Sample Code'], pbmc.s[as == s]), 'type']))
    v[s] <- length(table(meta[pbmc.s[as == s], 'type']))
  }
  v
  pbmc.s <- pbmc.s[as %in% names(v[v == 2])]
  print(str(pbmc.s))
  
  as <- as[as %in% names(v[v == 2])]
  table(as)
  mat <- d[[1]][, colnames(d[[1]]) %in% pbmc.s]
  
  fea.base <- meta[match(colnames(mat), meta[, 'Library Sample Code']), c('Age', 'Sex')]
  rownames(fea.base) <- pbmc.s
  colnames(fea.base) <- paste0('base:', colnames(fea.base))
  fea.base[,2] <- ifelse(fea.base[,2] == 'M', 0, 10)
  
  # ----------------------
  ### for all sample: select features for genes when leaving a study
  gene <- allgene[[paste0(arg[1], '_', arg[2])]]
  str(gene)
  # gene.clu <- gsub('cluster', 'c', gsub(';.*', '', gene))
  gene.clu <- gsub(';.*', '', gene)
  str(gene.clu)
  gene.g <- gsub('.*;','',gene)
  str(gene.g)
  mat <- lapply(unique(gene.clu), function(c){
    print(c)
    tmp <- d[[sub('cluster','c',c)]] ####
    tmp <- tmp[gene.g, colnames(tmp) %in% pbmc.s]
    str(tmp)
    tmpmat <- matrix(0, nrow=nrow(tmp), ncol = length(pbmc.s))
    dimnames(tmpmat) <- list(rownames(tmp), pbmc.s)
    tmpmat[rownames(tmp), colnames(tmp)] <- tmp
    rownames(tmpmat) <- paste0(c, ';', rownames(tmpmat)) ########
    tmpmat
    
  })
  fealist.gene  <- do.call(rbind, mat)
  saveRDS(fealist.gene, paste0('pred/xgboost_useDiffgene/res/leave_one_study_out_selected_features_gene_', arg[1], '_', arg[2],'.rds'))
  
  # ------------------
  ## select features for cytokine 
  mat <- lapply(unique(gene.clu), function(c){
    print(c)
    tmp <- d[[sub('cluster','c',c)]] ####
    tmp <- tmp[intersect(gene.g, cytokine.v), colnames(tmp) %in% pbmc.s]
    tmpmat <- matrix(0, nrow=nrow(tmp), ncol = length(pbmc.s))
    dimnames(tmpmat) <- list(rownames(tmp), pbmc.s)
    tmpmat[rownames(tmp), colnames(tmp)] <- tmp
    rownames(tmpmat) <- paste0(c, ';', rownames(tmpmat)) ########
    tmpmat
    
  })
  fealist.cyt  <- do.call(rbind, mat)
  saveRDS(fealist.cyt, paste0('pred/xgboost_useDiffgene/res/leave_one_study_out_selected_features_cyt_', arg[1], '_', arg[2],'.rds'))
  
  # --------------------
  ## select features for receptor
  mat <- lapply(unique(gene.clu), function(c){
    print(c)
    tmp <- d[[sub('cluster','c',c)]] ####
    tmp <- tmp[intersect(gene.g, receptor.v), colnames(tmp) %in% pbmc.s]
    tmpmat <- matrix(0, nrow=nrow(tmp), ncol = length(pbmc.s))
    dimnames(tmpmat) <- list(rownames(tmp), pbmc.s)
    tmpmat[rownames(tmp), colnames(tmp)] <- tmp
    rownames(tmpmat) <- paste0(c, ';', rownames(tmpmat)) ########
    tmpmat
    
  })
  fealist.rec  <- do.call(rbind, mat)
  saveRDS(fealist.rec, paste0('pred/xgboost_useDiffgene/res/leave_one_study_out_selected_features_receptor_', arg[1], '_', arg[2],'.rds'))
  
  
  ##  select features for tf
  mat <- lapply(unique(gene.clu), function(c){
    print(c)
    tmp <- d[[sub('cluster','c',c)]] ####
    tmp <- tmp[intersect(gene.g, tf.v), colnames(tmp) %in% pbmc.s]
    tmpmat <- matrix(0, nrow=nrow(tmp), ncol = length(pbmc.s))
    dimnames(tmpmat) <- list(rownames(tmp), pbmc.s)
    tmpmat[rownames(tmp), colnames(tmp)] <- tmp
    rownames(tmpmat) <- paste0(c, ';', rownames(tmpmat)) ########
    tmpmat
    
  })
  fealist.tf  <- do.call(rbind, mat)
  saveRDS(fealist.tf, paste0('pred/xgboost_useDiffgene/res/leave_one_study_out_selected_features_tf_', arg[1], '_', arg[2],'.rds'))
  
  # -----------------------------
  ## combine selected features of base, gene, cytokine, receptor, tf, cellprop
  
    m1 <- fealist.gene
    rownames(m1) <- paste0('gene:', rownames(m1))
    m2 <- fealist.cyt
    rownames(m2) <- paste0('cytokine:', rownames(m2))
    m3 <- fealist.rec
    rownames(m3) <- paste0('receptor:', rownames(m3))
    m4 <- fealist.tf
    rownames(m4) <- paste0('tf:', rownames(m4))          ################3
    m <- rbind(m1, m2, m3, m4, t(fea.base), cellProp.mat[, pbmc.s]) 
    rownames(m)[nrow(m)] <- 'cellProp:'
    allfealist <- m
  
  # -------------------------------
  ###  predict using selected features
  ## base = TRUE; gene = FALSE; cellProp = FALSE; cytokine = F; receptor = T
  getauc <- function(base = TRUE,  cytokine = FALSE, receptor = FALSE, tf = FALSE, cellProp = FALSE, gene = FALSE){
    print(c(base, cytokine, receptor, tf, cellProp, gene))
    auclist <- list()
    for (j in 1:length(unique(as))) {
      print(j)
      leaveid <- which(as == unique(as)[j]) 
      trainid <- setdiff(1:length(pbmc.s), leaveid)
      y <- meta[pbmc.s[trainid], 'type']
      fea.type <- NULL
      if (base == TRUE)  fea.type <- c(fea.type, 'base')
      if (cytokine == TRUE)  fea.type <- c(fea.type, 'cytokine')
      if (receptor == TRUE)  fea.type <- c(fea.type, 'receptor')
      if (tf == TRUE)  fea.type <- c(fea.type, 'tf')
      if (cellProp == TRUE)  fea.type <- c(fea.type, 'cellProp')
      if (gene == TRUE)  {
        fea.type <- c(fea.type, 'gene')
        fea.type <- setdiff(fea.type, c('cytokine', 'receptor','tf'))
      } 
      
      mat <- allfealist
      g <- sub('.*;', '', rownames(mat))
      type <- sub(':.*', '', rownames(mat))
      print(table(type))
      mat <- mat[type %in% fea.type, , drop = FALSE]
      str(mat)
      model <- xgboost(data = t(mat[, trainid, drop=F]), label = trans(y), nrounds = 5, objective='multi:softprob',num_class = 6) ## nrounds to be 20!!
      pred <- predict(model,t(mat[,leaveid, drop=F]))
      pred <- matrix(pred,nrow=6)
      str(pred)
      id <- apply(pred,2,which.max) - 1
      predy <- revtrans(id)
      truey <- meta[match(colnames(mat)[leaveid], meta[, 'Library Sample Code']), 'type']
      auclist[[j]] <- auc(multiclass.roc(trans(truey), trans(predy), levels = as.vector(trans(c(arg[1], arg[2])))))
    }
    names(auclist) <- unique(as)  
    return(unlist(auclist))
  }
  
  base = F
  cytokine = FALSE
  receptor = FALSE
  tf = FALSE
  cellProp = T
  gene = FALSE
  allauc <- data.frame(base = getauc(T), 
                       cyt = getauc(F,T),
                       base_cyt = getauc(T,T),
                       rec = getauc(F,F,T),
                       base_rec = getauc(T,F,T),
                       tf = getauc(F,F,F,T),
                       base_tf = getauc(T,F,F,T),
                       cellProp = getauc(F,F,F,F,T),
                       base_cellProp = getauc(T,F,F,F,T),
                       gene = getauc(F,F,F,F,F,T),
                       cyt_rec = getauc(F,T,T),
                       base_cyt_rec = getauc(T,T,T),
                       cyt_rec_cellProp = getauc(F,T,T,F,T),
                       base_cyt_rec_cellProp = getauc(T,T,T,F,T),
                       cyt_rec_tf = getauc(F,T,T,T),
                       base_cyt_rec_tf = getauc(T,T,T,T),
                       cyt_rec_tf_cellProp = getauc(F,T,T,T,T,F),
                       base_cyt_rec_tf_cellProp = getauc(T,T,T,T,T,F),
                       base_gene = getauc(T,F,F,T,F,T),
                       base_gene_cellProp = getauc(T,F,F,T,F,T),
                       base_all = getauc(T,T,T,T,T,T))
  
  saveRDS(allauc, paste0('pred/xgboost_useDiffgene/res/AUC_', arg[1],'_', arg[2], '.rds'))
  allauc <- round(allauc,2)
  pd <- reshape2::melt(as.matrix(allauc))
  colnames(pd) <- c('Study', 'Type', 'AUC')
  pd[,3] <- round(pd[,3],2)
  
  pd$Type <- factor(as.character(pd$Type), levels = names(sort(tapply(pd$AUC, list(as.character(pd$Type)), median))))
  pdf(paste0('pred/xgboost_useDiffgene/plot/AUC_', arg[1],'_',arg[2], '_boxplot.pdf'), width = 8, height = 4) 
  print(ggplot(data = pd) +
    geom_boxplot(aes(x = Type, y = AUC), outlier.shape = NA) +
    geom_jitter(aes(x = Type, y = AUC, color = Study), width = 0.3, cex = 1) +
    scale_color_brewer(palette = 'Set1')+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('Feature'))
  dev.off()
}



#### summarize the auc across studies
# pd <- sapply(list(c('HD','Mi'), c('Se', 'Mi'), c('HD', 'Se')), function(arg){
mat <- sapply(lapply(1:nrow(typecom), function(i) typecom[i, ]), function(arg){  
  allauc <- readRDS(paste0('pred/xgboost_useDiffgene/res/AUC_', arg[1], '_', arg[2],'.rds'))
  colMeans(allauc)
})
colnames(mat) <-   c('HD_Mi',  'Se_Mi', 'HD_Se') #####
pd <- reshape2::melt(mat)
pd[,3] <- round(pd[,3], 2)
colnames(pd) <- c('Type', 'Comparison', 'AUC')
pdf('pred/xgboost_useDiffgene/plot/AUC_combine.pdf', width = 25, height = 3) 
print(ggplot(data = pd, aes(x = Comparison, y = AUC, fill = Type)) +
        geom_bar(stat="identity", position=position_dodge())+
        geom_text(aes(label=AUC), vjust=1.6, color="black", angle = 45,
                  position = position_dodge(0.9), size=3.5)+
        scale_fill_manual(values = colorRampPalette(brewer.pal(8,'Set1'))(length(unique(pd$Type))))+
        theme_minimal())

dev.off()

library(pheatmap)
pdf('pred/xgboost_useDiffgene/plot/AUC_combine_hm.pdf', width = 4.5, height = 4) 
pheatmap(mat, 
         cluster_cols = FALSE,
         border_color = NA)
dev.off()



############### --------------------
## heatmap
# # ---------------------------------
library(pheatmap)
library(RColorBrewer)
pt.arg <- readRDS('/dcl02/hongkai/data/rli/covid/multi_sample/data/backbone_123/mds_tscan/backbone123_info_aug.rds')
pt <- pt.arg$Pseudotime
names(pt) <- pt.arg$Patient
pbmc.s <- intersect(colnames(d[[1]]),meta[, 'Library Sample Code'][meta$type %in% c('HD','Se','Mi')])
pbmc.s <- intersect(names(pt), pbmc.s)
pt <- pt[pbmc.s]
pbmc.s <- names(pt)

for (arg in list(c('HD','Mi'), c('Se', 'Mi'))){
  fealist.gene <- readRDS(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost_useDiffgene/res/leave_one_study_out_selected_features_gene_', arg[1], '_', arg[2],'.rds'))
  fealist.cyt <- readRDS(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost_useDiffgene/res/leave_one_study_out_selected_features_cytokine_', arg[1], '_', arg[2],'.rds'))
  fealist.rec <- readRDS(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost_useDiffgene/res/leave_one_study_out_selected_features_receptor_', arg[1], '_', arg[2],'.rds'))
  
  for (type in c('cytokine', 'receptor', 'gene', 'tf')){
    print(type)
    if (type == 'cytokine'){
      fealist.select <- fealist.cyt
    } else if (type == 'receptor'){
      fealist.select <- fealist.rec
    } else if (type == 'gene'){
      fealist.select <- fealist.gene
    }
    
    rn <- unlist(sapply(fealist.select, rownames)) ###
    if (type == 'gene'){
      fea.select <- names(which(table(rn) >= 4))  
    } else {
      fea.select <- names(which(table(rn) == 5))
    }
    
    print(str(fea.select))
    
    fea.mat <- t(sapply(fea.select, function(i){
      int <- intersect(pbmc.s, colnames(d[[as.numeric(gsub(';.*', '', sub('cluster', '', i)))]]))
      tmp <- d[[as.numeric(gsub(';.*', '', sub('cluster', '', i)))]][sub('.*;', '', i), int]
      if (length(int) < length(pbmc.s)) tmp <- c(tmp, rep(tmp[length(tmp)], length(pbmc.s) - length(int))  )
      
      tmp <- loess(tmp~seq(1,length(tmp)))$fitted
      return(tmp)
    }))
    
    
    fea.mat <- fea.mat[names(sort(apply(fea.mat, 1, cor, seq(1, ncol(fea.mat))))),]
    
    fea.mat.ori <- t(sapply(fea.select, function(i){
      int <- intersect(pbmc.s, colnames(d[[as.numeric(gsub(';.*', '', sub('cluster', '', i)))]]))
      tmp <- d[[as.numeric(gsub(';.*', '', sub('cluster', '', i)))]][sub('.*;', '', i), int]
      if (length(int) < length(pbmc.s)) tmp <- c(tmp, rep(tmp[length(tmp)], length(pbmc.s) - length(int))  )
      return(tmp)
    }))
    
    fea.mat.ori <- fea.mat.ori[rownames(fea.mat), ]
    colnames(fea.mat) <- colnames(fea.mat.ori)
    
    cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    cpl = c(rep(cpl[1], 100), cpl, rep(cpl[length(cpl)], 100)) 
    
    # pheatmap(fea.mat, cluster_cols = F, cluster_rows = F, color = cpl, scale='row',
    #          labels_col = as.character(pt), show_rownames = F, show_colnames = F)
    
    h <- ifelse(type == 'cytokine', 18, ifelse(type == 'receptor', 20, 20))
    pdf(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost_useDiffgene/plot/hm_',arg[1], '_',arg[2], '_',type,'_rownames.pdf'), width = 7, height = h)
    pheatmap(fea.mat, cluster_cols = F, cluster_rows = F, color = cpl, scale='row',
             labels_col = pt, show_rownames = T, show_colnames = T,
             fontsize = 3, border_color = FALSE)
    dev.off()
    
    pdf(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost_useDiffgene/plot/hm_',arg[1], '_',arg[2], '_', type,'.pdf'), width = 3.2, height = 5)
    pheatmap(fea.mat, cluster_cols = F, cluster_rows = F, color = cpl, scale='row',
             labels_col = as.character(pt), show_rownames = F, show_colnames = F,
             fontsize = 10, border_color = FALSE)
    dev.off()
  }
  
}

## check fitted values appropriate or not
# library(ggplot2)
# pd <- reshape2::melt(fea.mat.ori[100:250,])
# pd[,2] <- pt[as.character(pd[,2])]
# colnames(pd) <- c('feature','pseudotime','value')
# pdf(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost_useDiffgene/plot/hm_',arg[1], '_',arg[2], '_', type,'_fitting_100_250.pdf'), width = 19.2, height = 12)
# ggplot(data = pd, aes(x = pseudotime, y = value)) +
#   geom_point(size = 0.2, alpha = 0.5) +
#   geom_smooth(method = 'loess', alpha = .15, fill = 'blue', color = 'blue', size = 0.5)+
#   facet_wrap(~feature, ncol=16) +
#   theme_minimal()
# dev.off()
# 
# pdf(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost_useDiffgene/plot/hm_',arg[1], '_',arg[2], '_', type,'_100_250.pdf'), width = 5, height = 7.5)
# pheatmap(fea.mat[100:250,], cluster_cols = F, cluster_rows = F, color = cpl, scale='row',
#          show_rownames = F, show_colnames = F,
#          fontsize = 10, border_color = FALSE)
# dev.off()
# 
# ## check geom_smooth loess and the fittied loess: they are the same
# colnames(fea.mat) <- colnames(fea.mat.ori)
# ld <- reshape2::melt(fea.mat[76:78,])
# ld[,2] <- pt[as.character(ld[,2])]
# colnames(ld) <- c('feature','pseudotime','value')
# ggplot() +
#   geom_point(data = pd, aes(x = pseudotime, y = value)) +
#   geom_line(data = ld, aes(x = pseudotime, y = value), method = 'loess', alpha = .15, fill = 'red', color = 'red')+
#   facet_grid(~feature)



## number of features: cluster specific. e.g. cluster1;POU2F2 and cluster10;POU2F2 are two features
ddir <- '/dcl02/hongkai/data/whou/covid/pred/xgboost_useDiffgene/res/'
allf <- list.files("/dcl02/hongkai/data/whou/covid/pred/xgboost_useDiffgene/res")
num <- lapply((allf[grepl('leave_one_study_out_selected_features',allf)]), function(f){
  print(f)
  a = readRDS(paste0(ddir, f))
  data.frame(study = names(a), number = sapply(a, nrow), type = sub('.rds','',sub('.*features_', '', f)), stringsAsFactors = F)
})
pd <- do.call(rbind, num)
pd <- pd[!grepl('features',pd[,3]),]
pd <- pd[!pd[,3] %in% c('HD_Mi','Se_Mi', 'HD_Se'),]
pd[,1] <- as.factor(pd[,1])
pd[,3] <- as.factor(pd[,3])

pdf('/dcl02/hongkai/data/whou/covid/pred/xgboost_useDiffgene/plot/number_of_features.pdf', width = 7, height = 4.5)
ggplot(data = pd) + 
  geom_boxplot(aes(x = type, y = number)) + 
  geom_jitter(aes(x = type, y = number, color = study)) +
  theme_classic() +
  scale_color_brewer(palette = 'Set1')+
  ylab('number of features') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## number of features: remove cluster information and then unique: e.g. POU2F2 is one feature
num <- lapply(setdiff(list.files("/dcl02/hongkai/data/whou/covid/pred/data/data"),'pbmcnorm_combat.rds'), function(f){
  print(f)
  a = readRDS(paste0(ddir, f))
  v <- sapply(a, function(i){
    length(unique(sub('.*;', '', rownames(i))))
  })
  data.frame(study = names(a), number = v, type = sub('.rds','',sub('.*features_', '', f)), stringsAsFactors = F)
})
pd <- do.call(rbind, num)
pd <- pd[!grepl('features',pd[,3]),]
pd <- pd[!pd[,3] %in% c('HD_Mi','Se_Mi', 'HD_Se'),]
pd[,1] <- as.factor(pd[,1])
pd[,3] <- as.factor(pd[,3])

pdf('/dcl02/hongkai/data/whou/covid/pred/xgboost_useDiffgene/plot/number_of_unique_features.pdf', width = 7, height = 4.5)
ggplot(data = pd) + 
  geom_boxplot(aes(x = type, y = number)) + 
  geom_jitter(aes(x = type, y = number, color = study)) +
  theme_classic() +
  scale_color_brewer(palette = 'Set1')+
  ylab('number of unique features') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## check the intersect of the unique features of cytokines, receptors
fea.int <- lapply(setdiff(list.files("/dcl02/hongkai/data/whou/covid/pred/data/data"),'pbmcnorm_combat.rds'), function(f){
  print(f)
  a = readRDS(paste0(ddir, f))
  unique.feature <- lapply(a, function(i){
    unique(sub('.*;', '', rownames(i)))
  })
  names(which(table(unlist(unique.feature)) == length(unique.feature)))
})
pd <- do.call(rbind, num)
pd <- pd[!grepl('features',pd[,3]),]
pd <- pd[!pd[,3] %in% c('HD_Mi','Se_Mi', 'HD_Se'),]
pd[,1] <- as.factor(pd[,1])
pd[,3] <- as.factor(pd[,3])












