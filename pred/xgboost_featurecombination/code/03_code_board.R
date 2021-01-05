## conda_R/4.0.x
# within each cluster, combat, leave one study out
# use the training samples to calculate correlation between feasutures and severity
# xgboost predict 
# aim: see how much each cluster constribute to severity levels
require(xgboost)  
require(Matrix)
require(data.table)
if (!require('vcd')) install.packages('vcd')
library(sva)
library(pROC)
library(ggplot2)
library(RColorBrewer)
set.seed(12345)
## combat to rm batch effect
# d <- readRDS('/dcs01/gcode/zji4/covid/analysis/pbmc/diffgene/pb/res/pbmcnorm.rds')
# for (i in 1:length(d)){
#   d[[i]] <- ComBat(d[[i]],sub('.*-','',colnames(d[[i]])))
# }
# saveRDS(d,'/dcl02/hongkai/data/whou/covid/pred/data/data/pbmcnorm_combat.rds')
d <- readRDS('/dcl02/hongkai/data/whou/covid/pred/data/data/pbmcnorm_combat.rds')
meta <- readRDS('/dcl02/hongkai/data/covid/data/200916/meta.rds')
cellProp.mat <- readRDS('/dcs01/gcode/zji4/covid/old/analysis/pbmc/cluprop/fullprop/pbmc.rds')
cyto <- readRDS('/dcl02/hongkai/data/covid/data/200916/cytokine_receptor.rds')
cytokine.v <- unique(as.character(cyto[,1]))
receptor.v <- unique(as.vector(apply(cyto[,2:ncol(cyto)], 2, as.character)))

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


for (arg in list(c('HD','Mi'), c('HD', 'Se'), c('Se', 'Mi'))){
  pbmc.s <- intersect(colnames(d[[1]]),meta[, 'Library Sample Code'][meta$type %in% arg]) ###########
  as <- sub('.*-', '', pbmc.s)
  v <- NULL
  for (s in unique(as)){
    v[s] <- length(table(meta[match(pbmc.s[as == s], meta[, 'Library Sample Code']), 'type']))
  }
  pbmc.s <- pbmc.s[as %in% names(v[v == 2])]
  print(str(pbmc.s))
  
  as <- as[as %in% names(v[v == 2])]
  
  mat <- d[[1]][, colnames(d[[1]]) %in% pbmc.s]
  fea.base <- meta[match(colnames(mat), meta[, 'Library Sample Code']), c('Age', 'Sex')]
  rownames(fea.base) <- pbmc.s
  colnames(fea.base) <- paste0('base:', colnames(fea.base))
  fea.base[,2] <- ifelse(fea.base[,2] == 'M', 0, 10)
  
  # ---------------
  # for only HD, Se sample
  ### for all sample: select features for genes when leeaving a study
  get_fealist <- function(feature.type = 'gene'){
    fealist <- list()
    for (j in 1:length(unique(as))) {
      print(j)
      fea <- sapply(setdiff(1:length(d), c(28, 30)), function(cluid){
        print(cluid)
        mat <- d[[cluid]]
        if (feature.type == 'gene'){
          mat <- mat[, colnames(mat) %in% pbmc.s]
        } else if (feature.type == 'cytokine'){
          mat <- mat[rownames(mat) %in% cytokine.v, colnames(mat) %in% pbmc.s]
        } else if (feature.type == 'receptor'){
          mat <- mat[rownames(mat) %in% receptor.v, colnames(mat) %in% pbmc.s]
        }
        leaveid <- which(as == unique(as)[j]) 
        trainid <- setdiff(1:ncol(mat), leaveid)
        tryCatch({
          y <- meta[match(colnames(mat)[trainid], meta[, 'Library Sample Code']), 'type']
          model <- xgboost(data = t(mat[, trainid, drop=FALSE]), label = trans(y), nrounds = 20,objective='multi:softprob',num_class=3)  
          impt <- data.frame(xgb.importance(feature_names = rownames(mat), model = model))
          tmp <- mat[impt[,1],]
          tmpmat <- matrix(0, nrow=nrow(tmp), ncol = length(pbmc.s))
          dimnames(tmpmat) <- list(rownames(tmp), pbmc.s)
          tmpmat[rownames(tmp), colnames(tmp)] <- tmp
          rownames(tmpmat) <- paste0('cluster', cluid, ';', rownames(tmpmat))
          return(tmpmat)
        },warning=function(w){},error=function(e){})
      })
      fealist[[j]] <-  do.call(rbind, fea)
    }
    ##
    names(fealist) <- unique(as)  
    return(fealist)
  }
  
  fealist.gene <- get_fealist(feature.type = 'gene')
  saveRDS(fealist.gene, paste0('/dcl02/hongkai/data/whou/covid/pred/data/data/leave_one_study_out_selected_features_gene_', arg[1], '_', arg[2],'.rds'))
  
  # ------------------
  ## select features for cytokine 
  fealist.cyt <- get_fealist(feature.type = 'cytokine')
  saveRDS(fealist.cyt, paste0('/dcl02/hongkai/data/whou/covid/pred/data/data/leave_one_study_out_selected_features_cytokine_', arg[1], '_', arg[2],'.rds'))
  
  # --------------------
  ## select features for receptor
  fealist.rec <- get_fealist(feature.type = 'receptor')
  saveRDS(fealist.rec, paste0('/dcl02/hongkai/data/whou/covid/pred/data/data/leave_one_study_out_selected_features_receptor_', arg[1], '_', arg[2],'.rds'))
  
  
  # -----------------------------
  ## combine selected features of base, gene, cytokine, receptor
  allfealist <- sapply(1:length(fealist.gene), function(i){
    m1 <- fealist.gene[[i]]
    rownames(m1) <- paste0('gene:', rownames(m1))
    m2 <- fealist.cyt[[i]]
    rownames(m2) <- paste0('cytokine:', rownames(m2))
    m3 <- fealist.rec[[i]]
    rownames(m3) <- paste0('receptor:', rownames(m3))
    m <- rbind(m1, m2, m3, t(fea.base), cellProp.mat[i, pbmc.s]) 
    rownames(m)[nrow(m)] <- 'cellProp:'
    m
  })
  names(allfealist) <- names(fealist.gene)
  
  # -------------------------------
  ###  predict using selected features, make contigency table 
  ## base = TRUE; gene = FALSE; cellProp = FALSE; cytokine = F; receptor = T
  getauc <- function(base = TRUE,  cytokine = FALSE, receptor = FALSE,  cellProp = FALSE, gene = FALSE){
    auclist <- list()
    for (j in 1:length(unique(as))) {
      leaveid <- which(as == unique(as)[j]) 
      trainid <- setdiff(1:length(pbmc.s), leaveid)
      y <- meta[match(pbmc.s[trainid], meta[, 'Library Sample Code']), 'type']
      fea.type <- NULL
      if (base == TRUE)  fea.type <- c(fea.type, 'base')
      if (cytokine == TRUE)  fea.type <- c(fea.type, 'cytokine')
      if (receptor == TRUE)  fea.type <- c(fea.type, 'receptor')
      if (cellProp == TRUE)  fea.type <- c(fea.type, 'cellProp')
      if (gene == TRUE)  {
        fea.type <- c(fea.type, 'gene')
        fea.type <- setdiff(fea.type, c('cytokine', 'receptor'))
      } 
      
      mat <- allfealist[[j]]
      g <- sub('.*;', '', rownames(mat))
      type <- sub(':.*', '', rownames(mat))
      print(table(type))
      mat <- mat[type %in% fea.type, , drop = FALSE]
      invisible(capture.output(model <- xgboost(data = t(mat[, trainid]), label = trans(y), nrounds = 20, objective='multi:softprob',num_class = 3)))
      pred <- predict(model,t(mat[,leaveid]))
      pred <- matrix(pred,nrow=3)
      id <- apply(pred,2,which.max) - 1
      predy <- revtrans(id)
      truey <- meta[match(colnames(mat)[leaveid], meta[, 'Library Sample Code']), 'type']
      invisible(capture.output(auclist[[j]] <- auc(multiclass.roc(trans(truey), trans(predy), levels = as.vector(trans(c(arg[1], arg[2])))))))
    }
    names(auclist) <- unique(as)  
    return(unlist(auclist))
  }
  
  allauc <- data.frame(base = getauc(T), 
                       base_cyt = getauc(T,T),
                       base_rec = getauc(T,F,T),
                       base_cyt_rec = getauc(T,T,T),
                       base_cellProp = getauc(T,F,F,T),
                       base_cyt_rec_cellProp = getauc(T,T,T,T),
                       base_gene = getauc(T,F,F,T,T),
                       base_gene_cellProp = getauc(T,F,F,T,T),
                       base_all = getauc(T,T,T,T,T))
  allacu <- round(allauc,2)
  saveRDS(allauc, paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost/res/AUC_', arg[1],'_', arg[2], '.rds'))
  pd <- reshape2::melt(as.matrix(allauc))
  colnames(pd) <- c('Study', 'Type', 'AUC')
  pd[,3] <- round(pd[,3],2)
  
  pdf(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost/plot/AUC_', arg[1],'_',arg[2], '.pdf'), width = 18, height = 5) 
  print(ggplot(data = pd, aes(x = Study, y = AUC, fill = Type)) +
          geom_bar(stat="identity", position=position_dodge())+
          geom_text(aes(label=AUC), vjust=1.6, color="black",
                    position = position_dodge(0.9), size=3.5)+
          scale_fill_brewer(palette="Pastel1", direction = -1)+
          theme_minimal())
  dev.off()
}


#### summarize the auc across studies
pd <- sapply(list(c('HD','Mi'), c('Se', 'Mi'), c('HD', 'Se')), function(arg){
  allauc <- readRDS(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost/res/AUC_', arg[1], '_', arg[2],'.rds'))
  colMeans(allauc)
})
colnames(pd) <-   c('HD_Mi',  'Se_Mi', 'HD_Se')
pd <- reshape2::melt(pd)
pd[,3] <- round(pd[,3], 2)
colnames(pd) <- c('Type', 'Comparison', 'AUC')
pdf('/dcl02/hongkai/data/whou/covid/pred/xgboost/plot/AUC_combine.pdf', width = 18, height = 5) 
print(ggplot(data = pd, aes(x = Comparison, y = AUC, fill = Type)) +
        geom_bar(stat="identity", position=position_dodge())+
        geom_text(aes(label=AUC), vjust=1.6, color="black",
                  position = position_dodge(0.9), size=3.5)+
        scale_fill_brewer(palette="Pastel1", direction = -1)+
        theme_minimal())

dev.off()


############### --------------------
## heatmap
# # ---------------------------------
library(pheatmap)
library(RColorBrewer)

pt.arg <- readRDS('/dcl02/hongkai/data/rli/covid/multi_sample/data/backbone_123/mds_tscan/backbone123_info_aug.rds')
pt <- pt.arg$Pseudotime
names(pt) <- pt.arg$Patient
pbmc.s <- intersect(colnames(d[[1]]),meta[, 'Library Sample Code'][meta$type %in% c('HD','Se','Mi')]) pbmc.s <- intersect(names(pt), pbmc.s)
pt <- pt[pbmc.s]
pbmc.s <- names(pt)

for (arg in list(c('HD','Mi'), c('Se', 'Mi'))){
  fealist.gene <- readRDS(paste0('/dcl02/hongkai/data/whou/covid/pred/data/data/leave_one_study_out_selected_features_gene_', arg[1], '_', arg[2],'.rds'))
  fealist.cyt <- readRDS(paste0('/dcl02/hongkai/data/whou/covid/pred/data/data/leave_one_study_out_selected_features_cytokine_', arg[1], '_', arg[2],'.rds'))
  fealist.rec <- readRDS(paste0('/dcl02/hongkai/data/whou/covid/pred/data/data/leave_one_study_out_selected_features_receptor_', arg[1], '_', arg[2],'.rds'))
  
  for (type in c('cytokine', 'receptor', 'gene')){
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
    pdf(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost/plot/hm_',arg[1], '_',arg[2], '_',type,'_rownames.pdf'), width = 7, height = h)
    pheatmap(fea.mat, cluster_cols = F, cluster_rows = F, color = cpl, scale='row',
             labels_col = pt, show_rownames = T, show_colnames = T,
             fontsize = 3, border_color = FALSE)
    dev.off()
    
    pdf(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost/plot/hm_',arg[1], '_',arg[2], '_', type,'.pdf'), width = 3.2, height = 5)
    pheatmap(fea.mat, cluster_cols = F, cluster_rows = F, color = cpl, scale='row',
             labels_col = as.character(pt), show_rownames = F, show_colnames = F,
             fontsize = 10, border_color = FALSE)
    dev.off()
  }
  
}



## check fitted values appropriate or not
library(ggplot2)
pd <- reshape2::melt(fea.mat.ori[100:250,])
pd[,2] <- pt[as.character(pd[,2])]
colnames(pd) <- c('feature','pseudotime','value')
pdf(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost/plot/hm_',arg[1], '_',arg[2], '_', type,'_fitting_100_250.pdf'), width = 19.2, height = 12)
ggplot(data = pd, aes(x = pseudotime, y = value)) +
  geom_point(size = 0.2, alpha = 0.5) +
  geom_smooth(method = 'loess', alpha = .15, fill = 'blue', color = 'blue', size = 0.5)+
  facet_wrap(~feature, ncol=16) +
  theme_minimal()
dev.off()

pdf(paste0('/dcl02/hongkai/data/whou/covid/pred/xgboost/plot/hm_',arg[1], '_',arg[2], '_', type,'_100_250.pdf'), width = 5, height = 7.5)
pheatmap(fea.mat[100:250,], cluster_cols = F, cluster_rows = F, color = cpl, scale='row',
            show_rownames = F, show_colnames = F,
             fontsize = 10, border_color = FALSE)
dev.off()

## check geom_smooth loess and the fittied loess: they are the same
colnames(fea.mat) <- colnames(fea.mat.ori)
ld <- reshape2::melt(fea.mat[76:78,])
ld[,2] <- pt[as.character(ld[,2])]
colnames(ld) <- c('feature','pseudotime','value')
ggplot() +
  geom_point(data = pd, aes(x = pseudotime, y = value)) +
  geom_line(data = ld, aes(x = pseudotime, y = value), method = 'loess', alpha = .15, fill = 'red', color = 'red')+
  facet_grid(~feature)

