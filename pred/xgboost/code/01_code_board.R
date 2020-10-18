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
  ### for all sample: select features when leeaving a study
  fealist <- list()
  for (j in 1:length(unique(as))) {
    print(j)
    fea <- sapply(setdiff(1:length(d), c(28, 30)), function(cluid){
      print(cluid)
      mat <- d[[cluid]]
      mat <- mat[, colnames(mat) %in% pbmc.s]
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
  saveRDS(fealist, paste0('/dcl02/hongkai/data/whou/covid/pred/data/data/leave_one_study_out_selected_features_', arg[1], '_', arg[2],'.rds'))
  
  
  
  allfealist <- sapply(1:length(fealist), function(i){
    m1 <- fealist[[i]]
    rownames(m1) <- paste0('gene:', rownames(m1))
    m <- rbind(m1, t(fea.base))
    m <- rbind(m, cellProp.mat[i, pbmc.s])
    rownames(m)[nrow(m)] <- 'cellProp:'
    m
  })
  names(allfealist) <- names(fealist)
  
  # -----------------------------
  
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
      if ('cytokine' %in% fea.type & !'receptor' %in% fea.type){
        g.select <- intersect(g, cytokine.v)
      } else if (!'cytokine' %in% fea.type & 'receptor' %in% fea.type){
        g.select <- intersect(g, receptor.v)
      } else if (sum(c('cytokine', 'receptor') %in% fea.type) == 2){
        g.select <- intersect(g, unique(c(cytokine.v, receptor.v)))
      } else if ('gene' %in% fea.type){
        g.select <- g
      } else {
        g.select <- NULL
      }
      
      if (length(fea.type %in% c('cytokine', 'receptor')) > 0) fea.type <- setdiff(fea.type, c('cytokine', 'receptor'))
      mat <- mat[g %in% g.select | type %in% fea.type, ]
      
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







