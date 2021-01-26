rm(list=ls())
setwd('/Users/wenpinhou/Dropbox/covid/')
rdir <- 'pred/xgboost_survivalDynamic/result/'
pdir <- 'pred/xgboost_survivalDynamic/plot/'
d.pt <- readRDS('/Users/wenpinhou/Dropbox/covid/pred/data/data/sample_info_pb_opt.rds')
rownames(d.pt) <- d.pt$Patient
str(d.pt)
d.pt <- d.pt[d.pt[,'Clinical.Outcome'] %in% c('Recovered','Deceased'),]
d.pt <- d.pt[d.pt[,'Phenotype'] %in% c('Mi', 'Mod', 'Se'),]

stu <- d.pt$Study
sel <- names(which(sapply(unique(stu),function(s) {
  tab <- table(d.pt[stu==s,'Clinical.Outcome'])
  length(tab) == 2 & min(tab) > 1
})))
d.pt <- d.pt[stu %in% sel,]
pt <- d.pt$Pseudotime
qt <- quantile(pt,seq(0,1,length.out=21))

deceased.id = which(d.pt$Clinical.Outcome == 'Deceased')
windowcut <- pt[deceased.id[1:length(deceased.id)%%5 == 0]]
windowcut[length(windowcut)] <- max(pt)
windowcut <- c(deceased.id[1]-1, windowcut)



sp.split <- sapply(1:(length(windowcut)-2),function(i) {
 d.pt[pt > windowcut[i] & pt <= windowcut[i+2],'Patient']
  
})
names(sp.split) <- paste0('window', 1:length(sp.split))
cut <- d.pt$Pseudotime[which(d.pt$Clinical.Outcome=='Deceased')[3]]
sp.split[[1]] <- setdiff(sp.split[[1]],d.pt$Patient[d.pt$Pseudotime < cut])

type.split <- sapply(1:(length(sp.split)),function(i) {
  tab <- table(d.pt[sp.split[[i]],'Clinical.Outcome'])
  v <- rep(0,2)
  names(v) <- c('Deceased', 'Recovered')
  v[names(tab)] <- tab
  v
})
colnames(type.split) <- paste0('window', 1:ncol(type.split))
write.csv(type.split, paste0(pdir, 'window_sample_size.csv'))



# type.split <- sapply(1:(length(qt)-5),function(i) {
#   tab <- table(d.pt[pt > qt[i] & pt < qt[i+5],'Clinical.Outcome'])
#   v <- rep(0,2)
#   names(v) <- c('Deceased', 'Recovered')
#   v[names(tab)] <- tab
#   v
# })
# colnames(type.split) <- paste0('window', 1:ncol(type.split))
# write.csv(type.split, paste0(pdir, 'window_sample_size.csv'))
# 
# sp.split <- pat <- sapply(1:(length(qt)-5),function(i) {
#   d.pt[pt > qt[i] & pt < qt[i+5],'Patient']
# },simplify = F)

## pseudobulk data
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
# meta <- readRDS('pred/data/data/meta.rds')
cellProp.mat <- readRDS('pred/data/prop/prop.rds')
str(cellProp.mat)
cyto <- readRDS('pred/data/data/cytokine_receptor.rds')
str(cyto)
cytokine.v <- unique(as.character(cyto[,1]))
receptor.v <- unique(as.vector(apply(cyto[,2:ncol(cyto)], 2, as.character)))
str(cytokine.v)
str(receptor.v)
tf.v <- readRDS('pred/data/TF/TFonly.rds')
str(tf.v)
## ----------------
## define functions: Deceased = 0, Recovered = 1
## ----------------
trans <- function(k) {
  a <- rep(0,length(k))
  a[k=='Recovered'] <- 1
  a
}
revtrans <- function(k){
  a <- rep('Deceased', length(k))
  a[k==1] <- 'Recovered'
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


for (window in 1:length(sp.split)){
  print(window)
  sp.select <- sp.split[[window]]
  arg <- c('Recovered', 'Deceased')
  pbmc.s <- intersect(colnames(d[[1]]), sp.select)
  as <- sub('.*-', '', pbmc.s)
  table(as)
  pbmc.patient <- sub('\\..*-', '-', d.pt[pbmc.s, 'Patient'])
  length(unique(pbmc.patient))
  mat <- d[[1]][, pbmc.s]
  fea.base <- d.pt[pbmc.s, c('Age', 'Sex')]
  rownames(fea.base) <- pbmc.s
  colnames(fea.base) <- paste0('base:', colnames(fea.base))
  fea.base[,2] <- ifelse(fea.base[,2] == 'M', 0, 10)
  outcome <- d.pt[pbmc.s, 'Clinical.Outcome']      
  id1 <- which(outcome == 'Deceased')
  id1.patient <- pbmc.patient[id1]
  id1.patient <- id1[!duplicated(id1.patient)]
  str(id1.patient)
  id2 <- which(outcome == 'Recovered')
  id2.patient <- pbmc.patient[id2]
  id2.patient <- id2[!duplicated(id2.patient)]
  str(id2.patient)
  get_fealist <- function(feature.type = 'gene'){
    fealist  <- list()
    for (j in 1:5) {
       set.seed(j)
       if (length(id1.patient)<4){
         leaveid <- c(sample(pbmc.s[id1.patient], length(id1.patient)/2), sample(pbmc.s[id2.patient], length(id2.patient)/4))
       } else {
         leaveid <- c(sample(pbmc.s[id1.patient], length(id1.patient)/4), sample(pbmc.s[id2.patient], length(id2.patient)/4))
       }
       trainid <- setdiff(pbmc.s, leaveid)
       if (length(table(d.pt[leaveid,'Clinical.Outcome'])) >= 1 & length(table(d.pt[trainid,'Clinical.Outcome'])) == 2){
       data <- lapply(1:length(d), function(cluid){
          mat <- d[[cluid]]
          if (feature.type == 'gene'){
            mat <- mat[, colnames(mat) %in% pbmc.s, drop = F]
          } else if (feature.type == 'cytokine'){
            mat <- mat[rownames(mat) %in% cytokine.v, colnames(mat) %in% pbmc.s, drop=F]
          } else if (feature.type == 'receptor'){
            mat <- mat[rownames(mat) %in% receptor.v, colnames(mat) %in% pbmc.s, drop=F]
          } else if (feature.type == 'tf'){
            mat <- mat[rownames(mat) %in% tf.v, colnames(mat) %in% pbmc.s, drop=F]
          }
          s.int <- intersect(colnames(mat), pbmc.s)
          mat <- mat[, s.int]
          tmpmat <- matrix(NA, nrow=nrow(mat), ncol = length(pbmc.s))  ## missing as NA !!
          dimnames(tmpmat) <- list(rownames(mat), pbmc.s)
          tmpmat[rownames(mat), colnames(mat)] <- mat
          rownames(tmpmat) <- paste0('cluster', sub('c','',names(d)[cluid]), ';', rownames(tmpmat))
          tmpmat
       })
         mat <- do.call(rbind, data)
            tryCatch({
              y <- d.pt[trainid, 'Clinical.Outcome']
              model <- xgboost(data = t(mat[, trainid, drop=FALSE]), label = trans(y), nrounds = 20,objective='binary:logistic', verbose = F)  ###### nrounds to be 20 !!! numclu = 2
              impt <- data.frame(xgb.importance(feature_names = rownames(mat), model = model))
              saveRDS(impt, paste0(rdir, 'impt_window', window, '_', feature.type, '_seed', j, '.rds'))
              tmp <- mat[impt[,1], , drop=F]
            },warning=function(w){},error=function(e){})
          if (exists('tmp')){
            fealist[[j]] <-  tmp
          } else {
            fealist[[j]] <- NULL
          }
       }
    }
    return(fealist)
}

  getauc <- function(base = TRUE,  cytokine = FALSE, receptor = FALSE, tf = FALSE, cellProp = FALSE, gene = FALSE){
    print(c(base, cytokine, receptor, tf, cellProp, gene))
    auclist <- list()
    for (j in 1:5) {
      set.seed(j)
      if (length(id1.patient)<4){
         leaveid <- c(sample(pbmc.s[id1.patient], length(id1.patient)/2), sample(pbmc.s[id2.patient], length(id2.patient)/4))
       } else {
         leaveid <- c(sample(pbmc.s[id1.patient], length(id1.patient)/4), sample(pbmc.s[id2.patient], length(id2.patient)/4))
       }
      trainid <- setdiff(pbmc.s, leaveid)
      y <- d.pt[trainid, 'Clinical.Outcome']
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
      mat <- allfealist[[j]]
      g <- sub('.*;', '', rownames(mat))
      type <- sub(':.*', '', rownames(mat))
      mat <- mat[type %in% fea.type, , drop = FALSE]
      model <- xgboost(data = t(mat[, trainid, drop=F]), label = trans(y), nrounds = 20, objective='binary:logistic', verbose = F) ## nrounds to be 20!!
      pred <- predict(model,t(mat[,leaveid, drop=F]))
      # pred <- matrix(pred,nrow=2)
      # id <- apply(pred,2,which.max) - 1
      # predy <- revtrans(id)
      truey <- d.pt[leaveid, 'Clinical.Outcome']
      library(pROC)
      auclist[[j]] <- auc(roc(trans(truey), pred))
    }
    names(auclist) <- names(allfealist)
    return(unlist(auclist))
  }
  # ----------------------
  ### for all sample: select features for genes when leaving a study
  fealist.gene <- get_fealist(feature.type = 'gene')
  saveRDS(fealist.gene, paste0(rdir, 'selected_features_gene_window', window,'.rds'))
  
  # ------------------
  ## select features for cytokine 
  fealist.cyt <- get_fealist(feature.type = 'cytokine')
  saveRDS(fealist.cyt, paste0(rdir, 'selected_features_cytokine_window', window,'.rds'))
  
  # --------------------
  ## select features for receptor
  fealist.rec <- get_fealist(feature.type = 'receptor')
  saveRDS(fealist.rec, paste0(rdir, 'selected_features_receptor_window', window,'.rds'))
  
  fealist.tf <- get_fealist(feature.type = 'tf')
  saveRDS(fealist.tf, paste0(rdir, 'selected_features_tf_window', window,'.rds'))
  
  # -----------------------------
  ## combine selected features of base, gene, cytokine, receptor, tf, cellprop
  allfealist <- lapply(which(sapply(fealist.gene, length) > 0), function(i){
    m1 <- fealist.gene[[i]]
    rownames(m1) <- paste0('gene:', rownames(m1))
    m2 <- fealist.cyt[[i]]
    rownames(m2) <- paste0('cytokine:', rownames(m2))
    m3 <- fealist.rec[[i]]
    rownames(m3) <- paste0('receptor:', rownames(m3))
    m4 <- fealist.tf[[i]]
    rownames(m4) <- paste0('tf:', rownames(m4))          ################3
    m <- rbind(m1, m2, m3, m4, t(fea.base), cellProp.mat[i, pbmc.s]) 
    rownames(m)[nrow(m)] <- 'cellProp:'
    m
  })
  names(allfealist) <- names(fealist.gene)[which(sapply(fealist.gene, length) > 0)]
  
  # -------------------------------
  ###  predict using selected features
  ## base = TRUE; gene = FALSE; cellProp = FALSE; cytokine = F; receptor = T
  base = T
  cytokine = FALSE
  receptor = F
  tf = FALSE
  cellProp = F
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
  saveRDS(allauc, paste0(rdir, 'AUC_window',window, '.rds'))
  allauc <- round(allauc,2)
  pd <- reshape2::melt(as.matrix(allauc))
  colnames(pd) <- c('Seed', 'Type', 'AUC')
  pd[,3] <- round(pd[,3],2)
  pd[,1] <- as.factor(pd[,1])
  
  pd$Type <- factor(as.character(pd$Type), levels = names(sort(tapply(pd$AUC, list(as.character(pd$Type)), median))))
  pdf(paste0(pdir, 'AUC_window', window, '_boxplot.pdf'), width = 8, height = 4) 
  print(ggplot(data = pd) +
    geom_boxplot(aes(x = Type, y = AUC), outlier.shape = NA) +
    geom_jitter(aes(x = Type, y = AUC, color = Seed), width = 0.001, height = 0.001) +
    scale_color_brewer(palette = 'Set1')+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('Feature'))
  dev.off()
}
  
