num.correct.list <- readRDS('/dcl02/hongkai/data/whou/covid/pred/xgboost/res/num.correct.pred.rds')
meta <- readRDS('/dcl02/hongkai/data/covid/data/200916/meta.rds')

mat <- matrix(0, nrow = 30, ncol = length(unique(meta[,'Paper code'])))
rownames(mat) <- names(num.correct.list)
colnames(mat) <- unique(meta[,'Paper code'])
mat <- mat[, !colnames(mat) %in% c("Liao et al.Nature Medicine, 2020", "He et al. Protein&Cell, 2020", "Chua et al. Nature Biotechnology,2020", "Reyfman et al. Am J Respir Crit Care Med. 2019")]
for (i in 1:length(num.correct.list)){
  mat[i, names(num.correct.list[[i]])] <- num.correct.list[[i]]
}

ct <- readRDS('/dcl02/hongkai/data/covid/data/200916/celltype.rds')
rownames(mat) <- paste0(rownames(mat), ';', ct)

library(pheatmap)
pdf('/dcl02/hongkai/data/whou/covid/pred/xgboost/plot/num.correct.pred2.pdf', width = 7, height = 8)
pheatmap(mat, scale = 'none')
dev.off()


## plot cluster composition for each study
d <- readRDS('/dcl02/hongkai/data/whou/covid/pred/data/data/pbmcnorm_combat.rds')
num.s <- sapply(d, ncol)
names(num.s) <- paste0('cluster', 1:length(d))
num <- lapply(1:length(d), function(i){
  data.frame(cluster = paste0('cluster',i, ';', ct[i]), 
             study = meta[match(colnames(d[[i]]), meta[, 'Library Sample Code']), 'Paper code'],
             sample = colnames(d[[i]]))
})
num <- do.call(rbind, num)
ord <- as.numeric(gsub(';.*','',sub('cluster','',unique(num[,1]))))
num[,1] <- factor(num[,1], levels = rev(unique(num[,1])[ord]))

library(ggplot2)
pdf('/dcl02/hongkai/data/whou/covid/pred/xgboost/plot/cluster_study_proportion.pdf', width = 8, height = 6)
ggplot(data = num) + 
  geom_histogram(aes(x = cluster, fill= study), stat = 'count') +
  theme_classic() +
  scale_fill_brewer(palette = 'Set3') +
  theme(axis.text.x = element_text(angle = 90, color = 'black')) +
  ylab('Number of Samples') +
  coord_flip()
dev.off()




