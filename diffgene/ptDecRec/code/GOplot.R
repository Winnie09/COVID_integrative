library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(grid)

ctcv <- readRDS('/home-4/zji4@jhu.edu/work-zfs/covid/data/current/palette/celltype.rds')
ct <- readRDS('/home-4/zji4@jhu.edu/work-zfs/covid/data/current/pbmc/celltype.rds')
d <- fread('/home-4/zji4@jhu.edu/scratch/covid/analysis/pbmc/diffgene/withincluster/deceasept/fullGO/GO.csv',data.table = F)

d <- d[order(d$FDR,-d$FC),]
cd <- do.call(rbind,sapply(unique(d$Cluster),function(i) {
  tmp <- d[d$Cluster==i,]
  tmp <- tmp[tmp$FDR < 0.05 & tmp$FC > 2,,drop=F]
  if (nrow(tmp) > 0) {
    tmp[1:min(10,nrow(tmp)),]  
  } else {
    NULL
  }
},simplify = F))
ut <- unique(cd$Term)
d <- d[d$Term %in% ut,c('Cluster','Term','FDR','FC')]
term <- sort(unique(d$Term))
d <- d[d$FDR < 0.05 & d$FC > 2,]

ud <- unique(d$Cluster)
d$Cluster <- factor(d$Cluster,levels = ud[order(as.numeric(sub('c','',sub('_.*','',ud))))])

dmat <- dcast(d,Term~Cluster)
rownames(dmat) <- dmat[,1]
dmat <- as.matrix(dmat[,-1,drop=F])
dmat <- is.na(dmat)

ut <- rownames(dmat)[hclust(dist(dmat))$order]

d$Term <- factor(d$Term,levels=ut)
d$enr <- 'Enriched'

d2 <- expand.grid(unique(d$Cluster),unique(d$Term))
colnames(d2) <- c('Cluster','Term')
d2$FDR <- d2$FC <- d2$enr <- 'Non-enriched'

d <- rbind(d,d2)
d <- d[!duplicated(apply(d[,c('Cluster','Term')],1,paste0,collapse=':')),]

a <- ctcv[sub('.*:','',ct[sub('_.*','',levels(d$Cluster))])]

pdf(paste0('/home-4/zji4@jhu.edu/scratch/covid/analysis/pbmc/diffgene/withincluster/deceasept/fullplot/GO/GO.pdf'),width=9,height=max(3,0.18*length(unique(d$Term))))
print(ggplot(d,aes(x=Cluster,y=Term,fill=enr)) + geom_tile() + theme_classic() + scale_x_discrete(position='bottom') + theme(legend.position = 'bottom',legend.title=element_blank()) + ylab('') + xlab('') + scale_y_discrete(position = "right") + theme(strip.text.y =element_text(angle=180),strip.background = element_rect(size=0.3)) + scale_fill_manual(values=c('orange','dodgerblue3')) + theme(axis.text.x = element_text(size=10,angle = 45, hjust = 1, colour = a),legend.box.margin=margin(-20,-10,-10,-10)))
dev.off()

