suppressMessages(library(topGO))

setwd('/home-4/zji4@jhu.edu/scratch/covid/analysis/pbmc/diffgene/withincluster/deceasept/res')
allres <- NULL
for (f in list.files()) {
  print(f)
  d <- readRDS(f)
  back <- rownames(d)
  gl <- rownames(d)[d[,'fdr'] < 0.05]
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  
  suppressMessages({GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "Symbol")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  sigres <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classicFisher",numChar=1000)})
  sigres$classicFisher[sigres$classicFisher=="< 1e-30"] <- 0
  sigres <- sigres[sigres$Annotated >= 10,]
  sigres$FDR <- p.adjust(sigres$classicFisher,method="fdr")
  ptcount <- 1
  fc <- ((sigres[,"Significant"]+ptcount)/(sum(GOdata@allScores[GOdata@feasible]==1)+ptcount))/((sigres[,"Annotated"]+ptcount)/(sum(GOdata@feasible)+ptcount))
  sigres <- data.frame(sigres,FC=fc)
  sigres <- sigres[order(sigres$FDR,-sigres$FC),]
  if (nrow(sigres) > 0) {
    allres <- rbind(allres,data.frame(Cluster=sub('.rds','',f),sigres))
  }
}

write.csv(allres,file=paste0('/home-4/zji4@jhu.edu/scratch/covid/analysis/pbmc/diffgene/withincluster/deceasept/GO/GO.csv'))

