type <- readRDS('/home-4/zji4@jhu.edu/scratch/covid/analysis/pbmc/diffgene/withincluster/deceasept/fulltype/type.rds')

suppressMessages(library(topGO))

setwd('/home-4/zji4@jhu.edu/scratch/covid/analysis/pbmc/diffgene/withincluster/deceasept/fullres')
allres <- NULL
for (f in list.files()) {
  st <- type[[sub('.rds','',f)]]
  print(f)
  d <- readRDS(f)
  back <- rownames(d)
  for (ut in unique(st$type)) {
    gl <- rownames(st)[st$type==ut]
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
      allres <- rbind(allres,data.frame(Cluster=sub('.rds','',f),Type=ut,sigres))
    }  
  }
}

write.csv(allres,file=paste0('/home-4/zji4@jhu.edu/scratch/covid/analysis/pbmc/diffgene/withincluster/deceasept/fullGO/typeGO.csv'))


