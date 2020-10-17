d <- readRDS('/dcl02/hongkai/data/whou/covid/pred/data/data/pbmcnorm_combat.rds')
meta <- readRDS('/dcl02/hongkai/data/covid/data/200916/meta.rds')
ord1 <- readRDS('/dcl02/hongkai/data/rli/covid/multi_sample/data/jason_dist/backbone124_info.rds')
ord2 <- readRDS('/dcl02/hongkai/data/rli/covid/multi_sample/data/jason_dist/branch13_info.rds')

  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))
  i = 1
  expr <- d[[i]]
  
  
  design = matrix(rep(1,8), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), 'intercept')
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  
  ##
  pdt <- data.frame(curve1 = pseudotime, curve2 = pseudotime)
  rownames(pdt) <- names(pseudotime)
  pdt = pdt[colnames(expr), ]
  
  v <- (cellanno$sample %in% paste0('BM',seq(1,8)) + 0)
  v <- ifelse(v==1, 0.99, 0.01)
  cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
  rownames(cellWeights) <- colnames(counts)
  
  set.seed(12345)
  sce <- fitGAM(counts = expr, pseudotime = pdt, cellWeights = cellWeights,
                nknots = 6, verbose = FALSE,parallel=TRUE)
  saveRDS(sce, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'_sce.rds'))
  Final <- list()
  for (TestType in (c('startVsEndTest', 'associationTest'))){
    print(TestType)
    if (grepl('startVsEndTest', TestType)){
      Res <- startVsEndTest(sce)
    } else if (grepl('associationTest', TestType)){
      Res <- associationTest(sce)
    } 
    res <- data.frame(waldStat = Res[,'waldStat'], P.Value = Res[,'pvalue'] ,adj.P.Val = p.adjust(Res$pvalue, method='fdr'))
    row.names(res) <- row.names(Res)
    res <- res[order(res[,3], -res[,1]), ]
    Final[[TestType]] <- res
  }



