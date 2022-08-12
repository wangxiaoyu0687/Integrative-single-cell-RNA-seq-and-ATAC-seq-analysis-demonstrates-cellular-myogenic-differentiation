library(SCENIC)
library(tibble)
library(ggplot2)
library(Seurat)
library(cowplot)
library(dplyr)



### Initialize SCENIC settings
subset.merge <- readRDS("merge.Myogenic.rds")

cellInfo <- subset.merge@meta.data

scenicOptions <- initializeScenic(org="hgnc", dbDir="~/databases/TF", nCores = 48)
scenicOptions@settings$dbs['500bp'] <- "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
scenicOptions@settings$dbs['10kb'] <- "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"

saveRDS(cellInfo, file="int/cellInfo.Rds")
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


###  Co-expression network
exprMat <- as.matrix(GetAssayData(subset.merge, slot = "count"))
genesKept.count <- geneFiltering(as.matrix(exprMat), scenicOptions)
exprMat_filtered <- exprMat[genesKept.count, ]
runCorrelation(as.matrix(exprMat_filtered), scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(as.matrix(exprMat_filtered_log), scenicOptions, nParts = 10)

scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top50perTarget")) 
runSCENIC_3_scoreCells(scenicOptions, countData_log , skipHeatmap = T, skipTsne = T)
runSCENIC_4_aucell_binarize(scenicOptions,skipHeatmaps = T)


### Exploring output
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
aucellRankings <- loadInt(scenicOptions, "aucell_rankings")
regulons <- split(regulonTargetsInfo$gene, regulonTargetsInfo$TF)
regulonAUC <- AUCell::AUCell_calcAUC(regulons, aucellRankings,
                                     aucMaxRank=aucellRankings@nGenesDetected["1%"], nCores=2)
regulonActivity <- AUCell::getAUC(regulonAUC)[,rownames(cellInfo)]
write.csv(regulonActivity %>% as.data.frame() %>% rownames_to_column(var="TF"), row.names = F, file = "RegulonActivity.csv")

binary_thresholds <- AUCell::AUCell_exploreThresholds(regulonAUC,
                                               smallestPopPercent=getSettings(scenicOptions,"aucell/smallestPopPercent"),
                                               assignCells=TRUE, plotHist=FALSE,
                                               verbose=TRUE, nCores= 8)

### fixed when row count not match
binary_thresholds <- AUCell::getThresholdSelected(c(binary_thresholds, rep(1, nrow(regulonAUC) - length(binary_thresholds) )))
names(binary_thresholds) <- rownames(binary_thresholds)

regulonsCells <- setNames(lapply(names(binary_thresholds), function(x) {
  thres <- binary_thresholds[x]
  names(which(AUCell::getAUC(regulonAUC)[x, ] >= thres))
}), names(binary_thresholds))

regulonActivity_binary <- reshape2::melt(regulonsCells)
binaryRegulonActivity <- t(table(regulonActivity_binary[, 1], regulonActivity_binary[,2])) %>% as.data.frame.matrix()

write.csv(binaryRegulonActivity %>% as.data.frame() %>%  tibble::rownames_to_column(var="TF"), 
    row.names = F, file = "binaryRegulonActivity.csv")