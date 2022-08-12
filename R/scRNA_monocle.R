library(tibble)
library(ggplot2)
library(Seurat)
library(cowplot)
library(dplyr)
library(monocle)



subset.merge <- readRDS("merge.Myogenic.rds")


## Construct the cds object
cds <- newCellDataSet(cellData = subset.merge@assays$RNA@counts,
                         phenoData = new("AnnotatedDataFrame", 
                                         data = subset.merge@meta.data),
                         featureData = new("AnnotatedDataFrame", 
                                           data = data.frame(gene_short_name = rownames(subset.merge), 
                                                             stringsAsFactors=F, 
                                                             row.names = rownames(subset.merge)) ),
                         lowerDetectionLimit = 0.5,
                         expressionFamily =  VGAM::negbinomial.size()
)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


### Ordering cds
Idents(subset.merge) <- 'celltype'
DEGs <- FindAllMarkers(subset.merge, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
DEGs_top500 <- DEGs %>% group_by(cluster) %>% top_n(n = 500, wt = avg_logFC) %>% as.data.frame

cds<- setOrderingFilter(cds, ordering_genes = unique(DEGs_top500$gene))
cds <- estimateSizeFactors(cds)
cds <- reduceDimension(cds, max_components = 2,
                           norm_method = 'log', method = 'DDRTree', verbose = TRUE)
cds <- orderCells(cds)

saveRDSMC(cds,"pseudotime_Myogenic.rds")