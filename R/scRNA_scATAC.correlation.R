library(tibble)
library(ggplot2)
library(Seurat)
library(cowplot)
library(dplyr)
library(Signac)
library(GenomicRanges)


cluster_mean <- function(data, group.by, func){
      func = match.fun(func)
      group.by <- sort(unique(group.by))
      mean_exp <- sapply(group.by, 
                                    function(x) apply(data[,group.by == x, drop=F], 1, func)
                                 )
      colnames(mean_exp) <- group.by
      return(mean_exp)
}


obj_rna <- readRDS("merge.Myogenic.rds")
obj_atac <- readRDS("mergeATAC.Myogenic.rds")

### fix peaks style
peaks_counts <- obj_atac@assays$adjPeaks@counts
peaks_id <- StringToGRanges(rownames(peaks_counts), sep = c('-', '-'))
rownames(peaks_counts) <- GRangesToString(peaks_id, sep = c(':', '-'))

activity.matrix <-  CreateGeneActivityMatrix(peak.matrix = peaks_counts,
                         include.body = T,
                         downstream = 0,
                         annotation.file = "~/databases/refdata-cellranger-Sscrofa11-3.0.1/genes/genes.gtf",
                         seq.levels = paste0('',c(1:18)), upstream = 2000, verbose = TRUE
                       )

atac_activity <- cluster_mean(activity.matrix, group.by = obj_atac$celltype, sum)
atac_activity.norm <- NormalizeData(atac_activity,scale.factor = 10000, normalization.method = "RC")

rna_expression <- cluster_mean(obj_rna[["RNA"]]@data,  group.by = obj_rna$celltype, sum)
rna_expression.norm <- NormalizeData(rna_expression,scale.factor = 10000, normalization.method = "RC")


atac_activity.scaled <- ScaleData(atac_activity.norm)
rna_expression.scale <- ScaleData(rna_expression.norm)



Idents(obj_rna) <- 'celltype'
rna.DEGs <- FindAllMarkers(obj_rna, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
var.genes <- unique(rna.DEGs$gene)
var.genes2 = intersect(var.genes,row.names(activity.matrix))


corr_matrix <- apply(rna_expression.scale[var.genes2,], MARGIN = 2, function(rna){
  apply(atac_activity.scaled[var.genes2,], MARGIN = 2, function(atac){
    cor(rna, atac)
  })
})

write.csv(corr_matrix, file = "scRNA_scATAC.correlation.csv")