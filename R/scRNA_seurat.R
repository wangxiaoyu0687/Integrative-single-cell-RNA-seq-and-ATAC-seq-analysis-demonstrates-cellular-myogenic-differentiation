library(dplyr)
library(tibble)
library(cowplot)
library(ggplot2)
library(Seurat)


opt <- readRDS('seurat_process_parameters.rds')

sampleNames <- c("Y19SS133", "Y19SS134", "Y19SS143", "Y19SS151", "Y19SS167", "Y20SS001", "Y20SS00149", "Y20SS00155")
group_setting <- c("E18-ZZ", "E18-DZ", "E21-ZZ", "E28-ZZ", "E28-DZ", "E21-DZ", "E16-DZ", "E16-ZZ")

#  Load datasets and create Seurat objects 

Sus_scrofa.MT.genes = read.table('Sus_scrofa.MT.genes.txt', stringsAsFactors = F)$V1

seurat_obj.list <- lapply(seq_along(sampleNames), function(i){
          
          data_dir <- list.files(path = sprintf("%s/%s", opt$indir, sampleNames[i]), recursive = T, pattern = "matrix.mtx", full.names = T)
          data_dir <- data_dir[grepl("filtered_\\w+_bc_matr", data_dir, perl =T)]

          d = Read10X(dirname(data_dir[1]))
          obj <- CreateSeuratObject(counts = d, project = sampleNames[i] , min.cells = 3, min.features = 0)

          obj[["percent.mt"]] <- PercentageFeatureSet(d0, features = Sus_scrofa.MT.genes[Sus_scrofa.MT.genes %in% row.names(obj)] )
          obj[["percent.rp"]] <- PercentageFeatureSet(obj, pattern = "^R[Pp][LSls]")

          obj@meta.data$sampleName <- sampleNames[i]
          obj <- RenameCells(object = obj, add.cell.id = sampleNames[i] )
          obj <- NormalizeData(obj)
          obj <- FindVariableFeatures(obj, selection.method = "vst",nfeatures= 3000)
          obj <- ScaleData(obj)
          return(obj)
        })

names(seurat_obj.list) <- sampleNames


### filter low quality cells
seurat_obj.list <- lapply(seurat_obj.list, 
  function(obj) subset(obj, nFeature_RNA>200 & nCount_RNA < 30000 & percent.mt < 10) )


### Integration samples
selected_anchor.features <- SelectIntegrationFeatures(seurat_obj.list, nfeatures  = 3000)
seurat_anchors <- FindIntegrationAnchors(object.list = seurat_obj.list, dims = 1:60,
                                     anchor.features = selected_anchor.features,
                                     scale = FALSE)
merge.obj <- IntegrateData(anchorset = seurat_anchors, dims = 1:60)


DefaultAssay(merge.obj) <- "integrated"
merge.obj <- RunPCA(merge.obj, npcs = 60) %>% 
  ScaleData() %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = seq(0.1,0.9,0.1))

DefaultAssay(merge.obj) <- "RNA"


#proceed with analysis with cluster resolution 0.6
Idents(merge.obj) <- merge.obj@meta.data$RNA_snn_res.0.6
merge.obj$seurat_clusters <- merge.obj$RNA_snn_res.0.6
merge.obj$day_strain <- setNames(group_setting, sampleNames)[merge.obj$sampleNames]
saveRDS(merge.obj,"merge.rds")

#================= Finding differentially expressed features (cluster biomarkers) =================
DEGs.cluster <- FindAllMarkers(merge.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
