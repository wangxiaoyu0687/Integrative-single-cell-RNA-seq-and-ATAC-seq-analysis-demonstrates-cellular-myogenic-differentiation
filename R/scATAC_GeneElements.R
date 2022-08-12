library(tibble)
library(ggplot2)
library(Seurat)
library(cowplot)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(UpSetR)
library(RColorBrewer)

obj_atac <- readRDS("mergeATAC.Myogenic.rds")

### annotation peaks to gene elements, as follow annotated order setting the priority
default_elem <- "distal"
order_elem <- setNames(c("exon", "promoter", "5'-UTR", "3'-UTR", "intron"),
                       c("gene_annote/Exon.bed",
                         "gene_annote/promoter.bed",
                         "gene_annote/UTR5.bed",
                         "gene_annote/UTR3.bed",
                         "gene_annote/Intron.bed")
                         )
GRanges_list <- lapply(names(order_elem), function(bed_file){
      bed_data <- read.table(file = bed_file,header = F)
      colnames(bed_data) <- c("chr","start","end","ID","score","strand")
      gr <- makeGRangesFromDataFrame(bed_data)
      return(gr)
})
names(GRanges_list) <- order_elem

peaks <- StringToGRanges(rownames(obj_atac@assays$adjPeaks@counts), sep = c("-", "-"))

overlap_binary_list <- lapply(GRanges_list, function(gr){
      overlap_bin <- countOverlaps(query = peaks, subject = gr,type="any",ignore.strand=T) != 0
      return(overlap_bin)
})


peaks_annot <- apply(
      do.call(cbind,overlap_binary_list)[,order_elem],
      MARGIN = 1,
      function(x){
        if(!any(x)){return(default_elem)}
        order_elem[x][1]
      }
)
names(peaks_annot) <- rownames(obj_atac@assays$adjPeaks@counts)




### set peaks called in percentage 10% of cells in a sample are trustly real peaks
trust_peak_percentage_threshold = 0.1

pct_atac.group <- cluster_mean(obj_atac@assays$adjPeaks@counts>0, Cluster = obj_atac$day_strain)
peaks_trust <- apply(pct_atac.group, MARGIN = 2, function(x) rownames(pct_atac.group)[x>trust_peak_percentage_threshold] )

peaks_annot_count.mtx <- do.call(cbind, lapply(peaks_trust, function(peaks) table(peaks_annot[peaks])[c(order_elem, default_elem)] ))


### output plots
pdf("peaks_annotate_counts.samples.pdf", width = 10, height = 8)
    ggplot(data = reshape2::melt(peaks_annot_count.mtx,value.name="count")) +
      geom_bar(aes(x = Var2, y = count, fill = Var1),
               stat = 'identity', width = 0.7, position = 'stack') +
      scale_fill_manual(values = rev(brewer.pal(12,"Paired")[-11])[1:6], name = "Genomic Elements") +
      ylab("Number of peaks") + scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))) +
      theme_classic()+
      theme(axis.text.y = element_text(size = 14, face="bold"),
            axis.text.x = element_text(size = 14, face="bold"),
            axis.title.x= element_blank(),
            axis.title.y= element_text(size = 16, face="bold"),
            legend.text = element_text(size = 16, face="bold"),
            legend.title = element_text(size = 16, face="bold")
      ) + guides(fill = guide_legend(
        override.aes = list(size=6)
      ))

dev.off()




pct_atac.celltype <- cluster_mean(obj_atac@assays$adjPeaks@counts>0, Cluster = obj_atac$celltype)

pdf("common_peaks_count.upset.pdf", width = 15, height = 8)
  upset(apply(pct_atac.celltype,MARGIN = 2, function(x) ifelse(x>trust_peak_percentage_threshold,1,0)  ) %>% 
              as.data.frame() %>% 
              tibble::rownames_to_column(var = "peaks"),
        sets  = rev(levels(obj_atac$celltype)),
        nsets = nlevels(obj_atac$celltype), nintersects = NA, keep.order = T, order.by = "degree",decreasing = FALSE,
        mb.ratio = c(0.7, 0.3),point.size = 6,
        # scale.intersections = "log10",
        sets.x.label = "total peaks count", mainbar.y.label = "Number of Peaks Intersected",
        sets.bar.color = scales::hue_pal()(nlevels(obj_atac$celltype)), main.bar.color = "#4175b1",
        text.scale = setNames(c(2, 2, 2, 1.5, 2, 2.5),
                              c('intersection size title',
                                'intersection size tick labels',
                                'set size title',
                                'set size tick labels',
                                'set names',
                                'numbers above bars'))
        )
dev.off()
