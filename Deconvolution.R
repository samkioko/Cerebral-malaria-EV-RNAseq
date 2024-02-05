#load the packages
library(edgeR)
library(DESeq2)
library(limma)
library(org.Hs.eg.db)
library(tximport)
library(csaw)
library(readr)
library(ggplot2)
library(ggrepel)
library(WGCNA)
library(ggpubr)
library(bestNormalize)
library(ComplexHeatmap)
library(reshape2)
library(tidyverse)
library(caret)
library(MASS)
library(scatterplot3d)
library(pROC)
library(Boruta)
library(plotROC)
library(GenomicAlignments)
library(GenomicFeatures)
library(xbioc)
library(ComplexHeatmap)
library(variancePartition)
library(elasticnet)
library("FactoMineR")
library("factoextra")
library(Seurat)
library(patchwork)
library(cowplot)
library(mixOmics)
library(readxl)
library(tidyverse)
library(tximport)
library(csaw)
library(illuminaHumanv4.db)
library('annotate')
library('data.table')
library(ggallin)
library(ENIGMA)
library('ggplot2')
library(openxlsx)
library('limma')
library(ggpubr)
library(ggrepel)
library(Seurat)
library(bseqsc)

###########################################################################################################

metadata<-read.xlsx("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/metadata/Metadata_all_Retinopathy_samples.xlsx", 
                    rowNames=TRUE,check.names = FALSE)

load("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/pseudotime_analysisV3/Fig.4/Retinopathy_normalizedData.Rdata")
load("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/pseudotime_analysisV3/Fig7/Gtex1_f.Rdata")


metadata_gtex<-data.frame(sample=colnames(Gtex1_f), condition=colnames(Gtex1_f))
head(metadata_gtex)
rownames(metadata_gtex)<-metadata_gtex$sample


###################################################################################
sobj<-CreateSeuratObject(counts = Gtex1_f[,rownames(metadata_gtex)], 
                         meta.data = metadata_gtex)%>%SCTransform()

sobj
Idents(sobj)<-"condition"

markers_roc<-FindAllMarkers(sobj,test.use = "bimod",
                            min.cells.group = 1,assay = "RNA",
                            only.pos = TRUE)
head(markers_roc)
dim(markers_roc)
markers_roc_noDup<-as.data.frame(markers_roc %>% distinct %>%
                                   group_by(gene) %>% top_n(1, avg_log2FC))
rownames(markers_roc_noDup)<-markers_roc_noDup$gene
markers_roc_noDup<-markers_roc_noDup[order(markers_roc_noDup$gene),]
head(markers_roc_noDup)

table(markers_roc_noDup$cluster)
dim(markers_roc_noDup)

markers_roc_noDup_f<-markers_roc_noDup[markers_roc_noDup$gene%in%rownames(cpm),]

###################################################################################
markers_roc_noDupv2<-as.data.frame(markers_roc_noDup_f %>% distinct %>%
                                   group_by(cluster) %>% top_n(100, avg_log2FC))


table(markers_roc_noDupv2$cluster)
table(markers_roc_noDup$cluster)

#########################################################################
library(ComplexHeatmap)
jpeg(filename = "pseudotime_analysisV3/Fig7/signatureMatrix.jpeg",
     width = 7,height = 10,units = "in",res = 300)
htm<-Heatmap(t(scale(t(Gtex1_f[markers_roc_noDupv2$gene,]))), 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        show_row_names = FALSE,
        row_title_rot = 0,
        heatmap_legend_param = list(title = "expression", 
                                    direction="horizontal"),
        name = "expression",
        row_title = NULL,
        use_raster = FALSE,
        column_names_gp = gpar(fontsize = 20),
        row_gap = unit(0, "mm"),
        row_split = markers_roc_noDupv2$cluster)
htm_Pj<-draw(htm,merge_legend = TRUE, 
             heatmap_legend_side = "top", 
             annotation_legend_side = "top",
             show_heatmap_legend = TRUE,
             legend_title_gp = gpar(fontsize = 20))

dev.off()

################################################################################################
cpm_scaled <- (cpm - mean(cpm)) / sd(as.vector(cpm))
head(cpm_scaled)

Tissuematrix<-read.csv("Fig.4/Matrix_tissueEV_origin.csv",row.names = 1)
Tissuematrix


library(e1071)
deconEVorigin<-do.exLR_origin(avdata.m = as.matrix(log10(cpm+0.1)),
                              ref.m = as.matrix(Tissuematrix),
                              nu.v = c(0.2))


dim(deconEVorigin$est.ab.sum)

meta_f<-metadata[rownames(deconEVorigin$est.ab.sum),]
rownames(meta_f)
Heatmap(deconEVorigin$est.ab.sum[rownames(meta_f),],
        row_split = meta_f$condition)


######################################################################################
library(e1071)
library(scRecover)

scRecover(counts = cpm[,rownames(merge_meta_heatmap)],
          outputDir = "Fig7",
          labels = merge_meta_heatmap$condition,MAGIC = T)

cpm_i<-read.csv("Fig7tempFile/scimpute_count.csv", row.names = 1)
colnames(cpm_i)<-gsub(pattern = "X","",colnames(cpm_i))
colnames(cpm_i)<-gsub(pattern = "\\.","-",colnames(cpm_i))


deconEVorigin<-do.exLR_origin(avdata.m = as.matrix(t((t(log10(cpm_i+0.01))))),
                              ref.m = as.matrix(Tissuematrix),
                              nu.v = c(0.08))

merge_meta_heatmap<-read.xlsx("Fig.3/Heatmap_metadata.xlsx",rowNames = TRUE)

dim(deconEVorigin$est.ab.sum)

meta_f<-merge_meta_heatmap[rownames(deconEVorigin$est.ab.sum),]
rownames(meta_f)

jpeg(filename = "Fig.4/Heatmap_tissueProportion_supplFg1.jpeg",
     width = 6,height = 5,units = "in",res = 300)
Heatmap(deconEVorigin$est.ab.sum[rownames(meta_f),],
        row_split = meta_f$condition,
        cluster_columns = T,
        clustering_method_columns = "ward.D2",
        name = "absolute fraction",
        clustering_distance_columns = "spearman",
        cluster_rows = F, show_row_names = F)
dev.off()

save(deconEVorigin, file="Fig.4/deconEVorigin_tissue.RData")

####################################################################
cm_meta<-subset(meta_f, meta_f$condition!="CC")
cm_meta

decon_CM<-as.data.frame((deconEVorigin$est.ab.sum[rownames(cm_meta),]))
head(decon_CM)
class(decon_CM)
#decon_SVR$Leucocyte<-colSums(as.matrix(decon_SVRHem))
head(decon_CM)
decon_CM<-as.data.frame(t(decon_CM))

library(matrixStats)
decon_CM$mean_Fraction<-rowMeans(as.matrix(decon_CM))
decon_CM$sd<-rowSds(as.matrix(decon_CM))


decon_CM<-decon_CM[order(decon_CM$mean_Fraction),]

pct<-round((decon_CM$mean_Fraction/sum(decon_CM$mean_Fraction))*100,digits = 1)
decon_CM


jpeg(filename = "Fig.4/Pie_TissueProportion_scimpute.jpeg",
     width = 7,height = 5,units = "in",res = 300)
pie(decon_CM$mean_Fraction,radius = 1,init.angle = 150,
    labels = paste(rownames(decon_CM), sep = " ", pct, "%"),
    col = rainbow(length(decon_CM$mean_Fraction),end = 0.8))
dev.off()


########################################################################33
cellProp_meta<-merge(deconEVorigin$est.ab.sum, meta_f, by="row.names")
dim(cellProp_meta)

library(ggthemes)
cellProp_meta$condition<-factor(as.factor(cellProp_meta$condition),
                                levels = c("CC","CM-R⁻","CM-R⁺"))


jpeg(filename = "Fig.4/Pie_TissueProportion_boxplot.jpeg",
     width = 11,height = 7,units = "in",res = 300)
P1<-ggplot(cellProp_meta, aes(condition, Nerve))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size=2,position = position_jitter(0.1))+
  stat_compare_means(method = "t.test",size=7,
                     label = "p.signif",
                     ref.group = "CC")+
  theme_stata()+ ggtitle("Nerve")+
  xlab("")+ylab("")+
  theme(text = element_text(size = 25))

P2<-ggplot(cellProp_meta, aes(condition, Brain))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size=2,position = position_jitter(0.1))+
  stat_compare_means(method = "t.test",
                     label = "p.signif",
                     hide.ns =T, 
                     ref.group = "CC")+
  theme_stata()+ ggtitle("Brain")+
  xlab("")+ylab("absolute proportion")+
  theme(text = element_text(size = 25))

P2+P1
dev.off()

