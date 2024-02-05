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
library(GenomicAlignments)
library(GenomicFeatures)
library(org.Pf.plasmo.db)
library(variancePartition)
library(elasticnet)
library("FactoMineR")
library("factoextra")
library(Seurat)
library(patchwork)
library(cowplot)
library(mixOmics)
library(readxl)
library("DiscoRhythm")
library(tximport)
library(csaw)
library(openxlsx)
library(matrixTests)
library(ISOpureR)
library(tximport)
library(RUVSeq)
library(GEOquery)
library(lumi)
library(sva)
library(BRETIGEA)
library(MultiBaC)
library(noise)
library(ViSEAGO)
library(M3Drop)
library(readr)

########################################################################################################################
metadata<-read.xlsx("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/metadata/Metadata_all_Retinopathy_samples.xlsx", 
                    rowNames=TRUE,check.names = FALSE)
colnames(metadata)

load("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/hs_ret_all_kallisto_genelevel.Rdata")

hs_ret_all_kallisto_genelevel$counts
hs_ret_all_kallisto_genelevel$length


#start the process of normalizing the counts by gene length and lib size

#make the design
design<-model.matrix(~retinopathy, metadata)
dim(design)
head(design)

##################################################################################

expr_pos <- 
  noisyr::noisyr(approach.for.similarity.calculation = "counts", 
                 expression.matrix = hs_ret_all_kallisto_genelevel$counts)


cts<-hs_ret_all_kallisto_genelevel$counts[rownames(expr_pos),]
library(scRecover)
cpm_scRecover<-scRecover(counts = cts[,rownames(merge_meta_heatmap)], 
                         labels = merge_meta_heatmap$condition,
                         outputDir = "pseudotime_analysisV3/Fig7/")





normMat<-hs_ret_all_kallisto_genelevel$length[rownames(cts),]
dim(cts)

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.

eff.lib <- calcNormFactors(normCts,method = "RLE") * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
normMat
library(edgeR)
#Creating a DGEList object for use in edgeR.
colnames(design)<-make.names(colnames(design))
dgeObj <- DGEList(cts,
                  group =metadata$retinopathy,
                  samples = metadata)

dgeObj<-scaleOffset(dgeObj,normMat)
dgeObj$offset

###########################################################################################
dgeObj$counts

####################################################################
#calculate cpm
cpm_log<-edgeR::cpm(y = dgeObj,normalized.lib.sizes = TRUE,
                    log = TRUE, prior.count = 0.001)
head(cpm_log)
####################################################################

cpm<-edgeR::cpm(dgeObj,normalized.lib.sizes = TRUE, log=FALSE)
library(circlize)
col_fun = colorRamp2(c(-1,0, 2), c("Darkblue", "black","gold"))

metadata1<-metadata[order(metadata$bradno),]


#############################################################################################################
library(phenopath)
library(ggallin)
library(ggbeeswarm)
######################################################################
all(rownames(metadata)==colnames(cpm))
dim(cpm)
head(cpm)
dim(cpm_loess)

fit <- phenopath(exprs_obj = (t(sqrt(cpm[]))), 
                 x = as.factor(metadata[colnames(cpm),]$retinopathy),
                 elbo_tol = 1e-6)

# Extract the trajectory
z <- trajectory(fit)
z
metadata<-metadata[colnames(cpm),]
dim(metadata)
dim(cpm)
metadata$Trajectory<-z
metada_noDup<-metadata[,!duplicated(colnames(metadata))]
dim(metada_noDup)
colnames(metada_noDup)

##################################################################
dim(cts)
cpm
dim(cpm)

#save(cpm, metada_noDup, file = "cpm_and_metada_noDup")
library(ggallin)
table(metada_noDup$condition)
ggplot(metada_noDup, aes(x = condition, y=Trajectory))+
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means(ref.group = "negative",method = "wilcox")+
  scale_y_continuous(trans = ssqrt_trans)+
  geom_beeswarm(aes(color=study),dodge.width = 0.7)

head(metada_noDup)

######################################################################
#run loess with a low pass filter
metada_noDup<-metada_noDup[order(metada_noDup$Trajectory),]
metada_noDup$Trajectory

metada_noDup<-subset(metada_noDup, Trajectory>-1 & Trajectory<1)

cpm<-cpm[,rownames(metada_noDup)]
head(cpm)
dim(cpm)

cpm_t<-as.data.frame(t(sqrt(cpm)))
head(cpm_t[1:10])

##LOESS the whole data
vars <- colnames(cpm_t)
vars
## covariate
id <- c(1:nrow(cpm_t))
id

############################################################################################
## define a loess filter function (fitting loess regression line)
loess.filter <- function (x, span) loess(formula = paste(x, "id", 
                                                         sep = "~"),
                                         data = cpm_t,
                                         degree = 2,
                                         span = span)$fitted 
## apply filter column-by-column
new.dat_p <- as.data.frame(lapply(vars, loess.filter, span = 0.25),
                           col.names = colnames(cpm_t))

log2_loess_p<-t(new.dat_p)
head(log2_loess_p)
colnames(log2_loess_p)<-rownames(cpm_t)
head(log2_loess_p)

############################################################################################
cpm_loess<-as.data.frame((log2_loess_p))
cpm_loess[1:10]

pca1<-mixOmics::pca(t((cpm_loess)), scale = TRUE)

plotIndiv(pca1,ind.names = FALSE,group = metada_noDup$condition)

cosinor<-row_cosinor(x = cpm_loess,t = c(1:85),period = 85)

table(cosinor$pvalue<0.01)

#####################################################################################
metada_noDup$disease_stage <- ""
metada_noDup[which(metada_noDup$Trajectory >(0) ),"disease_stage"] <- "late"
metada_noDup[which(metada_noDup$Trajectory <(0) ),"disease_stage"] <- "early"
head(metada_noDup)
metada_noDup$disease_stage
rownames(metada_noDup)
colnames(cpm_loess)

pca1<-pca(X = t((cpm_loess)), scale = TRUE,ncomp = 20)
plotIndiv(pca1,ind.names = FALSE,group = metada_noDup$disease_stage)


comp1<-as.data.frame(pca1$variates$X)

set.seed(12345)
umap1<-uwot::umap(X = (comp1),n_components = 5,n_neighbors = 10)
umap1

umap1<-as.data.frame(umap1)


comp<-rename(.data = umap1,MDS1=V1, MDS2=V2)
head(comp)
comp$condition<-metada_noDup$condition
comp$condition

comp$condition<-factor(as.factor(comp$condition),
                       levels = c("CC","CM-R⁻","CM-R⁺"))
head(comp)

########################################################################################
library(ggthemes)
jpeg(filename = "pseudotime_analysisV3/Fig.2/UMAP_plotV2.jpeg",
     width = 6,height = 6,units = "in",res = 400)
ggplot(data = comp,aes(MDS1*-1,MDS2,shape=condition, 
                       fill=condition)) +
  #stat_ellipse(level = 0.9999, aes(color=Malaria), show.legend = FALSE)+
  #geom_mark_ellipse(aes(fill=Malaria),alpha=0.05, show.legend = FALSE,n=50)+
  geom_point(data = comp,aes(MDS1,MDS2,shape=condition, 
                             fill=condition),
             size=4,alpha=0.7)+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = c("gold","steelblue",
                               "firebrick"))+
  scale_color_manual(values = c("gold","steelblue",
                                "firebrick"))+
  theme_stata()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size = 20),
        legend.title=element_text(size = 20), 
        legend.text=element_text(size=20))+
  xlab("Dim 1")+
  ylab("Dim 2")
  
  #geom_vline(xintercept = 0, linetype="dashed")
#xlim(c(-.7, .5))+ylim(-.5,.5)

dev.off()

###############################################################################################3

library(slingshot)
comp$disease_stage<-metada_noDup$disease_stage
pseudotime<-slingshot(data = comp[,1:2],
                      clusterLabels = comp$disease_stage)
pseudotime@assays@data$pseudotime[,1]

metada_noDup$Trajectory_slingshot<-pseudotime@assays@data$pseudotime[,1]
metada_noDup$condition<-factor(as.factor(metada_noDup$condition),
                       levels = c("CC","CM-R⁻","CM-R⁺"))
head(comp)


metada_noDup$cohort<-metada_noDup$study
metada_noDup$cohort<-gsub("test", "one", metada_noDup$cohort)
metada_noDup$cohort<-gsub("validation", "two", metada_noDup$cohort)

comparemeans<-list(c("CC", "CM-R⁻"), c("CC","CM-R⁺"), c("CM-R⁺","CM-R⁻"))



jpeg(filename = "pseudotime_analysisV3/Fig.2/Box_plot_pseduotimeV2.jpeg",
     width = 4.5,height = 5,units = "in",res = 400)
ggplot(metada_noDup, aes(x = condition, y=Trajectory_slingshot))+
  geom_boxplot(aes(fill=condition),outlier.shape = NA, alpha=0.7, show.legend = FALSE)+
  stat_compare_means(ref.group = "negative",method = "wilcox")+
  #scale_y_continuous(trans = ssqrt_trans)+
  geom_point(position = position_jitter(width = 0.2),size=2)+
  #geom_beeswarm(aes(shape=condition),dodge.width = 0.7, size=3)+
  scale_fill_manual(values = c("gold","steelblue",
                               "firebrick"))+
  scale_color_manual(values = c("gold","steelblue",
                                "firebrick"))+
  theme_pubr(border = TRUE,legend = "right")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size = 20),
        legend.title=element_text(size = 20), 
        legend.text=element_text(size=20))+
  ylab("Pseudotime")+xlab("")+
  stat_compare_means(comparisons = comparemeans,
                     label = "p.signif",size=9,vjust = 0.5)
dev.off()

comp$pseudotime<-metada_noDup$Trajectory_slingshot
library(viridis)
##############################################################################
jpeg(filename = "pseudotime_analysisV3/Fig.2/UMAP_plot_pseudotime.jpeg",
     width = 6,height = 6.5,units = "in",res = 400)
ggplot(data = comp,aes(MDS1*-1,MDS2, shape=condition, 
                       fill=pseudotime)) +
  #stat_ellipse(level = 0.9999, aes(color=Malaria), show.legend = FALSE)+
  #geom_mark_ellipse(aes(fill=Malaria),alpha=0.05, show.legend = FALSE,n=50)+
  geom_point(data = comp,aes(MDS1,MDS2),
             size=4,alpha=0.7)+
  scale_shape_manual(values = c(21:25))+
  scale_fill_viridis(option = "magma")+
  theme_stata()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size = 20),
        legend.title=element_text(size = 20), 
        legend.text=element_text(size=20))+
  xlab("Dim 1")+
  ylab("Dim 2")+guides(shape="none")

dev.off()

###############################################################################################
library(readstata13)
mogeni<-read.dta13("~/Downloads/AnalysisDataset_SevMal_v1_1998_2019.dta")
head(mogeni)
colnames(mogeni)[1]<-c("Serialno_2")

mogeni
grep("ret_all",colnames(metada_noDup))
colnames(metada_noDup)

merge_metaMogeni<-merge(metada_noDup[,c(1,2,5,198,631)],mogeni,by="Serialno_2")
merge_metaMogeni

write.xlsx(merge_metaMogeni,file = "pseudotime_analysisV3/Fig.2/Forest_plot_mogeni_metadata.xlsx")


library(readstata13)
mogeni<-read.dta13("~/Desktop/AnalysisDataset_SevMal_v1_1998_2019.csv")

colnames(mogeni)


library(openxlsx)
df<-read.xlsx("pseudotime_analysisV3/Fig.2/Forest_plot_mogeni_metadata.xlsx",sheet = 2,
              rowNames = FALSE,check.names = FALSE)
head(df)
hrp2<-read.xlsx("pseudotime_analysisV3/HRP_Plasma_SK_06feb13.xlsx")
hrp2
colnames(hrp2)[1]<-"Serialno_2"
merge_df_hrp<-merge(df, hrp2[,c(1,3),drop=FALSE], by="Serialno_2", all.x=TRUE)
head(merge_df_hrp)
str(merge_df_hrp)

df$Parasitaemia<-scale(sqrt(as.numeric(df$Parasitaemia)))

df$`Age.(months)`<-scale(sqrt(as.numeric(df$`Age.(months)`)))
df$Respiratory.rate<-sqrt(as.numeric(df$Respiratory.rate))
df$Temperature<-sqrt(as.numeric(df$Temperature))
df$MUAC<-sqrt(as.numeric(df$MUAC))
df$Admission.days<-sqrt(as.numeric(df$Admission.days))

df$BCS.verbal<-(as.numeric(df$BCS.verbal))
df$BCS.motor<-(as.numeric(df$BCS.motor))
df$BCS.eyes<-(as.numeric(df$BCS.eyes))

head(merge_df_hrp)
head(df)
dim(merge_df_hrp)
class(merge_df_hrp)
merge_df_hrp$hrp2_plasma

library(brglm2)
fit_un <- lapply(merge_df_hrp[,c(4:ncol(merge_df_hrp))], function(x) coef(summary(glm((Pseudotime) ~ x+`Age.(months)`, 
                                                                     data =merge_df_hrp, 
                                                                     family = gaussian()
                                                                     ))))

fit_un


class(fit_un)
fit_un_df<-do.call(rbind.data.frame, fit_un) #convert the list to a data frame
fit_un_df

df_un<-fit_un_df[!grepl("Intercept|Age", rownames(fit_un_df)),] #remove rows with the intercept
rownames(df_un)<-gsub("\\.x.*", "", rownames(df_un))
df_un

merge_df_hrp
colnames(metadata)

########################################################################################
metadata$seq_run<-rownames(metadata)
colnames(metadata)

cm_pos_dichoto<-read.xlsx("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/metadata/retinopathy_for Kioko.xlsx")
cm_pos_dichoto

new_metadata<-merge(metadata[,c(5,630),drop=F], cm_pos_dichoto, by = "Serialno_2")
new_metadata

s1<-read.xlsx("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/supplementary_data_files/supplementary_data_1.xlsx",check.names = F)
s1

colnames(new_metadata)
colnames(merge_df_hrp)
new_s1<-merge(new_metadata, s1, by="seq_run",all.y=T)
new_s1<-merge(merge_df_hrp[,c(1,19)], new_s1, by="Serialno_2",all.y=T)

new_s1
dim(new_s1)
write.xlsx(new_s1,rowNames=TRUE, 
           file = "supplementary_data_file_revised/supplementary_data_1.xlsx")


df2<-read.xlsx(xlsxFile = "pseudotime_analysisV3/Fig.2/LinearRegression_predicting_pseudotimeV2_adjustedforAge.xlsx", 
               rowNames = FALSE)
head(df2)

library(ggplot2)
library(ggforestplot)
library(ggthemes)
library(Seurat)

jpeg(filename = "pseudotime_analysisV3/Fig.2/ForestplotV2_adjusted4Age.jpeg",
     width = 6.5,height = 6.3,units = "in",res = 400)
forestplot(df = df2,name = name,estimate = Estimate,se = SE,size=12,
           pvalue = Pvalue,colour = legend)+theme_stata()+
  theme(axis.text.y =  element_text(angle = 0),
        text = element_text(size = 20))+
  scale_color_manual(values = c("steelblue","black","firebrick"))+
  scale_x_continuous(breaks = c(-6,-3,0, 3, 6))+
  NoLegend()
dev.off() 

#######################################################################
#run cosinor to get 
all(colnames(cpm_loess)==rownames(metada_noDup))

cosinor_mod<-row_cosinor(x = cpm_loess,
                         t = metada_noDup$Trajectory_slingshot,
                         period = 15)
table(cosinor_mod$pvalue<0.05)

library(illuminaHumanv4.db)
cosinor_mod$symbol<-mapIds(x = illuminaHumanv4.db,
                           keys = rownames(cosinor_mod),
                           column = "SYMBOL",
                           keytype = "ENSEMBL")


cosinor_mod$geneName<-mapIds(x = illuminaHumanv4.db,
                           keys = rownames(cosinor_mod),
                           column = "GENENAME",
                           keytype = "ENSEMBL")
head(cosinor_mod)
cosinor_mod_sig<-subset(cosinor_mod, pvalue<0.05)
#
#######################################################################
write.xlsx(cosinor_mod_sig, rowNames=TRUE,
           "pseudotime_analysisV3/Fig.3/Cosinor_periodicGenes.xlsx")

library(clusterProfiler)
library(openxlsx)
library(illuminaHumanv4.db)
cosinor_mod_sig<-read.xlsx(xlsxFile = "pseudotime_analysisV3/Fig.3/Cosinor_periodicGenes.xlsx", rowNames = TRUE)

geneList<-scale(cosinor_mod_sig$acrophase)
names(geneList)<-rownames(cosinor_mod_sig)

geneList<-sort(geneList,decreasing = TRUE)


gsea<-gseGO(geneList = geneList,ont = "BP",
      OrgDb = illuminaHumanv4.db,
      keyType = "ENSEMBL",minGSSize = 10,
      maxGSSize = 500,pvalueCutoff = 0.05,
      pAdjustMethod = "none")

gsea<-setReadable(x = gsea,OrgDb = illuminaHumanv4.db,keyType = "ENSEMBL")
gsea<-as.data.frame(gsea)
gsea$NES
################################################################################################
library(rrvgo)

simMatrix<-calculateSimMatrix(gsea$ID,orgdb = "illuminaHumanv4.db",
                              ont ="BP" ,keytype = "ENSEMBL",method = "Rel")

scores <- setNames(-log10(gsea$pvalue), gsea$ID)
scores

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="illuminaHumanv4.db")

heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

dim(reducedTerms)
dim(gsea)

########################################################################################################################
colnames(reducedTerms)[1]<-"ID"

merge_GOterm<-merge(reducedTerms, gsea, by='ID')
head(merge_GOterm)

write.xlsx(merge_GOterm,file = "pseudotime_analysisV3/Fig.3/GSEA_sematic_relationship.xlsx")



colnames(metada_noDup)
ncol(metada_noDup)
metada_noDup$run_id<-rownames(metada_noDup)
colnames(df)

merge_meta_heatmap<-merge(metada_noDup[,c(633,5,1,2,632,631)], df, by="Serialno_2", all.x=TRUE)
colnames(merge_meta_heatmap)
rownames(merge_meta_heatmap)<-merge_meta_heatmap$run_id
head(merge_meta_heatmap)
str(merge_meta_heatmap)

write.xlsx(merge_meta_heatmap, file = "pseudotime_analysisV3/Fig.3/Heatmap_metadata.xlsx", rowNames=TRUE)

library(openxlsx)
merge_meta_heatmap<-read.xlsx("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/pseudotime_analysisV3/Fig.3/Heatmap_metadata.xlsx",rowNames = TRUE)


merge_meta_heatmap<-merge_meta_heatmap[order(merge_meta_heatmap$Trajectory_slingshot),]
merge_meta_heatmap
cosinor_mod_sig<-cosinor_mod_sig[order(cosinor_mod_sig$acrophase),]

cpm_loess_sig<-cpm_loess[rownames(cosinor_mod_sig),][,rownames(merge_meta_heatmap)]
cpm_loess
library(circlize)
col_funp = colorRamp2(c(-2, 0, 2), c("green", "black", "purple"))



Heatmap(t(scale(t(cpm_loess_sig))),cluster_columns = FALSE,
        cluster_rows = FALSE, show_column_names = FALSE,
        show_row_names = FALSE, use_raster = FALSE, 
        col = col_funp, row_split=cosinor_mod_sig$cluster
        )

library(circlize)
library(viridis)
metada_noDup$cohort

colret = list(condition = c("CM-R⁺"="firebrick","CC"="gold",
                            "CM-R⁻"="steelblue"),
              acidosis = c("0"="bisque",
                         "1"="darkseagreen"),
              died = c("0"="lightblue1",
                         "1"="tan1"),
              cohort = c("one"="lightblue",
                       "two"="tan1"),
              pseudotime = colorRamp2(c(0,4,8,12,16),viridis(n = 5,
                                                                     option = "B",
                                                                     begin = 0,end = 1)))
              
colret

merge_meta_heatmap$condition<-factor(as.factor(merge_meta_heatmap$condition),levels = c("CC","CM-R⁻","CM-R⁺"))

library(ComplexHeatmap)
ha1 = HeatmapAnnotation(pseudotime=merge_meta_heatmap$Trajectory_slingshot,
                        condition= merge_meta_heatmap$condition,
                        died=merge_meta_heatmap$OutcomeDied,
                        acidosis=merge_meta_heatmap$Metabolic.acidosis,
                        #cohort=merge_meta_heatmap$cohort,
                        annotation_name_side = "left",
                        show_legend =TRUE,col = colret,
                        show_annotation_name = TRUE,
                        na_col = "white",border = TRUE,
                        annotation_legend_param = list(legend_gp = gpar(fontsize = 20)),
                        simple_anno_size = unit(0.8, "cm"),
                        annotation_name_gp = gpar(fontsize = 20,fontface = "bold"))
ha1
graphics.off()
lab<-c("GFAP","NCAN","VEGFA","TGFB1","TGFBR1","AQ4",
       "SLC1A2","CREB1","CREBRF","GRIA2","GRIA1","ROBO2",
       "SEMA5B","TRPM3","TMEM196","GFRA2", "S100A9","S100A8",
       "ELANE", "TBX21","AZU1","DEFA3","CTSG","CTSA")
mark_at = which(res_sig$symbol %in% lab)
mark_at

har = rowAnnotation(foo = anno_mark(at = which(res_sig$symbol %in% lab),
                                    labels = res_sig$symbol[res_sig$symbol%in%lab]))


####################################################################################
jpeg(filename = "pseudotime_analysisV3/Fig.3/Pseudotime_htm_retinopathy_combinedV5.jpeg",
     units = "in",width = 8.5,height = 11,res = 300)
htm<-Heatmap(t(scale(t(cpm_loess_sig))),
             clustering_method_rows = "single",
             clustering_method_columns = "complete",
             clustering_distance_columns = "canberra",
             clustering_distance_rows = "canberra",
             heatmap_legend_param = list(title = "mean z score", 
                                         direction="vertical"),
             top_annotation = ha1,
             show_row_names = FALSE, 
             #row_labels = cosi_sig$symbol,
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             cluster_row_slices = FALSE, 
             show_row_dend = FALSE,
             #column_split = metada_noDup_noCC$study,
             col = col_funp,
             column_title = NULL,
             column_labels = c(1:85),
             row_title_rot = 0,
             column_names_rot = 0,column_names_side = "top",
             show_column_names = FALSE,
             row_gap = unit(0.3, "mm"),
             row_split = cosinor_mod_sig$cluster,
             row_names_gp = gpar(fontsize = 10),
             row_names_max_width = unit(5,"cm"),
             column_title_gp = gpar(fontsize = 20),
             row_title_gp = gpar(fontsize = 20),
             use_raster = FALSE)
htm
htm_Pj<-draw(htm,
             merge_legend = TRUE, 
             heatmap_legend_side = "right",
             legend_title_gp = gpar(fontsize = 20, 
                                    fontface = "bold"),
             annotation_legend_side = "top",
             show_heatmap_legend = TRUE)
dev.off()


#run GSVA BP followed by REVIGO

library(msigdbr)
library(fgsea)
library(stringr)
as.data.frame(msigdbr_collections())
cgp_gene_sets = msigdbr(species = "human", 
                        category = "C2", 
                        subcategory = "CP:REACTOME")
cgp_gene_sets
msigdbr_list = split(x = cgp_gene_sets$ensembl_gene, 
                     f = cgp_gene_sets$gs_name)
msigdbr_list

library("GSVA")
library(openxlsx)


cpm_loess_sig<-read.xlsx(xlsxFile = "CPM_loess_Retinopathy.xlsx", rowNames = TRUE)
dim(cpm_loess_sig)
head(cpm_loess_sig)

gsva_es<-gsva(gset.idx.list = msigdbr_list,
              method="gsva",kcdf="Gaussian",
              expr = as.matrix(cpm_loess_sig[,-c(86)]),
              min.sz=20,
              max.sz=500)

dim(gsva_es)
head(gsva_es)
gsva_es<-gsva_es[,rownames(merge_meta_heatmap)]
all(colnames(gsva_es)==rownames(merge_meta_heatmap))

library(matrixTests)
cosinor_BP<-row_cosinor(gsva_es,t = merge_meta_heatmap$Trajectory_slingshot,period = 15)
table(cosinor_BP$pvalue<0.05)
head(cosinor_BP)

##############################################################
write.xlsx(cosinor_BP, file = "pseudotime_analysisV3/Fig.3/Reactome_gsva.xlsx", rowNames=TRUE)

#use revigo to summarize gene Terms
#############################################################################################
library(rrvgo)

simMatrix<-calculateSimMatrix(rownames(cosinor_BP),orgdb = "illuminaHumanv4.db",
                              ont ="BP" ,keytype = "SYMBOL",method = "Rel")

scores <- setNames(-log10(cosinor_BP$pvalue), rownames(cosinor_BP))
scores

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="illuminaHumanv4.db")

heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

dim(reducedTerms)
dim(gsea)

########################################################################################################################
colnames(reducedTerms)[1]<-"ID"
cosinor_BP$ID<-rownames(cosinor_BP)

merge_GOterm<-merge(reducedTerms, cosinor_BP, by='ID')
head(merge_GOterm)

write.xlsx(merge_GOterm,file = "pseudotime_analysisV3/Fig.3/GSVA_sematic_relationship.xlsx")
########################################################################################
library(illuminaHumanv4.db)
cosinor_mod_sig$entrezid<-mapIds(x = illuminaHumanv4.db, 
                                 keys = rownames(cosinor_mod_sig),
                                 column = "ENTREZID",
                                 keytype = "ENSEMBL")


library(ViSEAGO)
#browseVignettes("ViSEAGO")
# connect to Bioconductor
Bioconductor<-ViSEAGO::Bioconductor2GO()
# load GO annotations from Bioconductor
myGENE2GO<-ViSEAGO::annotate(
  "org.Hs.eg.db",
  Bioconductor)
# create topGOdata for BP for each list of DE genes
cosinor_mod_sig<-merge_brain_genes
g1<-subset(cosinor_mod_sig, cluster=="c1")
head(g1)
dim(g1)

g1_BP<-ViSEAGO::create_topGOdata(
  geneSel=g1$entrezid,
  allGenes=cosinor_mod_sig$entrezid,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=10)

g1_BP

g2<-subset(cosinor_mod_sig, cluster=="c2")
head(g2)
dim(g2)

g2_BP<-ViSEAGO::create_topGOdata(
  geneSel=g2$entrezid,
  allGenes=cosinor_mod_sig$entrezid,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=10)

g2_BP

g3<-subset(cosinor_mod_sig, cluster=="c3")
head(g3)
dim(g3)

g3_BP<-ViSEAGO::create_topGOdata(
  geneSel=g3$entrezid,
  allGenes=cosinor_mod_sig$entrezid,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=10)

g3_BP

g4<-subset(cosinor_mod_sig, cluster=="c4")
head(g4)
dim(g4)

g4_BP<-ViSEAGO::create_topGOdata(
  geneSel=g4$entrezid,
  allGenes=cosinor_mod_sig$entrezid,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=10)

g4_BP

# perform topGO tests
elim_g1_BP<-topGO::runTest(
  g1_BP,
  algorithm ="elim",
  statistic = "fisher",
  cutOff=0.05)
elim_g1_BP

elim_g2_BP<-topGO::runTest(
  g2_BP,
  algorithm ="elim",
  statistic = "fisher",
  cutOff=0.05)
elim_g2_BP

elim_g3_BP<-topGO::runTest(
  g3_BP,
  algorithm ="elim",
  statistic = "fisher",
  cutOff=0.05)
elim_g3_BP

elim_g4_BP<-topGO::runTest(
  g4_BP,
  algorithm ="elim",
  statistic = "fisher",
  cutOff=0.05)
elim_g4_BP


###############################################################
# merge topGO cosi_sig
BP_scosi_sig<-ViSEAGO::merge_enrich_terms(
  cutoff=0.05,
  Input=list(
    cluster_1=c(
      "g1_BP",
      "elim_g1_BP"),
    cluster_2=c(
      "g2_BP",
      "elim_g2_BP"),
    cluster_3=c(
      "g3_BP",
      "elim_g3_BP"),
    cluster_4=c(
      "g4_BP",
      "elim_g4_BP")))


# barchart of significant (or not) GO terms by comparison
ViSEAGO::GOcount(BP_scosi_sig)


# create GO_SS-class object
myGOs<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_scosi_sig)

# compute Semantic Similarity (SS)
myGOs<-ViSEAGO::compute_SS_distances(
  myGOs,
  distance="Rel")

graphics.off()
# Create GOterms heatmap
Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=FALSE,
  showGOlabels =FALSE,
  GO.tree=list(
    tree=list(
      distance="Rel",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=4,
        minClusterSize =4
      )
    )
  ),
  samples.tree=NULL)

# display the heatmap
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms")

# colored MDSplot
ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOterms")

Wang_clusters_wardD2
ls("package:ViSEAGO")
show_table(Wang_clusters_wardD2)


# MDSplot
ViSEAGO::MDSplot(myGOs)


Viseago_results<-Wang_clusters_wardD2@enrich_GOs@data

head(Viseago_results)



write.xlsx(Viseago_results,rowNames=TRUE,
           file = "pseudotime_analysisV3/Fig.3/Pseudotime_Braingenes_elimFisher.xlsx",
           overwrite = TRUE)

g1

################################################################
library(TissueEnrich)

gs1<-GeneSet(geneIds=rownames(g1),geneIdType=ENSEMBLIdentifier(),
             organism="Homo Sapiens")
gs1

gs1output<-teEnrichment(inputGenes = gs1,rnaSeqDataset = 1,
                       multiHypoCorrection = FALSE,
                       #backgroundGenes = rownames(cosinor_mod_sig),
                       tissueSpecificGeneType = 2)
gs1output

gs1seEnrichmentOutput<-gs1output[[1]]
gs1enrichmentOutput<-setNames(data.frame(assay(gs1seEnrichmentOutput),
                                      row.names = rowData(gs1seEnrichmentOutput)[,1]), 
                           colData(gs1seEnrichmentOutput)[,1])
gs1enrichmentOutput$Tissue<-row.names(gs1enrichmentOutput)
head(gs1enrichmentOutput)
gs1enrichmentOutput$cluster<-"c1"
gs1enrichmentOutput

####################################################################################################
####################################################################################################
gs2<-GeneSet(geneIds=rownames(g2),geneIdType=ENSEMBLIdentifier(),
             organism="Homo Sapiens")
gs2

gs2output<-teEnrichment(inputGenes = gs2,rnaSeqDataset = 1,
                        multiHypoCorrection = FALSE,
                        #backgroundGenes = rownames(cosinor_mod_sig),
                        tissueSpecificGeneType = 2)
gs2output

gs2seEnrichmentOutput<-gs2output[[1]]
gs2enrichmentOutput<-setNames(data.frame(assay(gs2seEnrichmentOutput),
                                         row.names = rowData(gs2seEnrichmentOutput)[,1]), 
                              colData(gs2seEnrichmentOutput)[,1])
gs2enrichmentOutput$Tissue<-row.names(gs2enrichmentOutput)
head(gs2enrichmentOutput)
gs2enrichmentOutput$cluster<-"c2"
gs2enrichmentOutput


########################################################################################
graphics.off()
ggplot(gs2enrichmentOutput,aes(x=reorder(Tissue,
                                      -Log10PValue),
                            y=Log10PValue,
                            label = Tissue.Specific.Genes,
                            fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())


####################################################################################################
gs3<-GeneSet(geneIds=rownames(g3),geneIdType=ENSEMBLIdentifier(),
             organism="Homo Sapiens")
gs3

gs3output<-teEnrichment(inputGenes = gs3,rnaSeqDataset = 1,
                        multiHypoCorrection = FALSE,
                        #backgroundGenes = rownames(cosinor_mod_sig),
                        tissueSpecificGeneType = 2)
gs3output

gs3seEnrichmentOutput<-gs3output[[1]]
gs3enrichmentOutput<-setNames(data.frame(assay(gs3seEnrichmentOutput),
                                         row.names = rowData(gs3seEnrichmentOutput)[,1]), 
                              colData(gs3seEnrichmentOutput)[,1])
gs3enrichmentOutput$Tissue<-row.names(gs3enrichmentOutput)
head(gs3enrichmentOutput)
gs3enrichmentOutput$cluster<-"c3"
gs3enrichmentOutput


graphics.off()
ggplot(gs3enrichmentOutput,aes(x=reorder(Tissue,
                                         -Log10PValue),
                               y=Log10PValue,
                               label = Tissue.Specific.Genes,
                               fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())

############################################################################################

gs4<-GeneSet(geneIds=rownames(g4),geneIdType=ENSEMBLIdentifier(),
             organism="Homo Sapiens")
gs4

gs4output<-teEnrichment(inputGenes = gs4,rnaSeqDataset = 1,
                        multiHypoCorrection = FALSE,
                        #backgroundGenes = rownames(cosinor_mod_sig),
                        tissueSpecificGeneType = 2)
gs4output

gs4seEnrichmentOutput<-gs4output[[1]]
gs4enrichmentOutput<-setNames(data.frame(assay(gs4seEnrichmentOutput),
                                         row.names = rowData(gs4seEnrichmentOutput)[,1]), 
                              colData(gs4seEnrichmentOutput)[,1])
gs4enrichmentOutput$Tissue<-row.names(gs4enrichmentOutput)
head(gs4enrichmentOutput)
gs4enrichmentOutput$cluster<-"c4"
gs4enrichmentOutput


bind_enrichment<-rbind(gs1enrichmentOutput, gs2enrichmentOutput, gs3enrichmentOutput, gs4enrichmentOutput)
head(bind_enrichment)
library(ggthemes)
library(Seurat)
#####################################################################################################3

jpeg(filename = "pseudotime_analysisV3/Fig.3/HPA_dataset_enrichment.jpeg",
     units = "in",width = 10,height = 12,res = 300)
ggplot(bind_enrichment,aes(x=reorder(Tissue,
                                         -Log10PValue,decreasing=FALSE),
                               y=Log10PValue,
                               label = Tissue.Specific.Genes,
                               fill = Tissue))+
  geom_bar(stat = 'identity',color = "black")+
  labs(x='', y = '-log10(P-value)')+
  theme_bw()+
  theme(legend.position="none")+
  theme_stata()+
  theme(text = element_text(size = 20))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        #panel.grid.major= element_blank(),
       # panel.grid.minor = element_blank()
        )+
  facet_wrap(~cluster,ncol = 1)+#coord_flip()+
  xlab("Human Protein Atlas")+
  NoLegend()+
  geom_hline(yintercept = 1.3, linetype="dashed")
dev.off()

ls("package:ggthemes")

###############################################################################################################
cpm_loess<-read.xlsx("pseudotime_analysisV3/Fig.3/CPM_loess_Retinopathy.xlsx", rowNames = TRUE)
library(tidyverse)
hpa<-read_tsv("~/Documents/phd.upgrading/human_acute_malaria/after_yearone_viva/edgeR/rna_tissue_consensus.tsv")
head(hpa)



hpa_wide_format<-pivot_wider(data = hpa, id_cols = c("Gene","Gene name"),
                             values_from = "NX",
                             names_from = c("Tissue"))

hpa_wide_format<-as.data.frame(hpa_wide_format)
rownames(hpa_wide_format)<-hpa_wide_format$Gene
hpa_wide_format<-hpa_wide_format[,-1]
head(hpa_wide_format)


#####################################################################3
hpa_wide_format_sel<-hpa_wide_format[rownames(cosinor_mod_sig),]
head(hpa_wide_format_sel)
hpa_wide<-hpa_wide_format_sel[is.na(hpa_wide_format_sel)] <- 0
head(hpa_wide_format_sel)



annot<-read.xlsx("~/Documents/phd.upgrading/human_acute_malaria/robust_regression/HPA_consensus_wide_format.xlsx",
                 sheet = 2)
head(annot)
names(hpa_wide_format_sel)<-make.names(colnames(hpa_wide_format_sel))
dim(hpa_wide_format_sel)


scale<-t(scale(t(hpa_wide_format_sel[,-c(1)])))
head(scale)
scale<-na.omit(scale)
dim(scale)
merge<-merge(hpa_wide_format_sel[,1, drop=FALSE], 
             scale, by="row.names")
head(merge)
dim(merge)
annot
annot$legend<-relevel(as.factor(annot$legend),ref = "others")

designhpa<-model.matrix(~legend, annot)
colnames(designhpa)
library(limma)
fit<-lmFit(log2(hpa_wide_format_sel[,-c(1)]+0.1), 
           design = designhpa)
fit<-eBayes(fit, robust = TRUE)
resultshpa<-topTable(fit = fit, coef = 2, number = Inf)
head(resultshpa)
dim(resultshpa)
go<-subset(resultshpa,resultshpa$P.Value<0.05 & resultshpa$logFC>0)
dim(go)

merge_brain_genes<-merge(go,cosinor_mod_sig, by="row.names")
head(merge_brain_genes)


write.xlsx(merge_brain_genes,file = "pseudotime_analysisV3/Fig.3/Braingenes_Cluster3.xlsx")

library(illuminaHumanv4.db)
goenrich<-enrichGO(gene = rownames(go),
                   OrgDb = illuminaHumanv4.db,
                   keyType = "ENSEMBL",ont = "BP",
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05)
goenrich<-setReadable(x = goenrich,OrgDb = illuminaHumanv4.db,keyType = "ENSEMBL")

df<-as.data.frame(goenrich)
df



write.xlsx(df,file = "pseudotime_analysisV3/Fig.3/Goterm_brainGenes_cluster3.xlsx")

df

##############################################################################################################
merge_meta_heatmap<-read.xlsx("pseudotime_analysisV3/Fig.3/Heatmap_metadata.xlsx",rowNames = TRUE)


merge_meta_heatmap<-merge_meta_heatmap[order(merge_meta_heatmap$Trajectory_slingshot),]
merge_meta_heatmap
cosinor_mod_sig<-cosinor_mod_sig[order(cosinor_mod_sig$acrophase),]


brain_genesC3<-read.xlsx("pseudotime_analysisV3/Fig.3/Goterm_brainGenes_cluster3.xlsx",sheet = 2)
brain_genesC3
brain_genesC3<-brain_genesC3[!duplicated(brain_genesC3$genes),]

cosinor_mod_sig2<-cosinor_mod_sig[cosinor_mod_sig$symbol%in%brain_genesC3$genes,]

cpm_loess_sig<-cpm_loess[rownames(cosinor_mod_sig2),][,rownames(merge_meta_heatmap)]
cpm_loess_sig
dim(cpm_loess_sig)

library(circlize)
col_funp = colorRamp2(c(-2, 0, 2), c("green", "black", "purple"))



Heatmap(t(scale(t(cpm_loess_sig))),cluster_columns = FALSE,
        cluster_rows = FALSE, show_column_names = FALSE,
        show_row_names = TRUE, row_labels = cosinor_mod_sig2$symbol,
        use_raster = FALSE, 
        col = col_funp, row_split=brain_genesC3$phenotype)

library(circlize)
library(viridis)
metada_noDup$cohort

colret = list(condition = c("CM-R⁺"="firebrick","CC"="gold",
                            "CM-R⁻"="steelblue"),
              cohort = c("one"="bisque",
                         "two"="darkseagreen"),
              died = c("0"="lightblue1",
                       "1"="tan1"),
              acidosis = c("0"="lightblue",
                           "1"="tan1"),
              pseudotime = colorRamp2(c(0,4,8,12,16),viridis(n = 5,
                                                             option = "B",
                                                             begin = 0,end = 1)))

colret

merge_meta_heatmap$condition<-factor(as.factor(merge_meta_heatmap$condition),levels = c("CC","CM-R⁻","CM-R⁺"))

library(ComplexHeatmap)
ha1 = HeatmapAnnotation(pseudotime=merge_meta_heatmap$Trajectory_slingshot,
                        condition= merge_meta_heatmap$condition,
                        died=merge_meta_heatmap$OutcomeDied,
                        acidosis=merge_meta_heatmap$Metabolic.acidosis,
                        #cohort=merge_meta_heatmap$cohort,
                        annotation_name_side = "left",
                        show_legend =TRUE,col = colret,
                        show_annotation_name = TRUE,
                        na_col = "white",border = TRUE,
                        annotation_legend_param = list(legend_gp = gpar(fontsize = 20)),
                        simple_anno_size = unit(0.8, "cm"),
                        annotation_name_gp = gpar(fontsize = 20))
ha1

lab<-c("GFAP","NCAN","VEGFA","TGFB1","TGFBR1","AQP4","AQP9","GATA2",
       "BACE1","GDAP1","PLXNC1","HTR5A","GMFG","VCAN","PSEN1",
       "SLC1A2","CREB1","CREBRF","GRIA2","GRIA1","ROBO2",
       "HEPACAM",
       "SEMA5B","TRPM3","TMEM196","GFRA2", "S100A9","S100A8",
       "ELANE", "TBX21","AZU1","DEFA3","CTSG","CTSA")
mark_at = which(cosinor_mod_sig$symbol %in% lab)
mark_at

har = rowAnnotation(foo = anno_mark(at = which(cosinor_mod_sig$symbol %in% lab),
                                    labels = cosinor_mod_sig$symbol[cosinor_mod_sig$symbol%in%lab]))


library(ComplexHeatmap)
graphics.off()
####################################################################################
jpeg(filename = "pseudotime_analysisV3/Fig.3/Pseudotime_htm_retinopathy_combinedV6.jpeg",
     units = "in",width = 8.5,height = 8.5,res = 350)
htm<-Heatmap(t(scale(t(cpm_loess_sig))),
             clustering_method_rows = "single",
             clustering_method_columns = "complete",
             clustering_distance_columns = "canberra",
             clustering_distance_rows = "canberra",
             heatmap_legend_param = list(title = "mean z score", 
                                         direction="vertical"),
             top_annotation = ha1,
             show_row_names = F, 
             right_annotation = har,
             #row_labels = cosinor_mod_sig2$symbol,
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             cluster_row_slices = FALSE, 
             show_row_dend = FALSE,
             #column_split = metada_noDup_noCC$study,
             col = col_funp,
             column_title = NULL,
             column_labels = c(1:85),
             row_title_rot = 0,
             column_names_rot = 0,column_names_side = "top",
             show_column_names = FALSE,
             row_gap = unit(0.6, "mm"),
             row_split = cosinor_mod_sig$cluster,
             row_names_gp = gpar(fontsize = 10),
             row_names_max_width = unit(5,"cm"),
             column_title_gp = gpar(fontsize = 20),
             row_title_gp = gpar(fontsize = 20),
             use_raster = FALSE)
htm
htm_Pj<-draw(htm,
             merge_legend = TRUE, 
             heatmap_legend_side = "right",
             legend_title_gp = gpar(fontsize = 20),
             annotation_legend_side = "top",
             show_heatmap_legend = T)
dev.off()


cpm_loess

library(decoupleR)
library(progeny)
library(illuminaHumanv4.db)
ls("package:decoupleR")

net <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))
net

cpm_loess2<-cpm_loess[rownames(cosinor_mod_sig),]
cpm_loess2$symbol<-mapIds(x = illuminaHumanv4.db,
                          keys = rownames(cpm_loess2), 
                          column = "SYMBOL",keytype = "ENSEMBL")

library(clustermole)

markers <- clustermole_markers(species = "hs")
markers

markers_list <- split(x = markers$gene, f = markers$celltype_full)
g3<-subset(cosinor_mod_sig, cluster=="c2")
g3

my_overlaps <- clustermole_overlaps(genes = g3$symbol, species = "hs")
as.data.frame(my_overlaps)

my_overlaps_humantissue

library(dplyr)

my_overlaps_humantissue<-subset(my_overlaps, my_overlaps$db=="TISSUES" & my_overlaps$species=="Human")
as.data.frame(my_overlaps_humantissue)




library(rcellmarker)
res<-cells(g3$symbol,species = "human",keytype ="SYMBOL",pvalue = 0.25)
head(res)


library(clusterProfiler)
cosinor_mod_sig$entrezid<-mapIds(x = illuminaHumanv4.db,
                                 keys = rownames(cosinor_mod_sig),
                                 keytype = "ENSEMBL",column = "ENTREZID")

geneList<-scale(cosinor_mod_sig$acrophase)
names(geneList)<-cosinor_mod_sig$entrezid
geneList<-sort(geneList,decreasing = TRUE)

keggenrich<-gseKEGG(geneList = geneList,
                    organism = "hsa")


###################################################################
library(readr)

rna<-read_tsv(file = "pseudotime_analysisV3/Fig.4/rna_tissue_hpa.tsv")
head(rna)


metahpa<-read_tsv(file = "pseudotime_analysisV3/Fig.4/rna_tissue_hpa_description.tsv")

table(metahpa$Organ)
head(metahpa)
head(rna)

#metahpa<-as.data.frame(metahpa)
rownames(metahpa)<-metahpa$Tissue
head(metahpa)

merge_rna<-merge(metahpa,rna, by="Tissue")
head(merge_rna)

# dplyr
require(dplyr)
df<-merge_rna %>% group_by(Gene,Organ) %>% summarise(Value = max(nTPM))

df



library(tidyverse)
######################################################################
hpa_wide_format<-pivot_wider(data = df, id_cols = c("Gene"),
                             values_from = "Value",
                             names_from = c("Organ"))

hpa_wide_format<-as.data.frame(hpa_wide_format)
rownames(hpa_wide_format)<-hpa_wide_format$Gene
hpa_wide_format<-hpa_wide_format[,-c(1)]
head(hpa_wide_format)

#######################################################################
hpa_wide_format$max_column<-colnames(hpa_wide_format)[max.col(hpa_wide_format,
                                                              ties.method="first")]
head(hpa_wide_format)
colnames(hpa_wide_format)
dim(hpa_wide_format)

hpa_wide_format$SNR<-(rowMax(as.matrix(hpa_wide_format[,c(1:13)])))/(rowSds(as.matrix(hpa_wide_format[,c(1:13)])))
table(hpa_wide_format$SNR>3, hpa_wide_format$max_column)

##########################################################################
colnames(hpa_wide_format)
hpa_wide_format_f<-hpa_wide_format[,-c(3,9)]
head(hpa_wide_format_f)

hpa_wide_format_f2<-hpa_wide_format_f%>%dplyr::filter(!grepl(pattern = "reproductive",max_column))
dim(hpa_wide_format_f2)
head(hpa_wide_format_f2)
colnames(hpa_wide_format_f2)
library(ENIGMA)
hpa_wide_format_f2<-subset(hpa_wide_format_f2, SNR>5)
dim(hpa_wide_format_f2)

head(cpm_loess)
ls("package:ENIGMA")
cts
##########################################################################
egm<-create_ENIGMA(bulk = as.matrix(edgeR::cpm(hs_ret_all_kallisto_genelevel$counts)),
                   ref = as.matrix((hpa_wide_format_f2[,1:11])),
                   ref_type = "aggre")
egm<-batch_correct(egm)
egm<-get_cell_proportion(egm,method = "RLR",max.iter = 1000)

all(rownames(egm@result_cell_proportion)==rownames(metada_noDup))
Heatmap(egm@result_cell_proportion, 
        cluster_rows=FALSE
        )

#load the Gtex grand mean
load("pseudotime_analysisV3/Fig.4/GTex_cpmByTissue.Rdata")
head(GTex_cpmByTissue)
Gtex<-GTex_cpmByTissue[,-1]
head(Gtex)

metadata<-data.frame(sample=colnames(Gtex), condition=colnames(Gtex))
head(metadata)
rownames(metadata)<-metadata$sample

library(tidyverse)
head(Gtex)

colnames(Gtex)

#########################################################################################
sobj<-CreateSeuratObject(counts = Gtex, 
                         meta.data = metadata)%>%SCTransform()

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

#DoHeatmap(sobj,features = markers_roc_noDup$gene)

library(ComplexHeatmap)
Heatmap(t(scale(t(Gtex[rownames(markers_roc_noDup),]))), 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        show_row_names = FALSE,
        row_title_rot = 0,use_raster = FALSE,
        row_split = markers_roc_noDup$cluster,
        #column_split = metadata$condition,
        )

library(illuminaHumanv4.db)
rownames(markers_roc_noDup)<-gsub(pattern = "\\..*", 
                                  replacement = "",rownames(markers_roc_noDup))
markers_roc_noDup$symbol<-mapIds(keys = rownames(markers_roc_noDup),
                                 x=illuminaHumanv4.db,
                                 column = "SYMBOL",keytype = "ENSEMBL")
head(markers_roc_noDup)

markers_roc_noDup$geneName<-mapIds(keys = rownames(markers_roc_noDup),
                                 x=illuminaHumanv4.db,
                                 column = "GENENAME",keytype = "ENSEMBL")
head(markers_roc_noDup)

library(openxlsx)
write.xlsx(markers_roc_noDup, rowNames=TRUE,
           file = "pseudotime_analysisV3/Fig.4/GTEX_signature_FoldChange_BimodSeurat.xlsx")


markers_reload<-read.xlsx("pseudotime_analysisV3/Fig.4/GTEX_signature_FoldChange_BimodSeurat.xlsx", rowNames = TRUE)

#select the genes for the tissues I need from the set of markers
colnames(Gtex)
library(tidyverse)
markers_reload1<-markers_reload%>%dplyr::filter(grepl(pattern = "Blood|Brain",cluster))

dim(markers_reload1)

#####################################################################################################################
rownames(markers_reload)<-gsub(pattern = "\\..*", 
                                  replacement = "",rownames(markers_reload))
rownames(Gtex)<-gsub(pattern = "\\..*",
                     replacement = "",x = rownames(Gtex))
head(Gtex)
colnames(Gtex)


Gtex_signature<-Gtex[rownames(markers_reload1),][,c(7,4)]
head(Gtex_signature)
dim(Gtex_signature)

head(cpm)
colnames(Gtex_signature)

cpm_loess_sig[cpm_loess_sig<0]<-0

library(ENIGMA)
ls("package:ENIGMA")
ega<-create_ENIGMA(bulk = as.matrix((cpm)), 
                   ref = as.matrix((Gtex_signature)),
                   ref_type = "aggre")

#ega<-batch_correct(ega)
ega<-get_cell_proportion(ega, method = "CBS",
                         nu.v = c(0.1,0.2,0.4,0.6,0.8,1))

Heatmap((ega@result_cell_proportion[rownames(merge_meta_heatmap),]),
        row_split = merge_meta_heatmap$condition,
        cluster_row_slices = FALSE,
        cluster_rows = FALSE)



#############################################################################################################################
cor(ega@result_cell_proportion[rownames(merge_meta_heatmap),],
    merge_meta_heatmap$Trajectory_slingshot,
    method = "spearman")




library(EPIC)
library(MCPcounter)
library(immunedeconv)
library(ConsensusTME)
library(BRETIGEA)
library(illuminaHumanv4.db)

markers_df_human_brain$ensembl<-mapIds(x = illuminaHumanv4.db,
                                       keys = markers_df_human_brain$markers,
                                       keytype = "SYMBOL", column = "ENSEMBL")

msigdbr_list = split(x = markers_df_human_brain$ensembl, 
                     f = markers_df_human_brain$cell)
msigdbr_list

###########################################################################################

head(cpm_loess)
library(circlize)
col_funhtm = colorRamp2(c(0, 0.45, 0.8), c("blue", "black", "gold"))

deconTME<-deconvolute_consensus_tme_custom(gene_expression_matrix = as.matrix((cpm)),
                                 signature_genes = msigdbr_list,
                                 stat_method = "ssgsea")
Heatmap(t((t(deconTME[,rownames(merge_meta_heatmap)]))), 
        #column_split = merge_meta_heatmap$condition, 
        col = col_funhtm, name = "proportion",
        show_column_names = FALSE,
        cluster_columns = FALSE)





deconTME<-as.data.frame(deconTME)
cor(t(deconTME[,rownames(merge_meta_heatmap)]),
    merge_meta_heatmap$Trajectory_slingshot,
    method = "spearman")

head(cosinor_mod_sig)


##################################################################################################
brain_marker_noDup
term2gene<-data.frame(term=brain_marker_$cluster, gene=brain_marker_noDup$gene)
term2gene<-term2gene[!is.na(term2gene$gene),]
term2gene

library(clusterProfiler)
c1<-subset(cosinor_mod_sig, cluster=="c1")
head(c1)

c1go<-enricher(gene = rownames(c1),TERM2GENE = term2gene,
               pAdjustMethod = "none",
               pvalueCutoff = 2,
               minGSSize = 5,
               universe = term2gene$gene,
               maxGSSize = 1000)
c1go@result$cluster="c1"
c1go@result

##################################################################################################
c2<-subset(cosinor_mod_sig, cluster=="c2")
head(c2)

c2go<-enricher(gene = rownames(c2),TERM2GENE = term2gene,
               pAdjustMethod = "none",
               pvalueCutoff = 2,
               minGSSize = 5,
               universe = term2gene$gene,
               maxGSSize = 1000)
c2go@result$cluster="c2"
c2go@result

#######################################################################################
c3<-subset(cosinor_mod_sig, cluster=="c3")
head(c3)

c3go<-enricher(gene = rownames(c3),TERM2GENE = term2gene,
               pAdjustMethod = "none",
               pvalueCutoff = 2,
               minGSSize = 5,
               universe = term2gene$gene,
               maxGSSize = 1000)
c3go@result$cluster="c3"
c3go@result

#########################################################################################
c4<-subset(cosinor_mod_sig, cluster=="c4")
head(c4)

c4go<-enricher(gene = rownames(c4),TERM2GENE = term2gene,
               pAdjustMethod = "none",
               pvalueCutoff = 2,
               minGSSize = 5,
               universe = term2gene$gene,
               maxGSSize = 1000)
c4go@result$cluster="c4"
c4go@result


rbindgo<-rbind(c4go@result, c3go@result, c2go@result, c1go@result)
head(rbindgo)
colnames(rbindgo)

library(ggthemes)


jpeg(filename = "pseudotime_analysisV3/Fig.4/brainCell_enrichment.jpeg",
     height = 10,width = 6,units = "in",res = 300)
ggplot(rbindgo, aes(ID,-log10(pvalue), fill=-log10(pvalue)))+
  geom_histogram(stat = "identity", color="black")+
  facet_wrap(~cluster,nrow = 4,scales = "fixed")+
  theme_stata(scheme = "s2color")+
  theme(text = element_text(size = 25))+
  geom_hline(yintercept = 1.3, linetype="dashed", color="firebrick")
dev.off()


jpeg(filename = "pseudotime_analysisV3/Fig.4/Deconvoluted_RNA_proportion.jpeg",
     units = "in",width = 10,height = 4,res = 400)
htm<-Heatmap(t((t(deconTME))),
             clustering_method_rows = "ward.D2",
             clustering_method_columns = "complete",
             clustering_distance_columns = "canberra",
             clustering_distance_rows = "spearman",
             heatmap_legend_param = list(title = "mean z score", 
                                         direction="vertical"),
             top_annotation = ha1,
             show_row_names = TRUE, 
             #row_labels = cosinor_mod_sig2$symbol,
             cluster_columns = FALSE,
             cluster_rows = TRUE,
             cluster_row_slices = FALSE, 
             show_row_dend = TRUE,
             #column_split = metada_noDup_noCC$study,
             col = col_funhtm,
             column_title = NULL,
             column_labels = c(1:85),
             row_title_rot = 90,
             column_names_rot = 0,column_names_side = "top",
             show_column_names = FALSE,
             row_gap = unit(0.3, "mm"),
             #row_split = brain_genesC3$phenotype,
             row_names_gp = gpar(fontsize = 20),
             row_names_max_width = unit(5,"cm"),
             column_title_gp = gpar(fontsize = 20),
             row_title_gp = gpar(fontsize = 20),
             use_raster = FALSE)
htm
htm_Pj<-draw(htm,
             merge_legend = TRUE, 
             heatmap_legend_side = "right",
             legend_title_gp = gpar(fontsize = 20, 
                                    fontface = "bold"),
             annotation_legend_side = "top",
             show_heatmap_legend = TRUE)
dev.off()


save(rbindgo, deconTME, cosinor_mod_sig, file = "pseudotime_analysisV3/Fig.4/DeconvolutionResults_ssgea_BretigeaDermanis_markers.Rdata")

library(xCell)
ls("package:xCell")


do.exLR_origin <- function(avdata.m, ref.m, nu.v) {
  map.idx <- match(rownames(ref.m), rownames(avdata.m))
  rep.idx <- which(is.na(map.idx) == FALSE)
  data2.m <- avdata.m[map.idx[rep.idx], ]
  ref2.m <- ref.m[rep.idx, ]
  est.ab.lm <- list()
  est.lm <- list()
  nui <- 1
  for (nu in nu.v) {
    est.m <- matrix(NA, nrow = ncol(data2.m), ncol = ncol(ref2.m))
    colnames(est.m) <- colnames(ref2.m)
    rownames(est.m) <- colnames(data2.m)
    est.ab.m <- matrix(NA, nrow = ncol(data2.m), ncol = ncol(ref2.m))
    colnames(est.ab.m) <- colnames(ref2.m)
    rownames(est.ab.m) <- colnames(data2.m)
    
    for (s in seq_len(ncol(data2.m))) {
      svm.o <- svm(x = ref2.m, y = data2.m[, s], scale = TRUE, type = "nu-regression", kernel = "linear", nu = nu)
      coef.v <- t(svm.o$coefs) %*% svm.o$SV
      coef.v[which(coef.v < 0)] <- 1*10^-10
      est.ab.m[s,] <- coef.v
      total <- sum(coef.v)
      coef.v <- coef.v/total
      est.m[s, ] <- coef.v
    }
    est.lm[[nui]] <- est.m
    est.ab.lm[[nui]] <- est.ab.m
    nui <- nui + 1
  }
  
  #### select best nu using RMSE
  rmse.m <- matrix(NA, nrow = ncol(avdata.m), ncol = length(nu.v))
  for (nui in seq_along(nu.v)) {
    reconst.m <- ref2.m %*% t(est.lm[[nui]])
    s <- seq_len(ncol(avdata.m))
    rmse.m[s, nui] <- sqrt(colMeans((data2.m[, s] - reconst.m[, s])^2))
    message(nui)
  }
  colnames(rmse.m) <- nu.v
  nu.idx <- apply(rmse.m, 1, which.min)
  estF.m <- est.m
  for (s in seq_len(nrow(estF.m))) {
    estF.m[s, ] <- est.lm[[nu.idx[s]]][s, ]
  }
  estF.ab.m <- est.ab.m
  for (s in seq_len(nrow(estF.ab.m))) {
    estF.ab.m[s, ] <- est.ab.lm[[nu.idx[s]]][s, ]
  }
  #selecting min RMSE
  rmse.min.value <- as.data.frame(apply(rmse.m, 1, min))
  rownames(rmse.min.value) <- colnames(avdata.m)
  colnames(rmse.min.value)[1] <- 'RMSE'
  #caculating PCC
  pearson.corr.value <- c()
  for (i in 1:ncol(data2.m)) {
    cor.index <- cor.test(data2.m[, i], reconst.m[, i])
    cor.p <- as.numeric(cor.index$estimate)
    pearson.corr.value <- c(pearson.corr.value, cor.p)
  }
  estF.m <- cbind.data.frame(estF.m, pearson.corr.value, rmse.min.value)
  return(list(estF = estF.m, est.ab.sum = estF.ab.m, nu = nu.v[nu.idx]))
}

##########################################################################################
#deconvolute the origin of RNA
library(e1071)
deconEVorigin<-do.exLR_origin(avdata.m = as.matrix(HCC[,-c(1,23)]),
                              ref.m = as.matrix(Gtex1_f[markers_roc_noDupv2$gene,]),
                              nu.v = c(0.05))

library(ComplexHeatmap)
Heatmap((deconEVorigin$est.ab.sum),cluster_rows = FALSE)

save(cpm,Gtex_signature, file = "pseudotime_analysisV3/Fig.4/cpm&Gtexsiganture.Rdata")

library(granulator)
load("cpm&Gtexsiganture.Rdata")
decon_SVR<-deconvolute(m = cpm,
                       sigMatrix = Gtex_signature,
                       methods = "svr",use_cores = 30)

save(decon_SVR, file="decon_SVR.Rdata")


Tissuematrix<-read.csv("pseudotime_analysisV3/Fig.4/Matrix_tissueEV_origin.csv",row.names = 1)
Tissuematrix


Tissuematrix$TissueMax<-colnames(Tissuematrix)[max.col(Tissuematrix,ties.method="first")]
Tissuematrix
dim(Tissuematrix)

Heatmap(t(scale(t(Tissuematrix[,-17]))), 
        row_split = Tissuematrix$TissueMax, 
        cluster_rows = FALSE, cluster_columns = FALSE)

cpm1<-sqrt(cpm+0.01)
cpm_scaled <- (cpm1 - mean(cpm1)) / sd(as.vector(cpm1))
head(cpm_scaled)

cpm_scaled_loess <- (cpm_loess - mean(cpm_loess)) / sd(as.vector(cpm_loess))
head(cpm_scaled_loess)


GTex_scaled <- (Gtex_signature - mean(Gtex_signature)) / sd(as.vector(Gtex_signature))
head(GTex_scaled)

tpm<-as.data.frame(hs_ret_all_kallisto_genelevel$abundance[rownames(cpm),])
class(tpm)
head(tpm)
#################################################################################################
deconEVorigin<-do.exLR_origin(avdata.m = t((t(edgeR::cpm(hs_ret_all_kallisto_genelevel$abundance[,rownames(metadata)])))),
                              ref.m = as.matrix(Tissuematrix[,1:16]),
                              nu.v = c(0.2))

Heatmap(deconEVorigin$est.ab.sum,
        row_split =  metadata$condition,
        cluster_rows = FALSE)

decon_SVR<-as.data.frame((deconEVorigin$est.ab.sum))
head(decon_SVR)
class(decon_SVR)
#decon_SVR$Leucocyte<-colSums(as.matrix(decon_SVRHem))
head(decon_SVR)
decon_SVR<-as.data.frame(t(decon_SVR))
  
library(matrixStats)
decon_SVR$mean_Fraction<-rowMeans(as.matrix(decon_SVR[,1:16]))
decon_SVR$sd<-rowSds(as.matrix(decon_SVR[,1:16]))

cor(t(decon_SVR[,-c(17,18)]), merge_meta_heatmap$Trajectory_slingshot,method = "spearman")

cor(t(decon_SVR[,-c(17,18)]),method = "spearman")

decon_SVR<-decon_SVR[order(decon_SVR$mean_Fraction),]

pct<-round((decon_SVR$mean_Fraction/sum(decon_SVR$mean_Fraction))*100,digits = 1)
decon_SVR


jpeg(filename = "pseudotime_analysisV3/Fig.4/Pie_TissueProportionV1.jpeg",
     width = 7,height = 6,units = "in",res = 300)
pie(decon_SVR$mean_Fraction,radius = 1,init.angle = 160,
    labels = paste(rownames(decon_SVR), sep = " ", pct, "%"),
    col = rainbow(length(decon_SVR$mean_Fraction),end = 0.8))
dev.off()


library(ggpubr)

decon_SVR$Tissue<-rownames(decon_SVR)

library(Seurat)
ggpie(decon_SVR, "mean_Fraction", count_type = "full",
      lab.pos = "out",label_type = "vertical",label_threshold = 1,
      label = paste(rownames(decon_SVR), sep = " ", pct, "%"),
      fill = "Tissue",labe_size=4)+NoLegend()

decon_SVR$all<-paste(rownames(decon_SVR), sep = " ", pct, "%")

ggpie(decon_SVR, "mean_Fraction",group_key = "cut", count_type = "full",
      label = paste(rownames(decon_SVR), sep = " ", pct, "%"),
      label_type = "horizon", label_split = NULL,fill="Tissue",
      label_size = 2, label_pos = "out", label_threshold = 1)+NoLegend()
dim(decon_SVR)
decon_SVR$pos<-c(1:16)
library(ggthemes)

graphics.off()

ggplot(decon_SVR, aes(x = "" , y = mean_Fraction, fill = fct_inorder(Tissue))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Pastel1")+
  theme_void()+
  geom_label_repel(data =decon_SVR ,
                   aes(y = pos, label = paste0(round(mean_Fraction*100), "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE)
  


###########################################################################################################
hematoMatrix<-read.csv("pseudotime_analysisV3/Fig.4/Matrix_Hemopoietic.cell_EV_origin.csv",row.names = 1)
hematoMatrix

#################################################################################################

deconEVoriginHemato<-do.exLR_origin(avdata.m = as.matrix(t((t((tpm))))),
                              ref.m = as.matrix(t((t(hematoMatrix)))),
                              nu.v = c(0.1,0.25,0.5,0.75,1))

Heatmap(deconEVoriginHemato$est.ab.sum[rownames(merge_meta_heatmap),], cluster_rows = FALSE)
dim(decon_SVR)
cor(t(decon_SVR[,-89][,rownames(merge_meta_heatmap)]),merge_meta_heatmap$Trajectory_slingshot)

decon_SVRHem<-as.data.frame(t(deconEVoriginHemato$est.ab.sum))
head(decon_SVRHem)
class(decon_SVRHem)


library(matrixStats)
decon_SVRHem$mean_Fraction<-rowMeans(as.matrix(decon_SVRHem))
decon_SVRHem
dim(decon_SVRHem)

cor(t(decon_SVRHem[,-22]), merge_meta_heatmap$Trajectory_slingshot,method = "spearman")

Heatmap(cor(t(decon_SVRHem[,-22]),method = "spearman"))

decon_SVRHem<-decon_SVRHem[order(decon_SVRHem$mean_Fraction),]

pct<-round((decon_SVRHem$mean_Fraction/sum(decon_SVRHem$mean_Fraction))*100,digits = 1)
decon_SVRHem$mean_Fraction


jpeg(filename = "pseudotime_analysisV3/Fig.4/Pie_LM22ProportionV2.jpeg",
     width = 7,height = 6,units = "in",res = 300)
pie(decon_SVRHem$mean_Fraction,radius = 1,init.angle = 170,
    labels = paste(rownames(decon_SVRHem), sep = " ", pct, "%"),
    col = rainbow(length(decon_SVRHem$mean_Fraction),end = 0.9))
dev.off()



save(decon_SVR, decon_SVRHem, file = "pseudotime_analysisV3/Fig.4/DeconvlouteAbsoluteRNA_amounts.Rdata")


write.xlsx(rbindgo,file = "pseudotime_analysisV3/Fig.4/BrainEnrichment_fishers_test.xlsx",rowNames=TRUE)

######################################################################################################################
LM22$CellMax<-colnames(LM22)[max.col(LM22,ties.method="first")]
LM22

msigdbr_list = split(f = LM22$CellMax, 
                     x = rownames(LM22))
msigdbr_list


###########################################################################################
library(immunedeconv)

deconTME<-deconvolute_consensus_tme_custom(gene_expression_matrix = as.matrix((cts)),
                                           signature_genes = msigdbr_list,
                                           stat_method = "ssgsea")
Heatmap(t(scale(t(deconTME[,rownames(merge_meta_heatmap)]))), 
        #column_split = merge_meta_heatmap$condition, 
        #col = col_funhtm, 
        name = "proportion",
        show_column_names = FALSE,
        cluster_columns = FALSE)

#cor(t(deconTME[,rownames(merge_meta_heatmap)]), merge_meta_heatmap$Trajectory_slingshot,method = "spearman")

#####################################################################
library(scRNAseq)

DarmanisData<-DarmanisBrainData(ensembl = TRUE,
                                location = FALSE,
                                remove.htseq = TRUE)
class(DarmanisData)

assay(DarmanisData)
table(DarmanisData$cell.type)
library(Rmagic)
library(Seurat)
library(tidyverse)


pdata_darmanis<-as.data.frame(DarmanisData@colData)
head(pdata_darmanis)
pdata_darmanis<-pdata_darmanis%>%dplyr::filter(!grepl(pattern = "fetal|hybrid",cell.type))

pdata_darmanis$cell.type

dim(assay(DarmanisData)[,rownames(pdata_darmanis)])


sobj<-CreateSeuratObject(counts = assay(DarmanisData)[,rownames(pdata_darmanis)],
                         meta.data = pdata_darmanis)
sobj<-magic(sobj,t = 10)%>%SCTransform()
sobj@active.assay
sobj<-RunPCA(sobj,assay = "SCT")
DimPlot(sobj,group.by = "cell.type")
sobj<-RunUMAP(sobj,dims = 1:50,assay = "SCT")
DimPlot(sobj,group.by = "cell.type")
Idents(sobj)<-"cell.type"

brain_markers<-FindAllMarkers(sobj,only.pos = TRUE,
                              test.use = "bimod")

brain_marker_noDup<-as.data.frame(brain_markers %>% distinct %>%
                                   group_by(gene) %>% top_n(1, avg_log2FC))
rownames(brain_marker_noDup)<-brain_marker_noDup$gene
brain_marker_noDup<-brain_marker_noDup[order(brain_marker_noDup$gene),]
head(brain_marker_noDup)

table(brain_marker_noDup$cluster)
brain_marker_noDup$score<--log10(brain_marker_noDup$p_val)
brain_marker_top500<-as.data.frame(brain_marker_noDup %>% distinct %>%
                                    group_by(cluster) %>% top_n(500, avg_log2FC))
dim(brain_marker_top200)
head(brain_marker_top200)

signaturebrain<-AverageExpression(sobj,assays = "SCT",)$SCT[brain_marker_top200$gene,]
signaturebrain<-as.data.frame(signaturebrain)

library(ComplexHeatmap)
signaturebrain$CellMax<-colnames(signaturebrain)[max.col(signaturebrain,ties.method="first")]
signaturebrain$CellMax

Heatmap(t(scale(t(signaturebrain[,-7]))), 
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_title_rot = 0,
        cluster_row_slices = TRUE,
        row_split = signaturebrain$CellMax,
        show_row_names = FALSE)
#deconvolute the origin of RNA
library(e1071)

save(brain_marker_noDup, signaturebrain,brain_marker_top500 ,file = "pseudotime_analysisV3/Fig.4/DarmanisBrainSignature.Rdata")



deconEVorigin<-do.exLR_origin(avdata.m = (tpm),
                              ref.m = (as.matrix(signaturebrain[,-7])),
                              nu.v = c(0.1,0.25,0.5,0.75,1))

library(ComplexHeatmap)

jpeg(filename = "pseudotime_analysisV3/Fig.4/Deconvoluted_SVR_BrainV2.jpeg",
     units = "in",width = 15,height = 3,res = 300)
htm<-Heatmap(t(decon_SVRBrain[rownames(merge_meta_heatmap),][,1:6]),
        cluster_columns = FALSE,
        top_annotation = ha1,
        name = "proportion",
        #row_split = merge_meta_heatmap$condition,
        show_column_names = FALSE,
        column_title = NULL,
        column_labels = c(1:85),
        row_title_rot = 90,
        column_names_rot = 0,column_names_side = "top",
        row_gap = unit(0.3, "mm"),
        row_names_gp = gpar(fontsize = 20),
        row_names_max_width = unit(5,"cm"),
        column_title_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 20),
        use_raster = FALSE)
htm_Pj<-draw(htm,
             merge_legend = TRUE, 
             heatmap_legend_side = "right",
             legend_title_gp = gpar(fontsize = 20, 
                                    fontface = "bold"),
             annotation_legend_side = "top",
             show_heatmap_legend = TRUE)
dev.off()




cor(deconEVorigin$est.ab.sum[rownames(merge_meta_heatmap),],
    merge_meta_heatmap$Trajectory_slingshot,method = "spearman")

decon_SVRBrain<-as.data.frame(t(deconEVorigin$est.ab.sum))
head(decon_SVRBrain)
class(decon_SVRBrain)

save(decon_SVRBrain,file = "pseudotime_analysisV3/Fig.4/decon_SVRbrain.Rdata")

decon_SVRBrain<-decon_SVRBrain[,rownames(cm_meta)]

library(matrixStats)
decon_SVRBrain$mean_Fraction<-rowMeans(as.matrix(decon_SVRBrain))
decon_SVRBrain
dim(decon_SVRBrain)

graphics.off()
Heatmap(cor(t(decon_SVRBrain[,-89]),method = "spearman"))

decon_SVRBrain<-decon_SVRBrain[order(decon_SVRBrain$mean_Fraction),]

pct<-round((decon_SVRBrain$mean_Fraction/sum(decon_SVRBrain$mean_Fraction))*100,digits = 1)
decon_SVRBrain$mean_Fraction


jpeg(filename = "Fig.4/Pie_BrainProportion_CM.jpeg",
     width = 6,height = 4.5,units = "in",res = 400)
pie(decon_SVRBrain$mean_Fraction,radius = 1,init.angle = 190,
    labels = paste(rownames(decon_SVRBrain), sep = " ", pct, "%"),
    col = rainbow(length(decon_SVRBrain$mean_Fraction)))
dev.off()

decon_SVR_t<-as.data.frame(t(decon_SVR))
head(decon_SVR_t)

library(circlize)
col_funp = colorRamp2(c(0, 0.05, 0.1), c("blue", "white", "red"))


Heatmap(t((t(decon_SVRHem[,rownames(merge_meta_heatmap)]))),
        cluster_columns = FALSE, #col = col_funp,
        show_column_names = FALSE,cluster_rows = FALSE)


###############################################################################
#decon_SVRBrain

head(decon_SVR_t)
head(decon_SVRHem)
graphics.off()

htm1<-Heatmap(t(decon_SVRTissue[rownames(merge_meta_heatmap),]),
             cluster_columns = FALSE,
             #top_annotation = ha1,
             name = "absolute proportion",
             #row_split = merge_meta_heatmap$condition,
             show_column_names = FALSE,
             column_title = NULL,
             column_labels = c(1:85),
             row_title_rot = 90,
             column_names_rot = 0,column_names_side = "top",
             row_gap = unit(0.3, "mm"),
             row_names_gp = gpar(fontsize = 20),
             row_names_max_width = unit(10,"cm"),
             column_title_gp = gpar(fontsize = 20),
             row_title_gp = gpar(fontsize = 20),
             use_raster = FALSE)
htm1
htm<-Heatmap((decon_SVRHem[,rownames(merge_meta_heatmap)]),
              cluster_columns = FALSE,
              #top_annotation = ha1,
              name = "absolute proportion",
              #row_split = merge_meta_heatmap$condition,
              show_column_names = FALSE,
              column_title = NULL,
              column_labels = c(1:85),
              row_title_rot = 90,
              column_names_rot = 0,column_names_side = "top",
              row_gap = unit(0.3, "mm"),
              row_names_gp = gpar(fontsize = 20),
              row_names_max_width = unit(10,"cm"),
              column_title_gp = gpar(fontsize = 20),
              row_title_gp = gpar(fontsize = 20),
              use_raster = FALSE)
htm
htm2<-Heatmap((deconEVorigin2[,rownames(merge_meta_heatmap)]),
             cluster_columns = FALSE,
             #top_annotation = ha1,
             name = "absolute proportion",
             #row_split = merge_meta_heatmap$condition,
             show_column_names = FALSE,
             column_title = NULL,
             column_labels = c(1:85),
             row_title_rot = 90,
             column_names_rot = 0,column_names_side = "top",
             row_gap = unit(0.3, "mm"),
             row_names_gp = gpar(fontsize = 20),
             row_names_max_width = unit(10,"cm"),
             column_title_gp = gpar(fontsize = 20),
             row_title_gp = gpar(fontsize = 20),
             use_raster = FALSE)
htm2
htm3<-Heatmap((decon_SVRBrain[,rownames(merge_meta_heatmap)]),
              cluster_columns = FALSE,
              #top_annotation = ha1,
              name = "absolute proportion",
              #row_split = merge_meta_heatmap$condition,
              show_column_names = FALSE,
              column_title = NULL,
              column_labels = c(1:85),
              row_title_rot = 90,
              column_names_rot = 0,column_names_side = "top",
              row_gap = unit(0.3, "mm"),
              row_names_gp = gpar(fontsize = 20),
              row_names_max_width = unit(10,"cm"),
              column_title_gp = gpar(fontsize = 20),
              row_title_gp = gpar(fontsize = 20),
              use_raster = FALSE)




htm3


jpeg(filename = "pseudotime_analysisV3/Fig.4/Deconvoluted_SVR_heatmapALLV3.jpeg",
     units = "in",width = 17,height = 19,res = 300)
htm_Pj<-draw(htm2 %v%htm1 %v% htm3 %v% htm ,
             merge_legend = TRUE, 
             heatmap_legend_side = "right",
             legend_title_gp = gpar(fontsize = 20, 
                                    fontface = "bold"),
             annotation_legend_side = "top",
             show_heatmap_legend = TRUE)
dev.off()

decon_SVRBrain
dim(decon_SVRBrain)
dim(merge_meta_heatmap)

########################################################################################################################
library(metan)
decon_SVRBrain<-as.data.frame(t(decon_SVRBrain[,rownames(merge_meta_heatmap)]))
decon_SVRBrain$Pseudotime<-merge_meta_heatmap$Trajectory_slingshot
results_cor<-as.data.frame(corr_ci(.data = decon_SVRBrain,sel.var = "Pseudotime"))

library(ggforestplot)
forestplot(df = results_cor,estimate = Corr,name = V1,se = CI)

write.xlsx(results_cor, file = "pseudotime_analysisV3/Fig.4/BrainRNA_decon_cor_pseudotimeV2.xlsx")

results_cor<-read.xlsx("pseudotime_analysisV3/Fig.4/BrainRNA_decon_cor_pseudotimeV2.xlsx")
head(results_cor)
##################################################################################################################
jpeg(filename = "pseudotime_analysisV3/Fig.4/ForestbrainRNA_vs_pseudotimeV2.jpeg",
     width = 4.5,height = 3,units = "in",res = 600)
results_cor |> 
  ggplot(aes(x = reorder(V1, Corr, decreasing=TRUE), y = Corr,colour=legend,
             ymin = LL, ymax = UL)) +
  geom_vline(aes(xintercept = V1), col = "grey95", size = 5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(size=2) +
  geom_linerange(size=1) +
  ylab("Rho (95% CI)")+
  coord_flip() +
  labs(x = NULL) +
  theme_stata()+
  theme(plot.margin = margin(t = 5, r = -4, b = 5, l = 5),
        axis.text.y = element_text(angle=0),
        text = element_text(size=20))+
  scale_color_manual(values = c("steelblue", "black","firebrick"))+
  scale_y_continuous(breaks = c(-0.5,0,0.5))+
  NoLegend()
dev.off()  

#######################################################################################################
decon_SVR_t<-decon_SVR_t[rownames(merge_meta_heatmap),]

decon_SVR_t$Pseudotime<-merge_meta_heatmap$Trajectory_slingshot
results_cor<-as.data.frame(corr_ci(.data = decon_SVR_t,sel.var = "Pseudotime"))

forestplot(df = results_cor,estimate = Corr,name = V1,se = CI)

write.xlsx(results_cor, file = "pseudotime_analysisV3/Fig.4/Tissue_decon_cor_pseudotime.xlsx")

#results_cor<-read.xlsx("pseudotime_analysisV3/Fig.4/BrainRNA_decon_cor_pseudotime.xlsx")
#head(results_cor)
##################################################################################################################
jpeg(filename = "pseudotime_analysisV3/Fig.4/ForestTissueRNA_vs_pseudotime.jpeg",
     width = 5,height = 6,units = "in",res = 300)
results_cor |> 
  ggplot(aes(x = reorder(V1, Corr, decreasing=TRUE), y = Corr,
             ymin = LL, ymax = UL)) +
  geom_vline(aes(xintercept = V1), col = "grey95", size = 5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point() +
  geom_linerange() +  ylab("Rho (95% CI)")+
  coord_flip() +
  labs(x = NULL) +
  theme_stata()+
  theme(plot.margin = margin(t = 5, r = -4, b = 5, l = 5),
        axis.text.y = element_text(angle=0),
        text = element_text(size=20))+
  #scale_color_manual(values = c("steelblue", "black","firebrick"))+
  scale_y_continuous(breaks = c(-0.25,0,0.25))+
  NoLegend()
dev.off()  


#######################################################################################################
decon_SVR_t<-decon_SVR_t[rownames(merge_meta_heatmap),]

decon_SVR_t$Pseudotime<-merge_meta_heatmap$Trajectory_slingshot
results_cor<-as.data.frame(corr_ci(.data = decon_SVR_t,sel.var = "Pseudotime"))

forestplot(df = results_cor,estimate = Corr,name = V1,se = CI)

write.xlsx(results_cor, file = "pseudotime_analysisV3/Fig.4/BrainRNA_decon_cor_pseudotime.xlsx")

results_cor<-read.xlsx("pseudotime_analysisV3/Fig.4/BrainRNA_decon_cor_pseudotime.xlsx")
head(results_cor)
##################################################################################################################
jpeg(filename = "pseudotime_analysisV3/Fig.4/ForestTissueRNA_vs_pseudotime.jpeg",
     width = 5,height = 6,units = "in",res = 600)
results_cor |> 
  ggplot(aes(x = reorder(V1, Corr, decreasing=TRUE), y = Corr,
             ymin = LL, ymax = UL)) +
  geom_vline(aes(xintercept = V1), col = "grey95", size = 5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point() +
  geom_linerange() +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme_stata()+
  theme(plot.margin = margin(t = 5, r = -4, b = 5, l = 5),
        axis.text.y = element_text(angle=0),
        text = element_text(size=20))+
  #scale_color_manual(values = c("steelblue", "black","firebrick"))+
  scale_y_continuous(breaks = c(-0.25,0,0.25))+
  NoLegend()
dev.off()  


#######################################################################################################
head(decon_SVRHem)
decon_SVRHem<-as.data.frame(t(decon_SVRHem))
head(decon_SVRHem)
decon_SVRHem<-decon_SVRHem[rownames(merge_meta_heatmap),]

decon_SVRHem$Pseudotime<-merge_meta_heatmap$Trajectory_slingshot
results_cor<-as.data.frame(corr_ci(.data = decon_SVRHem,sel.var = "Pseudotime"))

forestplot(df = results_cor,estimate = Corr,name = V1,se = CI)

write.xlsx(results_cor, file = "pseudotime_analysisV3/Fig.4/ImmuneCell_decon_cor_pseudotime.xlsx")

#results_cor<-read.xlsx("pseudotime_analysisV3/Fig.4/BrainRNA_decon_cor_pseudotime.xlsx")
#head(results_cor)
##################################################################################################################
jpeg(filename = "pseudotime_analysisV3/Fig.4/ForestImmuneRNA_vs_pseudotime.jpeg",
     width = 5,height = 9,units = "in",res = 600)
results_cor |> 
  ggplot(aes(x = reorder(V1, Corr, decreasing=TRUE), y = Corr,
             ymin = LL, ymax = UL)) +
  geom_vline(aes(xintercept = V1), col = "grey95", size = 5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point() +
  geom_linerange() +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme_stata()+
  theme(plot.margin = margin(t = 5, r = -4, b = 5, l = 5),
        axis.text.y = element_text(angle=0),
        text = element_text(size=20))+
  #scale_color_manual(values = c("steelblue", "black","firebrick"))+
  scale_y_continuous(breaks = c(-0.3,0,0.3))+
  NoLegend()
dev.off()  

###########################################################################################
resultsFisher<-read.xlsx("pseudotime_analysisV3/Fig.4/BrainEnrichment_fishers_test.xlsx",sheet = 2,rowNames = TRUE)
resultsFisher

library(circlize)
col_funp = colorRamp2(c(0, 0.3, 3), c("steelblue", "white", "red"))


jpeg(filename = "pseudotime_analysisV3/Fig.4/DeconvolutedbrainFisher_pvalue.jpeg",
     units = "in",width = 4.5,height = 4,res = 300)
htm_enrich<-Heatmap(-log10(resultsFisher[,1:4]), 
        cluster_rows = FALSE, 
        #name = "-log10(P-value)",
        col = col_funp,border = TRUE,
        cluster_columns = FALSE,
        column_title = NULL,
        row_title_rot = 90,
        column_names_rot = 0,
        column_names_side = "top",
        row_gap = unit(0.3, "mm"),
        row_names_gp = gpar(fontsize = 20),
        column_names_gp = gpar(fontsize = 20),
        row_names_max_width = unit(10,"cm"),
        column_title_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 20),
        use_raster = FALSE,
        column_names_centered = TRUE,
        row_names_side="left",
        heatmap_legend_param = list(title = "-log10 (P-value)", 
                                    direction="horizontal"))
htm_enrich

htm_Pj<-draw(htm_enrich,
             merge_legend = TRUE, 
             heatmap_legend_side = "top",
             legend_title_gp = gpar(fontsize = 20),
             annotation_legend_side = "top",
             show_heatmap_legend = TRUE)
htm_Pj
dev.off()

head(cpm_loess)


##################################################################################################################
neuralgenes<-c("ENSG00000136040","ENSG00000168546","ENSG00000104381","ENSG00000173452","ENSG00000157219")
namesNeural<-c("PLXNC1","GFRA2", "GDAP1","TMEM196","HTR5A")

cpm_loess_neural<-cpm[neuralgenes,]
sqrt(cpm_loess_neural)
rownames(cpm_loess_neural)<-namesNeural
cpm_loess_neural_scale<-t(sqrt(t(cpm_loess_neural)))

merge_cpm_neural<-merge(t(cpm_loess_neural_scale), merge_meta_heatmap, by="row.names")
head(merge_cpm_neural)

merge_cpm_neural$cell.type="Neural"

############################################################################################################################
glialgenes<-c("ENSG00000130287","ENSG00000131095","ENSG00000171885","ENSG00000130755","ENSG00000038427")
namesglial<-c("NCAN","GFAP", "AQP4","GMFG","VCAN")

cpm_loess_g<-cpm[glialgenes,]
sqrt(cpm_loess_g)
rownames(cpm_loess_g)<-namesglial
cpm_loess_glial_scale<-t(sqrt(t(cpm_loess_g)))

merge_cpm_glial<-merge(t(cpm_loess_glial_scale), merge_meta_heatmap, by="row.names")
head(merge_cpm_glial)

merge_cpm_glial$cell.type="Glial"

###############################################################################

neuralMelt<-reshape2::melt(data = merge_cpm_neural,id.vars=-c(2:6))
head(neuralMelt)
glialMelt<-reshape2::melt(data = merge_cpm_glial,id.vars=-c(2:6))
head(glialMelt)
dim(neuralMelt)
dim(glialMelt)

melt_bind<-rbind(neuralMelt, glialMelt)
head(melt_bind)

library(ggforce)
library(ggthemes)
########################################################################################
glial<-ggplot(glialMelt, aes(Trajectory_slingshot, value))+
  facet_wrap(~variable,scales = "free_y",nrow = 1)+
  geom_smooth(se = FALSE,method = "loess",span=1, color='tan1')+
  theme_stata()+
  ggtitle("Glial cell markers")+
  ylab("Abundance")+xlab("Pseudotime")+
  theme(text = element_text(size = 20))+
  geom_vline(xintercept = 7,linetype="dashed")
  #scale_y_continuous(breaks = c(2,3.5,4.5,6,10,14,16,18,20))
glial
  
neural<-ggplot(neuralMelt, aes(Trajectory_slingshot, value))+
  facet_wrap(~variable,scales = "free_y",nrow = 1)+
  geom_smooth(se = FALSE,method = "loess",span=1, color='tan1')+
  theme_stata()+
  ggtitle("Neuron markers")+
  ylab("Abundance")+xlab("Pseudotime")+
  theme(text = element_text(size = 20))+
  geom_vline(xintercept = 7,linetype="dashed")
  #scale_y_continuous(breaks = c(2,3,4,6,8,10,12,14,16,18,19,20,21))
neural

jpeg(filename = "pseudotime_analysisV3/Fig.4/ExampleMarkers.jpeg",
     units = "in",width = 14,height = 10,res = 200)
neural/glial
dev.off()



max(merge_meta_heatmap$Trajectory_slingshot)

#####################################################################################################
library(openxlsx)
library(data.table)

ref_data<-fread("~/Documents/phd.upgrading/ENIGMA_results_retinopathyRNAseq/rna_single_cell_type.tsv")
ref_wide<-as.data.frame(pivot_wider(values_from = 4,
                                    names_from = 3,data =ref_data ))

meta_hpa<-read.xlsx("~/Documents/phd.upgrading/ENIGMA_results_retinopathyRNAseq/CellType_filtered_V2.xlsx",rowNames = TRUE,sheet = 2)
head(meta_hpa)
head(ref_wide)
rownames(ref_wide)<-ref_wide$Gene
meta_hpa$tissue
rownames(annot)<-annot$subtissue
sobjHpa<-CreateSeuratObject(ref_wide[,-c(1,2)],
                            meta.data = meta_hpa,min.cells = 1)

head(meta_hpa)
sobjHpa<-NormalizeData(sobjHpa)
sobjHpa<-ScaleData(sobjHpa)
sobjHpa
Idents(sobjHpa)<-"tissue4"
sobjHpa@active.assay
#######################################################################################################
#get markers from Larksson single cell HPA data
markers<-FindAllMarkers(sobjHpa,test.use = "bimod",only.pos = TRUE,
                        slot = "counts",
                        min.cells.group = 1,min.cells.feature = 1)
table(markers$cluster)
markers
marker_noDup<-as.data.frame(markers %>% distinct %>%
                                    group_by(gene) %>% top_n(1, avg_log2FC))
rownames(marker_noDup)<-marker_noDup$gene
marker_noDup<-marker_noDup[order(marker_noDup$gene),]
head(marker_noDup)

table(marker_noDup$cluster)
marker_noDup$score<--log10(marker_noDup$p_val)
marker_top500<-as.data.frame(marker_noDup %>% distinct %>%
                                     group_by(cluster) %>% top_n(200, avg_log2FC))
dim(marker_top500)
head(marker_top500)

sobjHpa@assays
###############################################################################################################
marker_top500<-marker_top500[order(marker_top500$cluster),]
signaturebrain<-AverageExpression(sobjHpa,slot = "counts")$RNA[marker_top500$gene,]
signaturebrain<-as.data.frame(signaturebrain)
head(signaturebrain)

Heatmap(t(scale(t(signaturebrain))),cluster_columns = FALSE,
        cluster_rows = FALSE,show_row_names = FALSE,use_raster = FALSE)

#tpm<-as.data.frame(hs_ret_all_kallisto_genelevel$abundance[rownames(cpm),])
htm1

class(tpm)
head(tpm)
deconEVorigin2<-do.exLR_origin(avdata.m = log2(HCC1[,-c(ncol(HCC1))]+0.01),
                              ref.m = (as.matrix(signaturebrain)),
                              nu.v = c(0.1,0.25,0.5,0.75,1))


deconEVorigin2$est.ab.sum[rownames(merge_meta_heatmap),]

##################################################################################################################
Heatmap(deconEVorigin2$est.ab.sum,#[,rownames(merge_meta_heatmap)],
        cluster_columns = FALSE,
        cluster_rows = FALSE, show_column_names = FALSE)

deconEVorigin2<-as.data.frame(t(deconEVorigin2$est.ab.sum))
head(deconEVorigin2)
class(deconEVorigin2)

deconEVorigin2<-decon_bloodTissue[,rownames(cm_meta)]
save(decon_bloodTissue,file = "pseudotime_analysisV3/Fig.4/decon_bloodTissue.Rdata")

library(matrixStats)
deconEVorigin2$mean_Fraction<-rowMeans(as.matrix(deconEVorigin2))
dim(deconEVorigin2)

graphics.off()
Heatmap(cor(t(deconEVorigin2[,-89]),method = "spearman"))

deconEVorigin2<-deconEVorigin2[order(deconEVorigin2$mean_Fraction),]

pct<-round((deconEVorigin2$mean_Fraction/sum(deconEVorigin2$mean_Fraction))*100,digits = 1)
deconEVorigin2$mean_Fraction

rownames(deconEVorigin2)<-factor(rownames(deconEVorigin2),levels = c("Tissue","Blood"))

jpeg(filename = "Fig.4/Pie_Blood_TissueProportion_CM.jpeg",
     width = 4.75,height = 4,units = "in",res = 400)
pie(deconEVorigin2$mean_Fraction,radius = 1,init.angle = 300,
    labels = paste(rownames(deconEVorigin2), sep = " ", pct, "%"),
    col = c("cyan", "firebrick"))
dev.off()


########################################################################################
decon_SVR_t<-as.data.frame(t(decon_SVR))
head(decon_SVR_t)

######################################################################33
meta.data = meta_hpa,min.cells = 1)

head(meta_hpa)
sobjHpa<-NormalizeData(sobjHpa)
sobjHpa<-ScaleData(sobjHpa)
sobjHpa
Idents(sobjHpa)<-"tissue3"
sobjHpa@active.assay
#######################################################################################################
#get markers from Larksson single cell HPA data
markers<-FindAllMarkers(sobjHpa,test.use = "bimod",only.pos = TRUE,
                        assay = "RNA",
                        min.cells.group = 1,min.cells.feature = 1)
table(markers$cluster)
markers
marker_noDup<-as.data.frame(markers %>% distinct %>%
                              group_by(gene) %>% top_n(1, avg_log2FC))
rownames(marker_noDup)<-marker_noDup$gene
marker_noDup<-marker_noDup[order(marker_noDup$gene),]
head(marker_noDup)

table(marker_noDup$cluster)
marker_noDup$score<--log10(marker_noDup$p_val)
marker_top500<-as.data.frame(marker_noDup %>% distinct %>%
                               group_by(cluster) %>% top_n(500, avg_log2FC))
dim(marker_top500)
head(marker_top500)

sobjHpa@assays
###############################################################################################################
marker_top500<-marker_top500[order(marker_top500$cluster),]
signaturebrain<-AverageExpression(sobjHpa)$RNA#[marker_top500$gene,]
signaturebrain<-as.data.frame(signaturebrain)
head(signaturebrain)

Heatmap(t((t(signaturebrain))),cluster_columns = FALSE,km=6,
        cluster_rows = FALSE,show_row_names = FALSE,use_raster = FALSE)

tpm<-as.data.frame(hs_ret_all_kallisto_genelevel$abundance[rownames(cpm),])


class(tpm)
head(tpm)
deconEVorigin2<-do.exLR_origin(avdata.m = sqrt(tpm),
                               ref.m = (as.matrix(signaturebrain)),
                               nu.v = c(0.1,0.25,0.5,0.75,1))


deconEVorigin2<-deconEVorigin2$est.ab.sum#[rownames(merge_meta_heatmap),]

##################################################################################################################
Heatmap(deconEVorigin2[,rownames(merge_meta_heatmap)],
        cluster_columns = FALSE,
        cluster_rows = FALSE, 
        show_column_names = FALSE)

deconEVorigin2<-as.data.frame(t(deconEVorigin2$est.ab.sum))
head(deconEVorigin2)
class(deconEVorigin2)

decon_bloodTissue<-deconEVorigin2
save(decon_bloodTissue,file = "pseudotime_analysisV3/Fig.4/decon_bloodTissue.Rdata")

library(matrixStats)
deconEVorigin2$mean_Fraction<-rowMeans(as.matrix(deconEVorigin2))
dim(deconEVorigin2)

graphics.off()
Heatmap(cor(t(deconEVorigin2[,-89]),method = "spearman"))

deconEVorigin2<-deconEVorigin2[order(deconEVorigin2$mean_Fraction),]

pct<-round((deconEVorigin2$mean_Fraction/sum(deconEVorigin2$mean_Fraction))*100,digits = 1)
deconEVorigin2$mean_Fraction

rownames(deconEVorigin2)<-factor(rownames(deconEVorigin2),levels = c("Tissue","Blood"))

jpeg(filename = "pseudotime_analysisV3/Fig.4/Pie_Blood_TissueProportion.jpeg",
     width = 4,height = 4,units = "in",res = 400)
pie(deconEVorigin2$mean_Fraction,radius = 1,init.angle = 170,
    labels = paste(rownames(deconEVorigin2), sep = " ", pct, "%"),
    col = c("cyan", "firebrick"))
dev.off()

########################################################################################
decon_SVR_t<-as.data.frame(t(decon_SVR))
head(decon_SVR_t)

####################################################################################
library(pathview)
library(illuminaHumanv4.db)
library(openxlsx)
results<-read.xlsx("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/pseudotime_analysisV3/Fig.3/Cosinor_periodicGenes.xlsx")
head(results)
class(results)

results$`pseudotime-zscore`<-scale(results$acrophase)
head(results)

write.csv(results, quote = FALSE,
          file = "~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/pseudotime_analysisV3/Fig.3/Cosinor_periodicGenesV3.csv")

keytypes(illuminaHumanv4.db)
results$entrezid<-mapIds(x = illuminaHumanv4.db,column = "ENTREZID",keytype = "ENSEMBL",
                         keys = as.character(results$ensembleIDs))
results

geneList<-scale(results$acrophase)
names(geneList)<-results$entrezid
geneList<-sort(geneList,decreasing = TRUE)
geneList


library(clusterProfiler)
ego<-gseKEGG(geneList = na.omit(geneList),
             keyType = "uniprot",
             pAdjustMethod = "none",pvalueCutoff = 1,
             organism = "hsa")

############################################################################################3
ego@result
ego<-setReadable(x = ego,OrgDb = illuminaHumanv4.db,keyType = "UNIPROT")
dim(ego@result)

write.xlsx(as.data.frame(ego@result), "~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/pseudotime_analysisV3/Fig6/Gseakegg.xlsx")

head(results)


results_sel<-read.xlsx("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/pseudotime_analysisV3/Fig.3/Cosinor_periodicGenes.xlsx",sheet = 2)
head(results_sel)
dim(results_sel)

setdiff(results_sel$gene, results$symbol)

#########################################################################################################
ensembleID<-results[results$symbol%in%results_sel$gene,]
head(ensembleID)

dim(ensembleID)
head(ensembleID)

cpm_sel<-cpm[ensembleID$ensembleIDs,]
rownames(cpm_sel)<-ensembleID$symbol
cpm_sel

duplicated(results_sel$gene)
rownames(results_sel)<-results_sel$gene

head(results_sel)

cosinosel<-results[results$symbol%in%results_sel$gene,]

###########################################################################################
cpm_sel_scale<-t(scale(t(cpm_sel)))
merge_category<-merge(t(cpm_sel_scale), merge_meta_heatmap, by="row.names")
colnames(merge_category)

cpm_sel
library(ggthemes)
merge_category_melt<-melt(merge_category,
                          id.vars=-c(2:42),variable.name = "gene" ,
                          value.name = "CPM")
merge_category_melt<-merge(merge_category_melt, results_sel, by="gene")
head(merge_category_melt)
merge_category_melt$categ<-merge_category_melt$`function`

library(ggallin)


ampa<-subset(merge_category_melt,categ=="AMPAR receptors", group=gene, color=gene)
class(ampa)

as.factor(results_sel$`function`)
merge_category_melt$gene<-as.factor(merge_category_melt$gene)
jpeg(filename = "selected_genes.jpeg",
     width = 22,height = 15,units = "in",res = 150)
ggplot(subset(merge_category_melt,categ=="NMDA receptors") ,
       aes(Trajectory_slingshot, CPM,
                                color=gene))+
  geom_smooth(method = "loess",se = FALSE,
              span=0.8)+
  #geom_dl(aes(label=gene), method=list("last.points",stat="smooth",span=0.8, rot=30))
  #facet_wrap(~categ, scales = "free_y",ncol = 3)+
  theme_stata(scheme = "s2color")+
  theme(axis.text.y = element_text(angle=0),
        legend.position = "right")+
  ylab("Abundance")+xlab("Pseudotime")
  NoLegend()
dev.off()

#graphics.off()

library(broom)

jpeg(filename = "Fig6/selected_genesV2.jpeg",
     width = 20,height = 20,units = "in",res = 200)
merge_category_melt %>%
  {ggplot(., aes(Trajectory_slingshot, CPM,
                 color=gene, group=gene,label=gene)) + 
      geom_smooth(method = "loess",se = FALSE,
                  span=0.8) + 
      guides(color = F)+
      geom_text_repel(data = group_by(., gene,categ) %>% 
                        do(augment(loess(CPM~Trajectory_slingshot, span = 0.8,.))) %>% 
                        filter( Trajectory_slingshot== max(Trajectory_slingshot)),
                      aes(Trajectory_slingshot, .fitted), nudge_x = 0.7,size=4)+
      theme_stata(scheme = "s2color")+
      theme(axis.text.y = element_text(angle=0),
            legend.position = "right",
            text = element_text(size = 25))+
      ylab("Abundance")+xlab("Pseudotime")+
      facet_wrap(~categ,scales = "free_y",nrow = 3)
  }
dev.off()


########################################################################
library(openxlsx)
library(illuminaHumanv4.db)
cosinor_mod_sig<-read.xlsx(xlsxFile = "~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/pseudotime_analysisV3/Fig.3/Cosinor_periodicGenes.xlsx", 
                           rowNames = TRUE)


cosinor_mod_sig$entrezid<-mapIds(x = illuminaHumanv4.db,
                                 column = "ENTREZID",
                                 keys = rownames(cosinor_mod_sig),
                                 keytype = "ENSEMBL")

library(rbioapi)
#############################################################################
c1<-subset(cosinor_mod_sig, cluster=="c1")
c1
c1_allenBrain <- rba_enrichr(gene_list = c1$symbol,organism = "human",
                                     gene_set_library = "PanglaoDB_Augmented_2021",
                             regex_library_name = FALSE)
head(c1_allenBrain)
c1_allenBrain$module="c1"
colnames(c1_allenBrain)<-paste(colnames(c1_allenBrain),
                               "c1", sep = "_")
colnames(c1_allenBrain)[1]<-"Term"
c1_allenBrain

#################################################################
c2<-subset(cosinor_mod_sig, cluster=="c2")
c2
c2_allenBrain <- rba_enrichr(gene_list = c2$symbol,
                             organism = "human",
                             gene_set_library = "PanglaoDB_Augmented_2021",
                             regex_library_name = FALSE)
head(c2_allenBrain)
colnames(c2_allenBrain)<-paste(colnames(c2_allenBrain),
                                    "c2", sep = "_")
colnames(c2_allenBrain)[1]<-"Term"



##########################################################
c3<-subset(cosinor_mod_sig, cluster=="c3")
c3
c3_allenBrain <- rba_enrichr(gene_list = c3$symbol,organism = "human",
                             gene_set_library = "PanglaoDB_Augmented_2021",
                             regex_library_name = FALSE)
head(c3_allenBrain)

colnames(c3_allenBrain)<-paste(colnames(c3_allenBrain),
                               "c3", sep = "_")
colnames(c3_allenBrain)[1]<-"Term"

####################################################################
c4<-subset(cosinor_mod_sig, cluster=="c4")
c4
c4_allenBrain <- rba_enrichr(gene_list = c4$symbol,organism = "human",
                             gene_set_library = "PanglaoDB_Augmented_2021",
                             regex_library_name = FALSE)
head(c4_allenBrain)

colnames(c4_allenBrain)<-paste(colnames(c4_allenBrain),
                               "c4", sep = "_")
colnames(c4_allenBrain)[1]<-"Term"
c4_allenBrain

###################################################################
library(tidyverse)
PanglaoB_enrichr<-list(c1_allenBrain, c2_allenBrain, c3_allenBrain, c4_allenBrain) %>% 
  reduce(left_join, by = "Term")
head(PanglaoB_enrichr)


write.xlsx(PanglaoB_enrichr,
           file = "PanglaoDB_enrichr.xlsx")

df<-read.xlsx("PanglaoDB_enrichr.xlsx",sheet = 2)
df_sig<-subset(df, sig=="sig")

library(circlize)
col_funp = colorRamp2(c(0, 1, 3), c("dodgerblue", "white", "red"))

anno = anno_mark(at = c(1:24), 
                 labels = df_sig$symbol_annot[1:24], 
                 which = "row",labels_gp = gpar(fontsize=25))

##################################################################
jpeg(filename = "PangaoDB_enrichment.jpeg",
     units = "in",width = 10,
     height = 16,res = 250)
htm<-Heatmap(-log10(df_sig[,2:5]),
             clustering_method_rows = "complete",
             clustering_method_columns = "complete",
             clustering_distance_columns = "canberra",
             clustering_distance_rows = "pearson",
        cluster_rows = T,
        cluster_columns = F,
        cluster_row_slices = F,
        col=col_funp,row_title = NULL,
        row_split = df_sig$group,
        row_labels = df_sig$Term,
        column_names_side = "top",
        column_names_rot = 0,
        heatmap_legend_param = list(title = "-log10(P-value)", 
                                    direction="horizontal"),
        #name = "-log10(P-value)",
        show_row_dend = FALSE,
        row_gap = unit(0, "mm"),
        column_gap = unit(0, "mm"),
        #row_labels = data$symbol_annot,
        cluster_column_slices = FALSE,
        column_names_gp = gpar(fontsize = 25),
        
        row_names_max_width = unit(5,"cm"),
        column_title_gp = gpar(fontsize = 15),
        row_title_gp = gpar(fontsize = 15),
        row_title_rot = 90)

htm_Pj<-draw(htm+rowAnnotation(mark = anno),
             merge_legend = TRUE, 
             heatmap_legend_side = "top", 
             annotation_legend_side = "top",
             show_heatmap_legend = TRUE,
             legend_title_gp = gpar(fontsize = 25))
dev.off()
################################################################

#kegg analysis
#########################################################################

library(rbioapi)
#############################################################################
c1<-subset(cosinor_mod_sig, cluster=="c1")
c1
c1_allenBrain <- rba_enrichr(gene_list = c1$symbol,organism = "human",
                             gene_set_library = "KEGG_2021_Human",
                             regex_library_name = FALSE)
head(c1_allenBrain)
c1_allenBrain$module="c1"
colnames(c1_allenBrain)<-paste(colnames(c1_allenBrain),
                               "c1", sep = "_")
colnames(c1_allenBrain)[1]<-"Term"
c1_allenBrain

#################################################################
c2<-subset(cosinor_mod_sig, cluster=="c2")
c2
c2_allenBrain <- rba_enrichr(gene_list = c2$symbol,
                             organism = "human",
                             gene_set_library = "KEGG_2021_Human",
                             regex_library_name = FALSE)
head(c2_allenBrain)
colnames(c2_allenBrain)<-paste(colnames(c2_allenBrain),
                               "c2", sep = "_")
colnames(c2_allenBrain)[1]<-"Term"



##########################################################
c3<-subset(cosinor_mod_sig, cluster=="c3")
c3
c3_allenBrain <- rba_enrichr(gene_list = c3$symbol,organism = "human",
                             gene_set_library = "KEGG_2021_Human",
                             regex_library_name = FALSE)
head(c3_allenBrain)

colnames(c3_allenBrain)<-paste(colnames(c3_allenBrain),
                               "c3", sep = "_")
colnames(c3_allenBrain)[1]<-"Term"

####################################################################
c4<-subset(cosinor_mod_sig, cluster=="c4")
c4
c4_allenBrain <- rba_enrichr(gene_list = c4$symbol,organism = "human",
                             gene_set_library = "KEGG_2021_Human",
                             regex_library_name = FALSE)
head(c4_allenBrain)

colnames(c4_allenBrain)<-paste(colnames(c4_allenBrain),
                               "c4", sep = "_")
colnames(c4_allenBrain)[1]<-"Term"
c4_allenBrain

#
###################################################################
library(tidyverse)
Kegg_enrichr<-list(c1_allenBrain, 
                   c2_allenBrain, 
                   c3_allenBrain, 
                   c4_allenBrain) %>% 
  reduce(left_join, by = "Term")
head(Kegg_enrichr)

write.xlsx(Kegg_enrichr, "Kegg_enrich.xlsx")


#############################################################################################
library(rbioapi)
rba_enrichr_libs()
cosinor_mod_sig<-read.xlsx(xlsxFile = "~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/pseudotime_analysisV3/Fig.3/Cosinor_periodicGenes.xlsx", 
                           rowNames = TRUE)


cosinor_mod_sig$entrezid<-mapIds(x = illuminaHumanv4.db,
                                 column = "ENTREZID",
                                 keys = rownames(cosinor_mod_sig),
                                 keytype = "ENSEMBL")

library(rbioapi)
#############################################################################
c1<-subset(cosinor_mod_sig, cluster=="c1")
c1
c1_allenBrain <- rba_enrichr(gene_list = c1$symbol,organism = "human",
                             gene_set_library = "WikiPathways_2019_Human",
                             regex_library_name = F)
#c1_allenBrain$WikiPathways_2019_Human
head(c1_allenBrain)
c1_allenBrain$module="c1"
colnames(c1_allenBrain)<-paste(colnames(c1_allenBrain),
                               "c1", sep = "_")
colnames(c1_allenBrain)[1]<-"Term"
c1_allenBrain

#################################################################
c2<-subset(cosinor_mod_sig, cluster=="c2")
c2
c2_allenBrain <- rba_enrichr(gene_list = c2$symbol,
                             organism = "human",
                             gene_set_library = "WikiPathways_2019_Human",
                             regex_library_name = FALSE)
head(c2_allenBrain)
colnames(c2_allenBrain)<-paste(colnames(c2_allenBrain),
                               "c2", sep = "_")
colnames(c2_allenBrain)[1]<-"Term"

library(rbioapi)
##########################################################
c3<-subset(cosinor_mod_sig, cluster=="c3")
c3
c3_allenBrain <- rba_enrichr(gene_list = c3$symbol,organism = "human",
                             gene_set_library = "WikiPathways_2019_Human",
                             regex_library_name = FALSE)
head(c3_allenBrain)

colnames(c3_allenBrain)<-paste(colnames(c3_allenBrain),
                               "c3", sep = "_")
colnames(c3_allenBrain)[1]<-"Term"

####################################################################
c4<-subset(cosinor_mod_sig, cluster=="c4")
c4
c4_allenBrain <- rba_enrichr(gene_list = c4$symbol,organism = "human",
                             gene_set_library = "WikiPathways_2019_Human",
                             regex_library_name = FALSE)
head(c4_allenBrain)

colnames(c4_allenBrain)<-paste(colnames(c4_allenBrain),
                               "c4", sep = "_")
colnames(c4_allenBrain)[1]<-"Term"
c4_allenBrain
c3_allenBrain

#################################################################################################################################
library(tidyverse)
PanglaoB_enrichr<-list(c1_allenBrain, c2_allenBrain, c3_allenBrain, c4_allenBrain) %>% 
  reduce(left_join, by = "Term")
head(PanglaoB_enrichr)


write.xlsx(PanglaoB_enrichr,
           file = "PanglaoDB_enrichr.xlsx")

###################################################################################################################################
library(clusterProfiler)
?enrichGO
ego1<-enrichWP(gene = c1$entrezid,organism = "Homo sapiens",pvalueCutoff = 1.5)
ego1<-setReadable(x = ego1,OrgDb = illuminaHumanv4.db,keyType = "ENTREZID")
ego1@result

colnames(ego1@result)<-paste(colnames(ego1@result),
                               "c1", sep = "_")
colnames(ego1@result)[1]<-"Term"
ego1@result


#######################################################################################################
ego2<-enrichWP(gene = c2$entrezid,organism = "Homo sapiens",pvalueCutoff = 1.5)
ego2<-setReadable(x = ego2,
                  OrgDb = illuminaHumanv4.db,
                  keyType = "ENTREZID")

colnames(ego2@result)<-paste(colnames(ego2@result),
                             "c2", sep = "_")
colnames(ego2@result)[1]<-"Term"
dim(ego2@result)


#########################################################################################################################
ego3<-enrichWP(gene = c3$entrezid,organism = "Homo sapiens",pvalueCutoff = 1.5)
ego3<-setReadable(x = ego3,OrgDb = illuminaHumanv4.db,keyType = "ENTREZID")
dim(ego3@result)

colnames(ego3@result)<-paste(colnames(ego3@result),
                             "c3", sep = "_")
colnames(ego3@result)[1]<-"Term"
ego3@result

##########################################################################################################################
ego4<-enrichWP(gene = c4$entrezid,organism = "Homo sapiens",pvalueCutoff = 1.5,
               pAdjustMethod = "BH")
ego4<-setReadable(x = ego4,OrgDb = illuminaHumanv4.db,keyType = "ENTREZID")
ego4@result
dim(ego4@result)
colnames(ego4@result)<-paste(colnames(ego4@result),
                             "c4", sep = "_")
colnames(ego4@result)[1]<-"Term"
ego4@result

############################################################################################################################
library(tidyverse)
Wiki_enrichr<-list(ego1@result, ego2@result, ego3@result, ego4@result) %>% 
  reduce(left_join, by = "Term")
head(Wiki_enrichr)

l<-list(Wiki_enrichr,PanglaoB_enrichr )
write.xlsx(l, rowNames=T,
           file = "WikiPathway_enrichr.xlsx")





#################################################################################################################################
setdiff(rownames(metadata),rownames(merge_meta_heatmap))

merge_meta_heatmap_one<-subset(merge_meta_heatmap, 
                               cohort=="one")
merge_meta_heatmap_one$fq1<-paste(rownames(merge_meta_heatmap_one),
                                  "R1.fastq.gz",sep = "")

merge_meta_heatmap_one$fq2<-paste(rownames(merge_meta_heatmap_one),
                                  "R2.fastq.gz",sep = "")
merge_meta_heatmap_one$title<-paste(rownames(merge_meta_heatmap_one),
                                    merge_meta_heatmap_one$condition,
                                    sep = "")

merge_meta_heatmap_one

#####################################################################
merge_meta_heatmap_two<-subset(merge_meta_heatmap, 
                               cohort=="two")
merge_meta_heatmap_two$fq1<-paste(rownames(merge_meta_heatmap_two),
                                  "_R1_001.fastq.gz",sep = "")

merge_meta_heatmap_two$fq2<-paste(rownames(merge_meta_heatmap_two),
                                  "_R2_001.fastq.gz",sep = "")
merge_meta_heatmap_two$condition
merge_meta_heatmap_two$title<-paste(rownames(merge_meta_heatmap_two),
                                    merge_meta_heatmap_one$condition,
                                    sep = "_")
merge_meta_heatmap_two
#####################################################################
bind_meta<-rbind(merge_meta_heatmap_one, merge_meta_heatmap_two)
bind_meta

write.xlsx(bind_meta,rowNames=T, 
           file = "~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/metadata_for_transferring2GEO_file.xlsx")
dim(cpm)
class(cpm)
cpm<-as.data.frame.matrix(cpm)
write.xlsx(cpm,rowNames=T, 
           file = "~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/Cerebral_malaria_EVRNAseq_CPM.xlsx")
table(bind_meta$condition)

library(ComplexHeatmap)
library(openxlsx)

library(circlize)
col_funp = colorRamp2(c(0, 1, 5), c("dodgerblue", "white", "red"))


df<-read.xlsx("wikipathway_supp_figure.xlsx",sheet = 1)
df

###################################################################
jpeg(filename = "pseudotime_analysisV3/Fig.4/wikipathway_supp_figure.jpeg",
     units = "in",width = 7,
     height = 2,res = 300)
htm<-Heatmap(-log10(df[,3:6][1:2,]),
             clustering_method_rows = "complete",
             clustering_method_columns = "complete",
             clustering_distance_columns = "canberra",
             clustering_distance_rows = "pearson",
             cluster_rows = T,
             cluster_columns = F,
             cluster_row_slices = F,
             col=col_funp,#row_title =df[1:19,]$function2 ,
             #row_split = df[1:2,]$function2,
             row_labels = df[1:2,]$Description,
             column_names_side = "top",
             column_names_rot = 0,
             heatmap_legend_param = list(title = "-log10(P-value)", 
                                         direction="horizontal"),
             #name = "-log10(P-value)",
             show_row_dend = FALSE,
             row_gap = unit(2, "mm"),
             column_gap = unit(0, "mm"),
             #row_labels = data$symbol_annot,
             cluster_column_slices = FALSE,
             column_names_gp = gpar(fontsize = 20),
             column_names_centered = TRUE,
             
             row_names_max_width = unit(18,"cm"),
             column_title_gp = gpar(fontsize = 20),
             #row_title_gp = gpar(fontsize = 25),
             row_title_gp = gpar(fontsize = 20),
             row_title_rot = 90)

htm_Pj<-draw(htm,#rowAnnotation(mark = anno),
             merge_legend = TRUE, 
             heatmap_legend_side = "top", 
             annotation_legend_side = "top",
             show_heatmap_legend = TRUE,
             legend_title_gp = gpar(fontsize = 20))
dev.off()



new_meta<-read.xlsx("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/supplementary_data_file_revised/supplementary_data_1.xlsx")

new_meta$haemorrhage_coded<-factor(as.factor(new_meta$haemorrhage_coded),
                                   levels = c("CC","CM-R⁻","CM-R⁺:haemorrhage-NO", 
                                              "CM-R⁺:haemorrhage-Yes"))



library(ggpubr)

jpeg("retinopathy_dichotomization.jpeg",width = 10,height = 5,units = "in",res = 300)
ggplot(new_meta, aes(haemorrhage_coded, pseudotime))+
  geom_boxplot()+
  geom_point(position = position_jitter(0.1),size=3)+
  xlab("")+theme_classic2()+
  theme(text = element_text(size = 15))
  
dev.off()


######################################################################################################
library(compareGroups)
library(openxlsx)

data_df<-read.xlsx("supplementary_data_file_revised/supplementary_data_1.xlsx")
data_df_nocc<-subset(data_df, condition!="CC")

data_df_nocc$Metabolic.acidosis

colnames(data_df_nocc)
dim(data_df_nocc)
str(data_df_nocc)
data_df_nocc$Parasitaemia<-as.numeric(data_df_nocc$Parasitaemia)
data_df_nocc$Respiratory.rate<-as.numeric(data_df_nocc$Respiratory.rate)
data_df_nocc$BCS.verbal<-as.numeric(data_df_nocc$BCS.verbal)
data_df_nocc$BCS.motor<-as.numeric(data_df_nocc$BCS.motor)
data_df_nocc$BCS.eyes<-as.numeric(data_df_nocc$BCS.eyes)
data_df_nocc$OutcomeDied<-as.factor(data_df_nocc$OutcomeDied)
data_df_nocc$Metabolic.acidosis<-as.factor(data_df_nocc$Metabolic.acidosis)

res <- compareGroups(condition ~ ., data = data_df_nocc[,c(7,10:24)])
res
table_res<-createTable(res,show.n = T)
table_res
export2xls(table_res,file = "supplementary_data_file_revised/Data_Table_S2.xls")


library(readstata13)
##########################################################
df<-read.dta13("prot_b_c_for_abdi_31072023.dta")
df

##################################################################
write.xlsx(df, "prot_b_c_for_abdi_31072023.xlsx")





















