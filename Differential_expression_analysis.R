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
library(tximport)
library(lumi)
library(sva)
library(BRETIGEA)
library(MultiBaC)
library(noise)
library(ViSEAGO)
library(readr)

########################################################################################################################
metadata<-read.xlsx("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/supplementary_data_file_revised/Data_Table_S1.xlsx", 
                    rowNames=TRUE,check.names = FALSE)
colnames(metadata)

load("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/hs_ret_all_kallisto_genelevel.Rdata")

hs_ret_all_kallisto_genelevel$counts
hs_ret_all_kallisto_genelevel$length


#start the process of normalizing the counts by gene length and lib size

#make the design
design<-model.matrix(~condition, metadata)
dim(design)
head(design)
colnames(design)[2:3]<-c("RP","RN")
##################################################################################

expr_pos <- 
 noisyr::noisyr(approach.for.similarity.calculation = "counts", 
                 expression.matrix = hs_ret_all_kallisto_genelevel$counts)


cts<-hs_ret_all_kallisto_genelevel$counts[rownames(expr_pos),][,rownames(metadata)]
dim(cts)

normMat<-hs_ret_all_kallisto_genelevel$length[rownames(cts),][,rownames(metadata)]
dim(normMat)

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
colnames(design)
dgeObj <- DGEList(cts,
                  group =metadata$condition,
                  samples = metadata)

dgeObj<-scaleOffset(dgeObj,normMat)
dgeObj$offset

###########################################################################################
dgeObj<-estimateDisp(dgeObj,design = design,trend.method = "loess",tagwise = TRUE)
contr<-makeContrasts(RP-RN, levels = design)
plotBCV(dgeObj)
fit<-glmFit(dgeObj,design = design,dispersion = dgeObj$common.dispersion)
lrt<-glmLRT(glmfit = fit,contrast = contr)
res<-topTags(lrt,n = Inf)$table
res

table(res$logFC>0, res$FDR<0.05)

####################################################################
library(illuminaHumanv4.db)
res$symbol<-mapIds(x = illuminaHumanv4.db,
                           keys = rownames(res),
                           column = "SYMBOL",
                           keytype = "ENSEMBL")


res$geneName<-mapIds(x = illuminaHumanv4.db,
                             keys = rownames(res),
                             column = "GENENAME",
                             keytype = "ENSEMBL")
head(res)

table(res$logFC>0, res$FDR<0.05)
# Classify genes into significantly up and down
res <- res %>% 
  mutate(status=factor(case_when(logFC>0 & FDR<0.05 ~ "519 genes",
                                 logFC<=0 & FDR<0.05 ~ "921 genes",
                                 TRUE ~ "not.signif"),
                       levels=c("519 genes", 
                                "not.signif", 
                                "921 genes")))




res_sig<-subset(res, FDR<0.05)
res_sig

#write.xlsx(list(res_sig, res), "Fig.1/Retinopathy_positive_vs_negative_edgeR.xlsx", rowNames=T)

####################################################################
#calculate cpm
cpm<-edgeR::cpm(y = dgeObj,normalized.lib.sizes = TRUE,
                    dispersion=dgeObj$common.dispersion,
                    log = F, prior.count = 0.001)
head(cpm)


library(circlize)
library(viridis)
metada_noDup$cohort
col_funp = colorRamp2(c(-1, 0, 2), c("dodgerblue", "black", "gold"))

colret = list(condition = c("CM-R⁺"="firebrick","CC"="gold",
                            "CM-R⁻"="steelblue"),
              acidosis = c("0"="bisque",
                           "1"="darkseagreen"),
              died = c("0"="lightblue1",
                       "1"="tan1"),
              cohort = c("one"="lightblue",
                         "two"="tan1"),
              `CM-R⁻ group`=c("3"="firebrick","1"="gold",
                              "2"="steelblue"),
              pseudotime = colorRamp2(c(0,4,8,12,16),viridis(n = 5,
                                                             option = "B",
                                                             begin = 0,end = 1)))

colret


metadata$condition<-factor(as.factor(metadata$condition),levels = c("CC","CM-R⁻","CM-R⁺"))

library(ComplexHeatmap)
metadata_nocc<-subset(metadata, condition!="CC")
ha1 = HeatmapAnnotation(#pseudotime=merge_meta_heatmap$Trajectory_slingshot,
                        #condition= factor(metadata_nocc$condition),
                        #died=merge_meta_heatmap$OutcomeDied,
                        #acidosis=merge_meta_heatmap$Metabolic.acidosis,
                        #cohort=merge_meta_heatmap$cohort,
                        `CM-R⁻ group`=metadata_RN$gaussian_cluster,
                        annotation_name_side = "left",
                        show_legend =TRUE,col = colret,
                        show_annotation_name = F,
                        na_col = "white",border = TRUE,
                        annotation_legend_param = list(legend_gp = gpar(fontsize = 17)),
                        simple_anno_size = unit(0.8, "cm"),
                        annotation_name_gp = gpar(fontsize = 17))
ha1
#graphics.off()
lab<-c("GNAT1","PAX5","VEGFA","ABCA4","GPR179","RGS11","CRB1",
       "CXCR2","CARD9","NFKB1","BSN","KCNA1","MYCBP2","OVOL2",
       "CD38", "CTLA4", "GZMK", "TLR7","CD209","RBP1", 
       "GLUL","GPX4","VEGFB","AICDA",
       "PDGFA", "PRG4","TIGIT","LY75")

res_sig<-res_sig[order(res_sig$geneName),]
mark_at = which(res_sig$symbol %in% lab)
mark_at

har = rowAnnotation(foo = anno_mark(at = which(res_sig$symbol %in% lab),
                                    labels = res_sig$symbol[res_sig$symbol%in%lab]))
har

####################################################################################
jpeg(filename = "pseudotime_analysisV3/Fig.1/RP_vs_RN_heatmap.jpeg",
     units = "in",width = 7,height = 7,res = 300)
htm<-Heatmap(t(scale(t(cpm[rownames(res_sig),][,rownames(metadata_nocc)]))),
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
             column_split = metadata_nocc$condition,
             col = col_funp,
             right_annotation = har,
             column_title = NULL,
             column_labels = c(1:85),
             row_title_rot = 90,
             column_names_rot = 0,
             column_names_side = "top",
             show_column_names = FALSE,
             row_gap = unit(0, "mm"),
             column_gap = unit(0, "mm"),
             row_split = res_sig$status,
             row_names_gp = gpar(fontsize = 10),
             row_names_max_width = unit(5,"cm"),
             column_title_gp = gpar(fontsize = 17),
             row_title_gp = gpar(fontsize = 17),
             use_raster = FALSE)
htm
htm_Pj<-draw(htm,
             merge_legend = TRUE, 
             heatmap_legend_side = "right",
             legend_title_gp = gpar(fontsize = 17),
             annotation_legend_side = "top",
             show_heatmap_legend = TRUE)
dev.off()


df<-read.xlsx("Enrichment_analysis.xlsx",sheet = 2)
head(df)
library(ggallin)
library(ggthemes)
library(ggforce)

jpeg(filename = "RP_vs_RN_enriched_cellTypesV2.jpeg",
     units = "in",width = 6.7,height = 7,res = 300)
ggplot(df, aes(x = reorder(Term, NES), y = NES)) +
  geom_bar(stat = "identity",alpha=0.7,
           show.legend = F,
           aes(fill = NES),  # Background color
           color = "gray30") + # Border color
  xlab("Group") +
  ylab("Value") +
  scale_fill_gradient2(low = "#6495ED",
                       mid = "aliceblue",
                       high = "firebrick")+
  coord_flip()+
  facet_col(~cellType,scales = "free_y")+
  theme_stata()+
  theme(axis.text.y = element_text(angle = 0),
        text = element_text(size = 18))+
  xlab("")+ylab("Net Enrichment Score (NES)")+
  scale_y_continuous(trans = ssqrt_trans)
dev.off()


########################################################################################
#
metadata_RN<-subset(metadata, condition=="CM-R⁻")
metadata_RN$pseudotime

ggplot(metadata_RN, aes(pseudotime))+
  geom_histogram(bins = 30)

#####################################################################################
library(mixtools)
library(plotGMM)
library(ggthemes)
mixmdl <- normalmixEM((metadata_RN$pseudotime), k = 3,lambda = 2)
p1<-data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..),bins = 30, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "gold", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "steelblue", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[3], mixmdl$sigma[3], lam = mixmdl$lambda[3]),
                colour = "firebrick", lwd = 1.5)+
  ylab("Density")+
  theme_stata()+
  theme(text = element_text(size = 17),
        axis.text.y = element_text(angle = 0))+
  xlab("Pseudotime")

##############################################################################
jpeg(filename = "pseudotime_analysisV3/Fig.1/RN_stratification.jpeg",
     units = "in",width = 5,height = 4,res = 600)
p1
dev.off()


post.df <- as.data.frame(cbind(x = mixmdl$x, mixmdl$posterior))
head(post.df, 10) 
################################################################################################
RNgroup <- post.df %>% 
  mutate(group=factor(case_when(comp.1>0.3  ~ "1",
                                 comp.2>0.3 ~ "2",
                                 TRUE ~ "3"),
                       levels=c("1", 
                                "2", 
                                "3")))
RNgroup
#################################################################################################
metadata_RN$gaussian_cluster<-RNgroup$group
metadata_RN$`Age.(months)`
df<-metadata_RN
df$Parasitaemia<-((as.numeric(df$Parasitaemia)))
df$hrp2_plasma<-as.numeric(df$hrp2_plasma)
df$`Age.(months)`<-((as.numeric(df$`Age.(months)`)))
df$Respiratory.rate<-sqrt(as.numeric(df$Respiratory.rate))
df$Temperature<-sqrt(as.numeric(df$Temperature))
df$MUAC<-sqrt(as.numeric(df$MUAC))
df$Admission.days<-sqrt(as.numeric(df$Admission.days))

df$BCS.verbal<-(as.numeric(df$BCS.verbal))
df$BCS.motor<-(as.numeric(df$BCS.motor))
df$BCS.eyes<-(as.numeric(df$BCS.eyes))


df2vs1<-subset(df, df$gaussian_cluster!="1")
colnames(df2vs1)

library(brglm2)
fit_un <- lapply(df2vs1[,c(9:23)], function(x) coef(summary(glm((gaussian_cluster) ~ x, 
                                                                                      data =df2vs1, 
                                                                                      family = binomial(),
                                                                type = "MPL_Jeffreys",method = "brglmFit"
                                                                ))))

fit_un


class(fit_un)
fit_un_df<-do.call(rbind.data.frame, fit_un) #convert the list to a data frame
fit_un_df

df_un<-fit_un_df[!grepl("Intercept", rownames(fit_un_df)),] #remove rows with the intercept
rownames(df_un)<-gsub("\\.x.*", "", rownames(df_un))
df_un

gp1_vs_gp2<-list(c("early","mid"), c("mid","late"),c("early","late"))

library(ggallin)
p1<-ggplot(df,aes(as.factor(gaussian_cluster),Parasitaemia))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(0.1),size=2)+
  stat_compare_means(comparisons = gp1_vs_gp2,size=6,
                     bracket.size = 0.5,
                     #vjust = 0.5,
                     tip.length = 0.03)+
  scale_y_continuous(trans = ssqrt_trans,
                     breaks = c(50000, 250000,500000, 1000000))+
  theme_stata()+
  theme(axis.text.y = element_text(angle = 0),
        text = element_text(size = 18))+
  ylab("Parasites/ul")+
  ggtitle("Parasitaemia")+
  xlab("")
p1

p2<-ggplot(df,aes(as.factor(gaussian_cluster),Haemoglobin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(0.1),size=2)+
  stat_compare_means(comparisons = gp1_vs_gp2,size=6,
                     bracket.size = 0.5,
                     #vjust = 0.5,
                     tip.length = 0.03)+
  #scale_y_continuous(trans = ssqrt_trans,
                     #breaks = c(50000, 250000,500000, 1000000))+
  theme_stata()+
  theme(axis.text.y = element_text(angle = 0),
        text = element_text(size = 20))+
  xlab("")+
  ylab("g/dL")+
  ggtitle("Haemoglobin")+
  xlab("CM-R⁻ subgroups")
p2

p3<-ggplot(df,aes(as.factor(gaussian_cluster),`Age.(months)`))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(0.1))+
  stat_compare_means(comparisons = gp1_vs_gp2,size=6,
                     bracket.size = 0.5,
                     #vjust = 0.5,
                     tip.length = 0.03)+
  #scale_y_continuous(trans = ssqrt_trans,
  #breaks = c(50000, 250000,500000, 1000000))+
  theme_stata()+
  theme(axis.text.y = element_text(angle = 0),
        text = element_text(size = 20))+
  ylab("Months")+
  ggtitle("Age")+
  xlab("CM-R⁻ subgroups")
p3


p4<-ggplot(df,aes(as.factor(gaussian_cluster),hrp2_plasma))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(0.1))+
  stat_compare_means(comparisons = gp1_vs_gp2,size=6,
                     bracket.size = 0.5,
                     #vjust = 0.5,
                     tip.length = 0.03)+
  #scale_y_continuous(trans = ssqrt_trans,
  #breaks = c(50000, 250000,500000, 1000000))+
  theme_stata()+
  theme(axis.text.y = element_text(angle = 0),
        text = element_text(size = 18))+
  ylab("ng/ul")+
  ggtitle("PfHRP2")+
  xlab("CM-R⁻ groups")
p4
p3
library(gridExtra)
plot<-list(p2,p3)

bottom<-textGrob("CM-R⁻ groups", gp = gpar(fontsize=18))

grid.arrange(grobs=plot,nrow=2, bottom=bottom)

jpeg(filename = "RN_stratification_boxplotsAgeHBV2.jpeg",
     units = "in",width = 8,height = 7,res = 600)
p3+p2
dev.off()

str(metadata_RN)



#Creating a DGEList object for use in edgeR.
design<-model.matrix(~gaussian_cluster, metadata_RN)
colnames(design)<-make.names(colnames(design))
colnames(design)
dgeObj <- DGEList(cts[,rownames(metadata_RN)],
                  group =metadata_RN$condition,
                  samples = metadata_RN)

dgeObj<-scaleOffset(dgeObj,normMat[,rownames(metadata_RN)])
dgeObj$offset

###########################################################################################
dgeObj<-estimateDisp(dgeObj,design = design,trend.method = "loess",tagwise = TRUE)
#contr<-makeContrasts(RP-RN, levels = design)
plotBCV(dgeObj)
fit<-glmFit(dgeObj,design = design,dispersion = dgeObj$common.dispersion)
lrt<-glmLRT(glmfit = fit,coef = 2:3)
res<-topTags(lrt,n = Inf)$table
res

table(res$FDR<0.05)

####################################################################
library(illuminaHumanv4.db)
res$symbol<-mapIds(x = illuminaHumanv4.db,
                   keys = rownames(res),
                   column = "SYMBOL",
                   keytype = "ENSEMBL")


res$geneName<-mapIds(x = illuminaHumanv4.db,
                     keys = rownames(res),
                     column = "GENENAME",
                     keytype = "ENSEMBL")
head(res)

table(res$FDR<0.05)


res_sig<-subset(res, FDR<0.05)
res_sig

cpm_rn_mean<-cpmByGroup(dgeObj,group = metadata_RN$gaussian_cluster)

#########################################################################
#find optimal cluster of genes using kmeans
# Elbow method
fviz_nbclust(t(scale(t(cpm_rn_mean[rownames(res_sig),]))), 
             kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(t(scale(t(cpm_rn_mean[rownames(res_sig),]))), 
             kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

######################################################
kmeans_cluster=kmeans(t(scale(t(cpm_rn_mean[rownames(res_sig),]))),
                      centers = 3)
kmeans_cluster
res_sig$kmeans_cluster=kmeans_cluster$cluster

Heatmap(t(scale(t(cpm_rn_mean[rownames(res_sig),]))), 
        cluster_rows = F, col = col_funp,
        cluster_columns = F, show_row_names = F, 
        show_column_names = T,use_raster = F,
        row_split = res_sig$kmeans_cluster,
        #column_split = metadata_RN$gaussian_cluster
        )


metadata_RN$gaussian_cluster<-gsub("1","early",metadata_RN$gaussian_cluster)
metadata_RN$gaussian_cluster<-gsub("2","mid",metadata_RN$gaussian_cluster)
metadata_RN$gaussian_cluster<-gsub("3","late",metadata_RN$gaussian_cluster)


colret = list(condition = c("CM-R⁺"="firebrick","CC"="gold",
                            "CM-R⁻"="steelblue"),
              acidosis = c("0"="bisque",
                           "1"="darkseagreen"),
              died = c("0"="lightblue1",
                       "1"="tan1"),
              cohort = c("one"="lightblue",
                         "two"="tan1"),
              `CM-R⁻ subgroups`=c("late"="firebrick","early"="gold",
                              "mid"="steelblue"),
              pseudotime = colorRamp2(c(0,4,8,12,16),viridis(n = 5,
                                                             option = "B",
                                                             begin = 0,end = 1)))

colret

metadata_RN$condition<-factor(as.factor(metadata_RN$condition),levels = c("CC","CM-R⁻","CM-R⁺"))

library(ComplexHeatmap)
metadata_nocc<-subset(metadata, condition!="CC")
ha1 = HeatmapAnnotation(#pseudotime=metadata_RN$pseudotime,
  #condition= factor(metadata_nocc$condition),
  #died=merge_meta_heatmap$OutcomeDied,
  #acidosis=merge_meta_heatmap$Metabolic.acidosis,
  #cohort=merge_meta_heatmap$cohort,
  `CM-R⁻ subgroups`=metadata_RN$gaussian_cluster,
  annotation_name_side = "left",
  show_legend =TRUE,col = colret,
  show_annotation_name = F,
  na_col = "white",border = TRUE,
  annotation_legend_param = list(legend_gp = gpar(fontsize = 17)),
  simple_anno_size = unit(0.5, "cm"),
  annotation_name_gp = gpar(fontsize = 17))
ha1

lab<-c("CTLA4","GZMK","VEGFB","PDGFA",
       "CHAT", "CASD1","STY6", #cholinergic neuronics
       "SPTBN4","NRCAM","FGF14",#Axon markers
       "BSN","PCLO","NRXN1", #presynaptic terminals
       "PRG4","SOX10","MYT1",#OPC
       "S100A9",
       "ZBTB7B",#thymocytes
       "PAX5",#B cell precursors
       "TBX21",# Th1 cells cells
       "THPO",#megakaryocytes
       "ELANE","CTSG","S100A8","NFKB1","AZU1","DEFA1",#neutrophils"
       "PF4","GMFG","GPX4","CD74","CTSA","DEFA3","DEFA1B")
res_sig<-res_sig[order(rownames(res_sig)),]
mark_at = which(res_sig$symbol %in% lab)
mark_at

har = rowAnnotation(foo = anno_mark(at = which(res_sig$symbol %in% lab),
                                    labels = res_sig$symbol[res_sig$symbol%in%lab]))
har

table(res_sig$kmeans_cluster)
#cluster 1 =913 genes, cluster2=635, cluster3=1838 genes
metadata_RN<-metadata_RN[order(metadata_RN$pseudotime),]
metadata_RN$gaussian_cluster<-factor(as.factor(metadata_RN$gaussian_cluster),
                                     levels = c("early", "mid", "late"))
#########################################################################################################

jpeg(filename = "RN_subgraoups_diffexpressedGenesV2.jpeg",
     units = "in",width = 6.5,height = 8.5,res = 400)
htm<-Heatmap(t(scale(t(cpm[rownames(res_sig),][,rownames(metadata_RN)]))), 
        cluster_rows = F, col = col_funp,
        cluster_columns = F, show_row_names = F, 
        show_column_names = F,use_raster = F,
        row_split = factor(as.factor(res_sig$kmeans_cluster),
                           levels = c("3","1","2")),
        column_split = metadata_RN$gaussian_cluster,
        row_gap = unit(0.3, "mm"),
        column_gap = unit(0.3, "mm"),
        row_title = NULL,column_title = NULL,
        row_names_gp = gpar(fontsize = 10),
        row_names_max_width = unit(5,"cm"),
        column_title_gp = gpar(fontsize = 17),
        row_title_gp = gpar(fontsize = 17),
        heatmap_legend_param = list(title = "mean z score", 
                                    direction="vertical"),
        top_annotation = ha1,
        right_annotation = har)
htm
htm_Pj<-draw(htm,
             merge_legend = TRUE, 
             heatmap_legend_side = "right",
             legend_title_gp = gpar(fontsize = 17),
             annotation_legend_side = "top",
             show_heatmap_legend = TRUE)
dev.off()


library(rbioapi)
c1<-subset(res_sig, res_sig$kmeans_cluster=="1")
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
c2<-subset(res_sig, res_sig$kmeans_cluster=="2")
c2
c2_allenBrain <- rba_enrichr(gene_list = c2$symbol,
                             organism = "human",
                             gene_set_library = "PanglaoDB_Augmented_2021",
                             regex_library_name = FALSE)
head(c2_allenBrain)
colnames(c2_allenBrain)<-paste(colnames(c2_allenBrain),
                               "c2", sep = "_")
colnames(c2_allenBrain)[1]<-"Term"
c2_allenBrain


##########################################################
c3<-subset(res_sig, res_sig$kmeans_cluster=="3")
c3
c3_allenBrain <- rba_enrichr(gene_list = c3$symbol,organism = "human",
                             gene_set_library = "PanglaoDB_Augmented_2021",
                             regex_library_name = FALSE)
head(c3_allenBrain)

colnames(c3_allenBrain)<-paste(colnames(c3_allenBrain),
                               "c3", sep = "_")
colnames(c3_allenBrain)[1]<-"Term"
c3_allenBrain

#KEGG_2021_Human
c1<-subset(res_sig, res_sig$kmeans_cluster=="1")
c1
c1_allenBrain <- rba_enrichr(gene_list = c1$symbol,organism = "human",
                             gene_set_library = "KEGG_2021_Human",
                             regex_library_name = FALSE)
head(c1_allenBrain)
c1_allenBrain$module="group1"
colnames(c1_allenBrain)<-paste(colnames(c1_allenBrain),
                               "group1", sep = "_")
colnames(c1_allenBrain)[1]<-"Term"
c1_allenBrain

#################################################################
c2<-subset(res_sig, res_sig$kmeans_cluster=="2")
c2
c2_allenBrain <- rba_enrichr(gene_list = c2$symbol,
                             organism = "human",
                             gene_set_library = "KEGG_2021_Human",
                             regex_library_name = FALSE)
head(c2_allenBrain)
colnames(c2_allenBrain)<-paste(colnames(c2_allenBrain),
                               "c2", sep = "_")
colnames(c2_allenBrain)[1]<-"Term"
c2_allenBrain


##########################################################
c3<-subset(res_sig, res_sig$kmeans_cluster=="3")
c3
c3_allenBrain <- rba_enrichr(gene_list = c3$symbol,organism = "human",
                             gene_set_library = "KEGG_2021_Human",
                             regex_library_name = FALSE)
head(c3_allenBrain)

colnames(c3_allenBrain)<-paste(colnames(c3_allenBrain),
                               "c3", sep = "_")
colnames(c3_allenBrain)[1]<-"Term"
c3_allenBrain


write.xlsx(res_sig,rowNames=T,"pseudotime_analysisV3/Fig.1/RN_subgroups_edgeR_results.xlsx")


############################################################

small_mat<-read.xlsx("EdgeR_literatureLab_Series Analysis Data..xlsx",sheet = 2)
colnames(small_mat)
str(small_mat)

jpeg(filename = "enrichment_analysis_Ret_negv2.jpeg",
     units = "in",width = 4.2,height = 5.5,res = 400)
heatMap<-Heatmap(as.matrix(small_mat[,7:9]), 
        column_names_rot = 0,cluster_rows = F,
        column_title = "CM-R⁻ subgroups",
        column_title_side = "bottom",
        column_names_side = "bottom",
        column_names_centered = T,
        row_labels = small_mat$Term,cluster_columns = F,
        #name = "enrichment score", 
        heatmap_legend_param = list(title = "score", 
                                    direction="horizontal"),
        #col = col_fun,
        layer_fun = function(j, i, x, y, width, height, fill) {
          # since grid.text can also be vectorized
          grid.text(sprintf("%.1f", pindex(as.matrix(small_mat[,7:9]), i, j)), x, y, 
                    gp = gpar(fontsize = 12))
        })
htm_smallmat
htm_smallmat<-draw(heatMap,
             merge_legend = TRUE, 
             heatmap_legend_side = "top",
             legend_title_gp = gpar(fontsize = 12),
             annotation_legend_side = "left",
             show_heatmap_legend = TRUE)
dev.off()


metadata1<-read.xlsx("~/Documents/PHD_THESIS_FIGURES/Retinopathy_EV_paper_final/metadata/Metadata_all_Retinopathy_samples.xlsx", 
                    rowNames=TRUE,check.names = FALSE)
colnames(metadata1)


metadata1$study<-gsub("test", "Hiseq", metadata1$study)
metadata1$study<-gsub("validation", "Nextseq", metadata1$study)

commonsample<-intersect(rownames(metadata1),colnames(cpm))

pca1<-pca(t(log2(cpm[,commonsample]+0.5)),scale = T)

jpeg(filename = "pseudotime_analysisV3/Fig.1/Batch_inspectionPCA.jpeg",
     units = "in",width = 5,height = 4.5,res = 400)
plotIndiv(pca1, group = metadata1[commonsample,]$study, 
          ind.names = F,legend = T,title = "PCA colored by Batch")
dev.off()

table(res_sig$kmeans_cluster)
635+1838+913
table(metadata_RN$gaussian_cluster)

#######################################################################################################
library(compareGroups)






