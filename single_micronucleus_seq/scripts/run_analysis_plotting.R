library(pheatmap)
library(ComplexHeatmap)
library(ggplot2)
library(ggsignif)
library(dplyr)
library(circlize) 
##script to plot figures

#load data
setwd("/Users/robinxu/Documents/Projects/Brueckner_et_al/sequencing_final_publication/ecMN_final/single_micronucleus_seq")

window_calls_counts_normalised_amplicon_anno_all_run <- readRDS("./data/processed/data_processed_r1_r2.rds")
sample_pass_vec <- readRDS("./data/processed/samples_pass.rds")
meta_qc_ext_combined <- readRDS("./data/processed/meta_r1_r2.rds")
index_MNsequences_on_ecDNA_r1 <- readRDS("./data/processed/index_MNsequences_on_ecDNA_r1.rds")
index_MNsequences_on_ecDNA_r2 <- readRDS("./data/processed/index_MNsequences_on_ecDNA_r2.rds")

MN_pass_vec_all_runs <- sample_pass_vec$MN
PN_pass_vec_all_runs <- sample_pass_vec$PN
pMN_pass_vec <- sample_pass_vec$pMN



#######
###Fig.3c/Ext. Fig. 9a
#######

#make seperate heatmaps for single micronuclei. pooled mironuclei and primary nuclei
#apply row clustering (hclust) and take the order for the combined heatmap

##single micronuclei

MN_mat_heatmap_df <- window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename %in% MN_pass_vec_all_runs,]

MN_mat_heatmap <- unstack(MN_mat_heatmap_df,log2CPM~samplename)
MN_mat_heatmap_sub <- MN_mat_heatmap[c(seq(10000,10700),seq(80000,80600)),]


col_anno <- data.frame(MN_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$chr)
row.names(col_anno) <- c(seq(10000,10700),seq(80000,80600))
colnames(col_anno) <- "Region"
col_anno$Amplicon <- MN_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$amplicon_anno
#count(col_anno$Region) #determine number of chr bins
#colour annotation
ann_colors = list(
  Amplicon = c(z = "white", MYCN = "#e41a1c", CDK4 = "#377eb8", ODC1 = "#4daf4a", MDM2 = "#984ea3", SMC6 = "#ff7f00"),
  Region = c(chr2 = "grey", chr12 = "black")
)

#pdf("./output/figures/figure3c_ext9/heatmapt_MN.pdf", width =10 , height = 8) #default height 3

heatmap_MN <- ComplexHeatmap::pheatmap(data.frame(t(MN_mat_heatmap_sub)), color=colorRampPalette(c("white","red","darkred"))(200), annotation_colors = ann_colors,cluster_rows = TRUE,gaps_col = 701 ,cluster_cols = FALSE, show_rownames = T, show_colnames = FALSE, annotation_col = col_anno[,c("Amplicon","Region")])
heatmap_MN
#dev.off()

ht_MN <- draw(heatmap_MN)

MN_heatmap_row_order <- rownames(data.frame(t(MN_mat_heatmap_sub))[row_order(ht_MN), ])

###single primary nuclei
PN_mat_heatmap_df <- window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename %in% PN_pass_vec_all_runs,]

PN_mat_heatmap <- unstack(PN_mat_heatmap_df,log2CPM~samplename)
PN_mat_heatmap_sub <- PN_mat_heatmap[c(seq(10000,10700),seq(80000,80600)),]


col_anno <- data.frame(PN_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$chr)
row.names(col_anno) <- c(seq(10000,10700),seq(80000,80600))
colnames(col_anno) <- "Region"
col_anno$Amplicon <- PN_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$amplicon_anno
#count(col_anno$Region) #determine number of chr bins
#colour annotation
ann_colors = list(
  Amplicon = c(z = "white", MYCN = "#e41a1c", CDK4 = "#377eb8", ODC1 = "#4daf4a", MDM2 = "#984ea3", SMC6 = "#ff7f00"),
  Region = c(chr2 = "grey", chr12 = "black")
)

#pdf("./output/figures/figure3c_ext9/heatmapt_PN.pdf", width =10 , height = 3)

heatmap_PN <- ComplexHeatmap::pheatmap(as.matrix(t(PN_mat_heatmap_sub)), color=colorRampPalette(c("white","red","darkred"))(200), annotation_colors = ann_colors,cluster_rows = TRUE,gaps_col = 701 ,cluster_cols = FALSE, show_rownames = T, show_colnames = FALSE, annotation_col = col_anno[,c("Amplicon","Region")])
heatmap_PN
#dev.off()
ht_PN <- draw(heatmap_PN)
PN_heatmap_row_order <- rownames(data.frame(t(PN_mat_heatmap_sub))[row_order(ht_PN),])



###pooled micronuclei
pMN_mat_heatmap_df <- window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename %in% pMN_pass_vec,]

pMN_mat_heatmap <- unstack(pMN_mat_heatmap_df,log2CPM~samplename)
pMN_mat_heatmap_sub <- pMN_mat_heatmap[c(seq(10000,10700),seq(80000,80600)),]


col_anno <- data.frame(pMN_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$chr)
row.names(col_anno) <- c(seq(10000,10700),seq(80000,80600))
colnames(col_anno) <- "Region"
col_anno$Amplicon <- pMN_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$amplicon_anno
#count(col_anno$Region) #determine number of chr bins
#colour annotation
ann_colors = list(
  Amplicon = c(z = "white", MYCN = "#e41a1c", CDK4 = "#377eb8", ODC1 = "#4daf4a", MDM2 = "#984ea3", SMC6 = "#ff7f00"),
  Region = c(chr2 = "grey", chr12 = "black")
)

#pdf("./output/figures/figure3c_ext9/heatmapt_pMN.pdf", width =10 , height = 3)

heatmap_pMN <- ComplexHeatmap::pheatmap(as.matrix(t(pMN_mat_heatmap_sub)), color=colorRampPalette(c("white","red","darkred"))(200), annotation_colors = ann_colors,cluster_rows = TRUE,gaps_col = 701 ,cluster_cols = FALSE, show_rownames = T, show_colnames = FALSE, annotation_col = col_anno[,c("Amplicon","Region")])
heatmap_pMN
#dev.off()
ht_pMN <- draw(heatmap_pMN)
pMN_heatmap_row_order <- rownames(data.frame(t(pMN_mat_heatmap_sub))[row_order(ht_pMN),])

##Master Heatmap - combine all 3 heatmaps into 1

#make amplicon specific annotations
#MYCN

amplicon_vec <- c("MYCN","MYCN","MYCN","MYCN","CDK4","CDK4","SMC6","ODC1","ODC1","MDM2")
amplicon_start <- c(10052,10095,10506,10592,80070,82170,10638,10366,80041,80509)
amplicon_end <-c(10056,10097,10508,10613,80076,82174,10677,10389,80045,80549)

amplicon_anno_df <- data.frame(amplicon = amplicon_vec, start = amplicon_start, end = amplicon_end )


#extend 8 bins
amplicon_anno_df$start <- amplicon_anno_df$start - 4
amplicon_anno_df$end <- amplicon_anno_df$end + 4
#make annotation

seq_vec <- c()
region_vec <- c()
sub_temp <- Master_mat_heatmap_df[Master_mat_heatmap_df$samplename=="A1",]
for (i in 1:nrow(amplicon_anno_df)){
  seq_vec <- c(seq_vec, seq(amplicon_anno_df[i,"start"],amplicon_anno_df[i,"end"]))
  region_vec <- c(region_vec,rep(paste(format(round(sub_temp[amplicon_anno_df[i,"start"],"start"] / 1e6, 1), trim = TRUE),"Mb",format(round(sub_temp[amplicon_anno_df[i,"end"],"end"] / 1e6, 1), trim = TRUE),"Mb"),length(seq(amplicon_anno_df[i,"start"],amplicon_anno_df[i,"end"]))))
}

amplicon_extend_anno_df <- data.frame(index = seq_vec, region_vec = region_vec)
rownames(amplicon_extend_anno_df) <- amplicon_extend_anno_df$index

##heatmap_reduced_ final

##Master Heatmap

heatmap_order <- c(pMN_heatmap_row_order,MN_heatmap_row_order,PN_heatmap_row_order)
breaks_rows <- c(length(pMN_heatmap_row_order),length(pMN_heatmap_row_order)+length(MN_heatmap_row_order))

Master_mat_heatmap_df <- window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename %in% heatmap_order,]

Master_mat_heatmap <- unstack(Master_mat_heatmap_df,log2CPM~samplename)
Master_mat_heatmap_sub <- Master_mat_heatmap[sort(amplicon_extend_anno_df$index),]


#row anno
row_anno <- meta_qc_ext_combined[meta_qc_ext_combined$sample %in% heatmap_order, ]
rownames(row_anno) <- row_anno$sample
row_anno$placeholder <- "a"
row_anno_sub <- row_anno[heatmap_order,c("class","placeholder")]

col_anno <- data.frame(Master_mat_heatmap_df[sort(amplicon_extend_anno_df$index),]$chr)
row.names(col_anno) <- c(sort(amplicon_extend_anno_df$index))
colnames(col_anno) <- "Region"
col_anno$Amplicon <- Master_mat_heatmap_df[sort(amplicon_extend_anno_df$index),]$amplicon_anno
col_anno$Region_coord <- amplicon_extend_anno_df[paste(sort(amplicon_extend_anno_df$index)),]$region_vec
#count(col_anno$Region) #determine number of chr bins
#colour annotation
ann_colors = list(
  Amplicon = c(z = "white", MYCN = "#e41a1c", CDK4 = "#377eb8", ODC1 = "#4daf4a", MDM2 = "#984ea3", SMC6 = "#ff7f00"),
  Region = c(chr2 = "grey", chr12 = "black")
)

dir.create("./output/figures/figure3c_ext9")
pdf("./output/figures/figure3c_ext9/heatmap_main_cleaned_new_v1_color_BYR_final.pdf", width =10 , height = 3)

heatmap_main <- pheatmap(as.matrix(t(Master_mat_heatmap_sub[,heatmap_order])), color=colorRampPalette(c("#2c7bb6","#ffffbf","#d7191c"))(200), annotation_colors = ann_colors,cluster_rows = F,gaps_row = breaks_rows,  gaps_col = c(21-8,40-(2*8),80-(3*8),99-(4*8),137-(5*8),193-(6*8),214-(7*8),237-(8*8),294-(9*8)) ,cluster_cols = FALSE, show_rownames = T, show_colnames = FALSE, annotation_row = row_anno_sub,annotation_col = col_anno[,c("Amplicon","Region","Region_coord")])
heatmap_main
dev.off()

#######
###Fig.3d/Ext. Fig
#######

#take the intersection of both ecDNA regions identified
index_MNsequences_on_ecDNA <- intersect(index_MNsequences_on_ecDNA_r1,index_MNsequences_on_ecDNA_r2)

#PN
PN_logCPM_called_region_mat <- unstack(window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename %in% PN_pass_vec_all_runs ,],CPM~samplename) #CPM matrix
PN_logCPM_called_region_mat_ecDNA <- t(PN_logCPM_called_region_mat)[,as.numeric(index_MNsequences_on_ecDNA)]
PN_logCPM_called_region_mat_ecDNA_mean <- rowMeans(PN_logCPM_called_region_mat_ecDNA)

PN_logCPM_lin_mat <- unstack(window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename %in% PN_pass_vec_all_runs ,],CPM~samplename) #CPM matrix
PN_logCPM_lin_mat <- t(PN_logCPM_lin_mat)[,-as.numeric(index_MNsequences_on_ecDNA)]
PN_logCPM_lin_mat_mean_wins <- rowMeans(t(apply(PN_logCPM_lin_mat, 1, function(x) Winsorize(x, val = quantile(x, probs = c(0.04, 0.96))))))



#MN
MN_logCPM_called_region_mat <- unstack(window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename %in% MN_pass_vec_all_runs ,],CPM~samplename) #CPM matrix
MN_logCPM_called_region_mat_ecDNA <- t(MN_logCPM_called_region_mat)[,as.numeric(index_MNsequences_on_ecDNA)]
MN_logCPM_called_region_mat_ecDNA_mean <- rowMeans(MN_logCPM_called_region_mat_ecDNA)

MN_logCPM_lin_mat <- unstack(window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename %in% MN_pass_vec_all_runs ,],CPM~samplename) #CPM matrix
MN_logCPM_lin_mat<- t(MN_logCPM_lin_mat)[,-as.numeric(index_MNsequences_on_ecDNA)]
MN_logCPM_lin_mat_mean_wins <- rowMeans(t(apply(MN_logCPM_lin_mat, 1, function(x) Winsorize(x, val = quantile(x, probs = c(0.04, 0.96))))))


#pooled MN
pMN_logCPM_called_region_mat <- unstack(window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename %in% pMN_pass_vec ,],CPM~samplename) #CPM matrix
pMN_logCPM_called_region_mat_ecDNA <- t(pMN_logCPM_called_region_mat)[,as.numeric(index_MNsequences_on_ecDNA)]
pMN_logCPM_called_region_mat_ecDNA_mean <- rowMeans(pMN_logCPM_called_region_mat_ecDNA)

pMN_logCPM_lin_mat <- unstack(window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename %in% pMN_pass_vec ,],CPM~samplename) #CPM matrix
pMN_logCPM_lin_mat <- t(pMN_logCPM_lin_mat)[,-as.numeric(index_MNsequences_on_ecDNA)]
pMN_logCPM_lin_mat_mean_wins <- rowMeans(t(apply(pMN_logCPM_lin_mat, 1, function(x) Winsorize(x, val = quantile(x, probs = c(0.04, 0.96))))))



#compute (circCPM/linCPM)
MN_CPM_ratio <- log2(as.vector(MN_logCPM_called_region_mat_ecDNA_mean)/MN_logCPM_lin_mat_mean_wins)
PN_CPM_ratio <- log2(as.vector(PN_logCPM_called_region_mat_ecDNA_mean)/PN_logCPM_lin_mat_mean_wins)
pMN_CPM_ratio <- log2(as.vector(pMN_logCPM_called_region_mat_ecDNA_mean)/pMN_logCPM_lin_mat_mean_wins)


shapiro.test(MN_CPM_ratio)
shapiro.test(PN_CPM_ratio)
shapiro.test(pMN_CPM_ratio)


##combine MN and PN for fig
MN_CPM_ratio_df <- data.frame(MN_CPM_ratio)
colnames(MN_CPM_ratio_df) <- "CPM_ratio"
MN_CPM_ratio_df$sample <- names(MN_CPM_ratio)
MN_CPM_ratio_df$grp <- "single micronucleus"


PN_CPM_ratio_df <- data.frame(PN_CPM_ratio)
colnames(PN_CPM_ratio_df) <- "CPM_ratio"
PN_CPM_ratio_df$sample <- names(PN_CPM_ratio)
PN_CPM_ratio_df$grp <- "Primary Nucleus"

pMN_CPM_ratio_df <- data.frame(pMN_CPM_ratio)
colnames(pMN_CPM_ratio_df) <- "CPM_ratio"
pMN_CPM_ratio_df$sample <- names(pMN_CPM_ratio)
pMN_CPM_ratio_df$grp <- "pooled micronuclei"

mean(MN_CPM_ratio_df$CPM_ratio)
mean(PN_CPM_ratio_df$CPM_ratio)
mean(pMN_CPM_ratio_df$CPM_ratio)

CPM_ratio_df<- rbind(pMN_CPM_ratio_df,MN_CPM_ratio_df,PN_CPM_ratio_df)


#testing
anno <- t.test(MN_CPM_ratio_df$CPM_ratio, PN_CPM_ratio_df$CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
anno1 <- t.test(MN_CPM_ratio_df$CPM_ratio, pMN_CPM_ratio_df$CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
anno2 <- t.test(pMN_CPM_ratio_df$CPM_ratio, PN_CPM_ratio_df$CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
test_vec <- c(anno,anno1,anno)

#plotting
CPM_ratio_df <- CPM_ratio_df %>%
  mutate(grp = factor(grp, levels = c("pooled micronuclei", "single micronucleus", "Primary Nucleus")))



ggplot(CPM_ratio_df, aes(x = grp, y = CPM_ratio, colour = grp, group = grp)) +
  geom_violin(position = position_dodge(width = 0.7), trim = FALSE) +
  geom_point(
    aes(colour = grp),
    position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.7),
    size = 1
  ) + 
  #  scale_color_okabeito(order=1:9)+
  theme_classic() +
  stat_summary(
    aes(group = grp),  
    fun = "mean",
    geom = "point",
    col = "black",
    size = 5,
    shape = 95,
    position = position_dodge(width = 0.7)  
  ) +coord_cartesian(ylim = c(0, NA)) +
  theme(legend.position="none")+
  geom_signif(
    comparisons = list(
      c("Primary Nucleus", "single micronucleus"),
      c("single micronucleus", "pooled micronuclei"),
      c("pooled micronuclei", "Primary Nucleus")
    ),
    annotation = formatC(test_vec, digits = 1),
    map_signif_level = TRUE,
    step_increase = 0.05
  ) #+
#theme(axis.title.x=element_blank(),
#      axis.text.x=element_blank(),
#      axis.title.y=element_blank(),
#     axis.text.y=element_blank())

dir.create("./output/figures/figure3d")
ggsave("./output/figures/figure3d/circ_lin_CPM_ratio_box_winsorised_Content_with_labels2.eps", width = 1.4,
       height = 1.4, device= "eps", dpi=450)


#######
###Fig.3f
#######

##perc amplicon occupancy

###PN
sample_vec <- c()
MYCN_frac_vec <- c()
MDM2_frac_vec <- c()
CDK4_frac_vec <- c()
ODC1_frac_vec <- c()
SMC6_frac_vec <- c()

for(sample in PN_pass_vec_all_runs){
  read_count_df <- window_calls_counts_normalised_amplicon_anno_all_run
  MYCN_length_norm <- sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MYCN",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MYCN",]$amplicon_anno)
  MDM2_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MDM2",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MDM2",]$amplicon_anno)
  CDK4_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "CDK4",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "CDK4",]$amplicon_anno)
  ODC1_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "ODC1",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "ODC1",]$amplicon_anno)
  SMC6_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "SMC6",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "SMC6",]$amplicon_anno)
  
  length_norm_sum <- sum(c(MYCN_length_norm,MDM2_length_norm,CDK4_length_norm,ODC1_length_norm,SMC6_length_norm))
  
  sample_vec <- c(sample_vec,sample)
  MYCN_frac_vec <- c(MYCN_frac_vec,(MYCN_length_norm / length_norm_sum))
  MDM2_frac_vec <- c(MDM2_frac_vec,(MDM2_length_norm / length_norm_sum))
  CDK4_frac_vec <- c(CDK4_frac_vec,(CDK4_length_norm / length_norm_sum))
  ODC1_frac_vec <- c(ODC1_frac_vec,(ODC1_length_norm / length_norm_sum))
  SMC6_frac_vec <- c(SMC6_frac_vec,(SMC6_length_norm / length_norm_sum))
  
}
gene <- c(rep("MYCN",length(sample_vec)),rep("MDM2",length(sample_vec)),rep("CDK4",length(sample_vec)),rep("ODC1",length(sample_vec)),rep("SMC6",length(sample_vec)))
fraction <- c(MYCN_frac_vec,MDM2_frac_vec,CDK4_frac_vec,ODC1_frac_vec,SMC6_frac_vec)
PN_fraction_df <- data.frame(gene,fraction)
PN_fraction_df$grp <- "Primary Nucleus"


###MN
sample_vec <- c()
MYCN_frac_vec <- c()
MDM2_frac_vec <- c()
CDK4_frac_vec <- c()
ODC1_frac_vec <- c()
SMC6_frac_vec <- c()

for(sample in MN_pass_vec_all_runs){
  read_count_df <- window_calls_counts_normalised_amplicon_anno_all_run
  MYCN_length_norm <- sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MYCN",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MYCN",]$amplicon_anno)
  MDM2_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MDM2",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MDM2",]$amplicon_anno)
  CDK4_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "CDK4",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "CDK4",]$amplicon_anno)
  ODC1_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "ODC1",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "ODC1",]$amplicon_anno)
  SMC6_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "SMC6",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "SMC6",]$amplicon_anno)
  
  length_norm_sum <- sum(c(MYCN_length_norm,MDM2_length_norm,CDK4_length_norm,ODC1_length_norm,SMC6_length_norm))
  
  sample_vec <- c(sample_vec,sample)
  MYCN_frac_vec <- c(MYCN_frac_vec,(MYCN_length_norm / length_norm_sum))
  MDM2_frac_vec <- c(MDM2_frac_vec,(MDM2_length_norm / length_norm_sum))
  CDK4_frac_vec <- c(CDK4_frac_vec,(CDK4_length_norm / length_norm_sum))
  ODC1_frac_vec <- c(ODC1_frac_vec,(ODC1_length_norm / length_norm_sum))
  SMC6_frac_vec <- c(SMC6_frac_vec,(SMC6_length_norm / length_norm_sum))
  
}

gene <- c(rep("MYCN",length(sample_vec)),rep("MDM2",length(sample_vec)),rep("CDK4",length(sample_vec)),rep("ODC1",length(sample_vec)),rep("SMC6",length(sample_vec)))
fraction <- c(MYCN_frac_vec,MDM2_frac_vec,CDK4_frac_vec,ODC1_frac_vec,SMC6_frac_vec)
MN_fraction_df <- data.frame(gene,fraction)
MN_fraction_df$grp <- "Micronucleus"


fraction_df <- rbind(MN_fraction_df,PN_fraction_df)

fraction_df$log2fraction <- log2(fraction_df$fraction+1)

dim(fraction_df[fraction_df$grp == "Micronucleus",])




#plot
fraction_df %>% 
  mutate(grp = factor(grp, levels = c("Primary Nucleus", "Micronucleus"))) %>% 
  ggplot(aes(x = gene, y = fraction, colour = grp, group = interaction(gene, grp))) +
  # Add the violin plot
  geom_violin(position = position_dodge(width = 0.7), trim = T) +  # Adjust alpha for transparency
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +  # Separate points with dodge
  # scale_color_okabeito(order=7:9)+
  theme_classic() +
  stat_summary(
    geom = "point",
    fun = "mean",  # Display the mean for each group
    col = "black",
    size = 5,
    shape = 95,
    position = position_dodge(width = 0.7)  # Dodge the means
  ) +
  labs(
    x = "Gene",
    y = "Fraction",
    colour = "Group"
  ) +
  theme(legend.position = "none")+ 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) 
dir.create("./output/figures/figure3f")
ggsave("./output/figures/figure3f/fraction_plotContent_stripped_violin.eps", width = 1.4,
       height = 1.4, device= "eps", dpi=450)

##permutation test

# Observed difference in means

results_permutation_test <- data.frame(
  amplicon = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (amplicon in c("MDM2","MYCN","ODC1","SMC6","CDK4")){
  fractions_ecDNA_MN <- fraction_df[fraction_df$grp == "Micronucleus" & fraction_df$gene == amplicon,]$fraction #test for all amplicons
  fractions_ecDNA_PN <- fraction_df[fraction_df$grp == "Primary Nucleus" & fraction_df$gene == amplicon,]$fraction
  
  obs_diff_mean <- mean(fractions_ecDNA_MN) - mean(fractions_ecDNA_PN)
  
  # Combine all data
  all_fractions_ecDNA <- c(fractions_ecDNA_MN,fractions_ecDNA_PN)
  
  # Number of permutations
  n_iter_perm <- 10000
  perm_diffs <- numeric(n_iter_perm)
  
  set.seed(43)
  for(i in 1:n_iter_perm){
    perm_labels <- sample(all_fractions_ecDNA)
    perm_diffs[i] <- mean(perm_labels[1:length(fractions_ecDNA_MN)]) -
      mean(perm_labels[(length(fractions_ecDNA_MN)+1):length(all_fractions_ecDNA)])
  }
  
  # Permutation p-value (two-sided)
  p_value <- mean(abs(perm_diffs) >= abs(obs_diff_mean))
  p_value
  results_permutation_test <- rbind(results_permutation_test, data.frame(
    amplicon = amplicon,
    p_value = p_value
  ))
  
}


#######
###Fig.3e
#######

##to identify linear fragments enriched in single MN, we compute the mean number of reads per chromosomes (nReads / chromosome length (nbin))
#we test for extreme outliers under the assumption, that all non enriched chromosomes will have a similar level of background mapping, with tru enrichment signal having a high above Background signal

sample_enrich_vec <- list()
chroms_enrich_vec <- list()

for(sample in MN_pass_vec_all_runs){
  #compute per chr read counts and n(bins), exclude circular regions
  chr_enrich <- window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename == sample & window_calls_counts_normalised_amplicon_anno_all_run$amplicon_anno == "z",] %>%
    group_by(chr) %>%
    summarise(
      total_nReads = sum(nReads, na.rm = TRUE),
      nBins = n(), norm = sum(nReads, na.rm = TRUE)/n(),
    )
  
  chr_enrich <- as.data.frame(chr_enrich)
  #remove chromosomes "chrM","NC_007605.1","chrY"
  chr_enrich_subset <- (chr_enrich[!(chr_enrich$chr %in% c("chrM","NC_007605.1","chrY")),])
  chr_enrich_scores <- chr_enrich_subset$norm
  names(chr_enrich_scores) <- chr_enrich_subset$chr
  #compute IQR
  Q1 <- quantile(chr_enrich_scores, 0.25)
  Q3 <- quantile(chr_enrich_scores, 0.75)
  IQR <- Q3 - Q1
  
  upper_bound <- Q3 + 10 * IQR
  
  # Identify outlier chromosomes
  outliers <- chr_enrich_scores[chr_enrich_scores > upper_bound]
  if (sum(outliers) > 0){
    print(sample)
    sample_enrich_vec <- c(sample_enrich_vec,sample)
    chroms_enrich_vec <- c(chroms_enrich_vec,list(names(outliers)))
  }
  if (sum(outliers) == 0){
    sample_enrich_vec <- c(sample_enrich_vec,sample)
    chroms_enrich_vec <- c(chroms_enrich_vec,list(names(outliers)))
  }
}



chr_enrich_df <- tibble(
  sample = sample_enrich_vec,
  chromosomes = chroms_enrich_vec  # No I() needed
)
chr_enrich_df <- as.data.frame(chr_enrich_df)


##now analysze each enriched chrom
sample_enrich_chr_len_sample <- c()
sample_enrich_chr_len_chr <- c()
sample_enrich_chr_len_nbin <- c()

for(sample in MN_pass_vec_all_runs){
  chrs_enrich_vec <- chr_enrich_df[chr_enrich_df$sample == sample,"chromosomes"]
  n_enriched_chr <- length(chrs_enrich_vec[[1]])
  
  if(n_enriched_chr > 0){
    exclude_chr_vec <- c(c("chrM","NC_007605.1","chrY"),chrs_enrich_vec[[1]])
    #now estimates background
    #90th quantile of non zero non enriched bins, filters out all zeros and disregards identified enriched chromosomes from above, to get a proper estimate of the background --> read threshold for enriched bin
    current_sample_enrich_df <- window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename == sample & window_calls_counts_normalised_amplicon_anno_all_run$amplicon_anno == "z" & !(window_calls_counts_normalised_amplicon_anno_all_run$chr %in% exclude_chr_vec),]
    filt__enrich_chr_nread <- current_sample_enrich_df$nReads
    read_thresh_enrich_bin <- quantile(filt__enrich_chr_nread[filt__enrich_chr_nread != 0], 0.90)
    
    #now subset each enriched chromosome from above and get the enriched chr length
    for (enriched_chromosome in (chrs_enrich_vec[[1]])){
      print(enriched_chromosome)
      current_sample_enrich_df_chr <-window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename == sample & window_calls_counts_normalised_amplicon_anno_all_run$chr == enriched_chromosome, ]$nReads
      enrich_chr_len <- length(current_sample_enrich_df_chr[current_sample_enrich_df_chr > read_thresh_enrich_bin])
      #print(enrich_chr_len)
      sample_enrich_chr_len_sample <- c(sample_enrich_chr_len_sample,sample)
      sample_enrich_chr_len_chr <- c(sample_enrich_chr_len_chr,enriched_chromosome)
      sample_enrich_chr_len_nbin <- c(sample_enrich_chr_len_nbin,enrich_chr_len)
      
      
    }
  }
  else{
    sample_enrich_chr_len_sample <- c(sample_enrich_chr_len_sample,sample)
    sample_enrich_chr_len_chr <- c(sample_enrich_chr_len_chr,"none")
    sample_enrich_chr_len_nbin <- c(sample_enrich_chr_len_nbin,0)
  }
  
  
}


per_ernrich_chr_df <- data.frame(sample = sample_enrich_chr_len_sample, chr = sample_enrich_chr_len_chr, len = sample_enrich_chr_len_nbin)


#now add the length of enriched ecDNA as additional score
ecDNA_enrich_len_vec <- c()
for (sample in per_ernrich_chr_df$sample) {
  ecDNA_enrich_len <- sum(window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename == sample & window_calls_counts_normalised_amplicon_anno_all_run$amplicon_anno != "z" & window_calls_counts_normalised_amplicon_anno_all_run$window_call_final == 1, ]$window_call_final)
  ecDNA_enrich_len_vec <- c(ecDNA_enrich_len_vec,ecDNA_enrich_len)
}

per_ernrich_chr_df_final <- cbind(per_ernrich_chr_df,ecDNA_enrich_len_vec)
#filter out samples with linChr < 1mb and no ecDNA
per_ernrich_chr_df_final_filtered <- per_ernrich_chr_df_final[!(per_ernrich_chr_df_final$len > 0 & per_ernrich_chr_df_final$len < 40 ), ]


per_ernrich_chr_df_final_filtered$linchr <- (per_ernrich_chr_df_final_filtered$len *25000) / 1000000 #convert to Mb
per_ernrich_chr_df_final_filtered$circchr <- (per_ernrich_chr_df_final_filtered$ecDNA_enrich_len_vec *25000) / 1000000 #convert to Mb



##add relative size of lin fragments to whole chromosome
current_chr_enrich_len_rel_vec <- c()
for (sample in per_ernrich_chr_df_final_filtered$sample) {
  current_chr_enrich <- per_ernrich_chr_df_final_filtered[per_ernrich_chr_df_final_filtered$sample == sample,"chr"]
  current_chr_enrich_len <- dim(window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$chr == current_chr_enrich &  window_calls_counts_normalised_amplicon_anno_all_run$amplicon_anno == "z" & window_calls_counts_normalised_amplicon_anno_all_run$samplename == sample,])[1]
  current_chr_enrich_len_rel <- per_ernrich_chr_df_final_filtered[per_ernrich_chr_df_final_filtered$sample == sample,"len"]/current_chr_enrich_len
  current_chr_enrich_len_rel_vec <- c(current_chr_enrich_len_rel_vec,current_chr_enrich_len_rel)
}

per_ernrich_chr_df_final_filtered$rel_len <- current_chr_enrich_len_rel_vec
per_ernrich_chr_df_final_filtered[is.na(per_ernrich_chr_df_final_filtered)] <- 0

#plot
ggplot(per_ernrich_chr_df_final_filtered, aes(x = (linchr), y = (circchr), colour = chr)) +
  geom_point(size = 3) +
  labs(
    x = "cumulative linear fragment length (Mb)",
    y = "cumulative ecDNA fragment length > 1Mb (Mb)"
  ) + theme_classic() +
  theme(legend.position="none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank())

ggplot() +
  # Points with rel_size == 0 → squares
  geom_point(
    data = per_ernrich_chr_df_final_filtered[per_ernrich_chr_df_final_filtered$rel_len == 0 & per_ernrich_chr_df_final_filtered$circchr > 0, ],
    aes(x = linchr, y = circchr, color = chr),
    shape = 18,    # Square
    size = 3       # Fixed size for visibility
  ) +
  
  # Points with rel_size > 0 → circles, sized by rel_size
  geom_point(
    data = per_ernrich_chr_df_final_filtered[per_ernrich_chr_df_final_filtered$rel_len > 0 & per_ernrich_chr_df_final_filtered$circchr > 0, ],
    aes(x = linchr, y = circchr, color = chr, size = rel_len),
    shape = 16      # Default: filled circle
  ) +
  
  # Points with rno lin no cird
  geom_point(
    data = per_ernrich_chr_df_final_filtered[per_ernrich_chr_df_final_filtered$linchr == 0 & per_ernrich_chr_df_final_filtered$circchr == 0, ],
    aes(x = linchr, y = circchr, color = chr),
    shape = 17,    # Square
    size = 3      
  )+
  
  scale_size(range = c(2, 10)) +  # Adjust min/max size
  theme_classic()#+
#  theme(legend.position="none")+
#  theme(axis.title.x=element_blank(),
#        axis.text.x=element_blank(),
#        axis.title.y=element_blank(),
#        axis.text.y=element_blank())

dir.create("./output/figures/figure3e")
ggsave("./output/figures/figure3e/chromosome_enrichment_stripped_withlabels.eps", width = 5.4,
       height = 5.4, device= "eps", dpi=450)


#pieplot

# Create data frame
df_ecdnapos <- data.frame(
  group = c("ecDNA", "ecDNA_lin"),
  value = c(24, 2)
)

# Create pie chart using coord_polar
ggplot(df_ecdnapos, aes(x = "", y = value, fill = group)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  theme_void() +
  theme(legend.position="none")
ggsave("./output/figures/figure3e/ecDNA_pie_chart.eps", width = 1.4,
       height = 1.4, device= "eps", dpi=450)


#######
###Ext Fig. 9c (signal plot)
#######


signal_mat_PN <- read.csv("/Users/robinxu/Documents/Projects/Brueckner_et_al/sMNseq/signalplot/PN_ecNDA_reg_25000bp.txt", sep="\t", header = 1)
signal_mat_MN <- read.csv("/Users/robinxu/Documents/Projects/Brueckner_et_al/sMNseq/signalplot/MN_ecNDA_reg_25000bp.txt", sep="\t", header = 1)
signal_mat_mMN <- read.csv("/Users/robinxu/Documents/Projects/Brueckner_et_al/sMNseq/signalplot/pMN_ecNDA_reg_25000bp.txt", sep="\t", header = 1)
signal_mat_Empty <- read.csv("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/Signal_plot_analysis/Empty_ecNDA_reg_25000bp.txt", sep="\t", header = 1)


chunksize <- dim(signal_mat_PN)[2]/18

get_region_mean <- function(vec, chunk_size=chunksize){
  
  return(rowMeans(data.frame(split(vec,ceiling(seq_along(vec) / chunk_size))), na.rm=T))
  
} 


PN_mean <- apply(signal_mat_PN,1,function(x) get_region_mean(as.numeric(as.vector(x)), chunksize))
MN_mean <- apply(signal_mat_MN,1,function(x) get_region_mean(as.numeric(as.vector(x)), chunksize))
mMN_mean <- apply(signal_mat_mMN,1,function(x) get_region_mean(as.numeric(as.vector(x)), chunksize))
Empty_mean <- apply(signal_mat_Empty,1,function(x) get_region_mean(as.numeric(as.vector(x)), chunksize))

#plot mean signal

relative_coord <- seq(1:length(rowMeans(PN_mean)))
rowMeans(MN_mean)
rowMeans(mMN_mean)

signal_plot_df <- rbind(data.frame(mean_CPM_Signal = rowMeans(PN_mean), grp = "PN", rel_coord = relative_coord),data.frame(mean_CPM_Signal = rowMeans(MN_mean), grp = "MN", rel_coord = relative_coord),data.frame(mean_CPM_Signal = rowMeans(mMN_mean), grp = "mMN", rel_coord = relative_coord))

dir.create("./output/figures/Ext_figure9c")
pdf("./output/figures/Ext_figure9c/signal_plot.pdf", width =4 , height = 4)

ggplot(signal_plot_df ,aes(x=rel_coord, y=(mean_CPM_Signal), group=grp, color=grp)) +
  geom_line() +
  ylab("mean CPM Signal") + theme_classic()

dev.off()




#######
###Ext Fig. 9b (Upset plot)
#######


#make a vectors for different amplicons, to which positve sample will be added
SMC6 <- c()
MYCN <- c()
ODC1 <- c()
MDM2 <- c()
CDK4 <- c()

for(x in PN_pass_vec_all_runs){
  
  sample_coords <- window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename == x & window_calls_counts_normalised_amplicon_anno_all_run$window_call_final > 0, ][1:3] #subset sample and positive windows
  #test for amplicon overlap
  sample_coord_grange <-  makeGRangesFromDataFrame(sample_coords)
  SMC6_sum <- sum(countOverlaps(SMC6_grange,sample_coord_grange))
  MYCN_sum <- sum(countOverlaps(MYCN_grange,sample_coord_grange))
  ODC1_sum <- sum(countOverlaps(ODC1_grange,sample_coord_grange))
  MDM2_sum <- sum(countOverlaps(MDM2_grange,sample_coord_grange))
  CDK4_sum <- sum(countOverlaps(CDK4_grange,sample_coord_grange))
  
  if(SMC6_sum > 0){
    SMC6 <- c(SMC6,x)}
  if(MYCN_sum > 0){
    MYCN <- c(MYCN,x)}
  if(ODC1_sum > 0){
    ODC1 <- c(ODC1,x)}
  if(MDM2_sum > 0){
    MDM2 <- c(MDM2,x)}
  if(CDK4_sum > 0){
    CDK4 <- c(MDM2,x)}
}

Amplicon_list <- list(MDM2 = MDM2, SMC6 = SMC6, MYCN = MYCN, ODC1 = ODC1, CDK4=CDK4 )

dir.create("./output/figures/Ext_figure9b")
pdf("./output/figures/Ext_figure9b/PN_upset.pdf", width = 5, height = 5 )

upset(fromList(Amplicon_list), order.by = "freq", mainbar.y.max = 20, empty.intersections = "on", sets.bar.color = "#8B0000", matrix.color = "#8B0000", main.bar.color = "#8B0000", nintersects = 5, keep.order = T, sets = c("MYCN","CDK4","MDM2","ODC1","SMC6"))
dev.off()

##MN Upset

SMC6 <- c()
MYCN <- c()
ODC1 <- c()
MDM2 <- c()
CDK4 <- c()

for(x in MN_pass_vec_all_runs){
  
  sample_coords <- window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename == x & window_calls_counts_normalised_amplicon_anno_all_run$window_call_final > 0, ][1:3] #subset sample and positive windows
  #test for amplicon overlap
  sample_coord_grange <-  makeGRangesFromDataFrame(sample_coords)
  SMC6_sum <- sum(countOverlaps(SMC6_grange,sample_coord_grange))
  MYCN_sum <- sum(countOverlaps(MYCN_grange,sample_coord_grange))
  ODC1_sum <- sum(countOverlaps(ODC1_grange,sample_coord_grange))
  MDM2_sum <- sum(countOverlaps(MDM2_grange,sample_coord_grange))
  CDK4_sum <- sum(countOverlaps(CDK4_grange,sample_coord_grange))
  
  if(SMC6_sum > 0){
    SMC6 <- c(SMC6,x)}
  if(MYCN_sum > 0){
    MYCN <- c(MYCN,x)}
  if(ODC1_sum > 0){
    ODC1 <- c(ODC1,x)}
  if(MDM2_sum > 0){
    MDM2 <- c(MDM2,x)}
  if(CDK4_sum > 0){
    CDK4 <- c(MDM2,x)}
}

Amplicon_list1 <- list(MDM2 = MDM2, SMC6 = SMC6, MYCN = MYCN, ODC1 = ODC1, CDK4=CDK4 )
pdf("./output/figures/Ext_figure9b/MN_upset.pdf", width = 5, height = 5 )

upset(fromList(Amplicon_list1), order.by = "freq", mainbar.y.max = 20,empty.intersections = "on", nintersects = 5, keep.order = T)
dev.off()


#######
###Ext Fig. 9d (paired MN/PN plot)
#######

#get a table where per sample occourences of amplicons are depicted
#takes data from upset plot above

#merge MN and PN
merged_Amplicon_list <- mapply(base::c, Amplicon_list, Amplicon_list1, SIMPLIFY = FALSE)

all_samples <- sort(unique(unlist(merged_Amplicon_list)))

# Create presence/absence matrix
presence_matrix <- sapply(merged_Amplicon_list, function(samples) all_samples %in% samples)

# Transpose to get genes as rows
presence_df <- as.data.frame(t(presence_matrix))
colnames(presence_df) <- all_samples
# Convert logical to numeric (TRUE/FALSE → 1/0)
presence_df[] <- lapply(presence_df, as.integer)

#now subset for pairs:

presence_df[, paste0(c("A3","B3","C2","D2","E1","F1","G4","H4"),"_R2")]
##plot


col_ann_class <- HeatmapAnnotation(
  class = c("MN","PN","MN","PN","MN","PN","MN","PN"),
  col = list(class = c("MN" = "#1f77b4", "PN" = "#ff7f0e")),
  annotation_name_side = "left"
)

dir.create("./output/figures/Ext_figure9d")
pdf("./output/figures/Ext_figure9d/pairs.pdf", width = 5, height = 5 )
Heatmap(
  as.matrix(presence_df[, paste0(c("A3","B3","C2","D2","E1","F1","G4","H4"),"_R2")]),
  name = "Amplicon_presence",
  col = c("0" = "white", "1" = "red"),
  top_annotation = col_ann_class,
  column_split = rep(1:(ncol(as.matrix(presence_df[, paste0(c("A3","B3","C2","D2","E1","F1","G4","H4"),"_R2")])) / 2), each = 2),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  rect_gp = gpar(col = "black", lwd = 1), 
)
dev.off()


sink("sessionInfo.txt")
sessionInfo()
sink()