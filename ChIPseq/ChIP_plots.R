##ChIP visualisation
##VIZ
library("Gviz")
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsignif)

setwd("/Users/robinxu/Documents/Projects/Brueckner_et_al/sequencing_final_publication/ecMN_final/ChIPseq/plots")

##visualization MYC
gtf <- import("/Volumes/Robin_Xu_APFS/brueckner_et_al_data_upload/ChIP/processed/hg38.ncbiRefSeq.gtf.gz")
# MYC anno
myc <- gtf[gtf$transcript_id == "NM_001354870.1" & gtf$type == "exon"]

myc_track <- GeneRegionTrack(myc,
                             name = "MYC",
                             showId = T,
                             fill = "orange",
                             stacking = "dense",
                             chromosome = "chr8")
axisTrack <- GenomeAxisTrack()

bw_track <- DataTrack(range = "/Volumes/Robin Xu/COLO320DM_H3K27ac_HU_sub_DMSO.bw", genome = "hg38", type = "l", 
                      chromosome = "chr8", name = "bigwig",ylim = c(-28,10))

displayPars(grTrack) <- list(stacking = "pack")
pdf("H3K27ac_MYC.pdf")
plotTracks(
  trackList = list(axisTrack,bw_track,myc_track), from = 127733090+1000 ,to = 127744607+1000, type = "histogram",sizes     = c(0.4,1,1)
)
dev.off()

#visualize normalised CPMs

#read in table from bigwigcompare
CPM_bw_input_norm <- read.csv("/Volumes/Robin_Xu_APFS/brueckner_et_al_data_upload/ChIP/processed/H3K27ac_hu_sub_dmso_Inputnorm.counts", sep="\t", header = T, col.names = c("chr","start","end","HU_Ac","UT_Ac"))



#subset for exons
bins_gr <- GRanges("chr8", IRanges(start = CPM_bw_input_norm$start, end = CPM_bw_input_norm$end))
hits <- overlapsAny(bins_gr, myc)

CPM_bw_myc <- CPM_bw_input_norm[hits, ]


CPM_long <- CPM_bw_myc %>%
  pivot_longer(cols = c(HU_Ac, UT_Ac),
               names_to = "treatment",
               values_to = "CPM") %>%
  mutate(treatment = gsub("_Ac", "", treatment))


pval <- wilcox.test(CPM_bw_myc$HU_Ac, CPM_bw_myc$UT_Ac)$p.value

#plot raw cpm values
ggplot(CPM_long, aes(x = treatment, y = CPM, fill = treatment)) +
  geom_violin(width = 0.5) +
  xlab("") + ylab("normalized CPM") +
  geom_signif(comparisons = list(c("HU", "UT")),
              map_signif_level = TRUE,
              annotation = formatC(pval, digits = 1)) +
  theme_classic() +
  stat_summary(fun = "mean", geom = "pointrange", colour = "black") + ylim(-40,40)


