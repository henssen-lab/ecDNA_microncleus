library(tidyr)
library(stringr)
library(ggplot2)
library(fitdistrplus)
library(MASS)
library(ggpubr)
library(ggsignif)
library(DescTools)
library(GenomicRanges)
library(regioneR)
library(UpSetR)


source("./R/data_handler.R")
source("./R/compute_thresholds.R")
source("./R/classify_windows.R")
source("./R/QC_tools.R")
source("./R/analysis_tools.R")


#0. load meta table and reconstructions
meta_r1 <- read.csv("./data/meta_r1.txt", sep="\t", header = F, col.names = c("path","sample","class"))
reconstrucion_bed <- read.csv("./data/TR14_hg38_reconstruction_decoil_with_CDK4_shasta_assembly.bed", header = 1, sep = "\t")
#1. read files
bed_file_list_r1 <- load_bed_files("./data/meta_r1.txt")
#2. compute thresholds
dir.create("./output/qc/r1/CMF", recursive = T)
sample_threshold_df_r1 <- generate_threshold_df(bed_file_list_r1,P=0.01,percentile=0.99,plot=T,out_dir_plots="./output/qc/r1/CMF")
#3. classify windows
dir.create("./output/qc/r1/plots", recursive = T)
window_calls_r1 <- classify_windows(bed_file_list_r1,sample_threshold_df_r1,plot = T, plot_outdir = "./output/qc/r1/plots",window_proximity_n=4)


#4. add read counts from external files and add window counts to the meta table
window_calls_counts_r1 <- count_positive_windows(window_calls_r1)

files <- list.files(path="./data/read_counts/r1", pattern="*.count", full.names=TRUE, recursive=FALSE)
read_counts <- data.frame()
for (file in files){
  sample <- read.table(file, header=F)
  sample$V2 <- str_split_1(str_split_1(file, "/")[5],"\\.")[1]
  print(sample)
  read_counts <- rbind(read_counts,sample)
}
colnames(read_counts) <- c("n_unique_reads","sample")

meta_qc_r1 <- merge(meta_r1, window_calls_counts_r1, by='sample')
meta_qc_r1 <- merge(meta_qc_r1,read_counts , by='sample')


###
## analysis
###



#5. CPM normalise read counts
window_calls_counts_normalised_r1 <- normalise_CPM(window_calls_r1,meta_qc_r1$n_unique_reads)

#make different granges for the different amplicons
SMC6_grange<- makeGRangesFromDataFrame(reconstrucion_bed[reconstrucion_bed$circ_id=="SMC6",],keep.extra.columns = TRUE)
MYCN_grange<- makeGRangesFromDataFrame(reconstrucion_bed[reconstrucion_bed$circ_id=="MYCN",],keep.extra.columns = TRUE)
ODC1_grange<- makeGRangesFromDataFrame(reconstrucion_bed[reconstrucion_bed$circ_id=="ODC1",],keep.extra.columns = TRUE)
MDM2_grange<- makeGRangesFromDataFrame(reconstrucion_bed[reconstrucion_bed$circ_id=="MDM2",],keep.extra.columns = TRUE)
CDK4_grange<- makeGRangesFromDataFrame(reconstrucion_bed[reconstrucion_bed$circ_id=="CDK4",],keep.extra.columns = TRUE)

#6. first filter
MN_pass_vec1_r1 <- meta_qc_r1[as.numeric(meta_qc_r1$n_unique_reads) >= 100000 & as.numeric(meta_qc_r1$n_positive_windows_vec) >=20 & meta_qc_r1$class == "MN",]

#7. subset window_calls_counts_normalised df and determine consensus sequence of circular DNA in MN
window_calls_counts_normalised_MN_r1 <- window_calls_counts_normalised_r1[window_calls_counts_normalised_r1$samplename %in% MN_pass_vec1_r1$sample,]
MN_consensus_calls_mat_r1 <- unstack(window_calls_counts_normalised_MN_r1,window_call_final~samplename)
MN_consensus_calls_mat_r1$rowsum <- rowSums(MN_consensus_calls_mat_r1)
MN_consensus_calls_vec_r1 <- rownames(MN_consensus_calls_mat_r1[MN_consensus_calls_mat_r1$rowsum>=5,])
MN_consensus_regions_r1<- window_calls_counts_normalised_MN_r1[window_calls_counts_normalised_MN_r1$samplename == "A1",][MN_consensus_calls_vec_r1,c(1,2,3)] ##########
MN_consensus_regions_r1 <- MN_consensus_regions_r1[MN_consensus_regions_r1$chr!="chrM",] #remove chrM

#8. overlap identified consensus regions with reference
reconstrucion_bed_grange <- makeGRangesFromDataFrame(reconstrucion_bed,keep.extra.columns = TRUE)
consensus_regions_grange_r1 <- makeGRangesFromDataFrame(MN_consensus_regions_r1) 

MNsequences_on_ecDNA_r1 <- subsetByOverlaps(consensus_regions_grange_r1,reconstrucion_bed_grange) #extract indices of all regions in the consensus, which are overlapping the ecDNA reconstruction
index_MNsequences_on_ecDNA_r1 <-rownames(as.data.frame(MNsequences_on_ecDNA_r1)) # for later quantitaive analysis

#9 add bin counts
#assign amplicons to each bin
amplicon_assignment_vec_r1 <- c()
for(i in seq_along(MNsequences_on_ecDNA_r1)) {
  grange_MN <- MNsequences_on_ecDNA_r1[i]
  current_gene <-  "n"
  SMC6_sum <- sum(countOverlaps(SMC6_grange,grange_MN))
  MYCN_sum <- sum(countOverlaps(MYCN_grange,grange_MN))
  ODC1_sum <- sum(countOverlaps(ODC1_grange,grange_MN))
  MDM2_sum <- sum(countOverlaps(MDM2_grange,grange_MN))
  CDK4_sum <- sum(countOverlaps(CDK4_grange,grange_MN))
  
  if(SMC6_sum != 0){current_gene <- "SMC6"}
  if(MYCN_sum != 0){current_gene <- "MYCN"}
  if(ODC1_sum != 0){current_gene <- "ODC1"}
  if(MDM2_sum != 0){current_gene <- "MDM2"}
  if(CDK4_sum != 0){current_gene <- "CDK4"}
  
  amplicon_assignment_vec_r1 <- c(amplicon_assignment_vec_r1,current_gene)
}

amplicon_index_df_r1 <- cbind(as.data.frame(MNsequences_on_ecDNA_r1)[,1:3],amplicon_assignment_vec_r1)
window_calls_counts_normalised_amplicon_anno_r1 <- window_calls_counts_normalised_r1
A1_df <- window_calls_counts_normalised_r1[window_calls_counts_normalised_r1$samplename == "A1",]
A1_df$amplicon <- "z"
A1_df$amplicon[match(row.names(amplicon_index_df_r1), row.names(A1_df))] <- amplicon_index_df_r1$amplicon_assignment_vec_r1
A1_df[row.names(A1_df) %in% row.names(amplicon_index_df_r1), ]

window_calls_counts_normalised_amplicon_anno_r1$amplicon_anno <- rep(A1_df$amplicon, length(unique(window_calls_counts_normalised_r1$samplename)))

#bin counts
sample_count_vec = c()
bp_covered_vec_vec = c()
n_short_fragments_vec = c()
n_long_fragments_vec = c()
for(sample in unique(window_calls_counts_normalised_amplicon_anno_r1$samplename)){
  window_calls_counts_normalised_amplicon_anno_sample_df_r1 <- window_calls_counts_normalised_amplicon_anno_r1[window_calls_counts_normalised_amplicon_anno_r1$samplename == sample & window_calls_counts_normalised_amplicon_anno_r1$window_call_final == 1,]
  fragment_i <- 1
  current_index <- 0
  prev_index <- 0
  fragment_vec <- c()
  if (dim(window_calls_counts_normalised_amplicon_anno_sample_df_r1)[1] != 0){
    for (i in seq(1:dim(window_calls_counts_normalised_amplicon_anno_sample_df_r1)[1])){
      sample_row <- window_calls_counts_normalised_amplicon_anno_sample_df_r1[i,]
      current_index <- as.numeric(row.names(sample_row))
      diff <- current_index - prev_index 
      if (i == 1){
        fragment_vec <- c(fragment_vec,fragment_i)
        prev_index <- current_index
      }
      if (diff <= 3 & i != 1){
        fragment_vec <- c(fragment_vec,fragment_i)
        prev_index <- current_index
      }
      if (diff > 3 & i != 1){
        fragment_i <- fragment_i + 1
        fragment_vec <- c(fragment_vec,fragment_i)
        prev_index <- current_index
      }
    }
    
    window_calls_counts_normalised_amplicon_anno_sample_df_r1$fragment <- fragment_vec
    bp_covered <- sum(table(window_calls_counts_normalised_amplicon_anno_sample_df_r1$fragment)) * 25000
    n_short_fragments <- sum(table(window_calls_counts_normalised_amplicon_anno_sample_df_r1$fragment) == 3)
    n_long_fragments <- dim(table(window_calls_counts_normalised_amplicon_anno_sample_df_r1$fragment)) - sum(table(window_calls_counts_normalised_amplicon_anno_sample_df_r1$fragment) == 3)
    
    sample_count_vec <- c(sample_count_vec,sample)
    bp_covered_vec_vec <- c(bp_covered_vec_vec,bp_covered)
    n_short_fragments_vec <- c(n_short_fragments_vec,n_short_fragments)
    n_long_fragments_vec <- c(n_long_fragments_vec,n_long_fragments)
    
  }
  else {
    sample_count_vec <- c(sample_count_vec,sample)
    bp_covered_vec_vec <- c(bp_covered_vec_vec,0)
    n_short_fragments_vec <- c(n_short_fragments_vec,0)
    n_long_fragments_vec <- c(n_long_fragments_vec,0)
  }
}
#long fragments >3 bins
bin_count_meta_r1 <- data.frame(sample = sample_count_vec, bp_covered = bp_covered_vec_vec, n_short_fragments = n_short_fragments_vec, n_long_fragments = n_long_fragments_vec)

meta_qc_ext_r1 <- merge(meta_qc_r1,bin_count_meta_r1, by  = "sample") #extended meta
meta_qc_ext_r1[meta_qc_ext_r1$class == "MN", ]

#9. second filter: filter samples based on exclusion cirteria (here: MN, <= 20 windows, <= 100000 reads, and exclude if only short fragments)
MN_pass_vec_r1 <- meta_qc_ext_r1[as.numeric(meta_qc_ext_r1$n_unique_reads) >= 100000 & as.numeric(meta_qc_ext_r1$n_positive_windows_vec) >=20 & meta_qc_ext_r1$class == "MN" & meta_qc_ext_r1$n_long_fragments != 0,]
pMN_pass_vec_r1 <- meta_qc_ext_r1[as.numeric(meta_qc_ext_r1$n_unique_reads) >= 100000 & as.numeric(meta_qc_ext_r1$n_positive_windows_vec) >=20 & meta_qc_ext_r1$class == "mMN" & meta_qc_ext_r1$n_long_fragments != 0,]
PN_pass_vec_r1 <- meta_qc_r1[meta_qc_r1$class == "PN" & as.numeric(meta_qc_r1$n_positive_windows_vec) >= 50,]



#############################################################################


#0. load meta table and reconstructions
meta_r2 <- read.csv("./data/meta_r2.txt", sep="\t", header = F, col.names = c("path","sample","class"))
#1. read files
bed_file_list_r2 <- load_bed_files("./data/meta_r2.txt")
#2
dir.create("./output/qc/r2/CMF", recursive = T)
sample_threshold_df_r2 <- generate_threshold_df(bed_file_list_r2,P=0.01,percentile=0.99,plot=T,out_dir_plots="./output/qc/r2/CMF")
#3. classify windows
dir.create("./output/qc/r2/plots", recursive = T)
window_calls_r2 <- classify_windows(bed_file_list_r2,sample_threshold_df_r2,plot = T, plot_outdir = "./output/qc/r2/plots",window_proximity_n=4)

#4. add read counts from external files and add window counts to the meta table
window_calls_counts_r2 <- count_positive_windows(window_calls_r2)

files_r2 <- list.files(path="./data/read_counts/r2", pattern="*.count", full.names=TRUE, recursive=FALSE)
read_counts_r2 <- data.frame()
for (file in files_r2){
  sample <- read.table(file, header=F)
  sample$V2 <- str_split_1(str_split_1(file, "/")[5],"\\.")[1]
  print(sample)
  read_counts_r2 <- rbind(read_counts_r2,sample)
}
colnames(read_counts_r2) <- c("n_unique_reads","sample")

meta_qc_r2 <- merge(meta_r2, window_calls_counts_r2, by='sample')
meta_qc_r2 <- merge(meta_qc_r2,read_counts_r2 , by='sample')


###
## analysis
###

#5. CPM normalise read counts
window_calls_counts_normalised_r2 <- normalise_CPM(window_calls_r2,meta_qc_r2$n_unique_reads)


#6. first filter
MN_pass_vec1_r2 <- meta_qc_r2[as.numeric(meta_qc_r2$n_unique_reads) >= 100000 & as.numeric(meta_qc_r2$n_positive_windows_vec) >=20 & meta_qc_r2$class == "MN",]

#7. subset window_calls_counts_normalised df and determine consensus sequence of circular DNA in MN
window_calls_counts_normalised_MN_r2 <- window_calls_counts_normalised_r2[window_calls_counts_normalised_r2$samplename %in% MN_pass_vec1_r2$sample,]
MN_consensus_calls_mat_r2 <- unstack(window_calls_counts_normalised_MN_r2,window_call_final~samplename)
MN_consensus_calls_mat_r2$rowsum <- rowSums(MN_consensus_calls_mat_r2)
MN_consensus_calls_vec_r2 <- rownames(MN_consensus_calls_mat_r2[MN_consensus_calls_mat_r2$rowsum>=5,])
MN_consensus_regions_r2 <- window_calls_counts_normalised_MN_r2[window_calls_counts_normalised_MN_r2$samplename == "A1_R2",][MN_consensus_calls_vec_r2,c(1,2,3)] ##########
MN_consensus_regions_r2 <- MN_consensus_regions_r2[MN_consensus_regions_r2$chr!="chrM",] #remove chrM

#8. overlap identified consensus regions with reference
reconstrucion_bed_grange <- makeGRangesFromDataFrame(reconstrucion_bed,keep.extra.columns = TRUE)
consensus_regions_grange_r2 <- makeGRangesFromDataFrame(MN_consensus_regions_r2) 

MNsequences_on_ecDNA_r2 <- subsetByOverlaps(consensus_regions_grange_r2,reconstrucion_bed_grange) #extract indices of all regions in the consensus, which are overlapping the ecDNA reconstruction
index_MNsequences_on_ecDNA_r2 <-rownames(as.data.frame(MNsequences_on_ecDNA_r2)) # for later quantitaive analysis



#assign amplicons to each bin
amplicon_assignment_vec_r2 <- c()
for(i in seq_along(MNsequences_on_ecDNA_r2)) {
  grange_MN <- MNsequences_on_ecDNA_r2[i]
  current_gene <-  "n"
  SMC6_sum <- sum(countOverlaps(SMC6_grange,grange_MN))
  MYCN_sum <- sum(countOverlaps(MYCN_grange,grange_MN))
  ODC1_sum <- sum(countOverlaps(ODC1_grange,grange_MN))
  MDM2_sum <- sum(countOverlaps(MDM2_grange,grange_MN))
  CDK4_sum <- sum(countOverlaps(CDK4_grange,grange_MN))
  
  if(SMC6_sum != 0){current_gene <- "SMC6"}
  if(MYCN_sum != 0){current_gene <- "MYCN"}
  if(ODC1_sum != 0){current_gene <- "ODC1"}
  if(MDM2_sum != 0){current_gene <- "MDM2"}
  if(CDK4_sum != 0){current_gene <- "CDK4"}
  
  amplicon_assignment_vec_r2 <- c(amplicon_assignment_vec_r2,current_gene)
}

amplicon_index_df_r2 <- cbind(as.data.frame(MNsequences_on_ecDNA_r2)[,1:3],amplicon_assignment_vec_r2)
window_calls_counts_normalised_amplicon_anno_r2 <- window_calls_counts_normalised_r2
A1_df_r2 <- window_calls_counts_normalised_r2[window_calls_counts_normalised_r2$samplename == "A1_R2",]
A1_df_r2$amplicon <- "z"
A1_df_r2$amplicon[match(row.names(amplicon_index_df_r2), row.names(A1_df_r2))] <- amplicon_index_df_r2$amplicon_assignment_vec
A1_df_r2[row.names(A1_df_r2) %in% row.names(amplicon_index_df_r2), ]

window_calls_counts_normalised_amplicon_anno_r2$amplicon_anno <- rep(A1_df_r2$amplicon, length(unique(window_calls_counts_normalised_r2$samplename)))


####
#bin counts


#bin counts
sample_count_vec = c()
bp_covered_vec_vec = c()
n_short_fragments_vec = c()
n_long_fragments_vec = c()
for(sample in unique(window_calls_counts_normalised_amplicon_anno_r2$samplename)){
  window_calls_counts_normalised_amplicon_anno_sample_df <- window_calls_counts_normalised_amplicon_anno_r2[window_calls_counts_normalised_amplicon_anno_r2$samplename == sample & window_calls_counts_normalised_amplicon_anno_r2$window_call_final == 1,]
  fragment_i <- 1
  current_index <- 0
  prev_index <- 0
  fragment_vec <- c()
  if (dim(window_calls_counts_normalised_amplicon_anno_sample_df)[1] != 0){
    for (i in seq(1:dim(window_calls_counts_normalised_amplicon_anno_sample_df)[1])){
      sample_row <- window_calls_counts_normalised_amplicon_anno_sample_df[i,]
      current_index <- as.numeric(row.names(sample_row))
      diff <- current_index - prev_index 
      if (i == 1){
        fragment_vec <- c(fragment_vec,fragment_i)
        prev_index <- current_index
      }
      if (diff <= 3 & i != 1){
        fragment_vec <- c(fragment_vec,fragment_i)
        prev_index <- current_index
      }
      if (diff > 3 & i != 1){
        fragment_i <- fragment_i + 1
        fragment_vec <- c(fragment_vec,fragment_i)
        prev_index <- current_index
      }
    }
    
    window_calls_counts_normalised_amplicon_anno_sample_df$fragment <- fragment_vec
    bp_covered <- sum(table(window_calls_counts_normalised_amplicon_anno_sample_df$fragment)) * 25000
    n_short_fragments <- sum(table(window_calls_counts_normalised_amplicon_anno_sample_df$fragment) == 3)
    n_long_fragments <- dim(table(window_calls_counts_normalised_amplicon_anno_sample_df$fragment)) - sum(table(window_calls_counts_normalised_amplicon_anno_sample_df$fragment) == 3)
    
    sample_count_vec <- c(sample_count_vec,sample)
    bp_covered_vec_vec <- c(bp_covered_vec_vec,bp_covered)
    n_short_fragments_vec <- c(n_short_fragments_vec,n_short_fragments)
    n_long_fragments_vec <- c(n_long_fragments_vec,n_long_fragments)
    
  }
  else {
    sample_count_vec <- c(sample_count_vec,sample)
    bp_covered_vec_vec <- c(bp_covered_vec_vec,0)
    n_short_fragments_vec <- c(n_short_fragments_vec,0)
    n_long_fragments_vec <- c(n_long_fragments_vec,0)
  }
}
#long fragments >3 bins
bin_count_meta_r2 <- data.frame(sample = sample_count_vec, bp_covered = bp_covered_vec_vec, n_short_fragments = n_short_fragments_vec, n_long_fragments = n_long_fragments_vec)

meta_qc_ext_r2 <- merge(meta_qc_r2,bin_count_meta_r2, by  = "sample") #extended meta


#9. second filter: filter samples based on exclusion cirteria (here: MN, <= 20 windows, <= 100000 reads, and exclude if only short fragments)
MN_pass_vec_r2 <- meta_qc_ext_r2[as.numeric(meta_qc_ext_r2$n_unique_reads) >= 100000 & as.numeric(meta_qc_ext_r2$n_positive_windows_vec) >=20 & meta_qc_ext_r2$class == "MN" & meta_qc_ext_r2$n_long_fragments != 0,]
PN_pass_vec_r2 <- meta_qc_ext_r2[meta_qc_ext_r2$class == "PN" & as.numeric(meta_qc_ext_r2$n_positive_windows_vec) >= 50,]

##merge all passes

MN_pass_vec_all_runs <- c(MN_pass_vec_r2$sample,MN_pass_vec_r1$sample)
PN_pass_vec_all_runs <- c(PN_pass_vec_r2$sample,PN_pass_vec_r1$sample)

pMN_pass_vec_all_runs <- pMN_pass_vec_r1$sample


length(MN_pass_vec_all_runs)
length(PN_pass_vec_all_runs)
length(pMN_pass_vec_all_runs)


##
####add _run2 to combine with pilot:

window_calls_counts_normalised_amplicon_anno_all_run <- rbind(window_calls_counts_normalised_amplicon_anno_r2,window_calls_counts_normalised_amplicon_anno_r1)

#subset table for all samples
window_calls_counts_normalised_amplicon_anno_all_run_sample_pass <- window_calls_counts_normalised_amplicon_anno_all_run[window_calls_counts_normalised_amplicon_anno_all_run$samplename %in% c(MN_pass_vec_all_runs,PN_pass_vec_all_runs,pMN_pass_vec_all_runs),]

saveRDS(window_calls_counts_normalised_amplicon_anno_all_run_sample_pass,"./data/processed/data_processed_r1_r2.rds")
saveRDS(list(MN = MN_pass_vec_all_runs,PN = PN_pass_vec_all_runs,pMN = pMN_pass_vec_all_runs),"./data/processed/samples_pass.rds")

meta_qc_ext_combined <- rbind(meta_qc_ext_r1,meta_qc_ext_r2)
saveRDS(meta_qc_ext_combined,"./data/processed/meta_r1_r2.rds")
saveRDS(index_MNsequences_on_ecDNA_r1,"./data/processed/index_MNsequences_on_ecDNA_r1.rds")
saveRDS(index_MNsequences_on_ecDNA_r2,"./data/processed/index_MNsequences_on_ecDNA_r2.rds")

