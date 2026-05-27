library(GenomicRanges)
library(regioneR)
library(plyranges)
library(ggplot2)
library(plyr)
library(ggridges)
library(stats)
library(ggpubr)


setwd("./ecDNA_microncleus/nanopore")
###24h_HU

#load CHP212 amplicon reconstruction and annotations

circle_regions <- toGRanges("./data/reconstruct.bed")

H3K27ac_regions <- toGRanges("./data/GSM2664355_CHP212.K27ac.rep3_peaks.liftoverhg38.narrowPeak")

#load bed files
left_fork <- read.csv("./data/leftForks_DNAscent_forkSense_stressSignatures_24h_HU.bed", header = F, sep = " ", dec = ".")
right_fork <- read.csv("./data/rightForks_DNAscent_forkSense_stressSignatures_24h_HU.bed", header = F, sep = " ", dec = ".")
forks_24h_HU <- rbind(left_fork,right_fork)


#convert Bed into Grange
column_names_forks <- c("chr","start","end","read_ID","read_start", "read_end", "strand_fr", "fork_length","EdU_length","BrdU_length","1","2","3","4","stall_score")
colnames(forks_24h_HU) <- column_names_forks
forks_24h_HU_grange <- GenomicRanges::makeGRangesFromDataFrame(forks_24h_HU, keep.extra.columns = TRUE)

#split into linear and circular DNA
circle_forks_24h_HU <- subsetByOverlaps(forks_24h_HU_grange,circle_regions)
linear_forks_24h_HU <- subsetByOverlaps(forks_24h_HU_grange,circle_regions, invert = T)
linear_forks_24h_HU_27ac <- subsetByOverlaps(linear_forks_24h_HU,H3K27ac_regions)#subset for active regions

##load control
left_fork_c <- read.csv("./data/leftForks_DNAscent_forkSense_stressSignatures_control.bed", header = F, sep = " ", dec = ".")
right_fork_c <- read.csv("./data/rightForks_DNAscent_forkSense_stressSignatures_control.bed", header = F, sep = " ", dec = ".")
forks_control <- rbind(left_fork_c,right_fork_c)

#convert Bed into Grange
column_names_forks <- c("chr","start","end","read_ID","read_start", "read_end", "strand_fr", "fork_length","EdU_length","BrdU_length","1","2","3","4","stall_score")
colnames(forks_control) <- column_names_forks
forks_control_grange <- GenomicRanges::makeGRangesFromDataFrame(forks_control, keep.extra.columns = TRUE)

#split into linear and circular DNA
circle_forks_control <- subsetByOverlaps(forks_control_grange,circle_regions)
linear_forks_control <- subsetByOverlaps(forks_control_grange,circle_regions, invert = T)
linear_forks_control_27ac <- subsetByOverlaps(linear_forks_control,H3K27ac_regions)#subset for active regions


###analysis and plotting


##analysis untreated
control_speed_lin_h3k27ac <- linear_forks_control_27ac[(elementMetadata(linear_forks_control_27ac)$stall_score >= -1)]$fork_length / 15
control_speed_circ <- circle_forks_control[(elementMetadata(circle_forks_control)$stall_score >= -1)]$fork_length / 15
fs_control_lin_27ac <- data.frame(fork_speed = control_speed_lin_h3k27ac)
fs_control_lin_27ac$treatment <- "control"
fs_control_lin_27ac$type <- "linear"
fs_control_circ_27ac <- data.frame(fork_speed = control_speed_circ)
fs_control_circ_27ac$treatment <- "control"
fs_control_circ_27ac$type <- "circular"


##analysis 24h HU
HU24_speed_lin_27ac <- linear_forks_24h_HU_27ac[(elementMetadata(linear_forks_24h_HU_27ac)$stall_score >= -1)]$fork_length / 16.4
HU24_speed_circ <- circle_forks_24h_HU[(elementMetadata(circle_forks_24h_HU)$stall_score >= -1)]$fork_length / 16.4
fs_HU24_lin_27ac <- data.frame(fork_speed = HU24_speed_lin_27ac)
fs_HU24_lin_27ac$treatment <- "HU_24h"
fs_HU24_lin_27ac$type <- "linear"
fs_HU24_circ_27ac <- data.frame(fork_speed = HU24_speed_circ)
fs_HU24_circ_27ac$treatment <- "HU_24h"
fs_HU24_circ_27ac$type <- "circular"


#data_forks <- rbind(fs_control_lin,fs_control_circ,fs_HU24_lin,fs_HU24_circ,fs_pulse_lin,fs_pulse_circ)
data_forks_27ac <- rbind(fs_control_lin_27ac,fs_control_circ_27ac,fs_HU24_lin_27ac,fs_HU24_circ_27ac)

data_forks_27ac$combo <- paste0(data_forks_27ac$treatment, "_",data_forks_27ac$type)

data_forks_27ac$combo <- factor(data_forks_27ac$combo, levels = c("control_linear", "HU_24h_linear", "control_circular","HU_24h_circular"))

comparisons <- combn(levels(data_forks_27ac$combo), 2, simplify = FALSE)

ggplot(data_forks_27ac, aes(x = combo, y = fork_speed, colour = combo, group = combo)) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.7)) +  # Separate points with dodge
  # scale_color_okabeito(order=7:9)+
  theme_classic()+
  geom_violin(trim = FALSE, alpha = 0.5) +
  stat_compare_means(comparisons = comparisons,
                     method = "t.test",
                     var.equal = FALSE,      # Welch's t-test
                     label = "p.signif") +
  theme_classic() +
  labs(
    x = NULL,
    y = "log2 (circular / linear DNA mean(CPM))",
    colour = NULL
  ) +
  theme(legend.position = "none") +
  #theme(
  #  legend.position = "none",
  #  axis.title.x = element_blank(),
  #  axis.text.x = element_blank(),
  #  axis.title.y = element_blank(),
  #  axis.text.y = element_blank()
  #)+ 
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "black", fatten = 2)
ggsave("./plots/forks_circ_linear_stripped_final_pulse.eps", width = 1.4,
       height = 1.4, device= "eps", dpi=450)




##compute differential HU effect as done in Jaworski et al 2025

#computes logFC for median fork speeds 
get_logFC <- function(df, n = 30) {
  ctrl <- sample(df$fork_speed[df$treatment == "control"], n, replace = FALSE)
  hu   <- sample(df$fork_speed[df$treatment == "HU_24h"],   n, replace = FALSE)
  
  log(median(hu) / median(ctrl))
}


set.seed(43)
n_iter <- 1000
logFC_ecDNA   <- replicate(n_iter, get_logFC(subset(data_forks_27ac, type=="circular"), n=30))
logFC_genomic <- replicate(n_iter, get_logFC(subset(data_forks_27ac, type=="linear"), n=30))

obs_diff <- mean(logFC_ecDNA) - mean(logFC_genomic)
obs_diff
1-exp(obs_diff)

#permutation test
#test mean diff as above
combined <- c(logFC_ecDNA, logFC_genomic)
labels   <- c(rep("circular", n_iter), rep("linear", n_iter))

n_perm <- 10000
perm_diff <- numeric(n_perm)

set.seed(43)
for (i in 1:n_perm) {
  perm_labels <- sample(labels)  # shuffle labels
  perm_diff[i] <- mean(combined[perm_labels=="circular"]) - mean(combined[perm_labels=="linear"])
}

# p-val(one sided)
p_val <- mean(perm_diff <= obs_diff)
p_val
1/10000


##plot
library(tibble)
library(ggplot2)

resampled_FC_df <- tibble(
  logFC = c(logFC_ecDNA, logFC_genomic),
  type  = rep(c("ecDNA", "genomic"), each = length(logFC_ecDNA))
)




ggplot(resampled_FC_df, aes(x = type, y = logFC, colour = type, group = type)) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.7)) +  # Separate points with dodge
  # scale_color_okabeito(order=7:9)+
  theme_classic() +
  labs(
    x = NULL,
    y = "log2 (circular / linear DNA mean(CPM))",
    colour = NULL
  ) +
  theme(legend.position = "none") +
  #theme(
  #  legend.position = "none",
  #  axis.title.x = element_blank(),
  #  axis.text.x = element_blank(),
  #  axis.title.y = element_blank(),
  #  axis.text.y = element_blank()
  #)+ 
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "black", fatten = 2)


ggplot(resampled_FC_df, aes(x = type, y = logFC, colour = type, group = type)) +
  geom_violin(trim = FALSE, alpha = 0.4) +  # <- this draws the violin
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "black", fatten = 2) +
  theme_classic() +
  labs(
    x = NULL,
    y = "logFC",
    colour = NULL
  ) +
  theme(legend.position = "none")
ggsave("./plots/subsampled_FC.eps", width = 1.4,
       height = 1.4, device= "eps", dpi=450)


ggplot(data.frame(perm_diff = perm_diff), aes(x = perm_diff)) +
  geom_density(fill = "white", alpha = 0.4) +  # density, not counts
  geom_vline(xintercept = obs_diff, color = "red", linetype = "dashed") +
  theme_classic() +
  labs(x = "Mean difference (ecDNA - genomic)", y = "Density") + xlim(-0.2,0.05)
ggsave("./plots/subsampled_FC_permutation.eps", width = 1.4,
       height = 1.4, device= "eps", dpi=450)

