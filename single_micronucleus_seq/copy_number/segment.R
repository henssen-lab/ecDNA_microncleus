####################
##segmentation single cell
####################
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(DNAcopy)
library(data.table)
library(lsa)
library(circlize)  
library(ComplexHeatmap)


#only run on one sample; ret
get_gc_vector <- function(input_path){
  # Import centromere and telomere annotation from UCSC
  session <- browserSession("UCSC")
  genome(session) <- "hg38"
  chr_sizes <- seqlengths(Hsapiens)
  telomere_table <- getTable(ucscTableQuery(session, "gap"))
  tel_ignore <- telomere_table[telomere_table$type %in% c("telomere"),]
  centromere_table <- getTable(ucscTableQuery(session, "cytoBand"))
  cent_ignore <- centromere_table[centromere_table$gieStain %in% c("acen"),]
  
  ignore_regions <- rbind(tel_ignore[,2:4],cent_ignore[,1:3])
  
  #load data
  per_bin_means_df <- read.csv(input_path, sep="\t", header = T)
  #subset autosomes + chrX+chrY
  autosomes <- paste0("chr",c(1:22))
  sex_chrom <- paste0("chr",c("X","Y"))
  
  per_bin_means_df <- per_bin_means_df[per_bin_means_df$chr %in% c(autosomes,sex_chrom),]
  #get GC_vector
  per_bin_means_df$gc <- NA
  #add gc perc per bin
  for (i in 1:dim(per_bin_means_df)[1]){
    print(paste0(i/dim(per_bin_means_df)[1]*100,"%"))
    current_bin <- per_bin_means_df[i,]
    binsize <- current_bin$end - current_bin$start
    #check if chr in chromsize
    if (current_bin$end <= chr_sizes[current_bin$chr]){per_bin_means_df[i,]$gc <- (letterFrequency(Views(Hsapiens, GRanges(seqnames = per_bin_means_df[i,]$chr, ranges = IRanges(start = per_bin_means_df[i,]$start + 1, end = per_bin_means_df[i,]$end +1 ))), "GC") / (binsize-letterFrequency(Views(Hsapiens, GRanges(seqnames = per_bin_means_df[i,]$chr, ranges = IRanges(start = per_bin_means_df[i,]$start + 1, end = per_bin_means_df[i,]$end +1 ))), "N")))}
  }
  
  gc_vec <- per_bin_means_df$gc
  return(gc_vec)
  
}


get_cn_profile <- function(path,gc_vec){#getsegmented
  
  per_bin_means_df <- read.csv(path, sep="\t", header = T)
  #subset autosomes + chrX+chrY
  autosomes <- paste0("chr",c(1:22))
  sex_chrom <- paste0("chr",c("X","Y"))
  
  per_bin_means_df <- per_bin_means_df[per_bin_means_df$chr %in% c(autosomes,sex_chrom),]
  per_bin_means_df$gc <- gc_vec
  #remove telomeric and centromeric regions
  ignore_regions_sample <- unique(queryHits(findOverlaps(GRanges(per_bin_means_df),GRanges(ignore_regions))))
  
  per_bin_means_df_subset_telcen_filtered <- per_bin_means_df[-ignore_regions_sample,] #subset for regions not centromeric or telomeric
  per_bin_means_df_subset_telcen_filtered <- per_bin_means_df_subset_telcen_filtered[per_bin_means_df_subset_telcen_filtered$log2_mean_cov != -Inf,] #subset for rows with counts
  per_bin_means_df_subset_telcen_filtered <- per_bin_means_df_subset_telcen_filtered[complete.cases(per_bin_means_df_subset_telcen_filtered), ] #remove rows with NA
  
  #GC curve correction by LOESS
  
  loess_model <- loess(log2_mean_cov ~ gc, data = per_bin_means_df_subset_telcen_filtered, span = 0.3)
  per_bin_means_df_subset_telcen_filtered$gc_bias <- predict(loess_model, per_bin_means_df_subset_telcen_filtered$gc)
  per_bin_means_df_subset_telcen_filtered$log2_mean_cov_corrected <- per_bin_means_df_subset_telcen_filtered$log2_mean_cov - per_bin_means_df_subset_telcen_filtered$gc_bias
  
  #median centering for corrected values (compute median only or autosomes)
  autosmome_median_corrected <- median(per_bin_means_df_subset_telcen_filtered[per_bin_means_df_subset_telcen_filtered$chr %in% autosomes,]$log2_mean_cov_corrected)
  per_bin_means_df_subset_telcen_filtered$log2_mean_cov_corrected_centered <- per_bin_means_df_subset_telcen_filtered$log2_mean_cov_corrected - autosmome_median_corrected
  #median centering for raw values (compute median only or autosomes)
  autosmome_median_raw <- median(per_bin_means_df_subset_telcen_filtered[per_bin_means_df_subset_telcen_filtered$chr %in% autosomes,]$log2_mean_cov)
  per_bin_means_df_subset_telcen_filtered$log2_mean_cov_raw_centered <- per_bin_means_df_subset_telcen_filtered$log2_mean_cov - autosmome_median_raw
  
  
  ##segmentation
  
  # Convert to CNA object
  cna_data <- CNA(
    genomdat = per_bin_means_df_subset_telcen_filtered$log2_mean_cov_corrected_centered,
    chrom = per_bin_means_df_subset_telcen_filtered$chr,
    maploc = per_bin_means_df_subset_telcen_filtered$start,
    data.type = "logratio"
  )
  
  # Smooth noisy data before segmentation
  cna_data <- smooth.CNA(cna_data)
  
  # Run CBS segmentation
  cna_segmented <- segment(cna_data, undo.splits = "sdundo", undo.SD = 1.5)
  
  # View results
  
  seg_df <- cna_segmented$output
  
  # Convert log2 data to data.table
  setDT(per_bin_means_df_subset_telcen_filtered)
  setDT(seg_df)
  
  # Ensure columns are named correctly for foverlaps() 
  setnames(seg_df, c("chrom","loc.start", "loc.end"), c("chr","start", "end"))  # Rename for clarity
  #setnames(seg_df, c("chrom"), c("chr"))  # Rename for clarity
  
  # Set keys for range-based matching (chromosome, start, end)
  setkey(seg_df, chr, start, end)
  setkey(per_bin_means_df_subset_telcen_filtered, chr, start, end)
  
  df_matched <- foverlaps(per_bin_means_df_subset_telcen_filtered, seg_df, by.x = c("chr", "start", "end"), by.y = c("chr", "start", "end"), nomatch = 0)
  
  
  df_matched$chr <- factor(df_matched$chr, levels = c(paste0("chr",1:22), "chrX", "chrY"))
  return(data.frame(df_matched))
  
}

##visualize copynumber profiles

plot_cn_profile <- function(CN_list,outdir_path){
  
  CN_list <- CN_profile_df_list
  for (sample in names(CN_list)){
    
    current_df <- CN_list[[sample]]
    
    pdf(paste0(outdir_path,"/",sample,"_CN_profile.pdf"), width = 10, height = 5)
    CN_plot <- ggplot() +
      geom_point(data = current_df, aes(x = i.start, y = log2_mean_cov_corrected_centered), color = "gray", alpha = 0.5, size = 0.5) +
      geom_segment(data = current_df, aes(x = start, xend = end, y = seg.mean, yend = seg.mean), 
                   color = "red", size = 1.2) +
      facet_wrap(~chr, scales = "free_x", nrow = 1) +
      theme_minimal() +
      labs(x = "Genomic Position", y = "Log2 Copy Ratio", title = sample)# + ylim(-10,10)
    print(CN_plot)
    dev.off()
  }
}

plot_cn_profile_anno <- function(CN_list,outdir_path){
  
  CN_list <- CN_profile_df_list
  for (sample in names(CN_list)){
    
    current_df <- CN_list[[sample]]
    
    if(sample %in% c(meta[meta$sample %in% c(qc_samples_vec),]$sample,"TR_14_bulk")){
      pdf(paste0(outdir_path,"/",meta[meta$sample == sample,]$class,"_",sample,"_CN_profile.pdf"), width = 10, height = 5)
      CN_plot <- ggplot() +
        geom_point(data = current_df, aes(x = i.start, y = log2_mean_cov_corrected_centered), color = "gray", alpha = 0.5, size = 0.5) +
        geom_segment(data = current_df, aes(x = start, xend = end, y = seg.mean, yend = seg.mean), 
                     color = "red", size = 1.2) +
        facet_wrap(~chr, scales = "free_x", nrow = 1) +
        theme_minimal() 
      labs(x = NULL, y = "Log2 Copy Ratio", title = sample)# + ylim(-10,10)
      print(CN_plot)
      dev.off()
      
      
    }
    
  }
}





#get gc vector once, to save memory and computation time
gc_vec <- get_gc_vector("/ecMN/copy_number/example_data/P3469_scMN_A1_S1_500kb.txt")


#initilize empty list for CN data frames
##make CN profiles and get segmentation per sample
CN_profile_df_list = list()
for(file in list.files(path = "/Users/robinxu/Documents/Projects/MNMseq/Copynumber_analysis/MNM_500kb/", full.names = F)){
  sample = paste0(strsplit(file, "_")[[1]][3],"_",strsplit(file, "_")[[1]][4])
  path = paste0("/Users/robinxu/Documents/Projects/MNMseq/Copynumber_analysis/MNM_500kb/",file)
  CN_profile <- get_cn_profile(path, gc_vec)
  CN_profile$region_ID <- paste0(CN_profile$chr,"_",CN_profile$i.start)
  CN_profile_df_list[[sample]] <- CN_profile
}

plot_cn_profile(CN_profile_df_list,"/Users/robinxu/Downloads/untitledfolder2")



