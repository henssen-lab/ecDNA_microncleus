library(DESeq2)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
#load DM data

in_dir <- "htseq_counts/"
files <- grep("DM",list.files(in_dir),value=TRUE)
sampleCondition <- c("HU","HU","HU","UT","UT","UT")
sampleName <- c("DM_HU_1","DM_HU_2","DM_HU_3","DM_UT_1","DM_UT_2","DM_UT_3")
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = files,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)



ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = in_dir,
                                       design= ~ condition)

dds <- DESeq(ddsHTSeq)
res <- results(dds, contrast = c("condition","HU","UT"))
res


sizeFactors(dds)
colSums(counts(dds))
colSums(counts(dds, normalized=T))
summary(res)


#match ENSEMBL Id with Gene Symbol and ENTREZ
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
#get diff MYC
res[res$symbol %in% c("MYC"),]



MYC <- plotCounts(dds, gene="ENSG00000136997.22", intgroup="condition", 
                 returnData=TRUE,normalized = T)


MYC$gene <- "MYC"

MYC$log2 <- log2(MYC$count + 0.1)

MYC$log2scaled <- scale(MYC$log2)


counts_to_plot <- rbind(MYC)

counts_to_plot$scaled_counts <- scale(counts_to_plot$log2)


pdf("./plots/RNA_zscore_counts.pdf", width =5 , height = 4)

ggplot(data=counts_to_plot, aes(x=factor(gene, level=c("MYC")), y=log2scaled, color=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge()) +
  facet_wrap(~gene, scale="free_x") + theme_classic() 
dev.off()




## select significant genes
res_filtered <- subset(res, pvalue<.05)
res_filtered <- subset(res_filtered, symbol != "NA")
res_filtered <- subset(res_filtered, entrez != "NA")


write.table(res_filtered, "./tables/DE_HU_vs_UT_sign.txt", sep = "\t", quote = F)

sign_genes_vec <-  (res_filtered$entrez)
FC_vec <- res_filtered$log2FoldChange

names(FC_vec) <-sign_genes_vec
gene_list = sort(FC_vec, decreasing = TRUE) #rank


#geneset enrichment

gene_sets = msigdbr(category = "H") #load all Hallmark genesets

map = gene_sets[, c("gs_name", "entrez_gene")] 
map$entrez_gene = as.character(map$entrez_gene)

geneset_enrichment_DM = GSEA(geneList = gene_list, TERM2GENE = map,nPermSimple = 10000)
GO <- gseGO(
  geneList = gene_list,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH", minGSSize = 3, maxGSSize = 800
)


require(DOSE)
pdf("./plots/GSE_dot.pdf", width = 8 , height = 8)
dotplot(geneset_enrichment_DM, showCategory=10, split=".sign") 
dev.off()
pdf("./plots/GSE_dot2.pdf", width = 15 , height = 8)
dotplot(geneset_enrichment_DM, showCategory=15, split=".sign") + facet_grid(.~.sign)
dev.off()

pdf("./plots/GSEA_Myc.pdf", width = 7 , height = 5)
gseaplot2(geneset_enrichment_DM, geneSetID = 1, title = geneset_enrichment_DM$Description[1], pvalue_table = F)
dev.off()
pdf("./plots/GSEA_Myc_pval.pdf", width = 12 , height = 5)
gseaplot2(geneset_enrichment_DM, geneSetID = 1, title = geneset_enrichment_DM$Description[1], pvalue_table = T)
dev.off()


sink("sessionInfo.txt")
sessionInfo()
sink()

