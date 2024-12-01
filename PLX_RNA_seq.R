# Install and load all libraries
library(AnnotationHub)
library(ensembldb)
library(tximport)
library(tidyverse)
library("DESeq2")
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(enrichR)
library(dplyr)
library("DEGreport")
library(org.Mm.eg.db)
library(DOSE)
library(clusterProfiler)
library(enrichplot)
library(mvtnorm)
library(rrcov)
library(stringr)
library("GeneOverlap")

ah <- AnnotationHub()
unique(ah$species) %>% View()
unique(ah$rdataclass) %>% View()

mouse_ens <- query(ah, c("Mus musculus", "EnsDb"))
# Extract annotations of interest (Latest version)
mouse_ens <- mouse_ens[["AH104895"]]
# Extract gene-level information
genes(mouse_ens, return.type = "data.frame") %>% View()
# Extract transcript-level information
transcripts(mouse_ens, return.type = "data.frame") %>% View()

# Create a transcript dataframe
txdb <- transcripts(mouse_ens, return.type = "data.frame") %>%
  dplyr::select(tx_id, gene_id)
txdb <- txdb[grep("ENSM", txdb$tx_id),]

# Create a gene-level dataframe
genedb <- genes(mouse_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, gene_name)

# Merge the two dataframes together
annotations <- inner_join(txdb, genedb)
annotations <- as.data.frame(annotations)
write.table(annotations,"Mouse_transcript_gene_anno.txt",quote = FALSE, row.names=FALSE, sep = "\t")
tx2gene <- read.delim("local directory") 
annotations

# Set working directory
setwd("local directory")
samples <- list.files(full.names = T, pattern="salmon")
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "./", "") %>% 
  str_replace("salmon_", "")
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("tx_id","gene_id")], countsFromAbundance="no", ignoreTxVersion = FALSE)
data <- txi$counts %>% 
  round() %>% 
  data.frame()
setwd("local directory")
write.table(data,"local directory/Tx_counts.txt",quote = FALSE, row.names=TRUE, sep = "\t")
saveRDS(txi,"txi")

# Load txi file and metadata to run DESeq2 
txi <- readRDS("txi")
metadata  <- read.csv("metadata.txt", header = TRUE, sep = "\t")
rownames(metadata) <- metadata$ID

# Running differential expression analysis
metadata$Sex = factor(metadata$Sex)
metadata$Treatment = factor(metadata$Treatment)

dds <-DESeqDataSetFromTximport(txi, colData = metadata, design = ~1)
dds <-DESeqDataSetFromMatrix(counts, colData = metadata, design = ~1)

dds$group <- factor(paste0(dds$Sex,"_",dds$Treatment))
design(dds) <- ~group
dds <-DESeq(dds)

# Get rid of the version .x at the end of gene name (Salmon)
rownames(dds) <- str_replace(rownames(dds) ,pattern = ".[0-9]+$", replacement = "")

# Perform log transformation
rld <- rlog(dds)
colnames(rld) <- rld$ID #Make sure the colnames is the ID

# Calculate euclidean distance between log-transformed counts across samples and save as a matrix
sampleDists <- dist(t(assay(rld)))
# SampleDists
sampleDistMatrix <- as.matrix(sampleDists)
# Plot differences using the pheatmap function. Change heatmap colors
pheatmap(sampleDistMatrix, col=colorRampPalette(c("blue", "white"))(100))


# Pull out relevant genes for enrichment, up and down-regulation gene lists
dev.new()
pdf("PCA.pdf")
p<-plotPCA(rld, intgroup = "group")
p + theme_bw()
  #geom_label(aes (label=metadata$ID)) +  ggtitle("PCA")
rld.sub <- rld[ , rld$Sex=="F"]
p<-plotPCA(rld.sub, intgroup = "group")
p + theme_bw() + geom_label(aes (label=metadata[metadata$Sex=="F","ID"]))
rld.sub <- rld[ , rld$Sex=="M"]
p<-plotPCA(rld.sub, intgroup = "group")
p + theme_bw() + geom_label(aes (label=metadata[metadata$Sex=="M","ID"]))
dev.off()
dev.new()
tiff("PCA.tiff", units="in",width=10, height =6, res=1000)
p<-plotPCA(rld, intgroup = "group")
p + theme_bw() +  geom_point(size = 5) + theme(axis.title=element_text(size=22), axis.text = element_text(size=18),legend.title = element_blank(), legend.text = element_text(size=14))
dev.off()

#to view results of F
results_F <-results(dds, contrast = c("group", "F_PLX", "F_Control"))
#to view results of M
results_M <-results(dds, contrast = c("group", "M_PLX", "M_Control"))
#to view results of FvsM control
results_F_M_Control <-results(dds, contrast = c("group", "F_Control", "M_Control"))
#to view results of FvsM PLX
results_F_M_PLX <-results(dds, contrast = c("group", "F_PLX", "M_PLX"))

F_sig_genes <-results_F[which(results_F$padj <0.05),]
F_sig_genes_up <-results_F[which(results_F$padj <0.05 & results_F$log2FoldChange>0),]
F_sig_genes_down <-results_F[which(results_F$padj <0.05 & results_F$log2FoldChange<0),]

M_sig_genes <-results_M[which(results_M$padj <0.05),]
M_sig_genes_up <-results_M[which(results_M$padj <0.05 & results_M$log2FoldChange>0),]
M_sig_genes_down <-results_M[which(results_M$padj <0.05 & results_M$log2FoldChange<0),]

F_M_Control_sig_genes <-results_F_M_Control[which(results_F_M_Control$padj <0.05),]
F_M_Control_sig_genes_up <-results_F_M_Control[which(results_F_M_Control$padj <0.05 & results_F_M_Control$log2FoldChange>0),]
F_M_Control_sig_genes_down <-results_F_M_Control[which(results_F_M_Control$padj <0.05 & results_F_M_Control$log2FoldChange<0),]

F_M_PLX_sig_genes <-results_F_M_PLX[which(results_F_M_PLX$padj <0.05),]
F_M_PLX_sig_genes_up <-results_M[which(results_F_M_PLX$padj <0.05 & results_F_M_PLX$log2FoldChange>0),]
F_M_PLX_sig_genes_down <-results_M[which(results_F_M_PLX$padj <0.05 & results_F_M_PLX$log2FoldChange<0),]

# Save norm counts 
norm_counts <-as.data.frame(counts(dds, normalized = TRUE))
write.table(norm_counts,"./Genelevel_deSeq2_normcounts.txt",quote = FALSE, sep = "\t", col.names=NA)

write.table(F_sig_genes,"./Genelevel_deSeq2_DEGs_F.txt",quote = FALSE, sep = "\t", col.names=NA)
write.table(M_sig_genes,"./Genelevel_deSeq2_DEGs_M.txt",quote = FALSE, sep = "\t", col.names=NA)
write.table(F_M_Control_sig_genes,"./Genelevel_deSeq2_DEGs_F_M_Control.txt",quote = FALSE, sep = "\t", col.names=NA)
write.table(F_M_PLX_sig_genes,"./Genelevel_deSeq2_DEGs_F_M_PLX.txt",quote = FALSE, sep = "\t", col.names=NA)


# Heatmap of Differentially expressed Genes
pheatmap(assay(rld)[rownames(F_sig_genes),], scale = "row")
# pheatmap(assay(rld)[rownames(C_vs_L_sig_genes),c(3,4,9,12,13,1,5,7,11,15,16)], scale = "row")
pdf("DEG_Heatmaps.pdf")
anno.df<-metadata[,2:3]
rownames(anno.df)<-metadata$ID
pheatmap(assay(rld)[rownames(F_sig_genes),rld$Sex=="F"], scale = "row",annotation_col  = anno.df, show_rownames=FALSE,main = "Female PLX vs Control")
pheatmap(assay(rld)[rownames(M_sig_genes),rld$Sex=="M"], scale = "row",annotation_col  = anno.df, show_rownames=FALSE,main = "Male PLX vs Control")
pheatmap(assay(rld)[rownames(F_M_Control_sig_genes),rld$Treatment=="Control"], scale = "row",annotation_col  = anno.df, show_rownames=TRUE,main = "Control Female vs Male")
pheatmap(assay(rld)[rownames(F_M_PLX_sig_genes),rld$Treatment=="PLX"], scale = "row",annotation_col  = anno.df, show_rownames=FALSE,main = "PLX Female vs Male")
dev.off()


# Enrichment analysis 
rld <- na.omit(rld)
all_genes <- as.character(rownames(rld))

# Enrichment analysis in Female
de_F_up <- rownames(F_sig_genes_up)
ego_F_up <- enrichGO(de_F_up,universe = all_genes,keyType = "ENSEMBL",OrgDb = "org.Mm.eg.db",ont="BP", readable=TRUE)
ego2_F_up <- simplify(ego_F_up, cutoff=0.7, by="p.adjust", select_fun=min)

de_F_up <- mapIds(org.Mm.eg.db, de_F_up, 'ENTREZID', 'ENSEMBL')
de_F_up <- de_F_up
kegg_F_up <- enrichKEGG(de_F_up,organism = "mmu", keyType = "kegg")

saveRDS(ego2_F_up, file = "enrichGO_res_F_up")
ego2_F_up<-readRDS("enrichGO_res_F_up")
write.csv(ego2_F_up@result,"enrichGO_res_F_up.csv")

de_F_down <- rownames(F_sig_genes_down)
ego_F_down <- enrichGO(de_F_down,universe = all_genes,keyType = "ENSEMBL",OrgDb = "org.Mm.eg.db",ont="BP", readable=TRUE)
ego2_F_down <- simplify(ego_F_down, cutoff=0.7, by="p.adjust", select_fun=min)

saveRDS(ego2_F_down, file = "enrichGO_res_F_down")
ego2_F_down<-readRDS("enrichGO_res_F_down")
write.csv(ego2_F_down@result,"enrichGO_res_F_down.csv")

tiff("Upregulated GO pathways_F.tiff",units = "in", width = 7, height = 7,res=1000)
dotplot(ego2_F_up, showCategory = ego2_F_up$Description[1:13][c(-2,-6,-12)]) + ggtitle("Upregulated GO pathways") +theme(plot.title = element_text(size = 16, face = "bold"), axis.text.y=element_text(size=15), axis.title.x =element_text(size=14))
dev.off()
tiff("Downregulated GO pathways_F.tiff",units = "in", width = 7, height = 7,res=1000)
dotplot(ego2_F_down, showCategory=10) + ggtitle("Downregulated GO pathways") +theme(plot.title = element_text(size = 16, face = "bold"),axis.text.y=element_text(size=15), axis.title.x =element_text(size=14))
dev.off()

pdf("GO_enrichment_F_figure.pdf",width = 8, height = 4)
barplot(ego2_F_up, order=T, x= "GeneRatio",showCategory=5,font.size=15) + ggtitle("GO Enrichment pathways_up DEGs")
barplot(ego2_F_down, order=T, x= "GeneRatio",showCategory=5,font.size=15) + ggtitle("GO Enrichment pathways_down DEGs")
dev.off()

# Enrichment analysis_Male
de_M_up <- rownames(M_sig_genes_up)
ego_M_up <- enrichGO(de_M_up,universe = all_genes,keyType = "ENSEMBL",OrgDb = "org.Mm.eg.db",ont="BP", readable=TRUE)
ego2_M_up <- simplify(ego_M_up, cutoff=0.7, by="p.adjust", select_fun=min)

saveRDS(ego2_M_up, file = "enrichGO_res_M_up")
ego2_M_up <- readRDS("enrichGO_res_M_up")
write.csv(ego2_M_up@result,"enrichGO_res_M_up.csv")

de_M_down <- rownames(M_sig_genes_down)
ego_M_down <- enrichGO(de_M_down,universe = all_genes,keyType = "ENSEMBL",OrgDb = "org.Mm.eg.db",ont="BP", readable=TRUE)
ego2_M_down <- simplify(ego_M_down, cutoff=0.7, by="p.adjust", select_fun=min)

saveRDS(ego2_M_down, file = "enrichGO_res_M_down")
ego2_M_down <- readRDS("enrichGO_res_M_down")
write.csv(ego2_M_down@result,"enrichGO_res_M_down.csv")

tiff("Upregulated GO pathways_M.tiff",units = "in", width = 7, height = 7,res=1000)
dotplot(ego2_M_up, showCategory = ego2_M_up$Description[1:11][-5]) + ggtitle("Upregulated GO pathways") +theme(plot.title = element_text(size = 16, face = "bold"),axis.text.y=element_text(size=15), axis.title.x =element_text(size=14))
dev.off()
tiff("Downregulated GO pathways_M.tiff",units = "in", width = 7, height = 7,res=1000)
dotplot(ego2_M_down, showCategory=10) + ggtitle("Downregulated GO pathways") +theme(plot.title = element_text(size = 16, face = "bold"),axis.text.y=element_text(size=15), axis.title.x =element_text(size=14))
dev.off()

pdf("GO_enrichment_M_figure.pdf",width = 8, height = 4)
barplot(ego2_M_up, order=T, x= "GeneRatio",showCategory=5,font.size=15) + ggtitle("GO Enrichment pathways_up DEGs")
barplot(ego2_M_down, order=T, x= "GeneRatio",showCategory=5,font.size=15) + ggtitle("GO Enrichment pathways_down DEGs")
dev.off()

# Enrichment analysis_Female vs Male PLX
de_F_M_PLX_up <- rownames(F_M_PLX_sig_genes_up)
ego_F_M_PLX_up <- enrichGO(de_F_M_PLX_up,universe = all_genes,keyType = "ENSEMBL",OrgDb = "org.Mm.eg.db",ont="BP", readable=TRUE)
ego2_F_M_PLX_up <- simplify(ego_F_M_PLX_up, cutoff=0.7, by="p.adjust", select_fun=min)

saveRDS(ego2_F_M_PLX_up, file = "enrichGO_res_F_M_PLX_up")
ego2_F_M_PLX_up <- readRDS("enrichGO_res_F_M_PLX_up")
write.csv(ego2_F_M_PLX_up@result,"enrichGO_res_F_M_PLX_up.csv")

de_F_M_PLX_down <- rownames(F_M_PLX_sig_genes_down)
ego_F_M_PLX_down <- enrichGO(de_F_M_PLX_down,universe = all_genes,keyType = "ENSEMBL",OrgDb = "org.Mm.eg.db",ont="BP", readable=TRUE)
ego2_F_M_PLX_down <- simplify(ego_F_M_PLX_down, cutoff=0.7, by="p.adjust", select_fun=min)

saveRDS(ego2_F_M_PLX_down, file = "enrichGO_res_F_M_PLX_down")
ego2_F_M_PLX_down <- readRDS("enrichGO_res_F_M_PLX_down")
write.csv(ego2_F_M_PLX_down@result,"enrichGO_res_F_M_PLX_down.csv")

tiff("Upregulated GO pathways_FvsM_PLX.tiff",units = "in", width = 7, height = 7,res=1000)
dotplot(ego2_F_M_PLX_up, showCategory = ego2_F_M_PLX_up$Description[c(1:3,5:8,13,18:19)]) + ggtitle("Upregulated GO pathways in Female on PLX") +theme(plot.title = element_text(size = 14, face = "bold"),axis.text.y=element_text(size=15), axis.title.x =element_text(size=14))
dev.off()
tiff("Downregulated GO pathways_FvsM_PLX.tiff",units = "in", width = 7, height = 7,res=1000)
dotplot(ego2_F_M_PLX_down, showCategory = ego2_F_M_PLX_down$Description[c(1:8,13:14)]) + ggtitle("Upregulated GO pathways in Male on PLX") +theme(plot.title = element_text(size = 14, face = "bold"),axis.text.y=element_text(size=15), axis.title.x =element_text(size=14))
dev.off()

pdf("GO_enrichment_F_M_PLX_figure.pdf",width = 8, height = 10)
dotplot(ego2_F_M_PLX_up,orderBy="GeneRatio", showCategory=20,font.size=14) + ggtitle("GO Enrichment pathways_up DEGs")
dotplot(ego2_F_M_PLX_down, showCategory=20,font.size=14) + ggtitle("GO Enrichment pathways_down DEGs")
dev.off()


# Enrichment analysis_Female vs Male Control
de_F_M_Control_up <- rownames(F_M_Control_sig_genes_up)
ego_F_M_Control_up <- enrichGO(de_F_M_Control_up,universe = all_genes,keyType = "ENSEMBL",OrgDb = "org.Mm.eg.db",ont="BP", readable=TRUE)
ego2_F_M_Control_up <- simplify(ego_F_M_Control_up, cutoff=0.7, by="p.adjust", select_fun=min)

saveRDS(ego2_F_M_Control_up, file = "enrichGO_res_F_M_Control_up")
write.csv(ego2_F_M_Control_up@result,"enrichGO_res_F_M_Control_up.csv")

de_F_M_Control_down <- rownames(F_M_Control_sig_genes_down)
ego_F_M_Control_down <- enrichGO(de_F_M_Control_down,universe = all_genes,keyType = "ENSEMBL",OrgDb = "org.Mm.eg.db",ont="BP", readable=TRUE)
ego2_F_M_Control_down <- simplify(ego_F_M_Control_down, cutoff=0.7, by="p.adjust", select_fun=min)

saveRDS(ego2_F_M_Control_down, file = "enrichGO_res_F_M_Control_down")
ego2_F_M_Control_down <- readRDS("enrichGO_res_F_M_Control_down")
write.csv(ego2_F_M_Control_down@result,"enrichGO_res_F_M_Control_down.csv")

tiff("Downregulated GO pathways_FvsM_Control.tiff",units = "in", width = 7, height = 7,res=1000)
dotplot(ego2_F_M_Control_down, showCategory = ego2_F_M_Control_down$Description[c(2,5:8,10:11,17,13,19)]) + ggtitle("Upregulated GO pathways in Male at Baseline") +theme(plot.title = element_text(size = 13, face = "bold"), axis.text.y=element_text(size=15), axis.title.x =element_text(size=14))
dev.off()


pdf("GO_enrichment_F_M_Control_figure.pdf",width = 8, height = 10)
dotplot(ego2_F_M_Control_down,  showCategory=20, font.size=14) + ggtitle("GO Enrichment pathways_down DEGs")
dev.off()

# Visualization 
# Female 
results_F <-results_F[order(results_F$padj),]
results_F <- results_F %>%
  data.frame() %>%
  tibble::rownames_to_column(var="GeneName")
ID2name <- mapIds(org.Mm.eg.db, results_F$GeneName, 'SYMBOL','ENSEMBL')
results_F$GeneName <- ID2name
results_F_tb <- results_F %>%
  as_tibble()
#Create TRUE/FALSE boleen for padj<0.05
results_F_tb <- results_F_tb %>% 
  mutate(threshold_OE = padj < 0.05)
results_F_tb$color <- rep("black",nrow(results_F_tb))
results_F_tb[which(results_F_tb$log2FoldChange < 0 & results_F_tb$padj < 0.05 ),"color"] <- "lightgreen"
results_F_tb[which(results_F_tb$log2FoldChange > 0 & results_F_tb$padj < 0.05 ),"color"] <- "red"
## Create an empty column to indicate which genes to label
results_F_tb <- results_F_tb %>% mutate(genelabels = "")
## Populate the genelabels column with contents of the genename column for the first 10 rows, i.e. the top 10 most significantly expressed genes
results_F_tb$genelabels[1:20] <- as.character(results_F_tb$GeneName[1:20])

# Male
results_M <-results_M[order(results_M$padj),]
results_M <- results_M %>%
  data.frame() %>%
  tibble::rownames_to_column(var="GeneName")
ID2name <- mapIds(org.Mm.eg.db, results_M$GeneName, 'SYMBOL','ENSEMBL')
results_M$GeneName <- ID2name
results_M_tb <- results_M %>%
  as_tibble()
#Create TRUE/FALSE boleen for padj<0.05
results_M_tb <- results_M_tb %>% 
  mutate(threshold_OE = padj < 0.05)
results_M_tb$color <- rep("black",nrow(results_M_tb))
results_M_tb[which(results_M_tb$log2FoldChange < 0 & results_M_tb$padj < 0.05 ),"color"] <- "lightgreen"
results_M_tb[which(results_M_tb$log2FoldChange > 0 & results_M_tb$padj < 0.05 ),"color"] <- "red"
## Create an empty column to indicate which genes to label
results_M_tb <- results_M_tb %>% mutate(genelabels = "")
## Populate the genelabels column with contents of the genename column for the first 10 rows, i.e. the top 10 most significantly expressed genes
results_M_tb$genelabels[1:20] <- as.character(results_M_tb$GeneName[1:20])

# Female vs Male PLX
results_F_M_PLX <-results_F_M_PLX[order(results_F_M_PLX$padj),]
results_F_M_PLX <- results_F_M_PLX %>%
  data.frame() %>%
  tibble::rownames_to_column(var="GeneName")
ID2name <- mapIds(org.Mm.eg.db, results_F_M_PLX$GeneName, 'SYMBOL','ENSEMBL')
results_F_M_PLX$GeneName <- ID2name
results_F_M_PLX_tb <- results_F_M_PLX %>%
  as_tibble()
#Create TRUE/FALSE boleen for padj<0.05
results_F_M_PLX_tb <- results_F_M_PLX_tb %>% 
  mutate(threshold_OE = padj < 0.05)
results_F_M_PLX_tb$color <- rep("black",nrow(results_F_M_PLX_tb))
results_F_M_PLX_tb[which(results_F_M_PLX_tb$log2FoldChange < 0 & results_F_M_PLX_tb$padj < 0.05 ),"color"] <- "lightgreen"
results_F_M_PLX_tb[which(results_F_M_PLX_tb$log2FoldChange > 0 & results_F_M_PLX_tb$padj < 0.05 ),"color"] <- "red"
## Create an empty column to indicate which genes to label
results_F_M_PLX_tb <- results_F_M_PLX_tb %>% mutate(genelabels = "")
## Populate the genelabels column with contents of the genename column for the first 10 rows, i.e. the top 10 most significantly expressed genes
results_F_M_PLX_tb$genelabels[1:20] <- as.character(results_F_M_PLX_tb$GeneName[1:20])

# Female vs Male Control
results_F_M_Control <-results_F_M_Control[order(results_F_M_Control$padj),]
results_F_M_Control <- results_F_M_Control %>%
  data.frame() %>%
  tibble::rownames_to_column(var="GeneName")
ID2name <- mapIds(org.Mm.eg.db, results_F_M_Control$GeneName, 'SYMBOL','ENSEMBL')
results_F_M_Control$GeneName <- ID2name
results_F_M_Control_tb <- results_F_M_Control %>%
  as_tibble()
#Create TRUE/FALSE boleen for padj<0.05
results_F_M_Control_tb <- results_F_M_Control_tb %>% 
  mutate(threshold_OE = padj < 0.05)
results_F_M_Control_tb$color <- rep("black",nrow(results_F_M_Control_tb))
results_F_M_Control_tb[which(results_F_M_Control_tb$log2FoldChange < 0 & results_F_M_Control_tb$padj < 0.05 ),"color"] <- "lightgreen"
results_F_M_Control_tb[which(results_F_M_Control_tb$log2FoldChange > 0 & results_F_M_Control_tb$padj < 0.05 ),"color"] <- "red"
## Create an empty column to indicate which genes to label
results_F_M_Control_tb <- results_F_M_Control_tb %>% mutate(genelabels = "")
## Populate the genelabels column with contents of the genename column for the first 10 rows, i.e. the top 10 most significantly expressed genes
results_F_M_Control_tb$genelabels[1:20] <- as.character(results_F_M_Control_tb$GeneName[1:20])

# Save volcano plots 
dev.new()
tiff("Female PLX vs. Control.tiff", units="in",width=6, height =6, res=1000)
ggplot(results_F_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = factor(color),size = factor(color, ordered = T)), shape=1)+
  geom_text_repel(aes(label = genelabels), max.overlaps = Inf) +
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", linewidth=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype="dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "gray")+
  scale_color_manual(breaks=c(),values=c("black","darkgreen","red")) +
  scale_size_manual(breaks=c(),values=c(1,2,2)) +
  ggtitle("Female PLX vs. Control") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()

dev.new()
tiff("Male PLX vs. Control.tiff", units="in",width=6, height =6, res=1000)
ggplot(results_M_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = factor(color),size = factor(color, ordered = T)), shape=1)+
  geom_text_repel(aes(label = genelabels), max.overlaps = Inf) +
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", linewidth=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype="dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "gray")+
  scale_color_manual(breaks=c(),values=c("black","darkgreen","red")) +
  scale_size_manual(breaks=c(),values=c(1,2,2)) +
  ggtitle("Male PLX vs Control") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()

dev.new()
tiff("Female vs. Male PLX.tiff", units="in",width=6, height =6, res=1000)
ggplot(results_F_M_PLX_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = factor(color),size = factor(color, ordered = T)), shape=1)+
  geom_text_repel(aes(label = genelabels), max.overlaps = Inf) +
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", linewidth=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype="dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "gray")+
  scale_color_manual(breaks=c(),values=c("black","darkgreen","red")) +
  scale_size_manual(breaks=c(),values=c(1,2,2)) +
  ggtitle("Female vs Male PLX") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()

dev.new()
tiff("Female vs. Male Control.tiff", units="in",width=6, height =6, res=1000)
ggplot(results_F_M_Control_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = factor(color),size = factor(color, ordered = T)), shape=1)+
  geom_text_repel(aes(label = genelabels), max.overlaps = Inf) +
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", linewidth=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype="dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "gray")+
  scale_color_manual(breaks=c(),values=c("black","darkgreen","red")) +
  scale_size_manual(breaks=c(),values=c(1,2,2)) +
  ggtitle("Female vs Male Control") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()