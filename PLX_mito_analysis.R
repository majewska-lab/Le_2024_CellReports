# Download MitoCarta3.0 database to use for analysis
library("GeneOverlap")
setwd("local directory")
ego2_F_M_PLX_up <- readRDS("enrichGO_res_F_M_PLX_up")
ego2_F_M_PLX_down <- readRDS("enrichGO_res_F_M_PLX_down")
mouse_mitogenes <- read.csv("mouse_mitogenes.csv",header = TRUE)
mouse_mitolist <- mouse_mitogenes$EnsemblGeneID

# Generate upregulated and downregulated gene lists 
allgenes_GO_up <- ego2_F_M_PLX_up$geneID
allgenes_GO_up <- unlist(strsplit(allgenes_GO_up, "/"))
allgenes_GO_up <- unique(allgenes_GO_up)
allgenes_GO_up <- mapIds(org.Mm.eg.db, allgenes_GO_up , 'ENSEMBL','SYMBOL')
allgenes_GO_up <- intersect(allgenes_GO_up, mouse_mitolist)
allgenes_GO_up <- mouse_mitogenes[mouse_mitogenes$EnsemblGeneID %in% allgenes_GO_up,]
allgenes_GO_up <- allgenes_GO_up$Symbol

allgenes_GO_down <- ego2_F_M_PLX_down$geneID
allgenes_GO_down <- unlist(strsplit(allgenes_GO_down, "/"))
allgenes_GO_down <- unique(allgenes_GO_down)
allgenes_GO_down <- mapIds(org.Mm.eg.db, allgenes_GO_down , 'ENSEMBL','SYMBOL')
allgenes_GO_down <- intersect(allgenes_GO_down, mouse_mitolist)
allgenes_GO_down <- mouse_mitogenes[mouse_mitogenes$EnsemblGeneID %in% allgenes_GO_down,]
allgenes_GO_down <- allgenes_GO_down$Symbol

# Make list of mitochondrial pathways 
mito_pathways <- read.csv("mito_pathways.csv", header=TRUE)
mito_pathways <- mito_pathways[,c(1,3)]
mito_pathways_names <- mito_pathways[,1]
mito_pathways_genes <- mito_pathways[,-1]
mito_pathways_genes <- strsplit(mito_pathways_genes,", ")
names(mito_pathways_genes) <- mito_pathways_names

my_list <- list(allgenes_GO_up,allgenes_GO_down)
names(my_list) <- c("ego2_F_M_PLX_up" , "ego2_F_M_PLX_down")

# Run GeneOverlap package
gom.obj <- newGOM(my_list, mito_pathways_genes, genome.size=NULL,spec=c('mm9.gene'))
pval_matrix <- getMatrix(gom.obj,name="pval")
pval_vec <- as.vector(pval_matrix)
pvaladj_matrix <- matrix(p.adjust(pval_vec,method = "BH"),nrow=2)
rownames(pvaladj_matrix) <- rownames(pval_matrix)
colnames(pvaladj_matrix) <- colnames(pval_matrix)
pvaladj_matrix <- t(pvaladj_matrix)
write.csv(pvaladj_matrix,"./mito_pathway_enrichment.csv",row.names = TRUE)


allgenes_GO_up <- as.data.frame(allgenes_GO_up)
allgenes_GO_up <- allgenes_GO_up[,1]
all_mito_GO_genes <- list(allgenes_GO_up,allgenes_GO_down)
names(all_mito_GO_genes) <- c("allgenes_GO_up","allgenes_GO_down")
all_mito_GO_genes <- plyr::ldply(all_mito_GO_genes, rbind)
write.csv(t(all_mito_GO_genes), "./all_mito_DEGs.csv",row.names = FALSE)