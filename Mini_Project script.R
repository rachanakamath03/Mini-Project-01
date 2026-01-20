#Load Libraries
library(readr)
library(dplyr)
library(ggplot2)
library(GEOquery)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap) 
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)

#Import RNA-sew count matrix 
counts <- read.csv(  
  "/Users/Desktop/Mini Project 1/DATASET/GSE86468_GEO.bulk.islet.processed.data.RSEM.raw.expected.counts.csv",
  row.names = 1,
  check.names = FALSE
)

dim(counts)
colnames(counts)[1:10]

#Retrieve sample Metadata from GEO
gse <- getGEO("GSE86468", GSEMatrix = TRUE)[[1]]
pheno <- pData(gse)
dim(pheno)
colnames(pheno)
table(pheno$'disease:ch1')

#Construct Sample Metadata Table
metadata <- data.frame(
  sample_id = pheno$title,
  condition = ifelse(
    pheno$'disease:ch1' == "Type 2 diabetic",
    "T2D",
    "Control"
  ),
  stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$sample_id
table(metadata$condition)

#Align RNA-seq counts with Metadata
all(colnames(counts) %in% rownames(metadata))
counts <- counts [, rownames(metadata)]
all(colnames(counts) == rownames(metadata))

#convert RSEM COUNTS ti integers
counts_round <- round(counts)

#Create DESeq2 object

dds <- DESeqDataSetFromMatrix(
  countData = counts_round,
  colData = metadata,
  design = ~ condition
)

#Filter Low-expression Genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

#Run Differential Expression Analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "T2D", "Control"))
summary(res)
head(res)

#Convert DESeq2 result into dataframe
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)

res_df <- res_df[, c("Gene",
                     "baseMean",
                     "log2FoldChange",
                     "lfcSE",
                     "stat",
                     "pvalue",
                     "padj")]

#Export DESeq2 results
write_csv(as.data.frame(res_df), "DESeq2_results.csv")



#Identify significant Differentially Expressed Genes
sig <- res_df[ which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
dim(sig)
#export DEG table
write.csv(as.data.frame(sig), "GSE86468_T2D_DEGs.csv")

                     
#Sepearte upregulated and downregulated genes
upregulated_genes <- subset(res_df,
                            padj < 0.05 & log2FoldChange > 1)
dim(upregulated_genes)
#Export upregulated genes
write.csv(upregulated_genes, "GSE86468_T2D_Upregulated.csv", row.names = FALSE)


downregulated_genes <- subset(res_df,
                              padj < 0.05 & log2FoldChange < -1)
dim(downregulated_genes)
#Export downregulated genes
write.csv(downregulated_genes, "GSE86468_T2D_DOWNREGULATED.csv", row.names = FALSE)


#MA and volcano plot
pdf("Type2D_vs_Control_MAplot.pdf", width = 7, height = 6)

plotMA(res, main = "Type 2 Diabetes vs Control - MA Plot")

dev.off()

#volcano plot
#remove NA FDR rows
res_df <- res_df[!is.na(res_df$padj), ]
sum(is.na(res_df$padj))
#define significance
res_df <- res_df %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Upregulated in T2D",
    padj < 0.05 & log2FoldChange < -1 ~ "Downregulated in T2D",
    TRUE ~ "Not Significant"
  ))

volcano_plot<- ggplot(res_df, aes(log2FoldChange, -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Upregulated in T2D" = "red",
    "Downregulated in T2D" = "blue",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  labs(
    title = "Volcano plot: Type 2 Diabetes vs Control",
    x = "log2 Fold Change (T2D / Control)",
    y = "-log10(padj)"
    ) +
  theme_minimal()
volcano_plot
#Export volcano plot
ggsave(
  filename = "Volcano Plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

#Enhanced volcano plot
EnhancedVolcano(res_df,
                lab = rownames(res_df),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.1,
                FCcutoff = 1,
                title = "Type 2 Diabetes vs Control",
                caption = "DESeq2, FDR < 0.1"
                )
ggsave("GSE86468_Volcano.png", width = 10, height = 8, dpi = 300)

#visualisation

vsd <- vst(dds, blind = FALSE)
norm_counts <- assay(vsd)

#Heatmap of significant differentially expressed genes 
sig_genes <- rownames(sig)
sig_norm_counts <- norm_counts[sig_genes, ]
#annotattiom
annotation_col <- data_frame(
  Condition = metadata$condition
)
rownames(annotation_col) <- rownames(metadata)
sig_norm_counts <- as.matrix(sig_norm_counts)
mode(sig_norm_counts) <- "numeric"

annotation_col <- as.data.frame(annotation_col)
rownames(annotation_col) <- rownames(metadata)
str(annotation_col)
pheatmap(sig_norm_counts,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
         main = "T2D vs Control – Significant DEGs")

#Principle component analysis(PCA)
pca <- pca <- prcomp(t(norm_counts))
percent_var <- (pca$sdev^2) / sum(pca$sdev^2) * 100

pca_data <- as.data.frame(pca$x)
pca_data$Condition <- metadata$condition

ggplot(pca_data, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3) +
  labs(
    title = "PCA of Human Pancreatic Islets (GSE86468)",
    x = sprintf("PC1 (%.1f%%)", percent_var[1]),
    y = sprintf("PC2 (%.1f%%)", percent_var[2])
  ) +
  theme_minimal()
dev.off()
ggsave("GSE86468_PCA.pdf", width = 7, height = 6)

#dispersion plot
plotDispEsts(dds)


#Sample-to-Sample Distance Heatmap
rld <- rlog(dds, blind = TRUE)

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = annotation_col,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Sample-to-Sample Distance Heatmap (GSE86468)")


#Functional Enrichment Analysis

#prep deg list
deg <- as.data.frame(sig)
deg$gene <- rownames(sig)
#map gene symbols to entrz id
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = deg$gene,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

deg$entrez <- entrez_ids
deg <- deg[!is.na(deg$entrez), ]
head(deg)

#KEGG pathway enrichment
kegg <- enrichKEGG(
  gene = deg$entrez,
  organism = "hsa",
  pvalueCutoff = 0.05
)

head(kegg)
dotplot(kegg, showCategory = 10) +
  ggtitle("KEGG pathways dysregulated in Type 2 Diabetes") +
  theme(
    axis.text.y = element_text(size = 8),   
    axis.text.x = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )
  

#GO Biological Process Enrichment 
go_bp <- enrichGO(
  gene = deg$entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE
)

head(go_bp)
dotplot(go_bp, showCategory = 10) +
  ggtitle("Biological Processes disrupted in T2D pancreatic islets") +
  theme(
    axis.text.y = element_text(size = 7),   # GO term labels
    axis.text.x = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )


