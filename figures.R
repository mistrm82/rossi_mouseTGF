load('geneSummaries.RData')

library(RColorBrewer)
library(pheatmap)
library(Biobase)

## Figure 1: Heatmap of genes at 1.5 FC

# Heatmap color palette and sample annotation
heatcolors.1 <- colorRampPalette(c("blue", "white", "red"))(6)
annotation <- data.frame(sampletype=pData(geneSummaries)[,'genotype'], row.names=colnames(geneSummaries))
ann_colors = list(sampletype = c(WT="black", KO="darkgrey"))

# Get expression data for significant genes
siggenes <- row.names(results)[which(results$threshold)]
sig_data <- exprs(geneSummaries)[siggenes,]

# Setup scaling
linear <- 2^sig_data
normalized <- function(x){ (x-min(x))/(max(x)-min(x))}
sig_data <- apply(linear, 1, function(x){ (x-min(x))/(max(x)-min(x))})

tiff('Figure1_GP.tiff', width=500, height=800)
pheatmap(t(sig_data), color = heatcolors.1, cluster_rows = T, border_color="black", show_rownames = F, 
         cluster_cols = T, show_colnames = T, clustering_distance_rows = "euclidean", annotation = annotation,
         clustering_distance_cols = "euclidean", annotation_colors =ann_colors, legend= FALSE, cellwidth =20, cellheight=1.5,
         fontsize = 10, fontsize_row = 10, height=40)
dev.off()

pdf('Figure1_GP.pdf', width=3, height=5, onefile=FALSE)
pheatmap(t(sig_data), color = heatcolors.1, cluster_rows = T, border_color="black", show_rownames = F, 
         cluster_cols = T, show_colnames = T, clustering_distance_rows = "euclidean", annotation = annotation,
         clustering_distance_cols = "euclidean", annotation_colors =ann_colors, legend= FALSE, 
         fontsize = 9, fontsize_row = 9, height=30)
dev.off()

postscript('Figure1_GP.ps', width=400, height =1000, onefile=FALSE,  horizontal=F)
pheatmap(t(sig_data), color = heatcolors.1, cluster_rows = T, border_color="black", show_rownames = F, 
         cluster_cols = T, show_colnames = T, clustering_distance_rows = "euclidean", annotation = annotation,
         clustering_distance_cols = "euclidean", annotation_colors =ann_colors, legend= FALSE, 
         fontsize = 9, fontsize_row = 9, height=30, cellwidth=16)
dev.off()


## Figure 3: TGFb1 target genes (down-regulated)

target.genes <- read.delim("TGF_targetgenes_all", header=T, sep="\t", row.names=1)

# Get expression data for those genes
select <- match(row.names(target.genes), row.names(geneSummaries))
expression <-exprs(geneSummaries)[select, ]
rownames(expression) <- fData(geneSummaries)$symbol[select]

linear <- 2^expression
normalized <- function(x){ (x-min(x))/(max(x)-min(x))}
norm.expr <- apply(linear, 1, function(x){ (x-min(x))/(max(x)-min(x))})
idx <- c(3,4,5,1,2,6)

# Plot heatmap
tiff('02_Figure_targetgenes.tiff', width=500, height=800)
pheatmap(t(norm.expr)[,idx], color = heatcolors.1, cluster_rows = F, border_color="black", show_rownames = T, 
         cluster_cols = F, show_colnames = T, clustering_distance_rows = "euclidean", annotation = annotation,
         clustering_distance_cols = "euclidean", annotation_colors =ann_colors, legend= FALSE, cellwidth =12,
         fontsize = 10, fontsize_row = 10, height=40)
dev.off()

## Figure 4: Heatmap of Butavsky et al Fig 8b genes

# Load the genelist
genelist <- read.delim("genelist.txt", header=T, sep="\t", row.names=1)
genelist <- genelist[which(genelist$X5.FoldChange == "X"),]

# Get expression data for those genes
select <- match(row.names(genelist), fData(geneSummaries)$symbol)
select <- select[which(!is.na(select))]
expression<-exprs(geneSummaries)[select, ]
rownames(expression) <- fData(geneSummaries)$symbol[select]

# Get trimmed genelist and associated results info
genelist <- genelist[match(row.names(expression), row.names(genelist)),]
genelist$fc <- results$logFC[match(row.names(genelist), results$symbol)]
genelist$sig <- results$threshold[match(row.names(genelist), results$symbol)]

# Generate comparison with Butavsky
comp.vec <- factor(levels=c("agree", "trend", "disagree"))
comp.vec[which(genelist$Butavsky == "+" & genelist$fc > 0 )] <- "trend"
comp.vec[which(genelist$Butavsky == "-" & genelist$fc < 0 )] <- "trend"
comp.vec[which(is.na(comp.vec))] <- "disagree"
genelist$compare <- comp.vec
genelist$compare[which(genelist$compare == "trend" & genelist$sig)] <- "agree"


# Setup for heatmap
ann_col <- annotation
ann_row <- data.frame(GeneFamily=genelist[,'cellFunction'], Concordance=genelist[,'compare'], row.names=row.names(genelist))
cbPalette.3 <- c(brewer.pal(11, "PuOr"), "deeppink3")
names(cbPalette.3) <- levels(genelist$cellFunction)

ann_colors = list(sampletype = c(WT="black", KO="darkgrey"), 
                  Concordance = c(agree="darkgreen", trend="forestgreen"),
                  GeneFamily = cbPalette.3)

### Optional other choices of palettes ###
cbPalette <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", 
               "#cab2d6", "#6a3d9a", "#f5deb3", "#b15928")
cbPalette.2 <- c("#660099", "#CC33FF", "#9966ff", "#C71585", "#FF6984", "#FF86C1",
                 "#ca2d6", "#f5deb3", "#ffff99", "#fdbf6f", "#ff7f00", "#b15928")


### Gene Pattern style scaling 
linear <- 2^expression
expression.2 <- apply(linear, 1, function(x){ (x-min(x))/(max(x)-min(x))})


tiff('03_Figure_Butavsky.tiff', width=900, height=1000)
pheatmap(t(expression.2)[order(genelist$Butavsky, genelist$cellFunction), idx], color = heatcolors.1, cluster_rows = F, 
         annotation_colors = ann_colors, legend = F, annotation_col = ann_col, annotation_row = ann_row, border_color="black", 
         cluster_cols = F, show_colnames = T, clustering_distance_cols = "euclidean", cellwidth=14,
         fontsize = 14, fontsize_row = 14, height=30)
dev.off()

pdf('03_Figure_Butavsky.pdf', width=12, height=12, onefile=FALSE)
pheatmap(t(expression.2)[order(genelist$Butavsky, genelist$cellFunction), idx], color = heatcolors.1, cluster_rows = F, 
         annotation_colors = ann_colors, legend = F, annotation_col = ann_col, annotation_row = ann_row, border_color="black", 
         cluster_cols = F, show_colnames = T, clustering_distance_cols = "euclidean", cellwidth=14,
         fontsize = 14, fontsize_row = 14, height=40)
dev.off()

postscript('03_Figure_Butavsky.ps', width=800, height=1400, onefile=FALSE,  horizontal=F)
pheatmap(t(expression.2)[order(genelist$Butavsky, genelist$cellFunction), idx], color = heatcolors.1, cluster_rows = F, 
         annotation_colors = ann_colors, legend = F, annotation_col = ann_col, annotation_row = ann_row, border_color="black", 
         cluster_cols = F, show_colnames = T, clustering_distance_cols = "euclidean", cellwidth=14,
         fontsize = 14, fontsize_row = 14, height=40)
dev.off()

## Figure 5: Heatmap of Microglia + GSEA sets

# Load the genelists
mg <- scan("microglia_genes.txt", what="character")
tgfb <- read.delim("additional data/HALLMARK_TGF_BETA_SIGNALING.txt", header=T, sep="\t", row.names=1)
g2m <- read.delim("additional data/HALLMARK_G2M_CHECKPOINT.txt", header=T, sep="\t", row.names=1)
interferon <- read.delim("additional data/HALLMARK_INTERFERON_RESPONSE.txt", header=T, sep="\t")

# Get expression
expression <- exprs(geneSummaries)
row.names(expression) <- fData(geneSummaries)$symbol
row.names(expression) <- toupper(row.names(expression))

# Figure A: MG + TGFbeta

# Get genes
# mg.caps <- toupper(mg)
gset <- unique(row.names(tgfb))

select <- match(gset, row.names(expression))
select <- select[which(!is.na(select))]

# Heatmap color palette and sample annotation
heatcolors.1 <- colorRampPalette(c("blue", "white", "red"))(6)
annotation <- data.frame(sampletype=pData(geneSummaries)[,'genotype'], row.names=colnames(geneSummaries))
ann_colors = list(sampletype = c(WT="black", KO="darkgrey"))

# Get KO order
toPlot <- exprs(geneSummaries)[select,]
row.names(toPlot) <- fData(geneSummaries)$symbol[select]

### Gene Pattern style scaling 
linear <- 2^toPlot
expression <- apply(linear, 1, function(x){ (x-min(x))/(max(x)-min(x))})
expression <- t(expression)
idx <- c(3,4,5,1,2,6)


tiff('05_Figure_TGF_GSEA.tiff', width=500, height=800)
pheatmap(expression[,idx], color = heatcolors.1, cluster_rows = F, annotation_colors = ann_colors,
         annotation=annotation, border_color="black", show_rownames = T, legend = FALSE, cellwidth=14,  
         cluster_cols = F, show_colnames = T, clustering_distance_cols = "euclidean", 
         fontsize = 10, fontsize_row = 10, height=80)
dev.off()


## Figure 6: Dendrogram

library(ggplot2)
library(ggdendro)

tiff('04_Figure_cluster_dendrogram.tiff')
ggplot(segment(ddata)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  theme_dendro() +
  geom_text(data=label(ddata), aes(x=x, y=y, label=label, hjust=-0.1), size=8) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none")

dev.off()

