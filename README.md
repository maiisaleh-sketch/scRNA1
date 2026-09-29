# sc# =========================================================
# Full scRNA-seq Pipeline — Seurat-based
# =========================================================

# --- Clear environment ---
rm(list = ls())
gc(full = TRUE)

library(data.table)
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)

# --- 1. Load raw count data ---
counts_dt <- fread("G:/scRNA/Project/GSE149614_HCC.scRNAseq.S71915.count.txt.gz")
 

gene_names <- counts_dt[[1]]
counts_dt[, V1 := NULL]

# --- 2. Convert to sparse matrix in chunks (memory-safe) ---
n <- ncol(counts_dt)
chunk_size <- 5000
n_chunks <- ceiling(n / chunk_size)
sparse_list <- vector("list", n_chunks)

for (i in seq_len(n_chunks)) {
  start <- (i - 1) * chunk_size + 1
  end <- min(i * chunk_size, n)
  chunk_mat <- as.matrix(counts_dt[, start:end, with = FALSE])
  sparse_list[[i]] <- Matrix(chunk_mat, sparse = TRUE)
  rm(chunk_mat)
  gc()
  cat("Processed columns", start, "to", end, "of", n, "\n")
}

counts_sparse <- do.call(cbind, sparse_list)
rownames(counts_sparse) <- gene_names
dim(counts_sparse)   # should be 25712 x 71915

rm(counts_dt, sparse_list)
gc(full = TRUE)

# =========================================================
# --- 3. Create Seurat object (do this immediately) ---
# min.cells = 3 replaces your manual gene-level filtering step
# min.features = 200 replaces part of your manual cell filtering
# =========================================================
seurat_obj <- CreateSeuratObject(
  counts = counts_sparse,
  project = "HCC01T",
  min.cells = 3,
  min.features = 200
)

seurat_obj
rm(counts_sparse)
gc(full = TRUE)

# --- 4. Compute % mitochondrial reads (Seurat built-in) ---
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

head(seurat_obj@meta.data)

# --- 5. Pre-QC visualization ---
VlnPlot(
  seurat_obj,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
####
# violin plot

VlnPlot(seurat_obj, features = 
          c("nFeature_RNA","nCount_RNA","percent.mt"),
        ncol = 3)

#Scatter Plot

FeatureScatter(seurat_obj, feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")


FeatureScatter(seurat_obj, feature1 = "nCount_RNA",
               feature2 = "percent.mt")


FeatureScatter(seurat_obj, feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")


# --- 6. Set thresholds and filter ---
min_counts   <- 500
max_mt_pct   <- 10

seurat_obj$pass_qc <- seurat_obj$nCount_RNA >= min_counts &
  seurat_obj$percent.mt <= max_mt_pct

table(seurat_obj$pass_qc)

seurat_obj_filtered <- subset(
  seurat_obj,
  subset = nCount_RNA >= min_counts & percent.mt <= max_mt_pct
)

seurat_obj_filtered

# --- 7. Post-QC visualization ---
VlnPlot(
  seurat_obj_filtered,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

# --- 8. Combined before/after QC figure ---
seurat_obj$stage <- "Before QC"
seurat_obj_filtered$stage <- "After QC"

merged_meta <- rbind(
  seurat_obj@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "stage")],
  seurat_obj_filtered@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "stage")]
)
merged_meta$stage <- factor(merged_meta$stage, levels = c("Before QC", "After QC"))

merged_long <- tidyr::pivot_longer(
  merged_meta,
  cols = c(nCount_RNA, nFeature_RNA, percent.mt),
  names_to = "metric", values_to = "value"
)

p1 <- ggplot(merged_long, aes(x = stage, y = value, fill = stage)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2, alpha = 0.5) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal() +
  scale_fill_manual(values = c("Before QC" = "grey70", "After QC" = "steelblue")) +
  labs(title = "QC Metrics: Before vs After Filtering", x = NULL, y = NULL) +
  theme(legend.position = "none")

p1

# =========================================================
# --- 9. Doublet detection (scDblFinder on the filtered Seurat object) ---
# =========================================================
library(scDblFinder)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(seurat_obj_filtered)
set.seed(123)
sce <- scDblFinder(sce)

table(sce$scDblFinder.class)

# Add doublet info back onto the Seurat object
seurat_obj_filtered$doublet_score <- sce$scDblFinder.score
seurat_obj_filtered$doublet_class <- sce$scDblFinder.class

VlnPlot(
  seurat_obj_filtered,
  features = "doublet_score",
  group.by = "doublet_class"
)

# --- 10. Remove doublets ---
seurat_obj_singlets <- subset(seurat_obj_filtered, subset = doublet_class == "singlet")
seurat_obj_singlets


# Are the cell barcodes in the same order?
identical(colnames(seurat_obj_filtered), colnames(sce))

# If FALSE, this is your bug — fix with:
sce <- sce[, colnames(seurat_obj_filtered)]  # reorder sce to match

# Then redo the assignment
seurat_obj_filtered$doublet_score <- sce$scDblFinder.score
seurat_obj_filtered$doublet_class <- sce$scDblFinder.class

# Re-verify the pattern makes sense
tapply(seurat_obj_filtered$doublet_score, seurat_obj_filtered$doublet_class, summary)

VlnPlot(
  seurat_obj_filtered,
  features = "doublet_score",
  group.by = "doublet_class"
)


# =========================================================
# --- 11. Normalization, HVG, scaling, PCA ---
# =========================================================
seurat_obj_singlets <- NormalizeData(seurat_obj_singlets, normalization.method = "LogNormalize", scale.factor = 1e4)
seurat_obj_singlets <- FindVariableFeatures(seurat_obj_singlets, selection.method = "vst", nfeatures = 2000)

# --- Visualize Variable Features ---

# Basic plot: all genes, HVGs highlighted in red
plot1 <- VariableFeaturePlot(seurat_obj_singlets)
plot1
dev.off()  # close any stuck graphics device
plot1      # re-print
# --- With top N genes labeled (very useful for identifying key drivers) ---
top10 <- head(VariableFeatures(seurat_obj_singlets), 10)

plot2 <- LabelPoints(
  plot = plot1,
  points = top10,
  repel = TRUE,
  xnudge = 0,
  ynudge = 0
)
plot2


seurat_obj_singlets <- ScaleData(seurat_obj_singlets, features = VariableFeatures(seurat_obj_singlets))
seurat_obj_singlets <- RunPCA(seurat_obj_singlets, features = VariableFeatures(seurat_obj_singlets), npcs = 50)

# --- 1. Basic PCA scatter (PC1 vs PC2) ---
DimPlot(seurat_obj_singlets, reduction = "pca") +
  ggtitle("PCA: PC1 vs PC2")

# --- 2. Color by a QC metric to check for confounding structure ---
FeaturePlot(seurat_obj_singlets, features = "percent.mt", reduction = "pca")
FeaturePlot(seurat_obj_singlets, features = "nCount_RNA", reduction = "pca")

# --- 3. Top genes driving each PC (loadings) ---
VizDimLoadings(seurat_obj_singlets, dims = 1:4, reduction = "pca")

# --- 4. Heatmap of top genes/cells for individual PCs ---
# Useful for judging whether a PC reflects real biological signal vs noise
DimHeatmap(seurat_obj_singlets, dims = 1, cells = 500, balanced = TRUE)

# Multiple PCs at once
DimHeatmap(seurat_obj_singlets, dims = 1:15, cells = 500, balanced = TRUE)


ElbowPlot(seurat_obj_singlets, ndims = 50)

# --- 12. Clustering and UMAP ---
seurat_obj_singlets <- FindNeighbors(seurat_obj_singlets, dims = 1:30)
seurat_obj_singlets <- FindClusters(seurat_obj_singlets, resolution = 0.5)
seurat_obj_singlets <- RunUMAP(seurat_obj_singlets, dims = 1:30)

DimPlot(seurat_obj_singlets, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP: Clusters")

# --- 1. Find markers for every cluster ---XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# This can take a while at ~67K cells — consider min.pct/logfc thresholds to speed it up
all_markers <- FindAllMarkers(
  seurat_obj_singlets,
  only.pos = TRUE,        # only genes upregulated in each cluster
  min.pct = 0.25,          # expressed in at least 25% of cells in the cluster
  logfc.threshold = 0.25
)
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

seurat_obj_singlets <- FindClusters(seurat_obj_singlets, resolution = 0.2)
table(Idents(seurat_obj_singlets))  # check how many clusters this gives you

seurat_obj_singlets <- RunUMAP(seurat_obj_singlets, dims = 1:30)
DimPlot(seurat_obj_singlets, reduction = "umap", label = TRUE)

####################

#### save a checkpoint immediately after clustering,
###before attempting FindAllMarkers:

saveRDS(seurat_obj_singlets, "seurat_obj_singlets.rds")

# Top 5 markers per cluster, ranked by fold change
library(dplyr)


# --- Find markers for every cluster (24 clusters now, should be faster) ---
all_markers <- FindAllMarkers(
  seurat_obj_singlets,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 0.5
)

# Then:
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

print(top_markers, n = 100)

# View all 120 rows to survey every cluster's top markers
print(top_markers, n = 120)

# Or as a scrollable data viewer (easier to read at this size)
View(top_markers)

liver_immune_markers <- c(
  "ALB", "APOA2", "APOC3", "TTR",
  "PECAM1", "VWF",
  "CD3D", "CD3E", "CD8A",
  "CD79A", "MS4A1",
  "CD68", "CD163", "LYZ",
  "NKG7", "GNLY",
  "COL1A1", "COL3A1",
  "EPCAM", "KRT19"
)
liver_immune_markers <- liver_immune_markers[liver_immune_markers %in% rownames(seurat_obj_singlets)]

DotPlot(seurat_obj_singlets, features = liver_immune_markers) + RotatedAxis()

# Try re-printing or forcing the plot to render
p <- DotPlot(seurat_obj_singlets, features = liver_immune_markers) +
  RotatedAxis()
p


# Cluster by cluster, easier to read and paste
for (cl in levels(Idents(seurat_obj_singlets))) {
  cat("\n=== Cluster", cl, "===\n")
  print(top_markers[top_markers$cluster == cl, c("gene", "avg_log2FC", "pct.1", "pct.2")])
}
View(top_markers)
# --- Full cluster annotation based on marker analysis ---

cluster_labels <- c(
  "0"  = "NK/Cytotoxic cells",
  "1"  = "TAM (SPP1+)",
  "2"  = "Malignant hepatocytes (Tumor 1)",
  "3"  = "Hepatocytes (periportal)",
  "4"  = "Stressed cells (heat-shock)",
  "5"  = "Liver sinusoidal endothelial cells",
  "6"  = "Kupffer cells",
  "7"  = "Hepatocytes (periportal, ALDH+)",
  "8"  = "Malignant hepatocytes (Tumor 2)",
  "9"  = "Macrophages (metallothionein-high)",
  "10" = "T cells (naive/memory)",
  "11" = "B cells",
  "12" = "T cells (activated)",
  "13" = "Malignant hepatocytes (Tumor 3)",
  "14" = "Hepatic stellate cells/CAFs",
  "15" = "Plasma cells",
  "16" = "Macrophages (TREM2+)",
  "17" = "Proliferating cells (cycling)",
  "18" = "Hepatocytes (pericentral)",
  "19" = "Hepatocytes (periportal)",
  "20" = "Hepatocytes (acute-phase)",
  "21" = "Tissue-resident macrophages",
  "22" = "TAM (SPP1+/TREM2+)",
  "23" = "Cholangiocytes"
)

# --- Apply labels as a new metadata column ---
seurat_obj_singlets$cell_type <- cluster_labels[as.character(Idents(seurat_obj_singlets))]

# --- Apply labels as a new metadata column (fixed) ---
seurat_obj_singlets$cell_type <- unname(cluster_labels[as.character(Idents(seurat_obj_singlets))])

# Confirm mapping worked correctly (no NAs, all clusters mapped)
table(seurat_obj_singlets$cell_type, useNA = "always")



# Confirm mapping worked correctly (no NAs, all clusters mapped)
table(seurat_obj_singlets$cell_type, useNA = "always")

# --- Visualize annotated UMAP ---
DimPlot(
  seurat_obj_singlets,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  label.size = 3
) +
  ggtitle("Annotated Cell Types") +
  theme(legend.position = "right")
# --- Save your fully annotated object ---
saveRDS(seurat_obj_singlets, "seurat_obj_annotated.rds")

# --- 1. Load both files ---
library(data.table)
library(dplyr)
library(tidyr)

cluster_markers <- fread("G:/scRNA/Project/data sheet 1.csv")
setnames(cluster_markers, "V1", "row_id", skip_absent = TRUE)  # first unnamed column is a row index

de_genes <- fread("G:/scRNA/Project/data sheet 2.csv")

# --- 2. Inspect structure ---
str(cluster_markers)
head(cluster_markers)
table(cluster_markers$cluster)  # should show 18 clusters, 20 genes each

str(de_genes)
head(de_genes)

# --- 3. Clean the DE gene lists ---
up_genes   <- de_genes$Upregulated_Genes[de_genes$Upregulated_Genes != ""]
down_genes <- de_genes$Downregulated_Genes[de_genes$Downregulated_Genes != ""]

length(up_genes)    # expect ~356
length(down_genes)  # expect ~168

up_genes <- unique(na.omit(up_genes))
down_genes <- unique(na.omit(down_genes))
# --- 4. Overlap analysis: which clusters' markers match the bulk signature? ---

overlap_summary <- cluster_markers %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    n_up_overlap = sum(gene %in% up_genes),
    n_down_overlap = sum(gene %in% down_genes),
    up_genes_matched = paste(gene[gene %in% up_genes], collapse = ", "),
    down_genes_matched = paste(gene[gene %in% down_genes], collapse = ", ")
  ) %>%
  arrange(desc(n_up_overlap))

print(overlap_summary, width = Inf)

# --- 5. Statistical significance of overlap (hypergeometric test) ---
# Background = total unique genes across all clusters' marker lists
total_genes_tested <- length(unique(cluster_markers$gene))

enrichment_test <- cluster_markers %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    n_up_overlap = sum(gene %in% up_genes)
  ) %>%
  rowwise() %>%
  mutate(
    p_value = phyper(
      q = n_up_overlap - 1,
      m = length(up_genes),
      n = total_genes_tested - length(up_genes),
      k = n_markers,
      lower.tail = FALSE
    )
  ) %>%
  arrange(p_value)

print(enrichment_test)

# --- 6. Visualize overlap ---
library(ggplot2)

overlap_long <- overlap_summary %>%
  select(cluster, n_up_overlap, n_down_overlap) %>%
  pivot_longer(cols = c(n_up_overlap, n_down_overlap),
               names_to = "direction", values_to = "n_genes")

ggplot(overlap_long, aes(x = reorder(cluster, -n_genes), y = n_genes, fill = direction)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  scale_fill_manual(values = c("n_up_overlap" = "firebrick", "n_down_overlap" = "steelblue"),
                    labels = c("Upregulated overlap", "Downregulated overlap")) +
  labs(title = "Cluster Marker Overlap with Bulk Tumor Signature",
       x = "Cluster", y = "Number of overlapping genes", fill = NULL)

# --- Diagnose: check the actual data going into the plot ---
str(overlap_long)
print(overlap_long)

# Check the numbers before pivoting
print(overlap_summary %>% select(cluster, n_up_overlap, n_down_overlap))

# Is n_genes actually numeric?
class(overlap_long$n_genes)
sum(overlap_long$n_genes, na.rm = TRUE)  # total genes matched across everything
# --- Check gene name formatting on both sides ---
head(cluster_markers$gene, 10)
head(up_genes, 10)
head(down_genes, 10)

# Check for hidden characters, quotes, or whitespace
sapply(head(cluster_markers$gene), function(x) nchar(x))
sapply(head(up_genes), function(x) nchar(x))

# Check exact class/type
class(cluster_markers$gene)
class(up_genes)


# --- Rebuild de_genes and the gene lists from scratch ---
library(data.table)

de_genes <- fread("G:/scRNA/Project/data sheet 2.csv")  # adjust path if needed

# Check it loaded correctly
str(de_genes)
head(de_genes)

# Extract and clean the two columns
up_genes <- de_genes$Upregulated_Genes[de_genes$Upregulated_Genes != ""]
down_genes <- de_genes$Downregulated_Genes[de_genes$Downregulated_Genes != ""]

up_genes <- trimws(unique(na.omit(up_genes)))
down_genes <- trimws(unique(na.omit(down_genes)))

# Confirm these are now proper character vectors with content
class(up_genes)
length(up_genes)
head(up_genes)

class(down_genes)
length(down_genes)
head(down_genes)

# --- Inspect de_genes directly ---
str(de_genes)
colnames(de_genes)
dim(de_genes)
head(de_genes, 10)

# Check what's actually in the Upregulated_Genes column before filtering
head(de_genes$Upregulated_Genes, 10)
class(de_genes$Upregulated_Genes)

# Check for NA vs empty string vs other blank representations
sum(de_genes$Upregulated_Genes == "", na.rm = TRUE)
sum(is.na(de_genes$Upregulated_Genes))
sum(de_genes$Upregulated_Genes == "" | is.na(de_genes$Upregulated_Genes))


# --- Get the EXACT column names R is actually using ---
colnames(de_genes)
names(de_genes)


# --- Bypass name matching entirely - use column position ---
up_genes <- de_genes[[1]]     # first column = upregulated genes
down_genes <- de_genes[[3]]   # third column = downregulated genes (2nd is blank spacer)

# Clean: remove NA and empty strings
up_genes <- trimws(up_genes)
up_genes <- unique(up_genes[up_genes != "" & !is.na(up_genes)])

down_genes <- trimws(down_genes)
down_genes <- unique(down_genes[down_genes != "" & !is.na(down_genes)])

# Verify
length(up_genes)
length(down_genes)
head(up_genes)
head(down_genes)
class(up_genes)

# Remove the header string that got included as data
up_genes <- up_genes[up_genes != "Upregulated_Genes"]
down_genes <- down_genes[down_genes != "Downregulated_Genes"]

length(up_genes)    # should now be 356
length(down_genes)  # should now be 168

head(up_genes)
head(down_genes)

sum(cluster_markers$gene %in% up_genes)
sum(cluster_markers$gene %in% down_genes)


library(dplyr)
library(tidyr)
library(ggplot2)

overlap_summary <- cluster_markers %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    n_up_overlap = sum(gene %in% up_genes),
    n_down_overlap = sum(gene %in% down_genes),
    up_genes_matched = paste(gene[gene %in% up_genes], collapse = ", "),
    down_genes_matched = paste(gene[gene %in% down_genes], collapse = ", ")
  ) %>%
  arrange(desc(n_up_overlap))

print(overlap_summary, width = Inf)

# Reshape for plotting
overlap_long <- overlap_summary %>%
  select(cluster, n_up_overlap, n_down_overlap) %>%
  pivot_longer(cols = c(n_up_overlap, n_down_overlap),
               names_to = "direction", values_to = "n_genes")

# Sanity check before plotting
sum(overlap_long$n_genes)  # should now be > 0

ggplot(overlap_long, aes(x = reorder(cluster, -n_genes), y = n_genes, fill = direction)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  scale_fill_manual(
    breaks = c("n_up_overlap", "n_down_overlap"),
    values = c("n_up_overlap" = "firebrick", "n_down_overlap" = "steelblue"),
    labels = c("Upregulated overlap", "Downregulated overlap")
  ) +
  labs(title = "Cluster Marker Overlap with Bulk Tumor Signature",
       x = "Cluster", y = "Number of overlapping genes", fill = NULL)

# --- Install if needed ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
 BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# --- 1. Convert gene symbols to Entrez IDs ---
up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Check mapping success rate
nrow(up_entrez)    # out of 356
nrow(down_entrez)  # out of 168
# --- 2. GO Biological Process enrichment ---
go_up <- enrichGO(
  gene = up_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

go_down <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

head(as.data.frame(go_up), 15)
head(as.data.frame(go_down), 15)

####XXXXXXXXXXXXX###########

# Reload your checkpoint rather than rebuilding from scratch
seurat_obj_singlets <- readRDS("seurat_obj_annotated.rds")

# Confirm it loaded correctly
seurat_obj_singlets
table(seurat_obj_singlets$cell_type)

library(data.table)

cluster_markers <- fread("G:/scRNA/Project/data sheet 1.csv")
de_genes <- fread("G:/scRNA/Project/data sheet 2.csv")

up_genes <- de_genes[[1]]
down_genes <- de_genes[[3]]

up_genes <- trimws(up_genes)
up_genes <- unique(up_genes[up_genes != "" & !is.na(up_genes) & up_genes != "Upregulated_Genes"])

down_genes <- trimws(down_genes)
down_genes <- unique(down_genes[down_genes != "" & !is.na(down_genes) & down_genes != "Downregulated_Genes"])

length(up_genes)   # should be 356
length(down_genes) # should be 168

library(clusterProfiler)
library(org.Hs.eg.db)

up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
nrow(up_entrez)


# --- Map down_genes too ---
down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
nrow(down_entrez)  # out of 168

# --- GO enrichment (local computation, should be fast and safe) ---
go_up <- enrichGO(
  gene = up_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

go_down <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

head(as.data.frame(go_up), 15)
head(as.data.frame(go_down), 15)

# View in RStudio's data viewer instead - much easier to read
View(as.data.frame(go_up))
View(as.data.frame(go_down))

# Or print with wider console width temporarily
options(width = 200)
head(as.data.frame(go_up)[, c("Description", "p.adjust", "Count")], 15)
head(as.data.frame(go_down)[, c("Description", "p.adjust", "Count")], 15)

library(enrichplot)
dotplot(go_up, showCategory = 15, title = "GO BP: Upregulated Genes")
dotplot(go_down, showCategory = 15, title = "GO BP: Downregulated Genes")


# --- Test network connectivity to KEGG first, before running the real call ---
# This should return quickly (a few seconds) if KEGG is reachable
test_conn <- tryCatch({
  readLines("https://rest.kegg.jp/info/kegg", n = 1, warn = FALSE)
}, error = function(e) {
  message("Cannot reach KEGG servers: ", e$message)
  NULL
})

test_conn  # if this prints something, KEGG is reachable; if NULL/error, that's your answer

# --- KEGG enrichment: upregulated genes only, first ---
kegg_up <- enrichKEGG(
  gene = up_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

head(as.data.frame(kegg_up), 15)

# --- KEGG enrichment: downregulated genes ---
kegg_down <- enrichKEGG(
  gene = down_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

head(as.data.frame(kegg_down), 15)

# --- View both with readable names ---
kegg_up_df <- as.data.frame(kegg_up)
kegg_down_df <- as.data.frame(kegg_down)

kegg_up_df[, c("Description", "p.adjust", "Count")]
kegg_down_df[, c("Description", "p.adjust", "Count")]

# --- Visualize ---
library(enrichplot)
dotplot(kegg_up, showCategory = 15, title = "KEGG: Upregulated Genes")
dotplot(kegg_down, showCategory = 15, title = "KEGG: Downregulated Genes")

# --- 1. View readable KEGG pathway names ---
kegg_up_df <- as.data.frame(kegg_up)
kegg_down_df <- as.data.frame(kegg_down)

kegg_up_df[, c("Description", "p.adjust", "Count")]
kegg_down_df[, c("Description", "p.adjust", "Count")]

# --- 2. Consolidated summary table across all analyses ---
library(dplyr)


# Top GO/KEGG terms
go_up_top <- as.data.frame(go_up) %>% slice_min(p.adjust, n = 10) %>%
  mutate(analysis = "GO_BP", direction = "Upregulated") %>%
  select(analysis, direction, Description, p.adjust, Count)

go_down_top <- as.data.frame(go_down) %>% slice_min(p.adjust, n = 10) %>%
  mutate(analysis = "GO_BP", direction = "Downregulated") %>%
  select(analysis, direction, Description, p.adjust, Count)

kegg_up_top <- as.data.frame(kegg_up) %>% slice_min(p.adjust, n = 10) %>%
  mutate(analysis = "KEGG", direction = "Upregulated") %>%
  select(analysis, direction, Description, p.adjust, Count)

kegg_down_top <- as.data.frame(kegg_down) %>% slice_min(p.adjust, n = 10) %>%
  mutate(analysis = "KEGG", direction = "Downregulated") %>%
  select(analysis, direction, Description, p.adjust, Count)

enrichment_summary <- bind_rows(go_up_top, go_down_top, kegg_up_top, kegg_down_top)
print(enrichment_summary, n = 40)


# Check what's currently set for na.print
options("na.print")

# Reset it to the default (NULL means use standard NA printing)
options(na.print = NULL)

# Now try printing again
print(enrichment_summary, n = 40)


View(enrichment_summary)
write.csv(enrichment_summary, "enrichment_summary_table.csv", row.names = FALSE)



# --- 3. Combined 4-panel summary figure ---
library(ggplot2)
library(patchwork)
library(enrichplot)

library(Seurat)

# Confirm it loaded
packageVersion("Seurat")


p_umap <- DimPlot(
  seurat_obj_singlets, reduction = "umap", group.by = "cell_type",
  label = TRUE, repel = TRUE, label.size = 2.5
) + ggtitle("Annotated Cell Types") + theme(legend.position = "none")

p_overlap <- ggplot(overlap_long, aes(x = reorder(cluster, -n_genes), y = n_genes, fill = direction)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  scale_fill_manual(
    breaks = c("n_up_overlap", "n_down_overlap"),
    values = c("n_up_overlap" = "firebrick", "n_down_overlap" = "steelblue"),
    labels = c("Up overlap", "Down overlap")
  ) +
  labs(title = "Cluster Marker Overlap with Bulk Signature", x = "Cluster", y = "N genes", fill = NULL) +
  theme(legend.position = "bottom")

#####

# Check what survived the restart
exists("cluster_markers")
exists("up_genes")
exists("down_genes")
exists("overlap_summary")


library(dplyr)
library(tidyr)

overlap_summary <- cluster_markers %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    n_up_overlap = sum(gene %in% up_genes),
    n_down_overlap = sum(gene %in% down_genes)
  ) %>%
  arrange(desc(n_up_overlap))

overlap_long <- overlap_summary %>%
  select(cluster, n_up_overlap, n_down_overlap) %>%
  pivot_longer(cols = c(n_up_overlap, n_down_overlap),
               names_to = "direction", values_to = "n_genes")

# Confirm it built correctly
exists("overlap_long")
print(overlap_summary)


p_go_up <- dotplot(go_up, showCategory = 8) + ggtitle("GO: Upregulated") +
  theme(axis.text.y = element_text(size = 7))

p_go_down <- dotplot(go_down, showCategory = 8) + ggtitle("GO: Downregulated") +
  theme(axis.text.y = element_text(size = 7))


p_overlap <- ggplot(overlap_long, aes(x = reorder(cluster, -n_genes), y = n_genes, fill = direction)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  scale_fill_manual(
    breaks = c("n_up_overlap", "n_down_overlap"),
    values = c("n_up_overlap" = "firebrick", "n_down_overlap" = "steelblue"),
    labels = c("Up overlap", "Down overlap")
  ) +
  labs(title = "Cluster Marker Overlap with Bulk Signature", x = "Cluster", y = "N genes", fill = NULL) +
  theme(legend.position = "bottom")

p_overlap  # confirm it renders before combining

# Combine into a 2x2 layout
combined_summary <- (p_umap | p_overlap) / (p_go_up | p_go_down)
combined_summary

# --- 4. Save the final figure ---
ggsave("HCC01T_analysis_summary.png", combined_summary, width = 16, height = 12, dpi = 300)

RNA1
