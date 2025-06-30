# Title: Blastoid Blastocyst Analysis —— in house toti cells （10X）
# Author: Gaozy
# Time: 2025-05-10
# Goal: 我们in house toti cells （10X）是否其他细胞相比，其整体表达谱更类似于小鼠2C胚胎

# comput172-r_env-R
# 0. Basic settings ----
setwd("/home/qcao02/gaozy/blastocyst/")

# 1. library ----
library(Matrix, lib.loc = "/usr/lib64/R/library") 
library(sva, lib.loc = "/home/qcao02/R/x86_64-conda-linux-gnu-library/4.4") 

.libPaths(c("/usr/lib64/R/library",
            "/usr/share/R/library",
            "/home/lushi02/.local/share/r-miniconda/envs/r-reticulate",
            "/home/lushi02/anaconda3/lib/R/library"))

library(Matrix, lib.loc = "/usr/lib64/R/library") 
library(Seurat, lib.loc = "/usr/lib64/R/library") 
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggcorrplot, lib.loc = "/home/qcao02/R/x86_64-conda-linux-gnu-library/4.4")
library(Hmisc)
library(dplyr, lib.loc = "/usr/lib64/R/library")
library(reshape2)
library(ggrepel)
library(patchwork)  # For combining plots
library(readr)

if (dir.exists("250513") == F) dir.create('250513')
setwd("250513")

on.exit({
  while(dev.cur() > 1) dev.off()
})

# 2. load data ----
load_data <- function(expr_path, meta_path, output_path, data_name){
  # ①.read data ----
  counts <- read.table(expr_path, sep  = " ", header = T, row.names = 1)
  meta <- read_delim(meta_path, delim = ",", col_names = TRUE) %>% 
    as.data.frame %>% tibble::column_to_rownames("Run")
  
  # ② create seurat obj ----
  seurat_data <- CreateSeuratObject(counts = counts, project = data_name, min.cells = 3)
  
  check_meta <- c("Developmental_stage", "develomental_stage", "embryo_derived_from")
  if ( all(check_meta %in% colnames(meta) == F) ) {
    seurat_data$cell_type <- meta[rownames(seurat_data@meta.data), "Cell_type"]
  } else {
    seurat_data$cell_type <- meta[rownames(seurat_data@meta.data), intersect(check_meta, colnames(meta))]
  }
  
  # ③ save ----
  if (dir.exists(output_path) == F) dir.create(output_path)
  
  saveRDS(seurat_data, paste0(output_path, data_name, ".rds"))
  print(paste("Save data as:", paste0(output_path, data_name,".rds -----------------------------------")))
}

## a) GSE45719 ----
load_data(expr_path = "/home/fengyan02/Project/YYYProject202306/processedData/scRNAData/GSE45719/kallistoResult/abundance_sum.tsv",
          meta_path = "/home/fengyan02/Project/YYYProject202306/rawdata/GSE45719/Metadata/SraRunTable.txt",
          output_path = "/home/qcao02/gaozy/blastocyst/rawdata/",
          data_name = "GSE45719")

GSE45719 <- readRDS("/home/qcao02/gaozy/blastocyst/rawdata/GSE45719.rds")
GSE45719$cell_type <- gsub("-cell stage blastomere.*", "C", GSE45719$cell_type)
GSE45719$cell_type <- gsub(".*blastocyst.*", "Pub Blastocyst", GSE45719$cell_type)

## b) GSE100597 ----
load_data(expr_path = "/home/fengyan02/Project/YYYProject202306/processedData/scRNAData/GSE100597/kallistoResult/abundance_sum.tsv",
          meta_path = "/home/fengyan02/Project/YYYProject202306/rawdata/GSE100597/Meatdata/SraRunTable.txt",
          output_path = "/home/qcao02/gaozy/blastocyst/rawdata/",
          data_name = "GSE100597")

GSE100597 <- readRDS("/home/qcao02/gaozy/blastocyst/rawdata/GSE100597.rds")

## c) GSE78140 ----
load_data(expr_path = "/home/fengyan02/Project/YYYProject202306/processedData/scRNAData/GSE78140/kallistoResult/abundance_sum.tsv",
          meta_path = "/home/fengyan02/Project/YYYProject202306/rawdata/GSE78140/Metadata/SraRunTable.txt",
          output_path = "/home/qcao02/gaozy/blastocyst/rawdata/",
          data_name = "GSE78140")

GSE78140 <- readRDS("/home/qcao02/gaozy/blastocyst/rawdata/GSE78140.rds")
GSE78140$cell_type <- gsub("-cell embryo.*", "C", GSE78140$cell_type)
GSE78140$cell_type <- gsub("Embryonic stem cell", "mESC", GSE78140$cell_type)
GSE78140 <- subset(GSE78140, cell_type %in% c("mESC", "Zygote", "2C", "4C", "8C"))

## d) GSE185000 ----
load_data(expr_path = "/home/fengyan02/Project/YYYProject202306/processedData/scRNAData/GSM5603033-42/kallistoResult/abundance_sum.tsv",
          meta_path = "/home/fengyan02/Project/YYYProject202306/rawdata/GSM5603033-42/Metadata/SraRunTable.txt",
          output_path = "/home/qcao02/gaozy/blastocyst/rawdata/",
          data_name = "GSE185000")

GSE185000 <- readRDS("/home/qcao02/gaozy/blastocyst/rawdata/GSE185000.rds")

## e) GSE145609 ----
expr_path = "/home/fengyan02/Project/YYYProject202306/processedData/scRNAData/GSE145609/kallistoResult/abundance_sum.tsv"
meta_path = "/home/fengyan02/Project/YYYProject202306/rawdata/GSE145609/Metadata/SraRunTable.txt"
output_path = "/home/qcao02/gaozy/blastocyst/rawdata/"
data_name = "GSE145609"

  counts <- read.table(expr_path, sep  = " ", header = T, row.names = 1)
  meta <- read_delim(meta_path, delim = ",", col_names = TRUE) %>% 
    as.data.frame %>% tibble::column_to_rownames("Run")
  
  seurat_data <- CreateSeuratObject(counts = counts, project = data_name, min.cells = 3)
  seurat_data$cell_type <- meta[colnames(seurat_data), "Cell_type"]
  
  if (dir.exists(output_path) == F) dir.create(output_path)
  
  saveRDS(seurat_data, paste0(output_path, data_name, ".rds"))
  print(paste("Save data as:", paste0(output_path, data_name,".rds -----------------------------------")))

GSE145609 <- readRDS("/home/qcao02/gaozy/blastocyst/rawdata/GSE145609.rds")
GSE145609$cell_type <- gsub(" embryo", "", GSE145609$cell_type)
GSE145609$cell_type <- gsub("embryonic stem cells", "mESC", GSE145609$cell_type)

## f) GSE183522 ----
load_data(expr_path = "/home/fengyan02/Project/YYYProject202306/processedData/scRNAData/GSM5558728-33/kallistoResult/abundance_sum.tsv",
          meta_path = "/home/fengyan02/Project/YYYProject202306/rawdata/GSM5558728-33/Metadata/SraRunTable.txt",
          output_path = "/home/qcao02/gaozy/blastocyst/rawdata/",
          data_name = "GSE183522")

GSE183522 <- readRDS("/home/qcao02/gaozy/blastocyst/rawdata/GSE183522.rds")
GSE183522$cell_type <- gsub("mouse totipotent potential stem cells", "TPS", GSE183522$cell_type)

## g) GSM5195025 ----
data <- Read10X(data.dir = "/home/fengyan02/Project/YYYProject202306/processedData/scRNAData/GSM5195025/run_SC-TBLC/outs/filtered_feature_bc_matrix")
GSM5195025 <- CreateSeuratObject(counts = data, assay = "RNA", project = "GSM5195025", min.cells = 3)

GSM5195025$cell_type <- "TBLC"

## h) GSM7798458 ----
data <- read.delim("/home/qcao02/gaozy/blastocyst/rawdata/GSM7798458/GSM7798458_mTBLC-2MYCP_counts.txt",
                   header = TRUE,
                   check.names = F)
meta <- read.table("/home/qcao02/gaozy/blastocyst/rawdata/GSM7798458/GSM7798458_mTBLC-2MYCP_celltype.txt", sep = "\t", header = TRUE)
rownames(meta) <- paste0("eTBLCs_", rownames(meta))
colnames(meta)[2] <- "cell_type"

#assay_data <- CreateAssayObject(counts = Matrix(as.matrix(data), sparse = TRUE), min.cells = 3, min.features = 200)
GSM7798458 <- CreateSeuratObject(counts = data, assay = "RNA", meta.data = meta,
                                 project = "GSM7798458", min.cells = 3, min.features = 200)
GSM7798458$barcodes <- NULL

GSM7798458 <- subset(GSM7798458, cell_type %nin% "feeder cells")

## i) GSM5603638_39 ----
load_data(expr_path = "/home/fengyan02/Project/YYYProject202306/processedData/scRNAData/GSM5603636-41/kallistoResult/abundance_sum.tsv",
          meta_path = "/home/fengyan02/Project/YYYProject202306/rawdata/GSM5603636-41/Metadata/SraRunTable.txt",
          output_path = "/home/qcao02/gaozy/blastocyst/rawdata/",
          data_name = "GSM5603636-41")

GSM5603636_41 <- readRDS("/home/qcao02/gaozy/blastocyst/rawdata/GSM5603636-41.rds")

GSM5603638_39 <- subset(GSM5603636_41, cells = c("SRR16119544", "SRR16119545"))
GSM5603638_39$cell_type <- "TLSC"

## j) Inhouse - toti cells ----
data <- Read10X(data.dir = "/home/qcao02/gaozy/blastocyst/250410_mapping/KSXY-QHDX-DO25021301-2/CellRanger_result/CellRanger_result/ED251379910/filtered_gene_bc_matrices/ref/")

inhouse <- CreateSeuratObject(counts = data, project = "Toti", assay = "RNA",
                              min.cells = 3, min.features = 300)

inhouse$cell_type <- "In-house Toti"

### ① calculate the percentage of mitochondria ----
inhouse[["percent.mt"]] <- PercentageFeatureSet(inhouse, pattern = "^mt-")

qc_plot <- VlnPlot(inhouse, 
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3,
                   pt.size = 0.1) +
  theme(plot.title = element_text(size=10))

ggsave("QC_metrics_before_filtering.pdf", qc_plot, width = 10, height = 4)

### ② filter mitochondria ----
inhouse_filtered <- subset(inhouse, subset = nFeature_RNA > 300 & percent.mt < 10)

# 过滤前后细胞数对比 - percent.mt < 10
cat("过滤前细胞数:", ncol(inhouse), "\n",
    "过滤后细胞数:", ncol(inhouse_filtered), "\n",
    "过滤比例:", round(1 - ncol(inhouse_filtered)/ncol(inhouse), 2)*100, "%\n")

# plot
qc_plot_filtered <- VlnPlot(inhouse_filtered, 
                            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                            ncol = 3,
                            pt.size = 0.1)

ggsave("QC_metrics_after_filtering.pdf", qc_plot_filtered, width = 10, height = 4)

inhouse$percent.mt <- NULL
saveRDS(inhouse_filtered, file = "inhouse_toti_filtered.rds") # 25493 features across 13936 samples

## j) merge data ----
seurat_list <- list( inhouse_filtered, GSE45719, GSE100597, GSE78140, GSE185000,
  GSE145609, GSE183522, GSM5195025, GSM7798458, GSM5603638_39)

merge_data <- Reduce(function(x, y) merge(x, y), seurat_list)

valid_cells <- rownames(merge_data@meta.data[!is.na(merge_data@meta.data$cell_type), ])
merge_data <- subset(merge_data, cells = valid_cells)

saveRDS(merge_data, "raw_merge_data.rds")

# 3. batch correction ----
## 1) ComBat_seq ----
method <- "ComBat_seq"
### a. batch correction ----
library(edgeR, lib.loc = "/home/qcao02/R/x86_64-conda-linux-gnu-library/4.4")
table(merge_data$orig.ident, merge_data$orig.ident) 

batch_list <- SplitObject(merge_data, split.by = "orig.ident")
#### Normalize for each dataset
for(i in 1:length(batch_list)){
  batch_list[[i]] <- NormalizeData(batch_list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}

combined_obj <- merge(batch_list[[1]], y=batch_list[2:length(batch_list)], add.cell.ids = names(batch_list), project = "Combined")
combined_obj <- subset(combined_obj, nCount_RNA != 0)

expr_mat <- as.matrix(combined_obj@assays$RNA@data)

batch_info <- combined_obj@meta.data$orig.ident

expr_corrected <- ComBat_seq(counts = expm1(expr_mat), batch = batch_info, group = NULL)

expr_corrected_log <- log1p(expr_corrected)

corrected_seurat_combat <- CreateSeuratObject(counts = expr_corrected_log, meta.data = combined_obj@meta.data)

### b. choose proper PCs ----
Idents(corrected_seurat_combat) <- "cell_type"
corrected_seurat_combat <- FindVariableFeatures(corrected_seurat_combat, selection.method = "vst", nfeatures = 2000)
corrected_seurat_combat <- ScaleData(corrected_seurat_combat, features = rownames(corrected_seurat_combat)) %>%
  RunPCA(npcs=50, features = rownames(corrected_seurat_combat))

### choose proper PCs
corrected_seurat_combat <- JackStraw(corrected_seurat_combat, num.replicate = 50, dims = 50)
corrected_seurat_combat <- ScoreJackStraw(corrected_seurat_combat, dims = 1:50) 

p1 <- JackStrawPlot(corrected_seurat_combat, dims = 1:50)
p2 <- ElbowPlot(corrected_seurat_combat, ndims = 50)

pdf(paste0("corrected_pca_elbow_", method, ".pdf"), width = 12, height = 8, onefile = F)
p1 / p2
dev.off()

## 2) Harmony ----
method <- "Harmony"

library(harmony)
Idents(merge_data) <- "orig.ident"

merge_data <- NormalizeData(merge_data, normalization.method = "LogNormalize", scale.factor = 10000)
merge_data <- FindVariableFeatures(merge_data, selection.method = "vst", nfeatures = 2000)
merge_data <- ScaleData(merge_data, features = rownames(merge_data)) %>%
  RunPCA(npcs=30, features = rownames(merge_data))

### a. run Harmony ----
corrected_seurat_harmony <- RunHarmony(object = merge_data, group.by.vars = "orig.ident", plot_convergence = TRUE)

### b. choose proper PCs ----
Idents(corrected_seurat_harmony) <- "cell_type"
corrected_seurat_harmony <- NormalizeData(corrected_seurat_harmony, normalization.method = "LogNormalize", scale.factor = 10000)
corrected_seurat_harmony <- FindVariableFeatures(corrected_seurat_harmony, selection.method = "vst", nfeatures = 2000)
corrected_seurat_harmony <- ScaleData(corrected_seurat_harmony, features = rownames(corrected_seurat_harmony)) %>%
  RunPCA(npcs=50, features = rownames(corrected_seurat_harmony))

### choose proper PCs
corrected_seurat_harmony <- JackStraw(corrected_seurat_harmony, num.replicate = 50, dims = 50)
corrected_seurat_harmony <- ScoreJackStraw(corrected_seurat_harmony, dims = 1:50) 

p1 <- JackStrawPlot(corrected_seurat_harmony, dims = 1:50)
p2 <- ElbowPlot(corrected_seurat_harmony, ndims = 50)

pdf(paste0("corrected_pca_elbow_", method, ".pdf"), width = 12, height = 12, onefile = F)
print(p1 / p2)
dev.off()

# 4. UMAP/tSNE/PCA分析，胚胎公共数据 ----
## 1) result - Combat_seq ----
### a. PCA ----
corrected_seurat_combat$orig.ident <- as.character(corrected_seurat_combat$orig.ident)
corrected_seurat_combat$orig.ident["GSM5603636-41_SRR16119544"] <- "GSM5603638"
corrected_seurat_combat$orig.ident["GSM5603636-41_SRR16119545"] <- "GSM5603639"
corrected_seurat_combat$orig.ident <- factor(corrected_seurat_combat$orig.ident)

corrected_seurat_combat <- ScaleData(corrected_seurat_combat, features = rownames(corrected_seurat_combat)) %>%
  RunPCA(npcs=30, features = rownames(corrected_seurat_combat))

corrected_seurat_combat$orig.ident <- factor(corrected_seurat_combat$orig.ident, levels = unique(corrected_seurat_combat$orig.ident))

## Get PCA coordinates and metadata
pca_var <- corrected_seurat_combat[["pca"]]@stdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)

p <- DimPlot( corrected_seurat_combat, reduction = "pca", dims = c(1, 2), group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 1) +
  labs( x = paste0("PC1 (", pca_var_perc[1], "%)"),
        y = paste0("PC2 (", pca_var_perc[2], "%)") )

pdf(paste0("corrected_pca_", method, ".pdf"), width=12, height=9.5, onefile = F)
p
dev.off()

### b. UMAP & tSNE ----
corrected_seurat_combat <- FindNeighbors(corrected_seurat_combat, dims = 1:30)
corrected_seurat_combat <- FindClusters(corrected_seurat_combat, resolution = 0.5)

corrected_seurat_combat <- RunTSNE(corrected_seurat_combat, dims = 1:30 )
corrected_seurat_combat <- RunUMAP(corrected_seurat_combat, dims = 1:30 )

# UMAP plot
p1 <- DimPlot(corrected_seurat_combat, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 1)  +
  ggtitle("UMAP") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# tSNE plot
p2 <- DimPlot(corrected_seurat_combat, reduction = "tsne", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 1)  +
  ggtitle("tSNE") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# Combine plots
final_plot <- (p1 | p2) + plot_layout(guides = "collect")

pdf(paste0("corrected_umap_tsne_", method, ".pdf"), width=25, height=12, onefile = F)
final_plot
dev.off()

saveRDS(corrected_seurat_combat, "combat.rds")

## 2) result - Harmony ----
### a. PCA ----
corrected_seurat_harmony$orig.ident <- as.character(corrected_seurat_harmony$orig.ident)
corrected_seurat_harmony$orig.ident["SRR16119544"] <- "GSM5603638"
corrected_seurat_harmony$orig.ident["SRR16119545"] <- "GSM5603639"
corrected_seurat_harmony$orig.ident <- factor(corrected_seurat_harmony$orig.ident)

corrected_seurat_harmony <- ScaleData(corrected_seurat_harmony, features = rownames(corrected_seurat_harmony)) %>%
  RunPCA(npcs=30, features = rownames(corrected_seurat_harmony))

corrected_seurat_harmony$cell_type <- factor(corrected_seurat_harmony$cell_type, levels = unique(corrected_seurat_harmony$cell_type))

## Get PCA coordinates and metadata
pca_var <- corrected_seurat_harmony[["pca"]]@stdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)

p <- DimPlot( corrected_seurat_harmony, reduction = "pca", dims = c(1, 2), group.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 1) +
  labs( x = paste0("PC1 (", pca_var_perc[1], "%)"),
        y = paste0("PC2 (", pca_var_perc[2], "%)") )

pdf(paste0("corrected_pca_", method, ".pdf"), width=8, height=6, onefile = F)
p
dev.off()

### b. UMAP & tSNE ----
corrected_seurat_harmony <- FindNeighbors(corrected_seurat_harmony, dims = 1:30)
corrected_seurat_harmony <- FindClusters(corrected_seurat_harmony, resolution = 0.5)

corrected_seurat_harmony <- RunTSNE(corrected_seurat_harmony, dims = 1:30 )
corrected_seurat_harmony <- RunUMAP(corrected_seurat_harmony, dims = 1:30 )

# UMAP plot
p1 <- DimPlot(corrected_seurat_harmony, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 1)  +
  ggtitle("UMAP") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# tSNE plot
p2 <- DimPlot(corrected_seurat_harmony, reduction = "tsne", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 1)  +
  ggtitle("tSNE") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# Combine plots
final_plot <- (p1 | p2) + plot_layout(guides = "collect")

pdf(paste0("corrected_umap_tsne_", method, ".pdf"), width=25, height=12, onefile = F)
final_plot
dev.off()

## 3) results - Raw merged data ----
method <- "raw_merged_data"

### a. PCA ----
Idents(merge_data) <- "cell_info"
merge_data <- NormalizeData(merge_data, normalization.method = "LogNormalize", scale.factor = 10000)
merge_data <- FindVariableFeatures(merge_data, selection.method = "vst", nfeatures = 2000)
merge_data <- ScaleData(merge_data, features = rownames(merge_data)) %>%
  RunPCA(npcs=30, features = rownames(merge_data))

merge_data$cell_info <- factor(merge_data$cell_info, levels = unique(merge_data$cell_info))

## Get PCA coordinates and metadata
pca_var <- merge_data[["pca"]]@stdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)

p <- DimPlot( merge_data, reduction = "pca", dims = c(1, 2), group.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 1) +
  labs( x = paste0("PC1 (", pca_var_perc[1], "%)"),
        y = paste0("PC2 (", pca_var_perc[2], "%)") )

pdf(paste0("corrected_pca_", method, ".pdf"), width = 12, height = 12, onefile = F)
p
dev.off()

### b. UMAP & tSNE ----
merge_data <- RunTSNE(merge_data, dims = 1:30 )
merge_data <- RunUMAP(merge_data, dims = 1:30 )

# UMAP plot
p1 <- DimPlot(merge_data, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 1)  +
  ggtitle("UMAP") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# tSNE plot
p2 <- DimPlot(merge_data, reduction = "tsne", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 1)  +
  ggtitle("tSNE") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# Combine plots
final_plot <- (p1 | p2) + plot_layout(guides = "collect")

pdf(paste0("corrected_umap_tsne_", method, ".pdf"), width=25, height=12, onefile = F)
final_plot
dev.off()

### c. save data ----
saveRDS(merge_data, "merge_data.rds")

# 5. Specific genes 表达 ----
Totipotent <- read.csv("Totipotent genes.csv", check.names = F)
Pluripotent <- read.csv("Pluripotent genes.csv", check.names = F)
Maternal <- read.csv("Maternal genes.csv", check.names = F)
Minor_ZGA <- read.csv("Minor ZGA genes.csv", check.names = F)
Major_ZGA <- read.csv("Major ZGA genes.csv", check.names = F)

marker_list <- list(
  "Totipotent genes" = Totipotent[,2],
  "Pluripotent genes" = Pluripotent[,2],
  "Maternal genes" = Maternal[,2],
  "Minor ZGA genes" = Minor_ZGA[,2],
  "Major ZGA genes" = Major_ZGA[,2]
)

marker_dt <- data.frame(
  gene = c(Totipotent[,2], Pluripotent[,2], Maternal[,2], Minor_ZGA[,2], Major_ZGA[,2]),
  type = c(
    rep("Totipotent genes", nrow(Totipotent)),
    rep("Pluripotent genes", nrow(Pluripotent)),
    rep("Maternal genes", nrow(Maternal)),
    rep("Minor ZGA genes", nrow(Minor_ZGA)),
    rep("Major ZGA genes", nrow(Major_ZGA))
  )
)
names(marker_dt$type) <-c(Totipotent[,2], Pluripotent[,2], Maternal[,2], Minor_ZGA[,2], Major_ZGA[,2])
gene_meta <- marker_dt$type
names(gene_meta) <-c(Totipotent[,2], Pluripotent[,2], Maternal[,2], Minor_ZGA[,2], Major_ZGA[,2])

### 1) heatmap - ComBat_seq ----
method <- "ComBat_seq"
Idents(corrected_seurat_combat) <- "cell_type"
options(future.globals.maxSize= 4000*1024^2)

pb_seurat <- Seurat:::PseudobulkExpression(corrected_seurat_combat, return.seurat = T)
pb_seurat@meta.data$cell_type <- rownames(pb_seurat@meta.data)

RhpcBLASctl::blas_set_num_threads(8)
Idents(pb_seurat) <- "cell_type"
pb_seurat <- NormalizeData(pb_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
pb_seurat <- FindVariableFeatures(pb_seurat, selection.method = "vst", nfeatures = 2000)
pb_seurat <- ScaleData(pb_seurat, features = VariableFeatures(pb_seurat)) %>%
  RunPCA(npcs=10, features = VariableFeatures(pb_seurat))

pb_seurat$cell_type <- factor(pb_seurat$cell_type, levels = unique(pb_seurat$cell_type))

exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data")
exp_data <- exp_data[intersect(rownames(exp_data), marker_dt$gene), ]
colnames(exp_data) <- pb_seurat$cell_type

meta <- pb_seurat@meta.data
exp_scale <- exp_data

pdf(paste0("corrected_marker_heatmap_", method, ".pdf"), width = 8, height = 63, onefile = F)
Heatmap(as.matrix(exp_scale),
        row_split = gene_meta[rownames(exp_scale)],
        column_split = meta$cell_type,
        show_column_names = TRUE,  # Show column names
        show_column_dend = TRUE,  # Show column dendrogram
        column_title = "Markers Expression Heatmap",
        col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
        heatmap_legend_param = list(title=''))
dev.off()

### 2) heatmap - Harmony ----
method <- "Harmony"
Idents(corrected_seurat_harmony) <- "cell_type"
options(future.globals.maxSize= 4000*1024^2)

pb_seurat <- Seurat:::PseudobulkExpression(corrected_seurat_harmony, return.seurat = T)
pb_seurat@meta.data$cell_type <- rownames(pb_seurat@meta.data)

RhpcBLASctl::blas_set_num_threads(8)
Idents(pb_seurat) <- "cell_type"
pb_seurat <- NormalizeData(pb_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
pb_seurat <- FindVariableFeatures(pb_seurat, selection.method = "vst", nfeatures = 2000)
pb_seurat <- ScaleData(pb_seurat, features = VariableFeatures(pb_seurat)) %>%
  RunPCA(npcs=10, features = VariableFeatures(pb_seurat))

pb_seurat$cell_type <- factor(pb_seurat$cell_type, levels = unique(pb_seurat$cell_type))

exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data")
exp_data <- exp_data[intersect(rownames(exp_data), marker_dt$gene), ]
colnames(exp_data) <- pb_seurat$cell_type

meta <- pb_seurat@meta.data
exp_scale <- exp_data

pdf(paste0("corrected_marker_heatmap_", method, ".pdf"), width = 8, height = 63, onefile = F)
Heatmap(as.matrix(exp_scale),
        row_split = gene_meta[rownames(exp_scale)],
        column_split = meta$cell_type,
        show_column_names = TRUE,  # Show column names
        show_column_dend = TRUE,  # Show column dendrogram
        column_title = "Markers Expression Heatmap",
        col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
        heatmap_legend_param = list(title=''))
dev.off()

### 3) heatmap - raw merge data ----
method <- "raw_merge_data"
Idents(merge_data) <- "cell_type"
options(future.globals.maxSize= 4000*1024^2)

pb_seurat <- Seurat:::PseudobulkExpression(merge_data, return.seurat = T)
pb_seurat@meta.data$cell_type <- rownames(pb_seurat@meta.data)

RhpcBLASctl::blas_set_num_threads(8)
Idents(pb_seurat) <- "cell_type"
pb_seurat <- NormalizeData(pb_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
pb_seurat <- FindVariableFeatures(pb_seurat, selection.method = "vst", nfeatures = 2000)
pb_seurat <- ScaleData(pb_seurat, features = VariableFeatures(pb_seurat)) %>%
  RunPCA(npcs=10, features = VariableFeatures(pb_seurat))

pb_seurat$cell_type <- factor(pb_seurat$cell_type, levels = unique(pb_seurat$cell_type))

exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data")
exp_data <- exp_data[intersect(rownames(exp_data), marker_dt$gene), ]
colnames(exp_data) <- pb_seurat$cell_type

meta <- pb_seurat@meta.data
exp_scale <- exp_data

pdf(paste0("corrected_marker_heatmap_", method, ".pdf"), width = 8, height = 63, onefile = F)
Heatmap(as.matrix(exp_scale),
        row_split = gene_meta[rownames(exp_scale)],
        column_split = meta$cell_type,
        show_column_names = TRUE,  # Show column names
        show_column_dend = TRUE,  # Show column dendrogram
        column_title = "Markers Expression Heatmap",
        col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
        heatmap_legend_param = list(title=''))
dev.off()

# 6. UMAP或者tSNE分析 totipotent genes 及pluripotent genes 表达情况 ----
## 1) ComBat_seq ----
method <- "ComBat_seq"
seurat_obj <- corrected_seurat_combat

seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Totipotent[,2], rownames(seurat_obj))), name = "Totipotent")
seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Pluripotent[,2], rownames(seurat_obj))), name = "Pluripotent")
seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Maternal[,2], rownames(seurat_obj))), name = "Maternal")
seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Minor_ZGA[,2], rownames(seurat_obj))), name = "Minor_ZGA")
seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Major_ZGA[,2], rownames(seurat_obj))), name = "Major_ZGA")

p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 1)  +
  ggtitle("UMAP") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

p2 <- FeaturePlot(seurat_obj, features = "Totipotent1") + ggtitle("Expression of Totipotent genes")
p3 <- FeaturePlot(seurat_obj, features = "Pluripotent1") + ggtitle("Expression of Pluripotent genes")
p4 <- FeaturePlot(seurat_obj, features = "Maternal1") + ggtitle("Expression of Maternal genes")
p5 <- FeaturePlot(seurat_obj, features = "Minor_ZGA1") + ggtitle("Expression of Minor ZGA genes")
p6 <- FeaturePlot(seurat_obj, features = "Major_ZGA1") + ggtitle("Expression of Major ZGA genes")

pdf(paste0("specific_feature_plot_", method, ".pdf"), width=28, height=14, onefile = F)
(p1 | p2 ) / (p3 | p4 | p5 | p6)
dev.off()

## 2) Harmony ----
method <- "Harmony"
seurat_obj <- corrected_seurat_harmony

seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Totipotent[,2], rownames(seurat_obj))), name = "Totipotent")
seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Pluripotent[,2], rownames(seurat_obj))), name = "Pluripotent")
seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Maternal[,2], rownames(seurat_obj))), name = "Maternal")
seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Minor_ZGA[,2], rownames(seurat_obj))), name = "Minor_ZGA")
seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Major_ZGA[,2], rownames(seurat_obj))), name = "Major_ZGA")

p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 1)  +
  ggtitle("UMAP") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

p2 <- FeaturePlot(seurat_obj, features = "Totipotent1") + ggtitle("Expression of Totipotent genes")
p3 <- FeaturePlot(seurat_obj, features = "Pluripotent1") + ggtitle("Expression of Pluripotent genes")
p4 <- FeaturePlot(seurat_obj, features = "Maternal1") + ggtitle("Expression of Maternal genes")
p5 <- FeaturePlot(seurat_obj, features = "Minor_ZGA1") + ggtitle("Expression of Minor ZGA genes")
p6 <- FeaturePlot(seurat_obj, features = "Major_ZGA1") + ggtitle("Expression of Major ZGA genes")

pdf(paste0("specific_feature_plot_", method, ".pdf"), width=28, height=14, onefile = F)
(p1 | p2 ) / (p3 | p4 | p5 | p6)
dev.off()

## 3) raw merge data ----
method <- "raw_merge_data"
seurat_obj <- merge_data

seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Totipotent[,2], rownames(seurat_obj))), name = "Totipotent")
seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Pluripotent[,2], rownames(seurat_obj))), name = "Pluripotent")
seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Maternal[,2], rownames(seurat_obj))), name = "Maternal")
seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Minor_ZGA[,2], rownames(seurat_obj))), name = "Minor_ZGA")
seurat_obj <- AddModuleScore(seurat_obj, features = list( intersect(Major_ZGA[,2], rownames(seurat_obj))), name = "Major_ZGA")

p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 1)  +
  ggtitle("UMAP") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

p2 <- FeaturePlot(seurat_obj, features = "Totipotent1") + ggtitle("Expression of Totipotent genes")
p3 <- FeaturePlot(seurat_obj, features = "Pluripotent1") + ggtitle("Expression of Pluripotent genes")
p4 <- FeaturePlot(seurat_obj, features = "Maternal1") + ggtitle("Expression of Maternal genes")
p5 <- FeaturePlot(seurat_obj, features = "Minor_ZGA1") + ggtitle("Expression of Minor ZGA genes")
p6 <- FeaturePlot(seurat_obj, features = "Major_ZGA1") + ggtitle("Expression of Major ZGA genes")

pdf(paste0("specific_feature_plot_", method, ".pdf"), width=28, height=14, onefile = F)
(p1 | p2 ) / (p3 | p4 | p5 | p6)
dev.off()
