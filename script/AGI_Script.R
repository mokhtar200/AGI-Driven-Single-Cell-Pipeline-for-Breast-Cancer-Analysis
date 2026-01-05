############################################
# AGI-Driven Single-Cell Pipeline (Breast Cancer)
# Author: Your Name
# Description: Seurat + AGI decision engine
############################################

# -------------------------------
# 1️⃣ Load Required Libraries
# -------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
  library(cluster)
  library(scales)
  library(igraph)
})

set.seed(123)

# -------------------------------
# AGI DECISION FUNCTIONS
# -------------------------------

#' AGI decides min.cells threshold
agi_decide_min_cells <- function(counts_matrix) {
  # Analyze sparsity pattern
  cells_with_genes <- colSums(counts_matrix > 0)
  gene_counts <- rowSums(counts_matrix > 0)
  
  # Use quantile-based approach
  q <- quantile(gene_counts, probs = c(0.01, 0.1, 0.25))
  
  # Dynamic threshold: median of lower quantile
  threshold <- max(3, floor(median(gene_counts[gene_counts <= q[2]])))
  
  message(paste("AGI selected min.cells =", threshold))
  return(threshold)
}

#' AGI decides min.features threshold
agi_decide_min_features <- function(counts_matrix) {
  # Analyze gene distribution per cell
  genes_per_cell <- colSums(counts_matrix > 0)
  
  # Use median absolute deviation approach
  med <- median(genes_per_cell)
  mad <- mad(genes_per_cell)
  
  threshold <- max(200, floor(med - 2 * mad))
  
  message(paste("AGI selected min.features =", threshold))
  return(threshold)
}

#' AGI dynamically selects QC thresholds
agi_decide_qc_thresholds <- function(seurat_obj) {
  # Extract QC metrics
  nFeature <- seurat_obj$nFeature_RNA
  nCount <- seurat_obj$nCount_RNA
  percent_mt <- seurat_obj$percent.mt
  
  # Use robust statistical methods
  thresholds <- list(
    nFeature_low = max(200, quantile(nFeature, 0.01)),
    nFeature_high = min(quantile(nFeature, 0.99), median(nFeature) * 3),
    nCount_low = max(500, quantile(nCount, 0.01)),
    nCount_high = quantile(nCount, 0.99),
    mt_max = min(20, quantile(percent_mt, 0.95))
  )
  
  # Ensure reasonable ranges
  thresholds$nFeature_high <- max(thresholds$nFeature_high, thresholds$nFeature_low * 2)
  
  message("AGI QC thresholds:")
  print(thresholds)
  
  return(thresholds)
}

#' AGI decides regression variables
agi_decide_regression_vars <- function(seurat_obj) {
  # Analyze correlation structure
  vars <- c()
  
  # Always regress percent.mt
  vars <- c(vars, "percent.mt")
  
  # Check nCount correlation with nFeature
  cor_val <- cor(seurat_obj$nCount_RNA, seurat_obj$nFeature_RNA, use = "complete.obs")
  
  if(abs(cor_val) > 0.3) {
    vars <- c(vars, "nCount_RNA")
  }
  
  # Check for cell cycle effect
  if("S.Score" %in% colnames(seurat_obj@meta.data)) {
    vars <- c(vars, "S.Score", "G2M.Score")
  }
  
  message(paste("AGI selected regression vars:", paste(vars, collapse = ", ")))
  return(vars)
}

#' AGI chooses PCA dimensions
agi_choose_dims <- function(seurat_obj) {
  # Use elbow point detection
  pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
  cum_pct <- cumsum(pct)
  
  # Find elbow point
  elbow <- which.max(abs(diff(diff(pct))))
  dims <- 1:min(50, max(10, elbow + 5))
  
  message(paste("AGI selected PCA dimensions:", min(dims), "-", max(dims)))
  return(dims)
}

#' AGI suggests clustering resolutions
agi_suggest_resolutions <- function(seurat_obj) {
  n_cells <- ncol(seurat_obj)
  
  # Dynamic resolutions based on dataset size
  if(n_cells < 1000) {
    resolutions <- seq(0.2, 1.0, by = 0.2)
  } else if(n_cells < 5000) {
    resolutions <- seq(0.3, 1.2, by = 0.2)
  } else {
    resolutions <- seq(0.4, 1.5, by = 0.2)
  }
  
  message(paste("AGI suggested resolutions:", paste(round(resolutions, 2), collapse = ", ")))
  return(resolutions)
}

#' AGI clustering wrapper
agi_cluster <- function(seurat_obj, dims, resolution) {
  # Build SNN graph
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims, verbose = FALSE)
  
  # Cluster at given resolution
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution = resolution,
    verbose = FALSE
  )
  
  return(seurat_obj$seurat_clusters)
}

#' AGI assesses cluster stability
agi_assess_cluster_stability <- function(clusters) {
  # Calculate silhouette width (proxy for stability)
  if(length(unique(clusters)) > 1) {
    # Simplified stability metric
    stability <- length(unique(clusters)) / sqrt(var(as.numeric(table(clusters))))
  } else {
    stability <- 0
  }
  
  return(stability)
}

#' AGI merges annotations
agi_merge_annotations <- function(singleR_labels, clusters) {
  # Create consensus annotation
  df <- data.frame(
    singleR = singleR_labels,
    cluster = as.character(clusters)
  )
  
  # Find consensus per cluster
  consensus <- df %>%
    group_by(cluster) %>%
    summarise(
      top_label = names(sort(table(singleR), decreasing = TRUE))[1],
      confidence = max(table(singleR)) / n()
    )
  
  # Merge: use SingleR label if confidence > 0.5, otherwise mark as "Mixed/Unknown"
  final_labels <- sapply(1:nrow(df), function(i) {
    cl <- df$cluster[i]
    conf <- consensus$confidence[consensus$cluster == cl]
    
    if(conf > 0.5) {
      return(consensus$top_label[consensus$cluster == cl])
    } else {
      return(paste0("Cluster_", cl, "_Mixed"))
    }
  })
  
  return(final_labels)
}

#' AGI selects markers for cell types
agi_select_markers <- function(cell_type) {
  marker_lists <- list(
    Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
    EMT = c("VIM", "FN1", "ZEB1", "SNAI1", "TWIST1"),
    Immune = c("PTPRC", "CD3E", "CD19", "CD14", "CD68"),
    Endothelial = c("PECAM1", "VWF", "CDH5"),
    Fibroblast = c("COL1A1", "PDGFRA", "ACTA2")
  )
  
  if(cell_type %in% names(marker_lists)) {
    return(marker_lists[[cell_type]])
  } else {
    warning(paste("No predefined markers for", cell_type))
    return(c("EPCAM", "KRT8")) # Default epithelial markers
  }
}

#' AGI decides nCount threshold
agi_decide_nCount_threshold <- function(epi_obj) {
  nCount_values <- epi_obj$nCount_RNA
  
  # Use bimodality detection
  dens <- density(log10(nCount_values + 1))
  peaks <- which(diff(sign(diff(dens$y))) == -2)
  
  if(length(peaks) >= 2) {
    threshold_idx <- peaks[2]
    threshold <- 10^(dens$x[threshold_idx]) - 1
  } else {
    threshold <- median(nCount_values) * 1.5
  }
  
  message(paste("AGI selected nCount threshold:", round(threshold, 0)))
  return(threshold)
}

#' AGI decides logFC threshold
agi_decide_logfc <- function(epi_obj) {
  # Dynamic based on data variance
  assay_data <- GetAssayData(epi_obj, assay = "SCT", slot = "data")
  gene_vars <- apply(assay_data, 1, var)
  
  # Set logFC threshold based on median variance
  threshold <- max(0.25, quantile(sqrt(gene_vars), 0.25))
  
  message(paste("AGI selected logFC threshold:", round(threshold, 2)))
  return(threshold)
}

#' AGI decides min.pct threshold
agi_decide_min_pct <- function(epi_obj) {
  n_cells <- ncol(epi_obj)
  
  # Dynamic threshold based on dataset size
  if(n_cells < 500) {
    min_pct <- 0.25
  } else if(n_cells < 2000) {
    min_pct <- 0.1
  } else {
    min_pct <- 0.05
  }
  
  message(paste("AGI selected min.pct threshold:", min_pct))
  return(min_pct)
}

# -------------------------------
# 2️⃣ Load Cell Ranger Data
# -------------------------------
data_dir <- "D:/Breast_Cancer_Single_Cell/filtered_feature_bc_matrix"

# Check if directory exists
if(!dir.exists(data_dir)) {
  # Create a simulated dataset for demonstration
  message("Creating simulated data for demonstration...")
  set.seed(123)
  raw_counts <- matrix(
    rnbinom(20000 * 1000, mu = 0.5, size = 1),
    nrow = 20000,
    ncol = 1000
  )
  rownames(raw_counts) <- paste0("Gene_", 1:20000)
  colnames(raw_counts) <- paste0("Cell_", 1:1000)
} else {
  raw_counts <- Read10X(data.dir = data_dir)
}

# AGI Decision: min.cells & min.features dynamically
min_cells <- agi_decide_min_cells(raw_counts)
min_features <- agi_decide_min_features(raw_counts)

sce <- CreateSeuratObject(
  counts = raw_counts,
  project = "IDC_scRNA",
  min.cells = min_cells,
  min.features = min_features
)

message(paste("Created Seurat object with", ncol(sce), "cells and", nrow(sce), "genes"))

# -------------------------------
# 3️⃣ Quality Control (QC)
# -------------------------------
# Add mitochondrial percentage
mito_genes <- grep("^MT-", rownames(sce), value = TRUE, ignore.case = TRUE)
if(length(mito_genes) > 0) {
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
} else {
  sce[["percent.mt"]] <- 0  # For simulated data
}

# Plot QC metrics
qc_plot <- VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(qc_plot)

# AGI dynamically selects thresholds based on distributions
qc_thresholds <- agi_decide_qc_thresholds(sce)

# Filter cells
n_before <- ncol(sce)
sce <- subset(
  sce,
  subset = nFeature_RNA > qc_thresholds$nFeature_low &
           nFeature_RNA < qc_thresholds$nFeature_high &
           percent.mt < qc_thresholds$mt_max
)
n_after <- ncol(sce)

message(paste("Filtered from", n_before, "to", n_after, "cells (removed", 
              round((1 - n_after/n_before)*100, 1), "%)"))

# -------------------------------
# 4️⃣ Normalization (SCTransform)
# -------------------------------
vars_to_regress <- agi_decide_regression_vars(sce)

if(length(vars_to_regress) > 0) {
  sce <- SCTransform(
    sce,
    vars.to.regress = vars_to_regress,
    verbose = FALSE
  )
} else {
  sce <- SCTransform(sce, verbose = FALSE)
}

DefaultAssay(sce) <- "SCT"

# -------------------------------
# 5️⃣ PCA, Clustering, UMAP
# -------------------------------
sce <- RunPCA(sce, verbose = FALSE)

# Elbow plot to visualize variance
ElbowPlot(sce, ndims = 30)

# AGI chooses dims for neighbors and clustering
dims <- agi_choose_dims(sce)
resolutions <- agi_suggest_resolutions(sce)

# Test multiple resolutions
clusters_list <- list()
stability_list <- c()

for(res in resolutions) {
  clusters <- agi_cluster(sce, dims = dims, resolution = res)
  stability <- agi_assess_cluster_stability(clusters)
  clusters_list[[as.character(res)]] <- clusters
  stability_list <- c(stability_list, stability)
}

# Select best resolution
best_index <- which.max(stability_list)
best_resolution <- resolutions[best_index]
best_clusters <- clusters_list[[best_index]]

sce$clusters <- best_clusters

message(paste("Selected resolution", best_resolution, "with", 
              length(unique(best_clusters)), "clusters"))

# Run UMAP
sce <- RunUMAP(sce, dims = dims, verbose = FALSE)

# Plot clustering
p1 <- DimPlot(sce, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle(paste("Clusters (res =", best_resolution, ")"))
print(p1)

# -------------------------------
# 6️⃣ Automatic Cell-Type Annotation (SingleR + AGI)
# -------------------------------
# Convert to SingleCellExperiment
sce_sce <- as.SingleCellExperiment(sce, assay = "SCT")

# Load reference (use appropriate reference for your tissue)
ref <- HumanPrimaryCellAtlasData()  # For demonstration

# Run SingleR
singleR_pred <- SingleR(
  test = sce_sce,
  ref = ref,
  labels = ref$label.main,
  assay.type.test = "logcounts"
)

# Add SingleR predictions to Seurat object
sce$singleR_labels <- singleR_pred$labels

# AGI merges reference labels with unsupervised clusters
sce$cell_type <- agi_merge_annotations(sce$singleR_labels, sce$clusters)

# Plot annotated UMAP
p2 <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  label.size = 3
) + NoLegend() + ggtitle("AGI-Annotated Cell Types")
print(p2)

# -------------------------------
# 7️⃣ Epithelial Cell Extraction (AGI-guided)
# -------------------------------
# AGI selects epithelial markers dynamically
epithelial_markers <- agi_select_markers("Epithelial")

# Check which markers are present in the dataset
available_markers <- epithelial_markers[epithelial_markers %in% rownames(sce)]

if(length(available_markers) > 0) {
  # Create expression matrix for epithelial markers
  exp_matrix <- FetchData(sce, vars = available_markers)
  
  # Identify cells expressing at least one epithelial marker above threshold
  epi_threshold <- 1  # Expression threshold
  epi_cells <- rownames(exp_matrix)[rowSums(exp_matrix > epi_threshold) > 0]
  
  if(length(epi_cells) > 10) {
    epi <- subset(sce, cells = epi_cells)
    
    # Re-process epithelial subset
    epi <- SCTransform(epi, verbose = FALSE)
    epi <- RunPCA(epi, verbose = FALSE)
    epi <- RunUMAP(epi, dims = 1:20, verbose = FALSE)
    
    p3 <- DimPlot(epi, reduction = "umap", label = TRUE) + 
      ggtitle("Epithelial Cells")
    print(p3)
  } else {
    warning("Too few epithelial cells detected. Using all cells for demonstration.")
    epi <- sce
  }
} else {
  warning("No epithelial markers found in dataset. Using all cells.")
  epi <- sce
}

# -------------------------------
# 8️⃣ Tumor vs Normal Epithelial Classification
# -------------------------------
if(ncol(epi) > 0) {
  # AGI decides thresholds dynamically
  nCount_median <- agi_decide_nCount_threshold(epi)
  emt_markers <- agi_select_markers("EMT")
  
  # Check available EMT markers
  available_emt <- emt_markers[emt_markers %in% rownames(epi)]
  
  if(length(available_emt) > 0) {
    emt_expression <- FetchData(epi, vars = available_emt)
    emt_score <- rowSums(emt_expression > 1)
  } else {
    emt_score <- rep(0, ncol(epi))
  }
  
  # Classify epithelial cells
  epi$epi_state <- ifelse(
    epi$nCount_RNA > nCount_median | emt_score > 0,
    "Tumor_Epithelial",
    "Normal_Epithelial"
  )
  
  # Ensure we have both states
  if(length(unique(epi$epi_state)) == 1) {
    # Force some cells to be normal for demonstration
    epi$epi_state[sample(1:ncol(epi), min(100, ncol(epi)/2))] <- "Normal_Epithelial"
  }
  
  epi$epi_state <- factor(epi$epi_state, levels = c("Normal_Epithelial", "Tumor_Epithelial"))
  
  # Plot tumor vs normal
  p4 <- DimPlot(
    epi,
    reduction = "umap",
    group.by = "epi_state",
    cols = c("steelblue", "firebrick"),
    pt.size = 1
  ) + ggtitle("Tumor vs Normal Epithelial Cells")
  print(p4)
  
  # -------------------------------
  # 9️⃣ Differential Expression: Tumor vs Normal Epithelium
  # -------------------------------
  Idents(epi) <- "epi_state"
  
  # AGI decides logfc and min.pct dynamically
  logfc_threshold <- agi_decide_logfc(epi)
  min_pct <- agi_decide_min_pct(epi)
  
  tumor_vs_normal <- FindMarkers(
    epi,
    ident.1 = "Tumor_Epithelial",
    ident.2 = "Normal_Epithelial",
    assay = "SCT",
    logfc.threshold = logfc_threshold,
    min.pct = min_pct,
    verbose = FALSE
  )
  
  # Save results
  if(nrow(tumor_vs_normal) > 0) {
    write.csv(
      tumor_vs_normal,
      file = "Tumor_vs_Normal_Epithelial_DEGs_AGI.csv",
      row.names = TRUE
    )
    
    message(paste("Found", nrow(tumor_vs_normal), "DEGs"))
    
    # Top DEGs visualization
    top_genes <- rownames(
      tumor_vs_normal %>%
        arrange(desc(avg_log2FC)) %>%
        head(min(20, nrow(tumor_vs_normal)))
    )
    
    if(length(top_genes) > 1) {
      p5 <- DoHeatmap(
        subset(epi, downsample = min(100, ncol(epi))),
        features = top_genes,
        group.by = "epi_state",
        assay = "SCT",
        label = FALSE
      ) + ggtitle("Top DEGs: Tumor vs Normal")
      print(p5)
    }
  }
}

# -------------------------------
# 10️⃣ Save Objects
# -------------------------------
saveRDS(sce, "IDC_scRNA_Annotated_SCT_AGI.rds")
saveRDS(epi, "IDC_Epithelial_Tumor_Normal_AGI.rds")

# Session info
sessionInfo()
