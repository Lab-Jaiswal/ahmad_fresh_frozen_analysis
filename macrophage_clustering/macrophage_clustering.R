library(AnnotationHub)
library(AnnotationDbi)
library(BiocSingular)
library(S4Vectors)
library(SummarizedExperiment)
library(SingleCellExperiment)

library(scran)
library(scuttle)
library(leidenAlg)
library(corral)
library(batchelor)
library(PCAtools)
library(uwot)
library(bluster)
library(RNOmni)
library(Seurat)
library(igraph)

library(Cairo)
library(dplyr)
library(ggplot2)
library(magrittr)
library(stringr)

add_celltype_labels <- function(sce_object, celltype_labels) {
    sce_metadata <- colData(sce_object) %>%
        as_tibble()
    celltype_labels_df <- tibble(label = factor(names(celltype_labels)), Celltype = factor(celltype_labels))
    sce_metadata_cell_labels <- left_join(sce_metadata, celltype_labels_df)
    colData(sce_object) <- DataFrame(sce_metadata_cell_labels)
    sce_object
}

filter_osca <- function(sce_object, gene_annot = NULL, filter_cells = TRUE, iterative_filter = FALSE) {
    if (!is.null(gene_annot)) {
        rowData(sce_object) %<>%
            as.data.frame() %>%
            left_join(gene_annot) %>%
            DataFrame()
    }

    if (filter_cells == TRUE) {
        if (iterative_filter == TRUE) {
            sce_object_filtered <- iterative_filter(sce_object)
        } else {
            mito_genes <- which(rowData(sce_object)$SEQNAME == "MT")

            # Compute QC metrics and exclude bad cells
            sce_object_qc <- addPerCellQC(sce_object, subset = list(percent_mt = mito_genes))
            colData(sce_object_qc)$log_sum <- log2(sce_object_qc$sum)
            colData(sce_object_qc)$log_detected <- log2(sce_object_qc$detected)

            low_counts <- isOutlier(sce_object_qc$sum, log = TRUE, type = "lower")
            low_features <- isOutlier(sce_object_qc$detected, log = TRUE, type = "lower")
            high_mt <- isOutlier(sce_object_qc$subsets_percent_mt_percent, type = "higher")
            excluded_cells <- low_counts | low_features | high_mt

            sce_object_filtered <- sce_object_qc[, !excluded_cells]
        }
    } else {
        sce_object_filtered <- sce_object
    }

    # Subset genes to only include ones which meet our variability cutoff
    set.seed(12345L)
    sce_object_gene_var <- modelGeneVarByPoisson(sce_object_filtered)
    sce_object_gene_var_top <- getTopHVGs(sce_object_gene_var)
    sce_object_gene_cv2 <- modelGeneCV2(sce_object_filtered)
    sce_object_gene_cv2_top <- getTopHVGs(sce_object_gene_cv2, var.field = "ratio", var.threshold = 1.0)
    sce_object_gene_shared <- intersect(sce_object_gene_var_top, sce_object_gene_cv2_top)

    rowSubset(sce_object_filtered, "hvg") <- sce_object_gene_shared
    sce_object_filtered
}

estimate_corral <- function(sce_object) {
    standardized_data <- counts(sce_object)[rowSubset(sce_object, "hvg"), ] %>% corral_preproc() %>% apply(1L, RankNorm) %>% t()
    pca_standardized <- pca(standardized_data)
    n_pcs <- chooseGavishDonoho(standardized_data, var.explained = pca_standardized$sdev^2.0, noise = 1.0)
    print(str_c("Chose ", n_pcs, " PCs"))
    reducedDim(sce_object, "corral") <- apply(pca_standardized$rotated[, 1L:n_pcs], 2L, RankNorm)
    sce_object
}

estimate_umap <- function(sce_object, reduced_dim = "corral") {
    # Run UMAP on latent space from corral.  All UMAP parameters are defaults.
    umap_data <- reducedDim(sce_object, reduced_dim) %>% umap()
    colnames(umap_data) <- c("UMAP_1", "UMAP_2")
    reducedDim(sce_object, "umap") <- umap_data
    sce_object
}

leiden_clustering <- function(sce_object, reduced_dim = "corral", k = 10L, resolution = 1.0) {
    # Infer clusters using Leiden algorith on SNN graph.  k = 10 is default parameter
    data_snn <- buildSNNGraph(sce_object, use.dimred = reduced_dim, k = k, type = "jaccard")
    set.seed(12345L)
    data_leiden <- find_partition(data_snn, E(data_snn)$weight, resolution = resolution) # 1.0 is default resolution
    colLabels(sce_object) <- factor(data_leiden)
    sce_object
}

umap_plot <- function(feature_color, prefix, sce_object,
                      data_type = "counts", var_type = "continuous",
                      facet_var = NULL, facet_ncol = NULL, suffix = "", no_legend = FALSE,
                      label_var = NULL, scale_name = NULL, file_name = NULL, plot_name = NULL,
                      width = 12.0, height = 10.0) {
    umap_data <- reducedDim(sce_object, type = "umap", withDimnames = TRUE)
    plot_data <- colData(sce_object) %>%
        as_tibble() %>%
        bind_cols(as_tibble(umap_data))

    if (data_type == "counts" || data_type == "reconstructed") {
        feature_plot <- "log_expr"
        scale_name_final <- expression(paste(log[2.0], "(counts)"))
        plot_name_final <- feature_color
        gene_annot <- rowData(sce_object)
        gene_row <- match(feature_color, gene_annot$Symbol)

        if (data_type == "counts") {
            sce_counts <- counts(sce_object) %>% as.matrix()
            plot_data$log_expr <- log2(sce_counts[gene_row, ])
        } else {
            assay_data <- assays(sce_object)
            reconstructed_data <- as.matrix(assay_data$reconstructed)
            plot_data$log_expr <- reconstructed_data[gene_row, ]
        }
    } else {
        feature_plot <- feature_color
        scale_name_final <- feature_color
        plot_name_final <- feature_color
    }

    if (var_type == "continuous") {
        plot_data %<>% arrange(!!sym(feature_plot)) %>% as.data.frame()
    } else {
        plot_data %<>% arrange(Barcode) %>% as.data.frame()
    }

    umap_plot <- ggplot(plot_data, aes(UMAP_1, UMAP_2, color = !!sym(feature_plot))) +
        geom_point(size = 0.1) +
        ggtitle(feature_color) +
        theme_classic() +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        theme(panel.border = element_rect(fill = NA),
              plot.title = element_text(hjust = 0.5))

    if (!is.null(scale_name)) {
        scale_name_final <- scale_name
    }

    if (var_type == "continuous") {
        umap_plot <- umap_plot + scale_color_gradient(low = "lightgrey", high = "blue", name = scale_name_final)
    }

    if (!is.null(label_var)) {
        umap_plot %<>% LabelClusters(id = label_var, color = "black")
    }

    if (no_legend) {
        umap_plot <- umap_plot + theme(legend.position = "none")
    }

    if (!is.null(facet_var) && !is.null(facet_ncol)) {
        facet_formula <- str_c("~ ", facet_var) %>% as.formula()
        umap_plot <- umap_plot + facet_wrap(facet_formula, ncol = facet_ncol)
    }

    if (!is.null(plot_name)) {
        plot_name_final <- plot_name
        umap_plot <- umap_plot + ggtitle(plot_name_final)
    }

    if (!is.null(file_name)) {
        file_name_final <- str_c(prefix, file_name, suffix)
    } else {
        file_name_final <- str_c(prefix, feature_color, suffix)
    }

    CairoPDF(file_name_final, width = width, height = height)
        print(umap_plot)
    dev.off()
}

qc_boxplot <- function(qc_var, prefix, sce_object, x_var = "label", fill_var = NULL,
                       suffix = "", scale_name = NULL, file_name = NULL, y_label = NULL,
                       angled_text = FALSE, height = 4.0, width = 12.0) {
    plot_data <- colData(sce_object) %>% as_tibble()
    if (!is.null(fill_var)) {
        box_plot_aes <- aes(!!sym(x_var), !!sym(qc_var), fill = !!sym(fill_var))
    } else {
        box_plot_aes <- aes(!!sym(x_var), !!sym(qc_var))
    }

    box_plot <- ggplot(plot_data, box_plot_aes) +
      geom_boxplot() +
      theme_classic() +
      theme(axis.title.x = element_blank(),
            panel.border = element_rect(color = "black", fill = NA)
      )

    if (angled_text) {
        box_plot <- box_plot + theme(axis.text.x = element_text(angle = 45.0, hjust = 1.0))
    }

    if (!is.null(y_label)) {
        box_plot <- box_plot + ylab(y_label)
    }

    if (!is.null(scale_name)) {
        scale_name_final <- scale_name
        box_plot <- box_plot + scale_fill_discrete(name = scale_name_final)
    }

    if (!is.null(file_name)) {
        file_name_final <- str_c(prefix, file_name, suffix, "_boxplot")
    } else {
        file_name_final <- str_c(prefix, qc_var, suffix, "_boxplot")
    }

    CairoPDF(file_name_final, height = height, width = width)
        print(box_plot)
    dev.off()
}


qc_umap_plots <- function(prefix, sce_object, qc_var, ...) {
    umap_plot(qc_var, prefix, sce_object, data_type = "metadata", ...)
    umap_plot(qc_var, prefix, sce_object, data_type = "metadata",
              width = 24.0, height = 15.0, facet_var = "Sample", facet_ncol = 4L, suffix = "_sample_names", ...)
    umap_plot(qc_var, prefix, sce_object, data_type = "metadata",
              width = 24.0, height = 15.0, facet_var = "Status", facet_ncol = 2L, suffix = "_status", ...)
}

qc_plots <- function(sce_object, prefix) {
    qc_umap_plots(prefix, sce_object, "label", var_type = "discrete", scale_name = "Cluster", file_name = "clusters", label_var = "label", plot_name = "", no_legend = TRUE)
    umap_plot("Sample", prefix, sce_object, data_type = "metadata", var_type = "discrete",
              width = 24.0, height = 15.0, facet_var = "Sample", facet_ncol = 4L,
              file_name = "sample_names", no_legend = TRUE, plot_name = "Sample")
    umap_plot("Status", prefix, sce_object, data_type = "metadata", var_type = "discrete",
              width = 24.0, height = 15.0, facet_var = "Status", facet_ncol = 2L,
              file_name = "status", no_legend = TRUE, plot_name = "Status")

    qc_umap_plots(prefix, sce_object, "doublet_scores", plot_name = "")
    qc_umap_plots(prefix, sce_object, "subsets_percent_mt_percent", scale_name = "% MT", file_name = "percent_mt", plot_name = "")
    qc_umap_plots(prefix, sce_object, "log_sum", plot_name = "")
    qc_umap_plots(prefix, sce_object, "log_detected", plot_name = "")

    qc_boxplot("doublet_scores", prefix, sce_object, y_label = "Doublet scores")
    qc_boxplot("doublet_scores", prefix, sce_object, fill_var = "Status", suffix = "_status", y_label = "Doublet scores")
}

npc_donoho <- function(original_data) {
    standardized_data <- apply(original_data, 1L, RankNorm) %>% t()
    pca_standardized <- pca(standardized_data)
    n_pcs <- chooseGavishDonoho(standardized_data, var.explained = pca_standardized$sdev^2.0, noise = 1.0)
    print(str_c("Chose ", n_pcs, " PCs"))
    n_pcs
}

merge_batches <- function(sce_object1, sce_object2) {
    sce_object1_hvgs <- rownames(sce_object1)[rowSubset(sce_object1, "hvg")]
    sce_object2_hvgs <- rownames(sce_object2)[rowSubset(sce_object2, "hvg")]
    shared_genes <- intersect(sce_object1_hvgs, sce_object2_hvgs)

    sce_object1 <- sce_object1[shared_genes, ]
    sce_object2 <- sce_object2[shared_genes, ]

    sce_batch_norm <- multiBatchNorm(sce_object1, sce_object2)

    npcs_batch1 <- ncol(reducedDim(sce_object1, type = "corral"))
    npcs_batch2 <- ncol(reducedDim(sce_object2, type = "corral"))

    npcs_merge <- max(npcs_batch1, npcs_batch2)

    sce_merge <- fastMNN(sce_batch_norm, d = npcs_merge, cos.norm = FALSE, BSPARAM = ExactParam())
    rowData(sce_merge) <- cbind(rowData(sce_object1), rowData(sce_merge))

    original_metadata <- bind_rows(select(as_tibble(colData(sce_object1)), Sample:doublet_scores),
                                   select(as_tibble(colData(sce_object2)), Sample:doublet_scores))
    merged_metadata <- colData(sce_merge) %>%
        as_tibble() %>%
        bind_cols(original_metadata) %>%
        DataFrame()
    colData(sce_merge) <- merged_metadata

    sce_merge
}

merge_counts <- function(sce_object1, sce_object2, sce_merged) {
    shared_genes <- intersect(rownames(sce_object1), rownames(sce_object2))

    sce_object1 <- sce_object1[shared_genes, is_in(colData(sce_object1)$Cell_ID, colData(sce_merged)$Cell_ID)]
    sce_object1_rowdata <- rowData(sce_object1) %>% as_tibble()
    sce_object1_rowdata_select <- select(sce_object1_rowdata, -hvg)
    rowData(sce_object1) <- DataFrame(sce_object1_rowdata_select)
    sce_object1_coldata <- colData(sce_object1) %>% as_tibble() %>% select(Cell_ID) %>% DataFrame()
    colData(sce_object1) <- sce_object1_coldata
    reducedDim(sce_object1, "corral") <- NULL

    sce_object2 <- sce_object2[shared_genes, is_in(colData(sce_object2)$Cell_ID, colData(sce_merged)$Cell_ID)]
    sce_object2_rowdata <- rowData(sce_object2) %>% as_tibble()
    sce_object2_rowdata_select <- select(sce_object2_rowdata, -hvg)
    rowData(sce_object2) <- DataFrame(sce_object2_rowdata_select)
    sce_object2_coldata <- colData(sce_object2) %>% as_tibble() %>% select(Cell_ID) %>% DataFrame()
    colData(sce_object2) <- sce_object2_coldata
    reducedDim(sce_object2, "corral") <- NULL

    shared_object <- cbind(sce_object1, sce_object2)
    shared_object_subset <- shared_object[, is_in(colData(shared_object)$Cell_ID, colData(sce_merged)$Cell_ID)]
    merged_metadata <- colData(sce_merged) %>% as_tibble()
    colData(shared_object_subset) %<>% as_tibble %>% left_join(merged_metadata) %>% DataFrame()
    shared_object
}

pseudobulk_counts <- function(sce_object, min_ncells = 10L) {
    celltype_metadata <- colData(sce_object) %>% as_tibble() %>% select(Sample, Celltype) %>% DataFrame()
    sce_pseudobulk <- aggregateAcrossCells(sce_object, id = celltype_metadata)
    sce_pseudobulk_filter <- sce_pseudobulk[, sce_pseudobulk$ncells >= min_ncells]
    sce_pseudobulk_filter
}

edger_bulk_dge <- function(pseudobulk_object) {
    edger_dge <- pseudoBulkDGE(pseudobulk_object,
        label = pseudobulk_object$Celltype,
        design = ~ batch + Status,
        coef = "StatusFrozen",
        condition = pseudobulk_object$Status,
        row.data = rowData(pseudobulk_object),
        method = "edgeR"
    )
    edger_dge
}

limma_dge <- function(pseudobulk_object) {
    limma_dge <- pseudoBulkDGE(pseudobulk_object,
        label = pseudobulk_object$Celltype,
        design = ~ batch + Status,
        coef = "StatusFrozen",
        row.data = rowData(pseudobulk_object),
        method = "voom"
    )
    limma_dge
}

# Read data and add gene annotation information
all_samples_sce_macros_v3 <- read_rds("../all_clustering/all_samples_sce_macrophages_v3.rda")
all_samples_sce_macros_v2 <- read_rds("../all_clustering/all_samples_sce_macrophages_v2.rda")

annotation_hub <- AnnotationHub()
ens_104 <- annotation_hub[["AH95744"]]
gene_annot <- AnnotationDbi::select(ens_104, keys = rownames(all_samples_sce_macros_v3), keytype = "GENEID", columns = c("GENEID", "GENEBIOTYPE", "DESCRIPTION", "SEQNAME", "GENESEQSTART", "GENESEQEND", "SEQSTRAND", "CANONICALTRANSCRIPT", "GENEIDVERSION"))
colnames(gene_annot)[1L] <- "ID"

all_samples_sce_filtered_v3 <- filter_osca(all_samples_sce_macros_v3, filter_cells = FALSE)
all_samples_sce_filtered_v2 <- filter_osca(all_samples_sce_macros_v2, filter_cells = FALSE)

all_samples_sce_corral_v3 <- estimate_corral(all_samples_sce_filtered_v3)
all_samples_sce_corral_v2 <- estimate_corral(all_samples_sce_filtered_v2)
write_rds(all_samples_sce_corral_v3, "all_samples_sce_corral_v3.rda")
write_rds(all_samples_sce_corral_v2, "all_samples_sce_corral_v2.rda")

all_samples_sce_umap_v3 <- estimate_umap(all_samples_sce_corral_v3)
all_samples_sce_umap_v2 <- estimate_umap(all_samples_sce_corral_v2)

all_samples_sce_leiden_v3 <- leiden_clustering(all_samples_sce_umap_v3)
all_samples_sce_leiden_v2 <- leiden_clustering(all_samples_sce_umap_v2)
write_rds(all_samples_sce_leiden_v3, "all_samples_sce_leiden_v3.rda")
write_rds(all_samples_sce_leiden_v2, "all_samples_sce_leiden_v2.rda")

all_samples_sce_macro_v3 <- all_samples_sce_leiden_v3[, !is_in(all_samples_sce_leiden_v3$label, c(11L, 13L, 14L, 15L))]
all_samples_sce_macro_v2 <- all_samples_sce_leiden_v2[, !is_in(all_samples_sce_leiden_v2$label, c(8L, 14L))]

qc_plots(all_samples_sce_leiden_v3, "all_samples_v3/")
qc_plots(all_samples_sce_leiden_v2, "all_samples_v2/")

# nolint start
# Plot marker genes
marker_genes <- c(
    "C3", # macrophages
    "CD68", # macrophages
    "CD14", # macrophages/monocytes
    "FCGR3A", # monocytes

    "S100A8", # inflammatory macrophages
    "IL1B", # inflammatory macrophages
    "CD300E", # inflammatory macrophages
    "CSF3R", # inflammatory macrophages
    "NFKB1", # inflammatory macrophages

    "APOE", # foam cells
    "TREM2", # foam cells
    "LPL", # foam cells
    "ABCA1", # foam cells
    "FABP5", # foam cells
    "FTL", # foam cells
    "PLTP", # foam cells

    "LYVE1", # TR macrophages
    "MRC1", # TR macrophages
    "FOLR2", # TR macrophages

    "FOXP3",
    "GZMB",
    "ACOT4",

    "HLA-DRB1", # MHC-hi macrophages
    "HLA-DQB1", # MHC-hi macrophages

    "FCER1A", # classical DCs
    "PLD4", # plasmacytoid DCs
    "XCR1", # XCR1+ DCs

    "FCGR3B", # neutrophils

    "MKI67", # Cycling cells
    "CD3E", # monocytes
    "CD79A" # monocytes
)
# nolint end

catch <- map(marker_genes, umap_plot, "all_samples_v3/", all_samples_sce_leiden_v3)
catch <- map(marker_genes, umap_plot, "all_samples_v2/", all_samples_sce_leiden_v2)

all_samples_sce_merge <- merge_batches(all_samples_sce_macro_v3, all_samples_sce_macro_v2)
write_rds(all_samples_sce_merge, "all_samples_sce_merge.rda")

all_samples_sce_merge_umap <- estimate_umap(all_samples_sce_merge, "corrected")

all_samples_sce_merge_leiden <- leiden_clustering(all_samples_sce_merge_umap, "corrected", resolution = 1.0)

qc_plots(all_samples_sce_merge_leiden, "all_samples_merge/")

catch <- map(marker_genes, umap_plot, "all_samples_merge/", all_samples_sce_merge_leiden, "reconstructed")
umap_plot("MORN4", "all_samples_merge/", all_samples_sce_merge_leiden, "reconstructed")

all_samples_sce_merge_macros <- all_samples_sce_merge[, !is_in(all_samples_sce_merge_leiden$label, c(10L, 11L, 13L))]
all_samples_sce_merge_macros_umap <- estimate_umap(all_samples_sce_merge_macros, "corrected")
all_samples_sce_merge_macros_leiden <- leiden_clustering(all_samples_sce_merge_macros_umap, "corrected", resolution = 1.0)
qc_plots(all_samples_sce_merge_macros_leiden, "all_samples_merge_macros/")
catch <- map(marker_genes, umap_plot, "all_samples_merge_macros/", all_samples_sce_merge_macros_leiden, "reconstructed")

# Find all markers
all_samples_markers <- findMarkers(all_samples_sce_merge_leiden, test = "t", direction = "up", pval.type = "some", assay.type = "reconstructed", row.data = rowData(all_samples_sce_merge_leiden))

# Assign rough cell types
all_samples_celltype_labels <- c(
    "0" = "LYVE1+ TR macrophages",
    "1" = "Inflammatory macrophages",
    "2" = "Inflammatory macrophages",
    "3" = "MHC-hi macrophages",
    "4" = "CD16+ monocytes",
    "5" = "cDCs",
    "6" = "Foam cells",
    "7" = "Foam cells",
    "8" = "Foam cells",
    "9" = "MHC-hi macrophages",
    "10" = "Inflammatory macrophages"
)

all_samples_sce_merge_macros_celltypes <- add_celltype_labels(all_samples_sce_merge_macros_leiden, all_samples_celltype_labels)

all_samples_silhouette <- reducedDim(all_samples_sce_merge_macros_celltypes, "corrected") %>%
    approxSilhouette(all_samples_sce_merge_macros_celltypes$Celltype) %>%
    as_tibble()
colData(all_samples_sce_merge_macros_celltypes)$silhouette_width <- all_samples_silhouette$width
write_rds(all_samples_sce_merge_macros_celltypes, "all_samples_sce_merge_celltypes.rda")

all_samples_sce_counts <- merge_counts(all_samples_sce_filtered_v3, all_samples_sce_filtered_v2, all_samples_sce_merge_celltypes)

all_samples_pseudobulk <- pseudobulk_counts(all_samples_sce_counts)

all_samples_edger_dge <- edger_bulk_dge(all_samples_pseudobulk)

all_samples_limma_dge <- limma_dge(all_samples_pseudobulk)
all_samples_limma_tibble <- as.list(all_samples_limma_dge) %>% map(as_tibble)
write_rds(all_samples_limma_tibble, "all_samples_limma_tibble.rda")
