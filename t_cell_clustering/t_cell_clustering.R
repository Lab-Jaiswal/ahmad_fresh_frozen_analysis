library(AnnotationDbi)
library(BiocSingular)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(S4Vectors)

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
library(magrittr)
library(dplyr)
library(ggplot2)
library(stringr)

add_celltype_labels <- function(sce_object, celltype_labels) {
    sce_metadata <- colData(sce_object) %>% as_tibble()
    celltype_labels_df <- tibble(label = factor(names(celltype_labels)), Celltype = factor(celltype_labels))
    sce_metadata_cell_labels <- left_join(sce_metadata, celltype_labels_df)
    colData(sce_object) <- DataFrame(sce_metadata_cell_labels)
    sce_object[, !is.na(sce_object$Celltype)]
}
filter_osca <- function(sce_object, gene_annot = NULL, filter_cells = TRUE, iterative_filter = FALSE) {
    if (!is.null(gene_annot)) {
        rowData(sce_object) %<>% as.data.frame %>% left_join(gene_annot) %>% DataFrame()
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

# Read data and add gene annotation information
all_samples_sce_tcells_v3 <- read_rds("../all_clustering/all_samples_sce_tcells_v3.rda")
all_samples_sce_tcells_v2 <- read_rds("../all_clustering/all_samples_sce_tcells_v2.rda")

all_samples_sce_filtered_v3 <- filter_osca(all_samples_sce_tcells_v3, filter_cells = FALSE)
all_samples_sce_filtered_v2 <- filter_osca(all_samples_sce_tcells_v2, filter_cells = FALSE)

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

qc_plots(all_samples_sce_leiden_v3, "all_samples_v3/")
qc_plots(all_samples_sce_leiden_v2, "all_samples_v2/")

# nolint start
marker_genes <- c(
    "CD3D", # T cells
    #"CD3E", # T cells
    "CD4", # CD4 T cells
    "CD8A", # CD8 T cells
    "TRDC", # Gamma-delta T-cells, NK cells
    "TRGC2", # Gamma-delta T-cells, NK cells
    "FOXP3", # T-regs
    "GZMB", # Effector T-cells
    "GZMK", # Memory T-cells
    "CX3CR1", # Memory T-cells
    "CD27", # Central Memory T-cells
    "CCR7", # Central Memory T-cells
    "CD44", # Central Memory T-cells
    "SELL", # Central Memory T-cells
    "IL7R", # Central Memory T-cells
    "GNLY", # Central Memory T-cells
    "S100A4", # Central Memory T-cells
    "GZMA", # Central Memory T-cells
    "HLX",
    "IL4",
    "IL5",
    "IL10",
    "IL13",
    "STAT6",
    "GATA3",
    "MAF",

    "NCAM1", # NK cells
    "XCL1" # active NK cells
)
#nolint end

catch <- map(marker_genes, umap_plot, "all_samples_v3/", all_samples_sce_leiden_v3)
catch <- map(marker_genes, umap_plot, "all_samples_v2/", all_samples_sce_leiden_v2)

all_samples_sce_merge <- merge_batches(all_samples_sce_leiden_v3, all_samples_sce_leiden_v2)
write_rds(all_samples_sce_merge, "all_samples_sce_merge.rda")

all_samples_sce_merge_umap <- estimate_umap(all_samples_sce_merge, "corrected")

all_samples_sce_merge_leiden <- leiden_clustering(all_samples_sce_merge_umap, "corrected", resolution = 1.0)
write_rds(all_samples_sce_merge_leiden, "all_samples_sce_merge_leiden.rda")

qc_plots(all_samples_sce_merge_leiden, "all_samples_merge/")

catch <- map(marker_genes, umap_plot, "all_samples_merge/", all_samples_sce_merge_leiden, "reconstructed")

# Find all markers
all_samples_markers <- findMarkers(all_samples_sce_merge_leiden, test = "t", direction = "up", pval.type = "some", assay.type = "reconstructed", row.data = rowData(all_samples_sce_merge_leiden))

# nolint start
# Assign rough cell types
all_samples_celltype_labels <- c(
    "0" = "CD8+ T-cells",
    "1" = "CD8+ T-cells",
    "2" = "CD4+ T-cells",
    "3" = "CD8+ T-cells",
    "4" = "CD4+ T-cells",
    "5" = "CD8+ T-cells",
    "6" = "Active NK cells",
    "7" = "CD4+ T-cells",
    "8" = "CD4+ T-cells",
    "9" = "CD8+ T-cells",
    "10" = "Resting NK cells",
    "11" = "CD4+ T-cells",
    "12" = "CD4+ T-cells",
    "13" = "CD4+ T-cells",
    "14" = "Active NK cells"
)
#nolint end

all_samples_sce_merge_celltypes <- add_celltype_labels(all_samples_sce_merge_leiden, all_samples_celltype_labels)

all_samples_silhouette <- reducedDim(all_samples_sce_merge_celltypes, "corrected") %>%
    approxSilhouette(all_samples_sce_merge_celltypes$Celltype) %>%
    as_tibble()
colData(all_samples_sce_merge_celltypes)$silhouette_width <- all_samples_silhouette$width
write_rds(all_samples_sce_merge_celltypes, "all_samples_sce_merge_celltypes.rda")
