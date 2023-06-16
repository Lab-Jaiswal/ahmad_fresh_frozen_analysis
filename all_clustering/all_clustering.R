library(AnnotationHub)
library(AnnotationDbi)
library(BiocSingular)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(S4Vectors)
library(ensembldb)

library(DropletUtils)
library(scran)
library(scuttle)
library(scDblFinder)
library(celldex)
library(SingleR)
library(leidenAlg)
library(igraph)
library(corral)
library(batchelor)
library(PCAtools)
library(uwot)
library(bluster)
library(RNOmni)
library(VGAM)
library(Seurat)
library(broom)
library(forestplot)

library(tictoc)
library(Cairo)
library(magrittr)
library(dplyr)
library(ggplot2)
library(rlang)
library(readr)
library(stringr)
library(tidyr)
library(tidyverse)

umap_plot <- function(feature_color, prefix, sce_object,
                      data_type = "counts", var_type = "continuous",
                      facet_var = NULL, facet_ncol = NULL, suffix = "", no_legend = FALSE,
                      label_var = NULL, scale_name = NULL, file_name = NULL, plot_name = NULL,
                      width = 12.0, height = 10.0) {
    umap_data <- reducedDim(sce_object, type = "umap", withDimnames = TRUE)
    plot_data <- colData(sce_object) %>% as_tibble() %>% bind_cols(as_tibble(umap_data))

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

filter_osca <- function(sce_object, gene_annot, filter_cells = TRUE) {
    rowData(sce_object) %<>% as.data.frame %>% left_join(gene_annot) %>% DataFrame()
    mito_genes <- which(rowData(sce_object)$SEQNAME == "MT")

    # Compute QC metrics and exclude bad cells
    sce_object_qc <- addPerCellQC(sce_object, subset = list(percent_mt = mito_genes))
    colData(sce_object_qc)$log_sum <- log2(sce_object_qc$sum)
    colData(sce_object_qc)$log_detected <- log2(sce_object_qc$detected)

    if (filter_cells == TRUE) {
        low_counts <- isOutlier(sce_object_qc$sum, log = TRUE, type = "lower")
        low_features <- isOutlier(sce_object_qc$detected, log = TRUE, type = "lower")
        high_mt <- isOutlier(sce_object_qc$subsets_percent_mt_percent, type = "higher")
        excluded_cells <- low_counts | low_features | high_mt

        sce_object_filtered <- sce_object_qc[, !excluded_cells]
    } else {
        sce_object_filtered <- sce_object_qc
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

normalize_osca <- function(sce_object) {
    size_factors <- librarySizeFactors(sce_object)
    size_factors_norm <- size_factors / mean(size_factors)
    norm_counts <- normalizeCounts(sce_object, size_factors_norm, center.size.factors = FALSE)
    normcounts(sce_object) <- norm_counts
    sce_object <- logNormCounts(sce_object)
    sce_object
}

npc_donoho <- function(original_data) {
    pca_scaled <- normcounts(original_data) %>% t() %>% scater::runPCA(center = TRUE, BSPARAM = ExactParam(), rank = nrow(original_data))
    n_pcs <- chooseGavishDonoho(original_data, var.explained = pca_scaled$sdev^2.0, noise = 1.0)
    n_pcs
}

# Compute doublet scores using number of PCs identified above.  k = 50 is default parameter
compute_doublets <- function(sce_object, n_neighbors = 50L) {
    sce_object_hvg <- sce_object[rowSubset(sce_object, "hvg"), ]
    npc_normcounts <- normalize_osca(sce_object_hvg) %>% npc_donoho()
    print(str_c("Chose ", npc_normcounts, " PCs"))
    doublet_scores <- computeDoubletDensity(sce_object_hvg, k = n_neighbors, dims = npc_normcounts, BSPARAM = ExactParam())
    colData(sce_object)$doublet_scores <- doublet_scores
    sce_object
}

estimate_corral <- function(sce_object) {
    standardized_data <- counts(sce_object)[rowSubset(sce_object, "hvg"), ] %>% corral_preproc() %>% apply(1L, RankNorm) %>% t()
    pca_standardized <- pca(standardized_data)
    n_pcs <- chooseGavishDonoho(standardized_data, var.explained = pca_standardized$sdev^2.0, noise = 1.0)
    print(str_c("Chose ", n_pcs, " PCs"))
    reducedDim(sce_object, "corral") <- apply(pca_standardized$rotated[, 1L:n_pcs], 2L, RankNorm)
    sce_object
}

estimate_umap <- function(sce_object, reduced_dim = "corral", normalize_dimred = FALSE) {
    # Run UMAP on latent space from corral.  All UMAP parameters are defaults.
    if (normalize_dimred) {
        normalized_dimred <- reducedDim(sce_object, reduced_dim) %>% apply(2L, RankNorm)
    } else {
        normalized_dimred <- reducedDim(sce_object, reduced_dim)
    }
    umap_data <- umap(normalized_dimred)
    set.seed(12345L)
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

get_hvgs <- function(sce_object) {
    sce_object_hvgs <- rownames(sce_object)[rowSubset(sce_object, "hvg")]
    sce_object_hvgs
}

subset_by_hvg <- function(sce_object, shared_genes) {
    sce_object <- sce_object[shared_genes, ]
    sce_object
}

extract_metadata <- function(sce_object) {
    select(as_tibble(colData(sce_object)), Sample:doublet_scores)
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

    set.seed(12345L)
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

merge_samples <- function(sce_object, batch_column) {
    tic("Batch normalization")
    sce_batch_norm <- multiBatchNorm(sce_object, batch = colData(sce_object)[[batch_column]], normalize.all = TRUE, subset.row = rowData(sce_object)$hvg)
    toc()

    set.seed(12345L)
    tic("Batch adjustment")
    sce_merge <- fastMNN(sce_batch_norm,
                         d = ncol(reducedDim(sce_object, "corral")),
                         cos.norm = FALSE, BSPARAM = ExactParam(),
                         batch = colData(sce_object)[[batch_column]],
                         subset.row = rowData(sce_object)$hvg,
                         correct.all = TRUE)
    rowData(sce_merge) <- cbind(rowData(sce_object), rowData(sce_merge))
    rm(sce_batch_norm)
    gc()
    toc()

    original_metadata <- colData(sce_object) %>% as_tibble() %>% select(Sample:doublet_scores)
    merged_metadata <- colData(sce_merge) %>%
        as_tibble() %>%
        bind_cols(original_metadata) %>%
        DataFrame()
    colData(sce_merge) <- merged_metadata
    counts(sce_merge) <- counts(sce_object)
    logcounts(sce_merge) <- logcounts(sce_object)

    sce_merge
}

merge_metadata <- function(sce_object1, sce_object2, sce_object_merge) {
    original_metadata <- bind_rows(select(as_tibble(colData(sce_object1)), Sample:doublet_scores),
                                   select(as_tibble(colData(sce_object2)), Sample:doublet_scores))
    umap_data <- reducedDim(sce_object_merge, "umap") %>% as_tibble()
    merged_metadata <- colData(sce_object_merge) %>%
        as_tibble() %>%
        bind_cols(umap_data) %>%
        bind_cols(original_metadata) %>%
        data.frame()
    merged_metadata
}

singler_annotation <- function(sce_object, sce_bulk, label_name, filename) {
    merge_assay_data <- assays(sce_object)
    reconstructed_data <- merge_assay_data$reconstructed
    rownames(reconstructed_data) <- rowData(sce_object)$Symbol

    singler_labels <- SingleR(reconstructed_data, sce_bulk, labels = sce_bulk$label.fine)
    singler_labels_df <- as_tibble(singler_labels, rownames = "Barcode")
    colData(sce_object)[[label_name]] <- factor(singler_labels_df$pruned.labels)
    umap_data <- reducedDim(sce_object, "umap") %>% as_tibble() %>% set_colnames(c("UMAP_1", "UMAP_2"))

    sce_metadata <- colData(sce_object) %>% as_tibble() %>% bind_cols(umap_data) %>% data.frame()
    singler_ggplot <- ggplot(sce_metadata, aes(UMAP_1, UMAP_2, color = !!sym(label_name))) +
        geom_point(size = 0.1) +
        theme_classic() +
        theme(legend.position = "none",
              panel.border = element_rect(fill = NA))

    singler_ggplot %<>% LabelClusters(id = label_name, color = "black")

    CairoPDF(filename, width = 12.0, height = 10.0)
        print(singler_ggplot)
    dev.off()

    sce_object
}

celltype_plots <- function(sce_object, prefix) {
    qc_umap_plots(prefix, sce_object, "Celltype", var_type = "discrete", label_var = "Celltype", no_legend = TRUE, file_name = "celltypes", plot_name = "")

    qc_boxplot("doublet_scores", prefix, sce_object, fill_var = "Status", suffix = "_celltypes", x_var = "Celltype", angled_text = TRUE, y_label = "Doublet scores")
    qc_boxplot("subsets_percent_mt_percent", prefix, sce_object, fill_var = "Status", suffix = "_celltypes", x_var = "Celltype", angled_text = TRUE, y_label = "% MT", file_name = "percent_mt")
    qc_boxplot("log_sum", prefix, sce_object, fill_var = "Status", suffix = "_celltypes", x_var = "Celltype", angled_text = TRUE, y_label = expression(paste(log[2.0], "(total counts)")))
    qc_boxplot("log_detected", prefix, sce_object, fill_var = "Status", suffix = "_celltypes", x_var = "Celltype", angled_text = TRUE, y_label = expression(paste(log[2.0], "(detected_genes)")))
    qc_boxplot("silhouette_width", prefix, sce_object, fill_var = "Status", suffix = "_celltypes", x_var = "Celltype", angled_text = TRUE, y_label = "Silhouette width")
}

add_celltype_labels <- function(sce_object, celltype_labels, celltype_column = "Celltype") {
    sce_metadata <- colData(sce_object) %>% as_tibble()
    celltype_labels_df <- tibble(label = factor(names(celltype_labels)), {{celltype_column}} := factor(celltype_labels))
    sce_metadata_cell_labels <- left_join(sce_metadata, celltype_labels_df)
    colData(sce_object) <- DataFrame(sce_metadata_cell_labels)
    sce_object
}

cluster_lm <- function(metadata_df, celltype_df) {
    metadata_df$norm_metric <- RankNorm(metadata_df$metric)

    metadata_lm <- lm(norm_metric ~ batch + Status, metadata_df)
    metadata_summary <- tidy(metadata_lm, conf.int = TRUE)
    metadata_summary
}

# Fit models for metrics
cluster_metrics <- function(sce_object) {
    metrics_long <- colData(sce_object) %>%
        as_tibble() %>%
        select(Barcode, batch, Celltype, Status, sum, detected, subsets_percent_mt_percent, silhouette_width) %>%
        pivot_longer(sum:silhouette_width, names_to = "metric_name", values_to = "metric")

    metrics_lm <- group_by(metrics_long, Celltype, metric_name) %>% group_modify(cluster_lm)
    metrics_lm
}

cluster_vglm <- function(metadata_df, celltype_df) {
    mt_percent_vglm <- vglm(cbind(subsets_percent_mt_sum, sum - subsets_percent_mt_sum) ~ batch + Status, family = betabinomial(lrho = "rhobitlink"), data = metadata_df, trace = TRUE)
    mt_percent_summary <- summary(mt_percent_vglm, HDEtest = FALSE)
    mt_percent_summary_df <- as_tibble(mt_percent_summary@coef3, rownames = "term")
    mt_percent_summary_df$metric_name <- "subsets_percent_mt_percent"

    silhouette_vglm <- vglm(silhouette_width ~ batch + Status, family = uninormal(lmean = fisherzlink()), data = metadata_df, trace = TRUE)
    silhouette_vglm_summary <- summary(silhouette_vglm, HDEtest = FALSE)
    silhouette_vglm_summary_df <- as_tibble(silhouette_vglm_summary@coef3, rownames = "term")
    silhouette_vglm_summary_df$metric_name <- "silhouette_width"

    sum_vglm <- vglm(sum ~ batch + Status, family = negbinomial(), data = metadata_df, trace = TRUE)
    sum_vglm_summary <- summary(sum_vglm, HDEtest = FALSE)
    sum_vglm_summary_df <- as_tibble(sum_vglm_summary@coef3, rownames = "term")
    sum_vglm_summary_df$metric_name <- "sum"

    detected_vglm <- vglm(detected ~ batch + Status, family = negbinomial(), data = metadata_df, trace = TRUE)
    detected_vglm_summary <- summary(detected_vglm, HDEtest = FALSE)
    detected_vglm_summary_df <- as_tibble(detected_vglm_summary@coef3, rownames = "term")
    detected_vglm_summary_df$metric_name <- "detected"

    metadata_summary <- bind_rows(mt_percent_summary_df, silhouette_vglm_summary_df, sum_vglm_summary_df, detected_vglm_summary_df)
    metadata_summary$CI95lo <- metadata_summary$Estimate - (1.96 * metadata_summary$`Std. Error`)
    metadata_summary$CI95hi <- metadata_summary$Estimate + (1.96 * metadata_summary$`Std. Error`)
    metadata_summary$OR <- exp(metadata_summary$Estimate)
    metadata_summary$OR_CI95hi <- exp(metadata_summary$CI95hi)
    metadata_summary$OR_CI95lo <- exp(metadata_summary$CI95lo)
    metadata_summary
}

cluster_metrics_vglm <- function(sce_object) {
    metrics_long <- colData(sce_object) %>%
        as_tibble() %>%
        select(Barcode, batch, Celltype, Status, sum, detected, subsets_percent_mt_sum, silhouette_width)
    metrics_vglm <- group_by(metrics_long, Celltype) %>% group_modify(cluster_vglm)
    metrics_vglm
}

cluster_lm_adj <- function(metadata_df, celltype_df) {
    metadata_lm <- lm(metric ~ Status, metadata_df)
    metadata_summary <- tidy(metadata_lm, conf.int = TRUE)
    metadata_summary
}

cluster_metrics_adj <- function(sce_object) {
    metrics <- colData(sce_object) %>% as_tibble()
    metrics$subsets_percent_mt_percent_adj <- (((metrics$subsets_percent_mt_percent / 100.0) * (nrow(metrics) - 1L)) + 0.5) / nrow(metrics)
    metrics$logit_percent_mt <- logitlink(metrics$subsets_percent_mt_percent_adj)
    metrics$logit_percent_mt_adj <- lm(logit_percent_mt ~ batch, metrics) %>% residuals()

    metrics$z_silhouette_width <- fisherzlink(metrics$silhouette_width)
    metrics$z_silhouette_width_adj <- lm(z_silhouette_width ~ batch, metrics) %>% residuals()

    metrics$log_sum_adj <- lm(log_sum ~ batch, metrics) %>% residuals()

    metrics$log_detected_adj <- lm(log_detected ~ batch, metrics) %>% residuals()

    metrics_long <- select(metrics, Barcode, Celltype, Status, log_sum_adj, log_detected_adj, logit_percent_mt_adj, z_silhouette_width_adj) %>%
        pivot_longer(log_sum_adj:z_silhouette_width_adj, names_to = "metric_name", values_to = "metric")

    metrics_lm <- group_by(metrics_long, Celltype, metric_name) %>% group_modify(cluster_lm_adj)
    metrics_lm$OR <- exp(metrics_lm$estimate)
    metrics_lm$OR_CI95hi <- exp(metrics_lm$conf.high)
    metrics_lm$OR_CI95lo <- exp(metrics_lm$conf.low)
    metrics_lm
}

merge_counts <- function(sce_object1, sce_object2, sce_merged) {
    shared_genes <- intersect(rownames(sce_object1), rownames(sce_object2))

    sce_object1 <- sce_object1[shared_genes, ]
    sce_object1_metadata <- rowData(sce_object1) %>% as_tibble()
    sce_object1_metadata_select <- select(sce_object1_metadata, -hvg)
    rowData(sce_object1) <- DataFrame(sce_object1_metadata_select)
    reducedDim(sce_object1, "corral") <- NULL
    sce_object2 <- sce_object2[shared_genes, ]
    sce_object2_metadata <- rowData(sce_object2) %>% as_tibble()
    sce_object2_metadata_select <- select(sce_object2_metadata, -hvg)
    rowData(sce_object2) <- DataFrame(sce_object2_metadata_select)
    reducedDim(sce_object2, "corral") <- NULL

    shared_object <- cbind(sce_object1, sce_object2)
    shared_object
}

subset_object <- function(shared_object, sce_merged) {
    shared_object_subset <- shared_object[, is_in(colData(shared_object)$Cell_ID, colData(sce_merged)$Cell_ID)]
    merged_metadata <- colData(sce_merged) %>% as_tibble()
    colData(shared_object_subset) %<>% as_tibble %>% select(Cell_ID) %>% left_join(merged_metadata) %>% DataFrame()
    shared_object_subset
}

avg_silhouette <- function(metadata_df, sample_df) {
    metadata_df$Z_silhouette_width <- fisherzlink(metadata_df$silhouette_width)
    mean_silhouette_width <- mean(metadata_df$Z_silhouette_width) %>% fisherzlink(inverse = TRUE)
    tibble(mean_silhouette_width = mean_silhouette_width)
}

combine_detected <- function(metadata_df, sample_df, sce_counts) {
    sample_counts <- counts(sce_counts[, is_in(sce_counts$Cell_ID, metadata_df$Cell_ID)]) %>% as.matrix()
    sample_nonzero <- apply(sample_counts > 0L, 1L, any) %>% which() %>% length()
    tibble(detected = sample_nonzero)
}

cluster_vglm_agg <- function(metadata_df, celltype_df, sce_object) {
    metadata_reduce <- select(metadata_df, Sample, batch, Status) %>% distinct()
    mt_percent_celltype <- group_by(metadata_df, Sample) %>% summarise(sum = sum(sum), percent_mt_sum = sum(subsets_percent_mt_sum)) %>% left_join(metadata_reduce)
    mt_percent_vglm <- vglm(cbind(percent_mt_sum, sum - percent_mt_sum) ~ batch + Status, family = betabinomial(lrho = "rhobitlink"), data = mt_percent_celltype, trace = TRUE)
    mt_percent_summary <- summary(mt_percent_vglm, HDEtest = FALSE)
    mt_percent_summary_df <- as_tibble(mt_percent_summary@coef3, rownames = "term")
    mt_percent_summary_df$metric_name <- "subsets_percent_mt_percent"

    silhouette_celltype <- group_by(metadata_df, Sample) %>% group_modify(avg_silhouette) %>% left_join(metadata_reduce)
    silhouette_vglm <- vglm(mean_silhouette_width ~ batch + Status, family = uninormal(lmean = fisherzlink()), data = silhouette_celltype, trace = TRUE)
    silhouette_vglm_summary <- summary(silhouette_vglm, HDEtest = FALSE)
    silhouette_vglm_summary_df <- as_tibble(silhouette_vglm_summary@coef3, rownames = "term")
    silhouette_vglm_summary_df$metric_name <- "silhouette_width"

    detected_celltype <- group_by(metadata_df, Sample) %>% group_modify(combine_detected, sce_object) %>% left_join(metadata_reduce)
    detected_vglm <- vglm(detected ~ batch + Status, family = negbinomial(), data = detected_celltype, trace = TRUE)
    detected_vglm_summary <- summary(detected_vglm, HDEtest = FALSE)
    detected_vglm_summary_df <- as_tibble(detected_vglm_summary@coef3, rownames = "term")
    detected_vglm_summary_df$metric_name <- "detected"

    sum_celltype <- group_by(metadata_df, Sample) %>% summarise(sum = sum(sum)) %>% left_join(metadata_reduce)
    sum_vglm <- vglm(sum ~ batch + Status, family = negbinomial(), data = sum_celltype, trace = TRUE)
    sum_vglm_summary <- summary(sum_vglm, HDEtest = FALSE)
    sum_vglm_summary_df <- as_tibble(sum_vglm_summary@coef3, rownames = "term")
    sum_vglm_summary_df$metric_name <- "sum"

    metadata_summary <- bind_rows(mt_percent_summary_df, silhouette_vglm_summary_df, detected_vglm_summary_df, sum_vglm_summary_df)
    metadata_summary$CI95lo <- metadata_summary$Estimate - (1.96 * metadata_summary$`Std. Error`)
    metadata_summary$CI95hi <- metadata_summary$Estimate + (1.96 * metadata_summary$`Std. Error`)
    metadata_summary$OR <- exp(metadata_summary$Estimate)
    metadata_summary$OR_CI95hi <- exp(metadata_summary$CI95hi)
    metadata_summary$OR_CI95lo <- exp(metadata_summary$CI95lo)
    metadata_summary
    metadata_summary
}

cluster_metrics_agg <- function(sce_object) {
    metadata_df <- colData(sce_object) %>% as_tibble()
    metrics_long <- select(metadata_df, Barcode, Sample, Cell_ID, batch, Celltype, Status, sum, detected, subsets_percent_mt_sum, silhouette_width)

    metrics_vglm <- group_by(metrics_long, Celltype) %>% group_modify(cluster_vglm_agg, sce_object)
    metrics_vglm
}

da_vglm <- function(abundances_df, sample_df) {
    n_samples <- group_by(abundances_df, Status) %>% summarise(Samples = n())
    replicates <- !any(n_samples$Samples < 2L)

    if (replicates) {
        abundances_vglm <- vglm(cbind(n_cells, total_cells - n_cells) ~ batch + Status, family = betabinomial(lrho = "rhobitlink"), data = abundances_df, trace = TRUE)
        abundances_summary <- summary(abundances_vglm, HDEtest = FALSE)
        abundances_summary_df <- as_tibble(abundances_summary@coef3, rownames = "term")
        abundances_summary_df$metric_name <- "abundances"
        abundances_summary_df
    } else {
        abundances_summary_df <- tibble(term = NA, Estimate = NA, "Std. Error" = NA, "z value" = NA, "Pr(>|z|)" = NA, metric_name = "abundances")
    }
    abundances_summary_df
}

sample_da_vglm <- function(sample_metadata, tissue_df) {
    metadata_reduce <- select(sample_metadata, batch, Sample, Status) %>% distinct()
    sample_abundances <- group_by(sample_metadata, Sample, Celltype) %>%
        summarise(n_cells = n()) %>%
        filter(Celltype != "?")
    total_cells <- group_by(sample_abundances, Sample) %>% summarise(total_cells = sum(n_cells))
    sample_abundances_df <- left_join(sample_abundances, total_cells) %>% left_join(metadata_reduce)

    abundance_summary_df <- group_by(sample_abundances_df, Celltype) %>% group_modify(da_vglm)
    abundance_summary_df$CI95lo <- abundance_summary_df$Estimate - (1.96 * abundance_summary_df$`Std. Error`)
    abundance_summary_df$CI95hi <- abundance_summary_df$Estimate + (1.96 * abundance_summary_df$`Std. Error`)
    abundance_summary_df$OR <- exp(abundance_summary_df$Estimate)
    abundance_summary_df$OR_CI95lo <- exp(abundance_summary_df$CI95lo)
    abundance_summary_df$OR_CI95hi <- exp(abundance_summary_df$CI95hi)
    abundance_summary_df$Tissue <- unique(tissue_df$Tissue)
    abundance_summary_df
}

sample_dirs_v3 <- c(
    # Fresh coronary samples
    Fresh_1667 = "/oak/stanford/groups/sjaiswal/single_cell_athero_human/all_data_redo/cellranger_output/1667_Fresh/outs/filtered_feature_bc_matrix",
    Fresh_1700 = "/oak/stanford/groups/sjaiswal/single_cell_athero_human/all_data_redo/cellranger_output/1700_Fresh/outs/filtered_feature_bc_matrix",

    # Frozen coronary samples
    Frozen_1667 = "/oak/stanford/groups/sjaiswal/single_cell_athero_human/all_data_redo/cellranger_output/1667_Frozen/outs/filtered_feature_bc_matrix",
    Frozen_1700 = "/oak/stanford/groups/sjaiswal/single_cell_athero_human/all_data_redo/cellranger_output/1700_Frozen/outs/filtered_feature_bc_matrix",

    # Frozen carotid samples
    Frozen_DTAN_4047 = "/oak/stanford/groups/sjaiswal/single_cell_athero_human/all_data_redo/cellranger_output/DTAN_4047_Frozen/outs/filtered_feature_bc_matrix"
)

sample_dirs_v2 <- c(
    # Fresh coronary samples
    Fresh_1495 = "/oak/stanford/groups/sjaiswal/single_cell_athero_human/all_data_redo/cellranger_output/1495_Fresh/outs/filtered_feature_bc_matrix",

    # Frozen coronary samples
    Frozen_1495 = "/oak/stanford/groups/sjaiswal/single_cell_athero_human/all_data_redo/cellranger_output/1495_Frozen/outs/filtered_feature_bc_matrix",

    # Fresh carotid samples
    Fresh_ROB_2026 = "/oak/stanford/groups/sjaiswal/single_cell_athero_human/all_data_redo/cellranger_output/ROB_2026_Fresh/outs/filtered_feature_bc_matrix",
    Fresh_DTAN_4047 = "/oak/stanford/groups/sjaiswal/single_cell_athero_human/all_data_redo/cellranger_output/DTAN_4047_Fresh/outs/filtered_feature_bc_matrix",

    # Frozen carotid samples
    Frozen_ROB_2026 = "/oak/stanford/groups/sjaiswal/single_cell_athero_human/all_data_redo/cellranger_output/ROB_2026_Frozen/outs/filtered_feature_bc_matrix"
)

# Read data and add gene annotation information
all_samples_sce_v3 <- read10xCounts(sample_dirs_v3)
all_samples_sce_v2 <- read10xCounts(sample_dirs_v2)

# Add phenotype information
all_samples_pheno <- tibble(
    Sample = c("Fresh_1495", "Frozen_1495", "Fresh_1667", "Frozen_1667", "Fresh_1700", "Frozen_1700", "Fresh_DTAN_4047", "Frozen_DTAN_4047", "Fresh_ROB_2026", "Frozen_ROB_2026"),
    Status = c("Fresh", "Frozen", "Fresh", "Frozen", "Fresh", "Frozen", "Fresh", "Frozen", "Fresh", "Frozen"),
    Tissue = c("Coronary", "Coronary", "Coronary", "Coronary", "Coronary", "Coronary", "Carotid", "Carotid", "Carotid", "Carotid")
)
colData(all_samples_sce_v3) %<>% as.data.frame %>% left_join(all_samples_pheno) %>% DataFrame()
colData(all_samples_sce_v2) %<>% as.data.frame %>% left_join(all_samples_pheno) %>% DataFrame()
colData(all_samples_sce_v3)$Cell_ID <- str_c(colData(all_samples_sce_v3)$Sample, colData(all_samples_sce_v3)$Barcode, sep = "_")
colData(all_samples_sce_v2)$Cell_ID <- str_c(colData(all_samples_sce_v2)$Sample, colData(all_samples_sce_v2)$Barcode, sep = "_")

annotation_hub <- AnnotationHub()
ens_104 <- annotation_hub[["AH95744"]]
gene_annot <- AnnotationDbi::select(ens_104, keys = rownames(all_samples_sce_v3), keytype = "GENEID", columns = c("GENEID", "GENEBIOTYPE", "DESCRIPTION", "SEQNAME", "GENESEQSTART", "GENESEQEND", "SEQSTRAND", "CANONICALTRANSCRIPT", "GENEIDVERSION"))
colnames(gene_annot)[1L] <- "ID"

all_samples_sce_filtered_v3 <- filter_osca(all_samples_sce_v3, gene_annot)
all_samples_sce_filtered_v2 <- filter_osca(all_samples_sce_v2, gene_annot)

all_samples_sce_doublets_v3 <- compute_doublets(all_samples_sce_filtered_v3)
all_samples_sce_doublets_v2 <- compute_doublets(all_samples_sce_filtered_v2)
write_rds(all_samples_sce_doublets_v3, "all_samples_sce_doublets_v3.rda")
write_rds(all_samples_sce_doublets_v2, "all_samples_sce_doublets_v2.rda")

all_samples_sce_corral_v3 <- estimate_corral(all_samples_sce_doublets_v3)
all_samples_sce_corral_v2 <- estimate_corral(all_samples_sce_doublets_v2)
write_rds(all_samples_sce_corral_v3, "all_samples_sce_corral_v3.rda")
write_rds(all_samples_sce_corral_v2, "all_samples_sce_corral_v2.rda")

all_samples_sce_umap_v3 <- estimate_umap(all_samples_sce_corral_v3)
all_samples_sce_umap_v2 <- estimate_umap(all_samples_sce_corral_v2)

all_samples_sce_leiden_v3 <- leiden_clustering(all_samples_sce_umap_v3)
all_samples_sce_leiden_v2 <- leiden_clustering(all_samples_sce_umap_v2)
write_rds(all_samples_sce_leiden_v3, "all_samples_sce_leiden_v3.rda")
write_rds(all_samples_sce_leiden_v2, "all_samples_sce_leiden_v2.rda")

all_samples_metadata_v3 <- colData(all_samples_sce_leiden_v3) %>%
    as_tibble() %>%
    bind_cols(as_tibble(reducedDim(all_samples_sce_leiden_v3, "umap"))) %>%
    as.data.frame()
all_samples_metadata_v2 <- colData(all_samples_sce_leiden_v2) %>%
    as_tibble() %>%
    bind_cols(as_tibble(reducedDim(all_samples_sce_leiden_v2, "umap"))) %>%
    as.data.frame()

qc_plots(all_samples_sce_leiden_v3, "all_samples_v3/")
qc_plots(all_samples_sce_leiden_v2, "all_samples_v2/")

umap_plot("CLDN5", "all_samples_v3/", all_samples_sce_leiden_v3)
umap_plot("CLDN5", "all_samples_v2/", all_samples_sce_leiden_v2)

umap_plot("CALD1", "all_samples_v3/", all_samples_sce_leiden_v3)
umap_plot("CALD1", "all_samples_v2/", all_samples_sce_leiden_v2)

umap_plot("CD68", "all_samples_v3/", all_samples_sce_leiden_v3)
umap_plot("CD68", "all_samples_v2/", all_samples_sce_leiden_v2)

umap_plot("CD3D", "all_samples_v3/", all_samples_sce_leiden_v3)
umap_plot("CD3D", "all_samples_v2/", all_samples_sce_leiden_v2)

all_samples_sce_immune_v3 <- all_samples_sce_leiden_v3[, !is_in(all_samples_sce_leiden_v3$label, c(17L, 18L))]
all_samples_sce_immune_v2 <- all_samples_sce_leiden_v2

# nolint start
# Plot marker genes
marker_genes <- c(
    "CD68", # macrophages
    "CD14", # macrophages/monocytes
    "FCGR3A", # monocytes

    "S100A8", # inflammatory macrophages
    "IL1B", # inflammatory macrophages
    "CD300E", # inflammatory macrophages
    "CSF3R", # inflammatory macrophages

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

    "HLA-DRB1", # MHC-hi macrophages
    "HLA-DQB1", # MHC-hi macrophages

    "FCER1A", # classical DCs
    "PLD4", # plasmacytoid DCs
    "XCR1", # XCR1+ DCs

    "FCGR3B", # neutrophils

    "MS4A2", # basophils/neutrophils
    "GATA2", # basophils/neutrophils
    "HPGD", # basophils

    "CD3D", # T cells
    #"CD3E", # T cells
    "CD4", # CD4 T cells
    "CD8A", # CD8 T cells
    "TRDC", # Gamma-delta T-cells, NK cells
    "TRGC2", # Gamma-delta T-cells, NK cells
    "FOXP3", # T-regs
    "GZMB", # Effector T-cells

    "NCAM1", # NK cells
    "XCL1", # active NK cells

    "CD19", # B cells
    "CD79A", # B cells

    "SDC1", # Plasma cells
    "JCHAIN", # Plasma cells

    "MKI67" # Cycling cells
)
# nolint end

catch <- map(marker_genes, umap_plot, "all_samples_v3/", all_samples_sce_leiden_v3)
catch <- map(marker_genes, umap_plot, "all_samples_v2/", all_samples_sce_leiden_v2)

#nolint start
all_samples_v2_labels <- c(
    "0" = "T-cells/NK cells",
    "1" = "T-cells/NK cells",
    "2" = "T-cells/NK cells",
    "3" = "B-cells",
    "4" = "T-cells/NK cells",
    "5" = "T-cells/NK cells",
    "6" = "Macrophages/cDCs",
    "7" = "T-cells/NK cells",
    "8" = "T-cells/NK cells",
    "9" = "T-cells/NK cells",
    "10" = "T-cells/NK cells",
    "11" = "T-cells/NK cells",
    "12" = "T-cells/NK cells",
    "13" = "Macrophages/cDCs",
    "14" = "Macrophages/cDCs",
    "15" = "Plasma cells",
    "16" = "T-cells/NK cells",
    "17" = "Plasmacytoid DCs"
)

all_samples_v3_labels <- c(
    "0" = "T-cells/NK cells",
    "1" = "B-cells",
    "2" = "Macrophages/cDCs",
    "3" = "T-cells/NK cells",
    "4" = "T-cells/NK cells",
    "5" = "Macrophages/cDCs",
    "6" = "B-cells",
    "7" = "Macrophages/cDCs",
    "8" = "T-cells/NK cells",
    "9" = "Macrophages/cDCs",
    "10" = "T-cells/NK cells",
    "11" = "T-cells/NK cells",
    "12" = "T-cells/NK cells",
    "13" = "Plasma cells",
    "14" = "Macrophages/cDCs",
    "15" = "Plasma cells",
    "16" = "Basophils",
    "17" = "Smooth muscle/endothelial cells",
    "18" = "Smooth muscle/endothelial cells",
    "19" = "Macrophages/cDCs",
    "20" = "Plasma cells"
)
# nolint end

all_samples_sce_celltypes_v3 <- add_celltype_labels(all_samples_sce_immune_v3, all_samples_v3_labels)
all_samples_sce_celltypes_v2 <- add_celltype_labels(all_samples_sce_immune_v2, all_samples_v2_labels)

all_samples_sce_tcells_v3 <- all_samples_sce_celltypes_v3[, str_detect(colData(all_samples_sce_celltypes_v3)$Celltype, fixed("T-cells"))]
all_samples_sce_tcells_v2 <- all_samples_sce_celltypes_v2[, str_detect(colData(all_samples_sce_celltypes_v2)$Celltype, fixed("T-cells"))]
write_rds(all_samples_sce_tcells_v3, "all_samples_sce_tcells_v3.rda")
write_rds(all_samples_sce_tcells_v2, "all_samples_sce_tcells_v2.rda")

all_samples_sce_macrophages_v3 <- all_samples_sce_celltypes_v3[, str_detect(colData(all_samples_sce_celltypes_v3)$Celltype, fixed("Macrophages"))]
all_samples_sce_macrophages_v2 <- all_samples_sce_celltypes_v2[, str_detect(colData(all_samples_sce_celltypes_v2)$Celltype, fixed("Macrophages"))]
write_rds(all_samples_sce_macrophages_v3, "all_samples_sce_macrophages_v3.rda")
write_rds(all_samples_sce_macrophages_v2, "all_samples_sce_macrophages_v2.rda")

other_cells_v3 <- !str_detect(colData(all_samples_sce_celltypes_v3)$Celltype, fixed("T-cells")) & !str_detect(colData(all_samples_sce_celltypes_v3)$Celltype, fixed("Macrophages"))
other_cells_v2 <- !str_detect(colData(all_samples_sce_celltypes_v2)$Celltype, fixed("T-cells")) & !str_detect(colData(all_samples_sce_celltypes_v2)$Celltype, fixed("Macrophages"))

all_samples_sce_merge_tcells <- read_rds("../t_cell_clustering/all_samples_sce_merge_celltypes.rda")
all_samples_sce_merge_macrophages <- read_rds("../macrophage_clustering/all_samples_sce_merge_celltypes.rda")
subclustered_cells <- union(all_samples_sce_merge_tcells$Cell_ID, all_samples_sce_merge_macrophages$Cell_ID)

all_samples_sce_final_v3 <- all_samples_sce_celltypes_v3[, other_cells_v3 | is_in(colData(all_samples_sce_celltypes_v3)$Cell_ID, subclustered_cells)]
all_samples_sce_final_v2 <- all_samples_sce_celltypes_v2[, other_cells_v2 | is_in(colData(all_samples_sce_celltypes_v2)$Cell_ID, subclustered_cells)]

all_samples_sce_leiden_Psaila <- read_rds("/oak/stanford/groups/sjaiswal/jk/052023_scrnaseq_reanalysis/scRNAseq/Psaila/all_samples_sce_leiden_Psaila.rda") %>% logNormCounts

find_markers <- findMarkers(all_samples_sce_leiden_Psaila, test.type = "t", pval.type = "any")
find_markers_df <- map(find_markers, as_tibble, rownames = "ID") %>% bind_rows(.id = "Cluster") 
find_markers_annot <- rowData(all_samples_sce_leiden_Psaila) %>% as_tibble() %>% select(Symbol, ID) %>% left_join(find_markers_df) 
filter(find_markers_annot, Cluster == 10) %>% arrange(p.value) %>% filter(summary.logFC > 0)

subset_samples_sce_Psaila <- all_samples_sce_leiden_Psaila[,is_in(colData(all_samples_sce_leiden_Psaila)$Sample, c("ID01", "ID02"))]
all_samples_sce_leiden_Psaila_merge <- merge_samples(subset_samples_sce_Psaila, "Sample")
all_samples_sce_merge <- merge_batches(all_samples_sce_final_v3, all_samples_sce_final_v2)
write_rds(all_samples_sce_merge, "all_samples_sce_merge.rda")

all_samples_sce_merge_umap <- estimate_umap(all_samples_sce_merge, "corrected")

all_samples_sce_merge_leiden <- leiden_clustering(all_samples_sce_merge_umap, "corrected", resolution = 1.0)
qc_plots(all_samples_sce_merge_leiden, "all_samples_merge/")

catch <- map(marker_genes, umap_plot, "all_samples_merge/", all_samples_sce_merge_leiden, "reconstructed")

# Load reference data and use SingleR to classify cells
monaco <- MonacoImmuneData()
lm22_se <- read_rds("../../new_code/lm22_se.rda")

all_samples_sce_merge_singler <- singler_annotation(all_samples_sce_merge_leiden, monaco, "monaco_celltypes", "all_samples_merge/monaco_labels") %>%
    singler_annotation(lm22_se, "lm22_celltypes", "all_samples_merge/lm22_labels")
write_rds(all_samples_sce_merge_singler, "all_samples_sce_merge_singler.rda")

# Find all markers
all_samples_markers <- findMarkers(all_samples_sce_merge_singler, test = "t", direction = "up", pval.type = "some", assay.type = "reconstructed", row.data = rowData(all_samples_sce_merge_singler))

#nolint start
# Assign rough cell types
all_samples_celltype_labels <- c(
    "0" = "T-cells/NK cells",
    "1" = "T-cells/NK cells",
    "2" = "B-cells",
    "3" = "T-cells/NK cells",
    "4" = "Macrophages/cDCs",
    "5" = "T-cells/NK cells",
    "6" = "Macrophages/cDCs",
    "7" = "T-cells/NK cells",
    "8" = "T-cells/NK cells",
    "9" = "Macrophages/cDCs",
    "10" = "T-cells/NK cells",
    "11" = "Plasma cells",
    "12" = "Plasmacytoid DCs",
    "13" = "B-cells",
    "14" = "Plasma cells",
    "15" = "Basophils"
)
#nolint end

all_samples_sce_merge_celltypes <- add_celltype_labels(all_samples_sce_merge_singler, all_samples_celltype_labels, "Rough_celltype")

# Run statistics
all_samples_silhouette <- reducedDim(all_samples_sce_merge_celltypes, "corrected") %>%
    approxSilhouette(all_samples_sce_merge_celltypes$Rough_celltype) %>%
    as_tibble()
colData(all_samples_sce_merge_celltypes)$silhouette_width <- all_samples_silhouette$width

# Add final fine grain clustering
tcells_celltypes <- colData(all_samples_sce_merge_tcells) %>% as_tibble() %>% select(Cell_ID, Celltype, silhouette_width)
macrophages_celltypes <- colData(all_samples_sce_merge_macrophages) %>% as_tibble() %>% select(Cell_ID, Celltype, silhouette_width)
subclustered_celltypes <- bind_rows(tcells_celltypes, macrophages_celltypes)

all_samples_sce_merge_tcells_orig <- all_samples_sce_merge_celltypes[, str_detect(all_samples_sce_merge_celltypes$Rough_celltype, fixed("T-cells/NK cells"))]
all_samples_sce_merge_tcells_metadata <- colData(all_samples_sce_merge_tcells_orig) %>% as_tibble() %>% select(-silhouette_width)
all_samples_sce_merge_tcells_labels <- left_join(all_samples_sce_merge_tcells_metadata, tcells_celltypes) %>% filter(!is.na(Celltype))
all_samples_sce_merge_tcells_filter <- all_samples_sce_merge_tcells_orig[, is_in(all_samples_sce_merge_tcells_orig$Cell_ID, all_samples_sce_merge_tcells_labels$Cell_ID)]
colData(all_samples_sce_merge_tcells_filter) <- DataFrame(all_samples_sce_merge_tcells_labels)

all_samples_sce_merge_macrophages_orig <- all_samples_sce_merge_celltypes[, str_detect(all_samples_sce_merge_celltypes$Rough_celltype, fixed("Macrophages/cDCs"))]
all_samples_sce_merge_macrophages_metadata <- colData(all_samples_sce_merge_macrophages_orig) %>% as_tibble() %>% select(-silhouette_width)
all_samples_sce_merge_macrophages_labels <- left_join(all_samples_sce_merge_macrophages_metadata, macrophages_celltypes) %>% filter(!is.na(Celltype))
all_samples_sce_merge_macrophages_filter <- all_samples_sce_merge_macrophages_orig[, is_in(all_samples_sce_merge_macrophages_orig$Cell_ID, all_samples_sce_merge_macrophages_labels$Cell_ID)]
colData(all_samples_sce_merge_macrophages_filter) <- DataFrame(all_samples_sce_merge_macrophages_labels)

all_samples_sce_merge_other <- all_samples_sce_merge_celltypes[, !str_detect(all_samples_sce_merge_celltypes$Rough_celltype, fixed("T-cells/NK cells")) & !str_detect(all_samples_sce_merge_celltypes$Rough_celltype, fixed("Macrophages/cDCs"))]
colData(all_samples_sce_merge_other)$Celltype <- colData(all_samples_sce_merge_other)$Rough_celltype

all_samples_sce_merge_celltypes2 <- cbind(all_samples_sce_merge_tcells_filter, all_samples_sce_merge_macrophages_filter, all_samples_sce_merge_other)
qc_umap_plots("all_samples_merge/", all_samples_sce_merge_celltypes2, "Celltype", var_type = "discrete", scale_name = "Cell type", file_name = "celltypes", label_var = "Celltype", plot_name = "", no_legend = TRUE)

write_rds(all_samples_sce_merge_celltypes2, "all_samples_sce_merge_celltypes2.rda")

# Run statistics
celltype_plots(all_samples_sce_merge_celltypes2, "all_samples_merge/")
all_samples_metrics <- cluster_metrics(all_samples_sce_merge_celltypes2)

all_samples_metrics_vglm <- cluster_metrics_vglm(all_samples_sce_merge_celltypes2)
write_rds(all_samples_metrics_vglm, "all_samples_metrics_vglm.rda")

all_samples_metrics_adj <- cluster_metrics_adj(all_samples_sce_merge_celltypes2)

# Merge counts
all_samples_sce_counts <- merge_counts(all_samples_sce_filtered_v3, all_samples_sce_filtered_v2)
all_samples_sce_counts_subset <- subset_counts(all_samples_sce_merge_celltypes2)

all_samples_metrics_agg <- cluster_metrics_agg(all_samples_sce_counts_subset)
write_rds(all_samples_metrics_agg, "all_samples_metrics_agg.rda")

v3_original_metadata <- colData(all_samples_sce_v3) %>% as_tibble()
v3_original_metadata$chemistry <- "v3"
v2_original_metadata <- colData(all_samples_sce_v2) %>% as_tibble()
v2_original_metadata$chemistry <- "v2"
all_samples_original_metadata <- bind_rows(v3_original_metadata, v2_original_metadata)
chemistry_df <- select(all_samples_original_metadata, Sample, chemistry) %>% distinct()

all_samples_original_ff <- group_by(all_samples_original_metadata, Status, Tissue, Sample) %>% tally()
all_samples_original_ff$Filtered <- "no"
all_samples_ff <- colData(all_samples_sce_merge_celltypes2) %>%
    as_tibble() %>%
    group_by(Status, Tissue, Sample) %>%
    tally()
all_samples_ff$Filtered <- "yes"
all_samples_ff_combined <- bind_rows(all_samples_original_ff, all_samples_ff) %>%
    pivot_wider(names_from = "Filtered", values_from = "n") %>%
    left_join(chemistry_df)
all_samples_ff_combined$unfiltered <- all_samples_ff_combined$no - all_samples_ff_combined$yes

ff_carotid <- filter(all_samples_ff_combined, Tissue == "Carotid") %>% arrange(desc(Status))
ff_carotid_vglm <- vglm(cbind(yes, unfiltered) ~ Status + chemistry, family = betabinomial(lrho = "rhobitlink"), data = ff_carotid)
ff_carotid_summary <- summary(ff_carotid_vglm, HDEtest = FALSE)
ff_carotid_summary_df <- as_tibble(ff_carotid_summary@coef3, rownames = "term")
ff_carotid_summary_df$OR <- exp(ff_carotid_summary_df$Estimate)

ff_coronary <- filter(all_samples_ff_combined, Tissue == "Coronary") %>% arrange(desc(Status))
ff_coronary_vglm <- vglm(cbind(yes, unfiltered) ~ Status + chemistry, family = betabinomial(lrho = "rhobitlink"), data = ff_coronary)
ff_coronary_summary <- summary(ff_coronary_vglm, HDEtest = FALSE)
ff_coronary_summary_df <- as_tibble(ff_coronary_summary@coef3, rownames = "term")
ff_coronary_summary_df$OR <- exp(ff_coronary_summary_df$Estimate)

ff_counts_all <- bind_rows(ff_carotid, ff_coronary)
ff_all_vglm <- vglm(cbind(yes, unfiltered) ~ Status + chemistry, family = betabinomial(lrho = "rhobitlink"), data = ff_counts_all)
ff_all_summary <- summary(ff_all_vglm, HDEtest = FALSE)
ff_all_summary_df <- as_tibble(ff_all_summary@coef3, rownames = "term")
ff_all_summary_df$OR <- exp(ff_all_summary_df$Estimate)
ff_counts_all_df <- select(ff_counts_all, Sample, Status, Tissue, no, yes, chemistry)
ff_counts_all_df$Sample %<>% str_remove_all("Fresh_") %>% str_remove_all("Frozen_")
ff_counts_all_df %<>% arrange(Sample, Status)
colnames(ff_counts_all_df)[4L:5L] <- c("Unfiltered", "Filtered")
write_tsv(ff_counts_all_df, "ff_counts_all_df.tsv")

t_cell_counts <- colData(all_samples_sce_merge_tcells) %>%
    as_tibble() %>%
    group_by(batch, Tissue, Status, Sample) %>%
    summarise(n_t_cells = n())

macrophage_cell_counts <- colData(all_samples_sce_merge_macrophages) %>%
    as_tibble() %>%
    group_by(batch, Tissue, Status, Sample) %>%
    summarise(n_macrophages = n())

all_cell_counts <- colData(all_samples_sce_merge_celltypes2) %>%
    as_tibble() %>%
    group_by(batch, Tissue, Status, Sample) %>%
    summarise(n_all_cells = n())

all_counts_merge <- left_join(all_cell_counts, t_cell_counts) %>%
    left_join(macrophage_cell_counts) %>%
    left_join(b_cell_counts)
all_counts_merge$not_t_cells <- all_counts_merge$n_all_cells - all_counts_merge$n_t_cells
all_counts_merge$not_macrophages <- all_counts_merge$n_all_cells - all_counts_merge$n_macrophages

all_counts_carotid <- filter(all_counts_merge, Tissue == "Carotid") %>% arrange(desc(Status))
t_cell_counts_carotid_vglm <- vglm(cbind(n_t_cells, not_t_cells) ~ Status + batch, family = betabinomial(lrho = "rhobitlink"), data = all_counts_carotid)
t_cell_counts_carotid_summary <- summary(t_cell_counts_carotid_vglm, HDEtest = FALSE)
t_cell_counts_carotid_summary_df <- as_tibble(t_cell_counts_carotid_summary@coef3, rownames = "term")
t_cell_counts_carotid_summary_df$OR <- exp(t_cell_counts_carotid_summary_df$Estimate)

macrophage_counts_carotid_vglm <- vglm(cbind(n_macrophages, not_macrophages) ~ Status + batch, family = betabinomial(lrho = "rhobitlink"), data = all_counts_carotid)
macrophage_counts_carotid_summary <- summary(macrophage_counts_carotid_vglm, HDEtest = FALSE)
macrophage_counts_carotid_summary_df <- as_tibble(macrophage_counts_carotid_summary@coef3, rownames = "term")
macrophage_counts_carotid_summary_df$OR <- exp(macrophage_counts_carotid_summary_df$Estimate)

all_counts_coronary <- filter(all_counts_merge, Tissue == "Coronary") %>% arrange(desc(Status))
t_cell_counts_coronary_vglm <- vglm(cbind(n_t_cells, not_t_cells) ~ Status + batch, family = betabinomial(lrho = "rhobitlink"), data = all_counts_coronary)
t_cell_counts_coronary_summary <- summary(t_cell_counts_coronary_vglm, HDEtest = FALSE)
t_cell_counts_coronary_summary_df <- as_tibble(t_cell_counts_coronary_summary@coef3, rownames = "term")
t_cell_counts_coronary_summary_df$OR <- exp(t_cell_counts_coronary_summary_df$Estimate)

macrophage_counts_coronary_vglm <- vglm(cbind(n_macrophages, not_macrophages) ~ Status + batch, family = betabinomial(lrho = "rhobitlink"), data = all_counts_coronary)
macrophage_counts_coronary_summary <- summary(macrophage_counts_coronary_vglm, HDEtest = FALSE)
macrophage_counts_coronary_summary_df <- as_tibble(macrophage_counts_coronary_summary@coef3, rownames = "term")
macrophage_counts_coronary_summary_df$OR <- exp(macrophage_counts_coronary_summary_df$Estimate)

all_samples_da_vglm <- colData(all_samples_sce_merge_celltypes2) %>%
    as_tibble() %>%
    group_by(Tissue) %>%
    group_map(sample_da_vglm) %>%
    bind_rows()
write_rds(all_samples_da_vglm, "all_samples_da_vglm.rda")

# GEO submission
colData(all_samples_sce_v2)$chemistry <- "v2"
colData(all_samples_sce_v3)$chemistry <- "v3"
all_samples_counts <- cbind(all_samples_sce_v2, all_samples_sce_v3)
rowData(all_samples_counts) %<>% as.data.frame %>% left_join(gene_annot) %>% DataFrame()
mito_genes <- which(rowData(all_samples_counts)$SEQNAME == "MT")

all_samples_counts_qc <- addPerCellQC(all_samples_counts, subset = list(mt = mito_genes))
colData(all_samples_counts_qc)$log_sum <- log2(all_samples_counts_qc$sum)
colData(all_samples_counts_qc)$log_detected <- log2(all_samples_counts_qc$detected)


merged_metadata <- colData(all_samples_sce_merge_celltypes2) %>%
    as_tibble() %>%
    select(Cell_ID, doublet_scores:silhouette_width)

colData(all_samples_counts_qc) %<>% as_tibble() %>% left_join(merged_metadata) %>% DataFrame()
colData(all_samples_counts_qc)$Celltype %<>% as.character %>% replace_na("Removed by QC") %>% factor()

colData(all_samples_counts_qc) %>% as_tibble() %>% write_tsv("all_samples_metadata.tsv")

save_sample_counts <- function(sample_name, sce_object) {
    sample_sce_object <- sce_object[, colData(sce_object)$Sample == sample_name]
    sample_counts_table <- counts(sample_sce_object)
    colnames(sample_counts_table) <- colData(sample_sce_object)$Cell_ID
    as.matrix(sample_counts_table) %>% as_tibble(rownames = "ensembl_id") %>% write_tsv(str_c(sample_name, ".tsv"))
}

unique(colData(all_samples_counts_qc)$Sample) %>%
    map(save_sample_counts, all_samples_counts_qc)

# Flow cytometry analysis
flow_cell_counts <- read_csv("./filtered_counts_with_flow_082222.csv") %>% filter(method == "flow")
flow_cell_counts$all_cells <- select(flow_cell_counts, `B-cells`:Other) %>% as.matrix() %>% rowSums()

flow_cell_counts_carotid <- filter(flow_cell_counts, str_detect(Sample, "ROB_2026|DTAN_4047"))
t_cell_counts_flow_carotid_vglm <- vglm(cbind(`T-cells`, all_cells - `T-cells`) ~ fresh_frozen, family = betabinomial(lrho = "rhobitlink"), data = flow_cell_counts_carotid)
t_cell_counts_flow_carotid_summary <- summary(t_cell_counts_flow_carotid_vglm, HDEtest = FALSE)
t_cell_counts_flow_carotid_summary_df <- as_tibble(t_cell_counts_flow_carotid_summary@coef3, rownames = "term")
t_cell_counts_flow_carotid_summary_df$OR <- exp(t_cell_counts_flow_carotid_summary_df$Estimate)

macrophage_counts_carotid_vglm <- vglm(cbind(Myeloid, all_cells - Myeloid) ~ fresh_frozen, family = betabinomial(lrho = "rhobitlink"), data = flow_cell_counts_carotid)
macrophage_counts_carotid_summary <- summary(macrophage_counts_carotid_vglm, HDEtest = FALSE)
macrophage_counts_carotid_summary_df <- as_tibble(macrophage_counts_carotid_summary@coef3, rownames = "term")
macrophage_counts_carotid_summary_df$OR <- exp(macrophage_counts_carotid_summary_df$Estimate)

flow_cell_counts_coronary <- filter(flow_cell_counts, !is_in(Sample, c("ROB2026", "DTAN4047")))
t_cell_counts_flow_coronary_vglm <- vglm(cbind(`T-cells`, all_cells - `T-cells`) ~ fresh_frozen, family = betabinomial(lrho = "rhobitlink"), data = flow_cell_counts_coronary)
t_cell_counts_flow_coronary_summary <- summary(t_cell_counts_flow_coronary_vglm, HDEtest = FALSE)
t_cell_counts_flow_coronary_summary_df <- as_tibble(t_cell_counts_flow_coronary_summary@coef3, rownames = "term")
t_cell_counts_flow_coronary_summary_df$OR <- exp(t_cell_counts_flow_coronary_summary_df$Estimate)

b_cell_counts_flow_coronary_vglm <- vglm(cbind(`B-cells`, all_cells - `B-cells`) ~ fresh_frozen, family = betabinomial(lrho = "rhobitlink"), data = flow_cell_counts_coronary)
b_cell_counts_flow_coronary_summary <- summary(b_cell_counts_flow_coronary_vglm, HDEtest = FALSE)
b_cell_counts_flow_coronary_summary_df <- as_tibble(b_cell_counts_flow_coronary_summary@coef3, rownames = "term")
b_cell_counts_flow_coronary_summary_df$OR <- exp(b_cell_counts_flow_coronary_summary_df$Estimate)

macrophage_counts_coronary_vglm <- vglm(cbind(Myeloid, all_cells - Myeloid) ~ fresh_frozen, family = betabinomial(lrho = "rhobitlink"), data = flow_cell_counts_coronary)
macrophage_counts_coronary_summary <- summary(macrophage_counts_coronary_vglm, HDEtest = FALSE)
macrophage_counts_coronary_summary_df <- as_tibble(macrophage_counts_coronary_summary@coef3, rownames = "term")
macrophage_counts_coronary_summary_df$OR <- exp(macrophage_counts_coronary_summary_df$Estimate)
