library(SingleCellExperiment)
library(SummarizedExperiment)
library(magick)
library(cowplot)
library(magrittr)
library(dplyr)
library(ggplot2)
library(readr)
library(rlang)
library(tidyverse)

umap_plot <- function(feature_color, sce_object, data_type = "counts", var_type = "continuous") {
    umap_data <- reducedDim(sce_object, type = "umap", withDimnames = TRUE)
    plot_data <- colData(sce_object) %>%
        as_tibble() %>%
        bind_cols(as_tibble(umap_data))

    feature_plot <- "log_expr"
    scale_name_final <- expression(paste(log[2.0], "(counts)"))
    gene_annot <- rowData(sce_object)
    gene_row <- match(feature_color, gene_annot$Symbol)

    assay_data <- assays(sce_object)
    reconstructed_data <- as.matrix(assay_data$reconstructed)
    plot_data$log_expr <- reconstructed_data[gene_row, ]

    plot_data %<>% arrange(!!sym(feature_plot)) %>% as.data.frame()

    umap_plot <- ggplot(plot_data, aes(UMAP_1, UMAP_2, color = !!sym(feature_plot))) +
        geom_point(size = 0.1) +
        ggtitle(feature_color) +
        theme_classic(base_family = "Noto Sans") +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        theme(panel.border = element_rect(fill = NA),
              plot.title = element_text(hjust = 0.5))

    umap_plot <- umap_plot + scale_color_gradient(low = "lightgrey", high = "blue", name = scale_name_final)

    umap_plot
}

all_samples_annotated <- read_rds("../figure2/all_samples_annotated.rda")

figure_a <- umap_plot("APOE", all_samples_annotated, data_type = "reconstructed")

plot_data <- colData(all_samples_annotated) %>% as_tibble()
gene_annot <- rowData(all_samples_annotated)
gene_row <- match("APOE", gene_annot$Symbol)
assay_data <- assays(all_samples_annotated)
reconstructed_data <- as.matrix(assay_data$reconstructed)
plot_data$log_expr <- reconstructed_data[gene_row, ]
plot_data_filter <- filter(plot_data, is_in(Celltype, c("Inflammatory Mφ", "Foam cells", "CD16+ monocytes", "cDCs", "LYVE1+ TR Mφ", "MHC-hi Mφ")))
plot_data_filter$Celltype %<>% str_replace_all("macrophages", fixed("Mφ")) %>%
    factor(levels = c(
        "CD16+ monocytes",
        "MHC-hi Mφ",
        "Inflammatory Mφ",
        "LYVE1+ TR Mφ",
        "Foam cells",
        "cDCs"
    ))

figure_b <- ggplot(plot_data_filter, aes(Celltype, log_expr)) +
  geom_boxplot() +
  ylab(expression(paste(log[2.0], "(counts)"))) +
  theme_classic(base_family = "Noto Sans") +
  theme(legend.text = element_text(size = 12.0),
        legend.title = element_text(size = 13.0),
        axis.text = element_text(size = 12.0),
        axis.title = element_text(size = 12.0),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 45.0, hjust = 1.0))

figure <- ggdraw() +
    draw_plot(figure_a, x = 0.01, y = 0.0, width = 0.64, height = 1.0)  +
    draw_plot(figure_b, x = 0.66, y = 0.2, width = 0.34, height = 0.60) +
    draw_label("A", size = 15.0, x = 0.0, y = 1.0, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("B", size = 15.0, x = 0.65, y = 0.8, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold")

figure_title <- ggdraw() +
    draw_label(label = "Figure S5", x = 0.5, hjust = 0.5, size = 18.0, fontfamily = "Noto Sans", fontface = "bold")

figure_grid <- plot_grid(figure_title, figure, ncol = 1L, rel_heights = c(0.1, 0.9))
figure_grid_notitle <- plot_grid(figure, ncol = 1L)
ggsave("figure_s5.pdf", figure_grid, width = 10.0, height = 6.0, units = "in", device = cairo_pdf)
ggsave("figure_s5_notitle.pdf", figure_grid_notitle, width = 10.0, height = 6.0, units = "in", device = cairo_pdf)
