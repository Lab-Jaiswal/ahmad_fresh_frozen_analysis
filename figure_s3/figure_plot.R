library(SingleCellExperiment)
library(Seurat)
library(cowplot)
library(magrittr)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyverse)

umap_plot <- function(feature_color, sce_object, plot_name = NULL) {
    umap_data <- reducedDim(sce_object, type = "umap", withDimnames = TRUE)
    plot_data <- colData(sce_object) %>%
        as_tibble() %>%
        bind_cols(as_tibble(umap_data))

    plot_data %<>% arrange(Cell_ID) %>% as.data.frame()
    plot_data$color_factor <- plot_data[[feature_color]]
    plot_data$color_factor[!plot_data$show_cells] <- NA
    plot_data[[feature_color]] %<>% str_replace_all("macrophages", fixed("Mφ")) %>% factor()

    umap_plot <- ggplot(data = plot_data, aes(UMAP_1, UMAP_2)) +
        geom_point(color = "lightgray", size = 0.1) +
        geom_point(aes(color = color_factor), size = 0.1) +
        ggtitle(plot_name) +
        theme_classic(base_family = "Noto Sans") +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        scale_color_discrete(na.value = "#00000000") +
        theme(panel.border = element_rect(fill = NA),
              plot.title = element_text(hjust = 0.5),
              legend.position = "none")

    umap_plot %<>% LabelClusters(id = feature_color, color = "black", seed = 12345L, family = "Noto Sans")

    umap_plot
}

all_samples_sce_merge_celltypes <- read_rds("../frozen_clustering/all_samples_sce_merge_celltypes2.rda")
colData(all_samples_sce_merge_celltypes)$Celltype %<>% str_replace_all("macrophages", factor("Mφ")) %>% factor()
write_rds(all_samples_sce_merge_celltypes, "all_samples_annotated.rda")

all_samples_frozen_coronary <- all_samples_sce_merge_celltypes
colData(all_samples_frozen_coronary)$show_cells <- all_samples_sce_merge_celltypes$Tissue == "Coronary" & all_samples_sce_merge_celltypes$Status == "Frozen"
all_samples_frozen_carotid <- all_samples_sce_merge_celltypes
colData(all_samples_frozen_carotid)$show_cells <- all_samples_sce_merge_celltypes$Tissue == "Carotid" & all_samples_sce_merge_celltypes$Status == "Frozen"

figure_a <- umap_plot("Celltype", all_samples_frozen_coronary, plot_name = "Frozen Coronary")
figure_b <- umap_plot("Celltype", all_samples_frozen_carotid, plot_name = "Frozen Carotid")

figure <- ggdraw() +
    draw_plot(figure_a, x = 0.01, y = 0.0, width = 0.49, height = 1.0) +
    draw_plot(figure_b, x = 0.51, y = 0.0, width = 0.49, height = 1.0) +
    draw_label("A", size = 15.0, x = 0.0, y = 1.0, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("B", size = 15.0, x = 0.5, y = 1.0, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold")

figure_title <- ggdraw() +
draw_label(label = "Figure S3", x = 0.5, hjust = 0.5, size = 16.0, fontfamily = "Noto Sans", fontface = "bold")

figure_grid <- plot_grid(figure_title, figure, ncol = 1L, rel_heights = c(0.1, 0.8))
figure_grid_notitle <- plot_grid(figure, ncol = 1L)
ggsave("figure_s3.pdf", figure_grid, width = 12.0, height = 6.0, units = "in", device = cairo_pdf)
ggsave("figure_s3_notitle.pdf", figure_grid_notitle, width = 12.0, height = 6.0, units = "in", device = cairo_pdf)
