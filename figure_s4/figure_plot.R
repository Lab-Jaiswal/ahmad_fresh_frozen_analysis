library(SingleCellExperiment)
library(Seurat)
library(cowplot)
library(magrittr)
library(dplyr)
library(ggplot2)
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

all_samples_annotated <- read_rds("../figure2/all_samples_annotated.rda")

all_samples_1495 <- all_samples_annotated
colData(all_samples_1495)$show_cells <- all_samples_annotated$Sample == "Fresh_1495" | all_samples_annotated$Sample == "Frozen_1495"
all_samples_1667 <- all_samples_annotated
colData(all_samples_1667)$show_cells <- all_samples_annotated$Sample == "Fresh_1667" | all_samples_annotated$Sample == "Frozen_1667"
all_samples_1700 <- all_samples_annotated
colData(all_samples_1700)$show_cells <- all_samples_annotated$Sample == "Fresh_1700" | all_samples_annotated$Sample == "Frozen_1700"
all_samples_DTAN_4047 <- all_samples_annotated
colData(all_samples_DTAN_4047)$show_cells <- all_samples_annotated$Sample == "Fresh_DTAN_4047" | all_samples_annotated$Sample == "Frozen_DTAN_4047"
all_samples_ROB_2026 <- all_samples_annotated
colData(all_samples_ROB_2026)$show_cells <- all_samples_annotated$Sample == "Fresh_ROB_2026" | all_samples_annotated$Sample == "Frozen_ROB_2026"

figure_a <- umap_plot("Celltype", all_samples_1495, plot_name = "Coronary 1")
figure_b <- umap_plot("Celltype", all_samples_1667, plot_name = "Coronary 2")
figure_c <- umap_plot("Celltype", all_samples_1700, plot_name = "Coronary 3")
figure_d <- umap_plot("Celltype", all_samples_DTAN_4047, plot_name = "Carotid 1")
figure_e <- umap_plot("Celltype", all_samples_ROB_2026, plot_name = "Carotid 2")

figure <- ggdraw() +
    draw_plot(figure_a, x = 0.01, y = 0.5, width = 0.32, height = 0.5) +
    draw_plot(figure_b, x = 0.34, y = 0.5, width = 0.32, height = 0.5) +
    draw_plot(figure_c, x = 0.67, y = 0.5, width = 0.32, height = 0.5) +
    draw_plot(figure_d, x = 0.17, y = 0.0, width = 0.32, height = 0.5) +
    draw_plot(figure_e, x = 0.52, y = 0.0, width = 0.32, height = 0.5) +
    draw_label("A", size = 15.0, x = 0.0, y = 1.0, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("B", size = 15.0, x = 0.33, y = 1.0, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("C", size = 15.0, x = 0.66, y = 1.0, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("D", size = 15.0, x = 0.16, y = 0.5, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("E", size = 15.0, x = 0.5, y = 0.5, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold")

figure_title <- ggdraw() +
draw_label(label = "Figure S4", x = 0.5, hjust = 0.5, size = 20.0, fontfamily = "Noto Sans", fontface = "bold")

figure_grid <- plot_grid(figure_title, figure, ncol = 1L, rel_heights = c(0.05, 1.0))
figure_grid_notitle <- plot_grid(figure, ncol = 1L)
ggsave("figure_s4.pdf", figure_grid, width = 18.0, height = 12.0, units = "in", device = cairo_pdf)
ggsave("figure_s4_notitle.pdf", figure_grid_notitle, width = 18.0, height = 12.0, units = "in", device = cairo_pdf)
