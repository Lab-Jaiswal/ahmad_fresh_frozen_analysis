library(SingleCellExperiment)
library(magick)
library(cowplot)
library(ggrepel)
library(magrittr)
library(dplyr)
library(ggplot2)
library(scales)
library(readr)
library(rlang)
library(tidyverse)

qc_boxplot <- function(qc_var, plot_data, x_var = "label", fill_var = NULL,
                       scale_name = NULL, y_label = NULL, angled_text = FALSE) {
   if (!is.null(fill_var)) {
        box_plot_aes <- aes(!!sym(x_var), !!sym(qc_var), fill = !!sym(fill_var))
    } else {
        box_plot_aes <- aes(!!sym(x_var), !!sym(qc_var))
    }

    box_plot <- ggplot(plot_data, box_plot_aes) +
      geom_boxplot() +
      theme_classic(base_family = "Noto Sans") +
      theme(legend.text = element_text(size = 12.0),
            legend.title = element_text(size = 13.0),
            axis.text = element_text(size = 12.0),
            axis.title = element_text(size = 12.0),
            axis.title.x = element_blank(),
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

    box_plot
}

volcano_plot <- function(top_table, plot_title, logfc_cutoff = NA, cutoff = 0.05, cutoff_column = "p_val_adj", log_column = "logFC", xlabel = "Log Fold Change") {
    top_table_filter <- filter(top_table, !is.na(P.Value)) %>%  arrange(P.Value)
    if (!is.na(logfc_cutoff)) {
        sig_logfc <- abs(top_table[[log_column]]) > logfc_cutoff
        top_table$Significant <- top_table$Significant & sig_logfc
    }
    top_table_filter$Log_Pvalue <- -log10(top_table_filter$P.Value)
    top_table_filter$Symbol_new <- NA
    top_table_filter$Symbol_new[1L:10L] <- top_table_filter$Symbol[1L:10L]

    p <- ggplot(top_table_filter, aes_string(x = log_column, y = "Log_Pvalue")) +
        geom_point(color = "gray") +
        theme_bw(base_family = "Noto Sans") +
        geom_hline(yintercept = -log10(0.05 / nrow(top_table_filter)), color = "red", linetype = "dashed") +
        geom_text_repel(label = top_table_filter$Symbol_new, color = muted("red"), segment.color = "black", box.padding = 0.5, min.segment.length = 0.0, size = 3.0, nudge_y = 0.25, seed = 12345L) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size = 11.0, face = "bold"),
              legend.title = element_blank(),
              panel.border = element_rect(size = 1.0, color = "black"),
              plot.background = element_blank(),
              plot.title = element_text(hjust = 0.5, face = "bold")) +
        xlab(xlabel) + ylab(expression(paste("-", log[10.0], " p-value"))) + ggtitle(plot_title)

    p
}

all_samples_limma_tibble_macro <- read_rds("../macrophage_clustering/all_samples_limma_tibble.rda")
all_samples_limma_tibble_tcell <- read_rds("../t_cell_clustering/all_samples_limma_tibble.rda")

all_samples_annotated <- read_rds("../figure2/all_samples_annotated.rda")

all_samples_metadata <- colData(all_samples_annotated) %>%
    as_tibble() %>%
    filter(is_in(Celltype, c("XCR1+ DCs", "Plasma cells", "Plasmacytoid DCs", "Basophils")))
all_samples_metadata$Celltype %<>% str_replace_all("macrophages", fixed("Mφ")) %>%
    factor(levels = c(
        "Plasma cells",
        "Plasmacytoid DCs",
        "XCR1+ DCs",
        "Basophils"
    ))

figure_a <- qc_boxplot("silhouette_width", all_samples_metadata, x_var = "Celltype", fill_var = "Status",
                       y_label = "Silhouette width", angled_text = TRUE)
figure_b <- volcano_plot(all_samples_limma_tibble_macro[["LYVE1+ TR macrophages"]], "LYVE1+ TR macrophages")
figure_c <- volcano_plot(all_samples_limma_tibble_macro[["CD16+ monocytes"]], "CD16+ monocytes")
figure_d <- volcano_plot(all_samples_limma_tibble_tcell[["Active NK cells"]], "Active NK cells")
figure_e <- volcano_plot(all_samples_limma_tibble_tcell[["Resting NK cells"]], "Resting NK cells")
figure_f <- volcano_plot(all_samples_limma_tibble_tcell[["CD4+ T-cells"]], "CD4+ T-cells")
figure_g <- volcano_plot(all_samples_limma_tibble_tcell[["CD8+ T-cells"]], "CD8+ T-cells")

figure <- ggdraw() +
    draw_plot(figure_a, x = 0.19, y = 0.75, width = 0.64, height = 0.25)  +
    draw_plot(figure_b, x = 0.01, y = 0.5, width = 0.49, height = 0.25) +
    draw_plot(figure_c, x = 0.51, y = 0.5, width = 0.49, height = 0.25) +
    draw_plot(figure_d, x = 0.01, y = 0.25, width = 0.49, height = 0.25) +
    draw_plot(figure_e, x = 0.51, y = 0.25, width = 0.49, height = 0.25) +
    draw_plot(figure_f, x = 0.01, y = 0.0, width = 0.49, height = 0.25) +
    draw_plot(figure_g, x = 0.51, y = 0.0, width = 0.49, height = 0.25) +
    draw_label("A", size = 15.0, x = 0.18, y = 1.0, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("B", size = 15.0, x = 0.0, y = 0.75, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("C", size = 15.0, x = 0.5, y = 0.75, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("D", size = 15.0, x = 0.0, y = 0.5, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("E", size = 15.0, x = 0.5, y = 0.5, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("F", size = 15.0, x = 0.0, y = 0.25, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("G", size = 15.0, x = 0.5, y = 0.25, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold")

figure_title <- ggdraw() +
    draw_label(label = "Figure S7", x = 0.5, hjust = 0.5, size = 20.0, fontfamily = "Noto Sans", fontface = "bold")

figure_grid <- plot_grid(figure_title, figure, ncol = 1L, rel_heights = c(0.05, 1.0))
figure_grid_notitle <- plot_grid(figure, ncol = 1L)
ggsave("figure_s7.pdf", figure_grid, width = 12.0, height = 24.0, units = "in", device = cairo_pdf)
ggsave("figure_s7_notitle.pdf", figure_grid_notitle, width = 12.0, height = 24.0, units = "in", device = cairo_pdf)
