library(SingleCellExperiment)
library(magick)
library(cowplot)
library(magrittr)
library(dplyr)
library(ggplot2)
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
            axis.text = element_text(size = 12.0, hjust = 0.0),
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

all_samples_annotated <- read_rds("../figure2/all_samples_annotated.rda")

all_samples_metadata <- colData(all_samples_annotated) %>%
    as_tibble() %>%
    filter(!is_in(Celltype, c("XCR1+ DCs", "Plasma cells", "Plasmacytoid DCs", "Basophils")))
all_samples_metadata$Celltype %<>% str_replace_all("macrophages", fixed("Mφ")) %>%
    factor(levels = c(
        "CD16+ monocytes",
        "MHC-hi Mφ",
        "Inflammatory Mφ",
        "LYVE1+ TR Mφ",
        "Foam cells",
        "cDCs",
        "CD4+ T-cells",
        "CD8+ T-cells",
        "Active NK cells",
        "Resting NK cells",
        "B-cells"
    ))

figure_a <- ggplot(all_samples_metadata, aes(Status, fill = `Celltype`)) +
    geom_bar(position = "fill", color = "black") +
    facet_wrap(~ Tissue, ncol = 2L) +
    theme_classic(base_family = "Noto Sans") +
    theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(size = 12.0),
        legend.title = element_text(size = 13.0),
        axis.text = element_text(size = 12.0),
        axis.title = element_text(size = 12.0),
        strip.text = element_text(face = "bold", size = 14.0),
        panel.background = element_rect(color = "black", fill = NA)) +
    ylab("Relative Proportion")

figure_c <- qc_boxplot("subsets_percent_mt_percent", all_samples_metadata, x_var = "Celltype", fill_var = "Status",
                       y_label = "% reads in MT genes", angled_text = TRUE)
figure_d <- qc_boxplot("log_detected", all_samples_metadata, x_var = "Celltype", fill_var = "Status",
                       y_label = expression(paste(log[2.0], "(detected genes)")), angled_text = TRUE)

ggsave("figure3D.pdf", figure_d, width = 12.0, height = 4.0, units = "in", device = cairo_pdf)

figure <- ggdraw() +
    draw_plot(figure_a, x = 0.01, y = 0.66, width = 0.55, height = 0.33)  +
    draw_image("../../new_code/figure2_panels/figure_2B.png", x = 0.57, y = 0.66, width = 0.4, height = 0.33) +
    draw_plot(figure_c, x = 0.01, y = 0.33, width = 0.99, height = 0.33) +
    draw_plot(figure_d, x = 0.01, y = 0.0, width = 0.99, height = 0.33) +
    draw_label("A", size = 15.0, x = 0.0, y = 1.0, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("B", size = 15.0, x = 0.56, y = 1.0, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("C", size = 15.0, x = 0.0, y = 0.66, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold") +
    draw_label("D", size = 15.0, x = 0.0, y = 0.33, hjust = 0.0, vjust = 1.0, fontfamily = "Noto Sans", fontface = "bold")

figure_title <- ggdraw() +
    draw_label(label = "Figure 3", x = 0.5, hjust = 0.5, size = 20.0, fontfamily = "Noto Sans", fontface = "bold")
figure_grid <- plot_grid(figure_title, figure, ncol = 1.0, rel_heights = c(0.05, 1.0))

ggsave("figure3.pdf", figure_grid, width = 12.0, height = 12.0, units = "in", device = cairo_pdf)
