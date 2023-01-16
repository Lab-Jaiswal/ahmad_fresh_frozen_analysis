library(Cairo)
library(cowplot)
library(tidyverse)

figure <- ggdraw() +
    draw_image("../../new_code/figure1_panels/figure_1A_new.png", x = 0.0, y = 0.0, width = 0.5, height = 1.0) +
    draw_image("../../new_code/figure1_panels/figure_1D_new.png", x = 0.51, y = 0.0, width = 0.5, height = 1.0) +
    draw_label("A", size = 15.0, x = 0.0, y = 1.0, hjust = 0.0, vjust = 1.0, fontface = "bold") +
    draw_label("B", size = 15.0, x = 0.5, y = 1.0, hjust = 0.0, vjust = 1.0, fontface = "bold")

figure_title <- ggdraw() +
draw_label(label = "Figure 1", x = 0.5, hjust = 0.5, size = 18.0, fontface = "bold")

figure_grid <- plot_grid(figure_title, figure, ncol = 1L, rel_heights = c(0.1, 0.9))

CairoPDF("figure1", width = 12.0, height = 6.0)
    print(figure_grid)
dev.off()
