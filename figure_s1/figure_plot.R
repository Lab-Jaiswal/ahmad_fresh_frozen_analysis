library(cowplot)
library(tidyverse)

figure <- ggdraw() +
    draw_image("./figure_s1.png", x = 0.0, y = 0.0, width = 1.0, height = 1.0)

figure_title <- ggdraw() +
    draw_label(label = "Figure S1", x = 0.5, hjust = 0.5, size = 18.0, fontface = "bold", fontfamily = "Noto Sans")

figure_grid <- plot_grid(figure_title, figure, ncol = 1L, rel_heights = c(0.1, 0.9))
figure_grid_notitle <- plot_grid(figure, ncol = 1L)

ggsave("figure_s1.pdf", figure_grid, width = 11.0, height = 6.0, units = "in", device = cairo_pdf)
ggsave("figure_s1_notitle.pdf", figure_grid_notitle, width = 11.0, height = 6.0, units = "in", device = cairo_pdf)
