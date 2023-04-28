library(openxlsx)
library(magrittr)
library(dplyr)
library(readr)
library(stringr)
library(tidyverse)

add_de_workbook <- function(index, table_list, wb) {
    de_table <- filter(table_list[[index]], !is.na(P.Value))
    de_table %<>% mutate(across(logFC:B, signif, digits = 3L)) %>% arrange(P.Value)
    pval_cols <- colnames(de_table) %>% str_detect("P.Value") %>% which()
    logfc_cols <- colnames(de_table) %>% str_detect("logFC") %>% which()

    workbook_name <- names(table_list)[index]

    addWorksheet(wb, sheetName = workbook_name)
    writeDataTable(wb, sheet = index, x = de_table)

    sig_style <- createStyle(fontColour = "red")
    conditionalFormatting(wb, index, cols = pval_cols,
                          rows = seq_len(nrow(de_table)), rule = "<0.05", style = sig_style)
    conditionalFormatting(wb, index, cols = logfc_cols,
                          rows = seq_len(nrow(de_table)), style = c("#63BE7B", "white", "red"),
                          type = "colourScale")

    setColWidths(wb, index, cols = 1L:4L, widths = "auto")
    setColWidths(wb, index, cols = 5L, widths = 45L)
    setColWidths(wb, index, cols = 6L:ncol(de_table), widths = "auto")

    pageSetup(wb, index, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, index, firstRow = TRUE)
}

# Table S1
celltype_da <- read_rds("../all_samples_da_vglm.rda") %>% filter(term == "StatusFrozen")
colnames(celltype_da)[1L] <- "Cell type"
colnames(celltype_da)[4L] <- "SE"
colnames(celltype_da)[5L] <- "Z"
colnames(celltype_da)[6L] <- "P"

celltype_da_select <- select(celltype_da, Tissue, `Cell type`, Estimate, SE, CI95lo:OR, OR_CI95lo, OR_CI95hi, Z, P) %>%
    mutate(across(where(is.numeric), signif, 3L))
celltype_da_coronary <- filter(celltype_da_select, Tissue == "Coronary") %>% select(-Tissue)
celltype_da_coronary$P_adjust <- p.adjust(celltype_da_coronary$P, method = "bonferroni")
celltype_da_coronary_filter <- filter(celltype_da_coronary, abs(Z) < 100.0)
celltype_da_carotid <- filter(celltype_da_select, Tissue == "Carotid") %>% select(-Tissue)
celltype_da_carotid$P_adjust <- p.adjust(celltype_da_carotid$P, method = "bonferroni")
celltype_da_carotid_filter <- filter(celltype_da_carotid, abs(Z) < 100.0)

wb_s1 <- createWorkbook()
addWorksheet(wb = wb_s1, sheetName = "Coronary pseudobulk")
writeDataTable(wb = wb_s1, sheet = 1L, x = celltype_da_coronary_filter, withFilter = FALSE)
setColWidths(wb = wb_s1, sheet = 1L, cols = seq_len(ncol(celltype_da_coronary)), widths = "auto")

addWorksheet(wb = wb_s1, sheetName = "Carotid pseudobulk")
writeDataTable(wb = wb_s1, sheet = 2L, x = celltype_da_carotid_filter, withFilter = FALSE)
setColWidths(wb = wb_s1, sheet = 2L, cols = seq_len(ncol(celltype_da_carotid)), widths = "auto")

saveWorkbook(wb = wb_s1, "table_s1.xlsx", overwrite = TRUE)

# Table S2
cluster_metrics <- read_rds("../all_samples_metrics_vglm.rda") %>% filter(term == "StatusFrozen")
colnames(cluster_metrics)[1L] <- "Cell type"
colnames(cluster_metrics)[4L] <- "SE"
colnames(cluster_metrics)[5L] <- "Z"
colnames(cluster_metrics)[6L] <- "P"

cluster_metrics_select <- select(cluster_metrics, metric_name, `Cell type`, Estimate, SE, CI95lo:OR, OR_CI95lo, OR_CI95hi, Z, P) %>%
    mutate(across(where(is.numeric), signif, 3L))
cluster_metrics_select$metric_name %<>% str_replace_all("subsets_percent_mt_percent", fixed("percent_mt"))
percent_mt <- filter(cluster_metrics_select, metric_name == "percent_mt") %>% select(-metric_name)
percent_mt$P_adjust <- p.adjust(percent_mt$P, method = "bonferroni")
detected <- filter(cluster_metrics_select, metric_name == "detected") %>% select(-metric_name)
detected$P_adjust <- p.adjust(detected$P, method = "bonferroni")
sum_umi <- filter(cluster_metrics_select, metric_name == "sum") %>% select(-metric_name)
sum_umi$P_adjust <- p.adjust(sum_umi$P, method = "bonferroni")

cluster_metrics_agg <- read_rds("../all_samples_metrics_agg.rda") %>% filter(term == "StatusFrozen")
colnames(cluster_metrics_agg)[1L] <- "Cell type"
colnames(cluster_metrics_agg)[4L] <- "SE"
colnames(cluster_metrics_agg)[5L] <- "Z"
colnames(cluster_metrics_agg)[6L] <- "P"

cluster_metrics_agg_select <- select(cluster_metrics_agg, metric_name, `Cell type`, Estimate, SE, CI95lo:OR, OR_CI95lo, OR_CI95hi, Z, P) %>%
    mutate(across(where(is.numeric), signif, 3L))
cluster_metrics_agg_select$metric_name %<>% str_replace_all("subsets_percent_mt_percent", fixed("percent_mt"))
percent_mt_agg <- filter(cluster_metrics_agg_select, metric_name == "percent_mt") %>% select(-metric_name)
percent_mt_agg$P_adjust <- p.adjust(percent_mt_agg$P, method = "bonferroni")
detected_agg <- filter(cluster_metrics_agg_select, metric_name == "detected") %>% select(-metric_name)
detected_agg$P_adjust <- p.adjust(detected_agg$P, method = "bonferroni")
sum_umi_agg <- filter(cluster_metrics_agg_select, metric_name == "sum") %>% select(-metric_name)
sum_umi_agg$P_adjust <- p.adjust(sum_umi_agg$P, method = "bonferroni")

wb_s2 <- createWorkbook()
addWorksheet(wb = wb_s2, sheetName = "percent mt")
writeDataTable(wb = wb_s2, sheet = 1L, x = percent_mt, withFilter = FALSE)
setColWidths(wb = wb_s2, sheet = 1L, cols = seq_len(ncol(percent_mt)), widths = "auto")

addWorksheet(wb = wb_s2, sheetName = "percent mt pseudobulk")
writeDataTable(wb = wb_s2, sheet = 2L, x = percent_mt_agg, withFilter = FALSE)
setColWidths(wb = wb_s2, sheet = 2L, cols = seq_len(ncol(percent_mt_agg)), widths = "auto")

addWorksheet(wb = wb_s2, sheetName = "unique genes")
writeDataTable(wb = wb_s2, sheet = 3L, x = detected, withFilter = FALSE)
setColWidths(wb = wb_s2, sheet = 3L, cols = seq_len(ncol(detected)), widths = "auto")

addWorksheet(wb = wb_s2, sheetName = "unique genes pseudobulk")
writeDataTable(wb = wb_s2, sheet = 4L, x = detected_agg, withFilter = FALSE)
setColWidths(wb = wb_s2, sheet = 4L, cols = seq_len(ncol(detected_agg)), widths = "auto")

addWorksheet(wb = wb_s2, sheetName = "UMIs")
writeDataTable(wb = wb_s2, sheet = 5L, x = sum_umi, withFilter = FALSE)
setColWidths(wb = wb_s2, sheet = 5L, cols = seq_len(ncol(sum_umi)), widths = "auto")

addWorksheet(wb = wb_s2, sheetName = "UMIs pseudobulk")
writeDataTable(wb = wb_s2, sheet = 6L, x = sum_umi_agg, withFilter = FALSE)
setColWidths(wb = wb_s2, sheet = 6L, cols = seq_len(ncol(sum_umi_agg)), widths = "auto")

saveWorkbook(wb = wb_s2, "table_s2.xlsx", overwrite = TRUE)

# Table S3
silhouette_width <- filter(cluster_metrics_select, metric_name == "silhouette_width") %>% select(-metric_name)
silhouette_width$P_adjust <- p.adjust(silhouette_width$P, method = "bonferroni")

silhouette_width_agg <- filter(cluster_metrics_agg_select, metric_name == "silhouette_width") %>% select(-metric_name)
silhouette_width_agg$P_adjust <- p.adjust(silhouette_width_agg$P, method = "bonferroni")

wb_s3 <- createWorkbook()
addWorksheet(wb = wb_s3, sheetName = "silhouette width")
writeDataTable(wb = wb_s3, sheet = 1L, x = silhouette_width, withFilter = FALSE)
setColWidths(wb = wb_s3, sheet = 1L, cols = seq_len(ncol(silhouette_width)), widths = "auto")

addWorksheet(wb = wb_s3, sheetName = "silhouette width pseudobulk")
writeDataTable(wb = wb_s3, sheet = 2L, x = silhouette_width_agg, withFilter = FALSE)
setColWidths(wb = wb_s3, sheet = 2L, cols = seq_len(ncol(silhouette_width_agg)), widths = "auto")
saveWorkbook(wb = wb_s3, "table_s3.xlsx", overwrite = TRUE)

# Table S4
de_table_list_mac <- read_rds("../macrophage_clustering/all_samples_limma_tibble.rda")
names(de_table_list_mac) %<>% str_replace_all("macrophages", fixed("mac"))

de_table_list_ordered_mac <- de_table_list[c("Foam cells",
    "Inflammatory mac",
    "MHC-hi mac",
    "cDCs",
    "LYVE1+ TR mac",
    "CD16+ monocytes")
]

de_table_list_t_cell <- read_rds("../t_cell_clustering/all_samples_limma_tibble.rda")
names(de_table_list_mac) %<>% str_replace_all("macrophages", fixed("mac"))

de_table_list_ordered_t_cell <- de_table_list_t_cell[c("Active NK cells",
    "Resting NK cells",
    "CD4+ T-cells",
    "CD8+ T-cells")
]

de_table_list_ordered <- c(de_table_list_ordered_mac, de_table_list_ordered_t_cell)

wb_s4 <- createWorkbook()

seq_along(de_table_list_ordered) %>% map(add_de_workbook, de_table_list_ordered, wb_s4)

saveWorkbook(wb_s4, "table_s4.xlsx", overwrite = TRUE)
