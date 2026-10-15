#!/usr/bin/env Rscript
# ============================================================
# 02_Normalize_Batch 图：逐样本文库大小 + 细胞数总览
#
# 输入：01_QC/data/gse150290_qc_meta.csv + gse150290_doublet_scores.csv
# 产出：02_Normalize_Batch/figures/library_and_cells_overview.png
#       + 02_Normalize_Batch/data/library_summary.json
#
# 关键设计要点：
#   - 细胞口径 = 02 归一化的同一批细胞：QC 保留细胞剔除 Scrublet
#     predicted_doublet（与 normalize_batch.R 第 1b 节同款过滤，
#     不重算任何值，只读 01_QC 的原始 meta）
#   - 中位文库仅作信息参考线（scale.factor 固定 10000）
#   - 顶面板 y 轴锚定真 0 基线；无标题（y 轴标签已说明面板内容）
#   - 底面板无 x 轴标题：逐样本身份无信息量，52 个标签纯噪音
#
# Usage (Windows native R 4.5.1)：
#   "E:\R-4.5.1\bin\Rscript.exe" 02_Normalize_Batch\scripts\viz_library_scaling.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

# ── 路径 ─────────────────────────────────────────────────────
BASE     <- "02_Normalize_Batch"
META_CSV <- "01_QC/data/gse150290_qc_meta.csv"
FIG_DIR  <- file.path(BASE, "figures")
DATA_DIR <- file.path(BASE, "data")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

COL_A <- "#4C72B0"   # 癌旁正常
COL_B <- "#C44E52"   # 肿瘤

# ── 保留细胞 meta 的逐样本统计 ────────────────────────────────
meta <- fread(META_CSV)   # 列：barcode, sample, patient, tissue, n_umis, n_genes, pct_mt

# 与 02 归一化同口径：剔除 Scrublet predicted_doublet（threshold 0.25
# 已在 01_QC 打分时判定），只读原始 meta 不重算——保证本图与
# normalize_batch.R 的细胞集合一致
DB_CSV <- "01_QC/data/gse150290_doublet_scores.csv"
if (file.exists(DB_CSV)) {
  db <- fread(DB_CSV)
  meta <- meta[!barcode %in% db$barcode[db$predicted_doublet == TRUE]]
} else {
  cat("WARNING: doublet scores not found; showing all QC-pass cells\n")
}

M <- median(meta$n_umis)  # 仅信息参考线
cat(sprintf("Kept cells: %d, median library = %.0f\n", nrow(meta), M))

df <- meta[, .(n_cells    = .N,
               median_lib = as.numeric(median(n_umis)),
               p10_lib    = as.numeric(quantile(n_umis, 0.10)),
               p90_lib    = as.numeric(quantile(n_umis, 0.90)),
               p25_lib    = as.numeric(quantile(n_umis, 0.25)),
               p75_lib    = as.numeric(quantile(n_umis, 0.75)),
               tissue     = tissue[1]), by = sample]
df <- df[order(median_lib)]
df[, `:=`(short_name    = sub("^GSM[0-9]+_", "", sample),
          idx           = .I - 1,
          tissue_label  = factor(tissue, levels = c("A", "B"),
                                 labels = c("Adjacent normal (A)", "Tumor (B)")))]

cat(sprintf("Samples: %d; median library range %.0f - %.0f\n",
            nrow(df), min(df$median_lib), max(df$median_lib)))

TISSUE_COLS <- c("Adjacent normal (A)" = COL_A, "Tumor (B)" = COL_B)

theme_v <- theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", colour = NA))

# ── 图：文库分布 + 细胞数合并 ───────────────────────────────
p3_top <- ggplot(df) +
  geom_rect(aes(xmin = idx - 0.4, xmax = idx + 0.4, ymin = p10_lib, ymax = p90_lib,
                fill = tissue_label), alpha = 0.25) +
  geom_rect(aes(xmin = idx - 0.25, xmax = idx + 0.25, ymin = p25_lib, ymax = p75_lib,
                fill = tissue_label), alpha = 0.55) +
  geom_point(aes(idx, median_lib, colour = tissue_label), size = 1.2,
             shape = 21, fill = "white", stroke = 0.3) +
  geom_hline(yintercept = M, colour = "#333333", linewidth = 0.9,
             linetype = "dashed", alpha = 0.8) +
  annotate("text", x = nrow(df) - 0.5, y = M,
           label = sprintf("median = %s", format(M, big.mark = ",")),
           hjust = 1, vjust = -0.6, size = 3.2, fontface = "bold",
           colour = "#333333") +
  scale_fill_manual(values = TISSUE_COLS, name = NULL) +
  scale_colour_manual(values = TISSUE_COLS, name = NULL) +
  # y 锚定真 0 基线；无标题（y 轴标签已说明面板内容）
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Library size (UMI)") +
  theme_v

# 无 x 轴标题：逐样本身份无信息量，52 个标签纯噪音
p3_bot <- ggplot(df, aes(idx, n_cells)) +
  geom_bar(aes(fill = tissue_label), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_manual(values = TISSUE_COLS, name = NULL) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Cells") +
  theme_v

p3 <- p3_top / p3_bot + plot_layout(heights = c(2, 1))
ggsave(file.path(FIG_DIR, "library_and_cells_overview.png"), p3,
       width = 16, height = 8, dpi = 300, device = ragg::agg_png)
cat("Saved: library_and_cells_overview.png\n")

# ── 汇总 JSON ─────────────────────────────────────────────────
writeLines(toJSON(list(
  median_library = as.numeric(M), n_samples = nrow(df),
  library_median_range = as.numeric(range(df$median_lib)),
  tissue_breakdown = list(
    A = list(n_samples = sum(df$tissue == "A"),
             mean_median_lib = round(mean(df$median_lib[df$tissue == "A"]), 0)),
    B = list(n_samples = sum(df$tissue == "B"),
             mean_median_lib = round(mean(df$median_lib[df$tissue == "B"]), 0)))
), pretty = TRUE, auto_unbox = TRUE), file.path(DATA_DIR, "library_summary.json"))
cat("Saved: data/library_summary.json\nDone.\n")
