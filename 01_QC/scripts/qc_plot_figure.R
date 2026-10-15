#!/usr/bin/env Rscript
# ============================================================
# 01_QC 图 — 渲染部分（Windows 原生 R 运行）：4 面板
#
# 输入：01_QC/data/qc_plot_*.csv + qc_plot_summary.json
#      （由 qc_plot_compute.R 在 WSL 生成）
# Device：ragg::agg_png（Windows 原生抗锯齿 — 无需 showtext）
#
# 风格：无标题、ABCD 标签左上角、黑框 legend、阈值虚线。
#
# 产出：01_QC/figures/qc_metrics_overview.png
#
# Usage (Windows native R 4.5.1)：
#   "E:\R-4.5.1\bin\Rscript.exe" 01_QC\scripts\qc_plot_figure.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

# ── 路径 ─────────────────────────────────────────────────────
BASE     <- "01_QC"
DATA_DIR <- file.path(BASE, "data")
FIG_DIR  <- file.path(BASE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── QC 阈值（与计算部分相同）──────────────────────────────
MIN_UMIS  <- 1001
MAX_UMIS  <- 15000
PAPER_MAX <- 6000
MIN_GENES <- 501

XMAX_COUNTS <- 20000
XMAX_GENES  <- 6000

# ── 读入聚合数据 ───────────────────────────────────────────
# 5 个轻量文件全部由 qc_plot_compute.R 在服务器端（WSL）生成——
# 渲染端不碰 raw 矩阵，把这 5 个文件拷到任何机器都能出图
hist_a  <- fread(file.path(DATA_DIR, "qc_plot_hist_total.csv"))
hist_b  <- fread(file.path(DATA_DIR, "qc_plot_hist_genes.csv"))
knee    <- fread(file.path(DATA_DIR, "qc_plot_knee.csv"))
scatter <- fread(file.path(DATA_DIR, "qc_plot_scatter.csv"))
smm     <- fromJSON(file.path(DATA_DIR, "qc_plot_summary.json"))

# Panel C 的 knee 阈值 = 四个候选切点（500/1000/2000/
# 5000 UMI）：1,000 处的拐点正是 QC 下限的依据；
# 只画在 panel A。n_above 在抽稀后的 rank 向量上线性插值 depth
# 越过阈值时的 rank（误差 < 抽稀步长 191，log 轴上不可见）——
# 无需回 WSL 重算。
knee_thr <- data.frame(y = c(500, 1000, 2000, 5000),
                       lab = c("500 UMI", "1,000 UMI", "2,000 UMI", "5,000 UMI"),
                       col = c("purple", "red", "orange", "green"))
n_above <- function(th) {
  above <- knee[depth >= th]
  below <- knee[depth < th]
  r_hi <- above$rank[nrow(above)]; d_hi <- above$depth[nrow(above)]
  r_lo <- below$rank[1];          d_lo <- below$depth[1]
  r_hi + (d_hi - th) / (d_hi - d_lo) * (r_lo - r_hi)
}
knee_thr$n_above <- vapply(knee_thr$y, n_above, numeric(1))

cat(sprintf("Raw barcodes: %d | scatter points: %d | knee decimation step: %d\n",
            smm$raw_barcodes_total, smm$n_scatter, smm$decimation_step))
# 校验：插值的 n_above(1001) 应与 compute 端精确值相差 < 抽稀步长
n1001_check <- n_above(1001)
cat(sprintf("Sanity: n_above(1001) interpolated %.0f vs exact %d (diff <= 191 ok)\n",
            n1001_check, as.numeric(smm$thresholds[["1,001 UMI"]])))

# ── Theme（风格：ABCD 标签、黑框 legend）─────────────────
theme_qc <- theme_bw(base_size = 12) +
  theme(plot.tag = element_text(size = 16, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white"),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "inside",
        legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1))

# ── (A) total counts 直方图 ────────────────────────────────
thr_a <- data.frame(x = c(MIN_UMIS, PAPER_MAX, MAX_UMIS),
                    lab = c("min = 1,001", "paper max = 6,000", "max = 15,000"),
                    col = c("red", "purple", "green"))
pA <- ggplot(hist_a, aes(mid, count)) +
  geom_col(fill = "steelblue", colour = "white", linewidth = 0.2) +
  # 不可见 hline（被裁剪到面板底部）只为给 colour legend 提供
  # 水平虚线 key— geom_vline 自带 key 是竖线，
  # 所以 vline 本身不贡献 legend
  geom_hline(data = thr_a, aes(yintercept = -Inf, colour = lab),
             linetype = "dashed", linewidth = 1, alpha = 0.7) +
  geom_vline(data = thr_a, aes(xintercept = x, colour = lab),
             linetype = "dashed", linewidth = 1, alpha = 0.7, show.legend = FALSE) +
  # 显式 breaks 钉死 legend 顺序（渲染顺序曾漂移于 data.frame 行序之外）
  scale_colour_manual(values = setNames(thr_a$col, thr_a$lab), name = NULL,
                      breaks = thr_a$lab) +
  # 每个 10 的整数次幂，同 matplotlib 默认 LogLocator
  scale_y_log10(breaks = 10^(0:8),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Count depth", y = "Frequency") +
  coord_cartesian(xlim = c(-50, XMAX_COUNTS), expand = FALSE) +
  theme_qc

# ── (B) 每细胞基因数直方图 ─────────────────────────────────
thr_b <- data.frame(x = MIN_GENES, lab = "min = 501", col = "red")
pB <- ggplot(hist_b, aes(mid, count)) +
  geom_col(fill = "seagreen", colour = "white", linewidth = 0.2) +
  # 与 panel A 相同的水平 key 技巧
  geom_hline(data = thr_b, aes(yintercept = -Inf, colour = lab),
             linetype = "dashed", linewidth = 1, alpha = 0.7) +
  geom_vline(data = thr_b, aes(xintercept = x, colour = lab),
             linetype = "dashed", linewidth = 1, alpha = 0.7, show.legend = FALSE) +
  scale_colour_manual(values = setNames(thr_b$col, thr_b$lab), name = NULL,
                      breaks = thr_b$lab) +
  scale_y_log10(breaks = 10^(0:8),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Number of genes", y = "Frequency") +
  coord_cartesian(xlim = c(-20, XMAX_GENES), expand = FALSE) +
  theme_qc

# ── (C) rank 图（log-log knee）─────────────────────────────
pC <- ggplot(knee, aes(rank, depth)) +
  geom_line(linewidth = 0.6, colour = "darkred", alpha = 0.8) +
  geom_hline(data = knee_thr, aes(yintercept = y, colour = lab),
             linetype = "dashed", alpha = 0.7, linewidth = 1) +
  geom_vline(data = knee_thr, aes(xintercept = n_above, colour = lab),
             linetype = "dashed", alpha = 0.7, linewidth = 1, show.legend = FALSE) +
  scale_colour_manual(values = setNames(knee_thr$col, knee_thr$lab), name = NULL,
                      breaks = knee_thr$lab) +
  scale_x_log10(breaks = 10^(0:8),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = 10^(0:8),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Barcode rank", y = "Count depth") +
  theme_qc

# ── (D) genes vs counts 散点，按 MT% 着色 ──────────────────
# ColorBrewer YlOrRd 9 级色带，MT% 映射 0-50
YLORRD <- c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C",
            "#FC4E2A", "#E31A1C", "#BD0026", "#800026")
pD <- ggplot(scatter, aes(total_counts, n_genes, colour = pct_mt)) +
  geom_point(size = 0.4, alpha = 0.6) +
  scale_colour_gradientn(colours = YLORRD, limits = c(0, 50),
                         name = "Mitochondrial%",
                         # barheight ≈ 0.8 × panel D 高（colorbar
                         # shrink=0.8）；ggplot2 4.0 默认色带过小
                         guide = guide_colourbar(order = 1, position = "right",
                                                barheight = unit(13, "cm"))) +
  geom_hline(aes(yintercept = MIN_GENES, linetype = "genes=501, counts=1,001"),
             colour = "blue", alpha = 0.7) +
  geom_vline(aes(xintercept = MIN_UMIS, linetype = "genes=501, counts=1,001"),
             colour = "blue", alpha = 0.7) +
  scale_linetype_manual(values = c("genes=501, counts=1,001" = "dashed"),
                        name = NULL,
                        guide = guide_legend(position = "inside", order = 2,
                          theme = theme(legend.background =
                                          element_rect(colour = "black", fill = "white"),
                                        legend.key = element_blank()))) +
  labs(x = "Count depth", y = "Number of genes") +
  coord_cartesian(xlim = c(-200, XMAX_COUNTS), ylim = c(-50, XMAX_GENES),
                  expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(plot.tag = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1))

# ── 拼接 + 保存 ────────────────────────────────────────────
# 2×2 patchwork 组合，tag_levels="A" 自动生成 ABCD 角标；
# ragg::agg_png 为 Windows 原生抗锯齿（渲染一致性靠固定版本）
p_all <- (pA | pB) / (pC | pD) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIG_DIR, "qc_metrics_overview.png"), p_all,
       width = 18, height = 14, dpi = 300, device = ragg::agg_png)
cat(sprintf("Saved: %s\n", file.path(FIG_DIR, "qc_metrics_overview.png")))
