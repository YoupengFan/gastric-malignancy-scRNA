#!/usr/bin/env Rscript
# ============================================================
# 03_Cluster_Annotate 注释 QC — 图件部分（Windows R）
# 读 annotation_qc_export.R（WSL）产出的 annotation_qc_*.csv；
# 零重算。所有图用 ragg 渲染。
#
# 产出（03_Cluster_Annotate/figures/）：
#   qc_marker_violin.png            代表 marker 逐类型 violin
#                                   （Zhang 2021 Fig D 风格）
#   qc_marker_dotplot.png           逐基因检出率 × 表达 DotPlot
#                                   （全部面板基因，Seurat DotPlot 风格）
#   qc_umap_feature.png             UMAP feature 图，12 面板（3×4），
#                                   Zhang 2021 Fig E 风格
#   umap_annotation.png             按最终注释着色的 UMAP，violin 配色
#   qc_suspected_contam.png         补充：疑似跨谱系污染 flag（灰+红）
#   qc_umap_clean.png               补充配套：剔除疑似细胞的注释总览
#                                   （带簇名标签）
#   qc_umap_clean_nolabel.png       同视图不带簇标签
#   qc_confusion_heatmap.png        投票 × 簇 mode，行归一化
#   qc_confusion_heatmap_clean.png  同上但剔除疑似污染细胞
#                                   （蓝色 ramp，与红色全数据版区分）
#
# Usage (Windows 原生 R 4.5.1)：
#   "E:\R-4.5.1\bin\Rscript.exe" 03_Cluster_Annotate\scripts\annotation_qc_figures.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE     <- "03_Cluster_Annotate"
DATA_DIR <- file.path(BASE, "data")
FIG_DIR  <- file.path(BASE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

stats <- fread(file.path(DATA_DIR, "annotation_qc_stats.csv"))
type_ord <- stats$annotation   # 丰度序（大者在前）

theme_qc <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          plot.caption = element_text(size = 7, hjust = 0))
}

# ══ 1) Violin 数据（供代表 marker violin 使用）═════════════
v <- fread(file.path(DATA_DIR, "annotation_qc_violin.csv"))

# ══ 2) 代表 marker violin（每类型一个基因）═══════════════
# Zhang 2021 Fig D 风格：每个注释类型一个规范 marker，跨 12 类
# 展示。行 = marker 基因（丰度序），strip 标签仅基因名（省略 Y
# 刻度/标题，log1p 由论文图注解释）。每个 violin 按所在列的细
# 胞类型着色——每行一色一类——用共享的 VIOLIN_COLS 配色，保证
# violin 填充与 UMAP 点的类型颜色始终一致。
REP_MARKER <- c(T_cells = "CD3D", B_cells = "MS4A1", NK_cells = "NKG7",
                Myeloid = "CD14", Mast_cells = "TPSAB1",
                Plasma_cells = "MZB1", Fibroblasts = "COL1A1",
                Endothelial = "VWF", Epithelial = "EPCAM",
                Enteroendocrine = "CHGA", Parietal = "ATP4A",
                Chief = "PGC")
# VIOLIN_COLS — violin、注释、clean UMAP 三图共享的 12 色分类配色。
VIOLIN_COLS <- c(Plasma_cells = "#00BCD8", Epithelial = "#eda100",
                 T_cells = "#1baf7a", Chief = "#e87ba4",
                 B_cells = "#eb6834", Fibroblasts = "#008300",
                 Myeloid = "#4a3aa7", Endothelial = "#e34948",
                 Enteroendocrine = "#8b5a2b", Mast_cells = "#17a2b8",
                 NK_cells = "#9467bd", Parietal = "#F0E442")
rep_dt <- data.table(ctype = names(REP_MARKER), gene = unname(REP_MARKER))
vr <- merge(v[gene %in% unname(REP_MARKER)], rep_dt, by = "gene")
vr[, annotation := factor(annotation, levels = type_ord)]
vr[, gene := factor(gene, levels = unname(REP_MARKER[type_ord]))]
p_rep <- ggplot(vr, aes(annotation, expr)) +
  geom_violin(aes(fill = annotation), scale = "width",
              linewidth = 0.15, colour = "#454440") +
  scale_fill_manual(values = VIOLIN_COLS, guide = "none") +
  facet_grid(gene ~ ., scales = "free_y", switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_qc() +
  theme(strip.text.y.left = element_text(angle = 0, size = 8,
                                         face = "italic"),
        strip.background = element_blank(),
        panel.spacing = unit(0.4, "lines"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(file.path(FIG_DIR, "qc_marker_violin.png"), p_rep,
       width = 7, height = 7.5, dpi = 300, device = ragg::agg_png)
cat("Saved: qc_marker_violin.png\n")

# ══ 3) marker 审计 DotPlot（全部面板基因 × 12 类）═════════
# 规范 "marker 是否有代表性" 图（Seurat DotPlot 风格）：行 = 每
# 个面板基因（按面板分组、字典序），点大小 = raw counts 检出比
# 例，颜色 = mean log1p 表达。有效 marker = 自身类型上的大深色
# 点；对角线外的小/淡点 = 高表达基因的 ambient RNA（如 PGC）
# 或跨型泄漏。无标题/图注：解读交给论文图注。
ms <- fread(file.path(DATA_DIR, "annotation_qc_marker_stats.csv"))
ms[, annotation := factor(annotation, levels = type_ord)]
ms[, gene := factor(gene, levels = unique(gene))]     # 字典序
ms[, panel := factor(panel, levels = unique(panel))]
p_mk <- ggplot(ms, aes(annotation, gene)) +
  geom_point(aes(size = det_pct, colour = mean_log1p), alpha = 0.9) +
  scale_size_area(max_size = 4.5, name = "cells\ndetected (%)") +
  scale_colour_gradient(low = "#f0efe9", high = "#a31621",
                        name = "mean\nlog1p") +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_qc(8) +
  theme(axis.text.y = element_text(size = 7, face = "italic"),
        strip.text.y.left = element_text(angle = 0, size = 7),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))
ggsave(file.path(FIG_DIR, "qc_marker_dotplot.png"), p_mk,
       width = 7.5, height = 10.5, dpi = 300, device = ragg::agg_png)
cat("Saved: qc_marker_dotplot.png\n")

# 控制台审计：最弱 marker（自身检出 ≤ 1.5× 最佳其他）
mk_sum <- fread(file.path(DATA_DIR, "annotation_qc_marker_summary.csv"))
cat("\n=== marker audit: genes with detection specificity < 1.5 ===\n")
print(mk_sum[spec_det < 1.5][order(spec_det),
      .(panel, gene, det_own, det_max_other, spec_det, expr_ratio)])

# ══ 4) UMAP feature 图（Zhang 2021 Fig E 风格）═══════════
# 12 面板（3 行 × 4 列），每注释类型一个：浅灰细胞云 + 该类型
# 代表 marker（REP_MARKER——与 qc_marker_violin.png 同集）的
# 红色表达梯度。面板 strip = 细胞类型；基因名斜体放面
# 板右上角（左上角在 10K 布局中被簇遮住）。无轴、无刻度、无图
# 例：视觉焦点只留红色信号。
TYPE_PRETTY <- c(T_cells = "T cells", B_cells = "B cells", NK_cells = "NK cells",
                 Myeloid = "Myeloid", Mast_cells = "Mast cells",
                 Plasma_cells = "Plasma cells", Fibroblasts = "Fibroblasts",
                 Endothelial = "Endothelial", Epithelial = "Epithelial",
                 Enteroendocrine = "Enteroendocrine", Parietal = "Parietal",
                 Chief = "Chief")
um <- fread(file.path(DATA_DIR, "annotation_qc_umap.csv"),
            select = c("umap_1", "umap_2", unname(REP_MARKER)))
um_long <- melt(um, id.vars = c("umap_1", "umap_2"),
                measure.vars = unname(REP_MARKER),
                variable.name = "gene", value.name = "expr")
um_long <- merge(um_long, rep_dt, by = "gene")   # 基因 → 类型
um_long[, type := factor(ctype, levels = type_ord)]
um_long <- um_long[order(expr)]                  # 零点先画（垫底）
# 逐基因着色截断：非零表达的 80 分位（0.90 太严：着色细胞过
# 少）。低于截断的细胞保持灰色——只有真实表达
# 区域拿到粉/红渐变（无全局淡粉洗色；灰/红对比干净）。分位即
# 调谐旋钮：调高着色细胞变少，调低变多。
cut_dt <- um_long[expr > 0, .(cut = quantile(expr, 0.80)), by = gene]
um_long <- merge(um_long, cut_dt, by = "gene")
lab_dt <- um_long[, .(x = max(umap_1), y = max(umap_2)), by = type]
lab_dt <- merge(lab_dt, unique(rep_dt[, .(ctype, gene)]),
                by.x = "type", by.y = "ctype")
p_um <- ggplot(um_long, aes(umap_1, umap_2)) +
  geom_point(colour = "#d9d7d1", size = 0.08, alpha = 0.5) +
  geom_point(data = um_long[expr >= cut], aes(colour = expr), size = 0.14) +
  scale_colour_gradient(low = "#ffe3de", high = "#ff1e1e", guide = "none") +
  geom_text(data = lab_dt, aes(x = x, y = y, label = gene),
            hjust = 1, vjust = 1, fontface = "italic", size = 3,
            inherit.aes = FALSE) +
  facet_wrap(~ type, ncol = 4, labeller = labeller(type = TYPE_PRETTY)) +
  theme_bw(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 9.5, face = "bold"),
        panel.spacing = unit(0.5, "lines"))
ggsave(file.path(FIG_DIR, "qc_umap_feature.png"), p_um,
       width = 10, height = 7.5, dpi = 300, device = ragg::agg_png)
cat("Saved: qc_umap_feature.png\n")

# ══ 5) UMAP 注释总览 ═════════════════════════════════════
# ragg 300 dpi、系统字体，VIOLIN_COLS 着色——每类一色，与
# qc_marker_violin.png 一致。
# 绘制顺序 = 丰度序（type_ord，大者在前），稀有类型后画置顶、
# 保持可见。混合细胞（投票 ≠ 簇 mode）如实绘制：总览展示的是
# 逐细胞投票，非簇纯度。
um_a <- fread(file.path(DATA_DIR, "annotation_qc_umap.csv"),
              select = c("umap_1", "umap_2", "annotation"))
um_a[, annotation := factor(annotation, levels = type_ord)]
p_ann <- ggplot(um_a, aes(umap_1, umap_2, colour = annotation)) +
  geom_point(size = 0.08, alpha = 0.7) +
  scale_colour_manual(values = VIOLIN_COLS, name = NULL,
                      guide = guide_legend(override.aes = list(size = 2.5,
                                                               alpha = 1))) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(3, "mm"))
ggsave(file.path(FIG_DIR, "umap_annotation.png"), p_ann,
       width = 8, height = 6, dpi = 300, device = ragg::agg_png)
cat("Saved: umap_annotation.png\n")

# ══ 6) 疑似跨谱系污染（补充图）════════════════════════════
# 灰细胞云 + suspected_contam 细胞（投票与簇 mode 在不同谱系：
# epi/immune/stroma）红色 flag。与 qc_umap_feature.png 同一色彩
# 语言（灰云、红信号）；补充图：展示 flag 本身；qc_umap_clean.png（下
# 节）是剔除这些细胞后的配套视图。
um_s <- fread(file.path(DATA_DIR, "annotation_qc_umap.csv"),
              select = c("umap_1", "umap_2", "suspected_contam"))
um_s <- um_s[order(suspected_contam)]   # 灰先画，flag 置顶
p_cont <- ggplot(um_s, aes(umap_1, umap_2)) +
  geom_point(data = um_s[suspected_contam == FALSE],
             colour = "#d9d7d1", size = 0.08, alpha = 0.5) +
  geom_point(data = um_s[suspected_contam == TRUE],
             colour = "#ff1e1e", size = 0.14, alpha = 0.8) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank())
ggsave(file.path(FIG_DIR, "qc_suspected_contam.png"), p_cont,
       width = 8, height = 6, dpi = 300, device = ragg::agg_png)
cat("Saved: qc_suspected_contam.png\n")

# ══ 7) Clean UMAP：剔除疑似污染细胞 ═══════════════════════
# 与注释总览同视图，但不绘制 suspected_contam 细胞——只展示投票与
# 簇身份跨谱系一致的细胞。qc_suspected_contam.png（flag 图）的
# 补充配套；非裁决，被移除细胞仍留在数据与所有下游分析中。
um_c <- fread(file.path(DATA_DIR, "annotation_qc_umap.csv"),
              select = c("umap_1", "umap_2", "annotation",
                         "suspected_contam"))
um_c <- um_c[suspected_contam == FALSE]
um_c[, annotation := factor(annotation, levels = type_ord)]
# 簇标签：每注释的中位 UMAP 位置（对残余散点稳健），英文名来自
# TYPE_PRETTY。纯墨色文字（无光晕）。
lab_c <- um_c[, .(x = median(umap_1), y = median(umap_2)), by = annotation]
lab_c[, label := TYPE_PRETTY[as.character(annotation)]]
# Myeloid 簇色（#4a3aa7）与墨色（#454440）太近，不适合中心标
# 签；按固定 UMAP y 偏移（数据尺度单位）把标签上移到簇心上方。
MYELOID_LABEL_DY <- 2.0
lab_c[annotation == "Myeloid",
      y := um_c[annotation == "Myeloid", median(umap_2)] + MYELOID_LABEL_DY]
p_clean_base <- ggplot(um_c, aes(umap_1, umap_2, colour = annotation)) +
  geom_point(size = 0.08, alpha = 0.7) +
  scale_colour_manual(values = VIOLIN_COLS, name = NULL,
                      guide = guide_legend(override.aes = list(size = 2.5,
                                                               alpha = 1))) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(3, "mm"))
# 两个导出：带标签版存为 qc_umap_clean.png；无标签视图同样保留。
p_clean <- p_clean_base +
  geom_text(data = lab_c, aes(x = x, y = y, label = label),
            colour = "#454440", size = 3.2, fontface = "bold",
            inherit.aes = FALSE)
ggsave(file.path(FIG_DIR, "qc_umap_clean.png"), p_clean,
       width = 8, height = 6, dpi = 300, device = ragg::agg_png)
cat("Saved: qc_umap_clean.png (labelled)\n")
ggsave(file.path(FIG_DIR, "qc_umap_clean_nolabel.png"), p_clean_base,
       width = 8, height = 6, dpi = 300, device = ragg::agg_png)
cat("Saved: qc_umap_clean_nolabel.png\n")

# ══ 8) 混淆热图（投票 × 簇 mode）══════════════════════════
# 完整 12×12 网格（含零计数格）。标签 = 行比例，填充白（0.00）
# → 深红（1.00），论文风格参考图。
conf <- fread(file.path(DATA_DIR, "annotation_qc_confusion.csv"))
conf[, vote := factor(vote, levels = rev(type_ord))]
conf[, cluster_mode := factor(cluster_mode, levels = type_ord)]
grid <- CJ(vote = levels(conf$vote), cluster_mode = levels(conf$cluster_mode))
conf <- merge(grid, conf, by = c("vote", "cluster_mode"), all.x = TRUE)
conf[is.na(N), N := 0L]
conf[, prop := N / sum(N), by = vote]
conf[, txt_col := ifelse(prop > 0.45, "white", "#454440")]
p_conf <- ggplot(conf, aes(cluster_mode, vote)) +
  geom_tile(aes(fill = prop), colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", prop), colour = txt_col), size = 2.1) +
  scale_fill_gradient(low = "white", high = "#8b0000",
                      limits = c(0, 1), name = "Row proportion") +
  scale_colour_identity() +
  labs(x = "cluster-mode annotation (cluster plurality)",
       y = "per-cell vote (final annotation)") +
  theme_qc(8) +
  theme(axis.text.y = element_text(size = 7),
        legend.position = "right")
ggsave(file.path(FIG_DIR, "qc_confusion_heatmap.png"), p_conf,
       width = 7.5, height = 6.5, dpi = 300, device = ragg::agg_png)
cat("Saved: qc_confusion_heatmap.png\n")

# ══ 9) 同热图，剔除疑似污染细胞 ══════════════════════════
# 跨谱系疑似细胞剔除后的残余投票 × 簇 mode 混合（与
# qc_umap_clean.png 配套）。蓝色 ramp 与红色全数据版视觉区分。
# 由逐细胞导出（annotation + cluster mode + flag）计算（与全集
# 的 qc_confusion.csv 口径逐位一致），无需单独导出。
um_cc <- fread(file.path(DATA_DIR, "annotation_qc_umap.csv"),
               select = c("annotation", "cluster_annotation",
                          "suspected_contam"))
um_cc <- um_cc[suspected_contam == FALSE]
conf_cl <- um_cc[, .N, by = .(vote = annotation,
                              cluster_mode = cluster_annotation)]
conf_cl[, vote := factor(vote, levels = rev(type_ord))]
conf_cl[, cluster_mode := factor(cluster_mode, levels = type_ord)]
grid_cl <- CJ(vote = levels(conf_cl$vote),
              cluster_mode = levels(conf_cl$cluster_mode))
conf_cl <- merge(grid_cl, conf_cl, by = c("vote", "cluster_mode"),
                 all.x = TRUE)
conf_cl[is.na(N), N := 0L]
conf_cl[, prop := N / sum(N), by = vote]
conf_cl[, txt_col := ifelse(prop > 0.45, "white", "#454440")]
p_conf_cl <- ggplot(conf_cl, aes(cluster_mode, vote)) +
  geom_tile(aes(fill = prop), colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", prop), colour = txt_col), size = 2.1) +
  scale_fill_gradient(low = "white", high = "#0b3d91",
                      limits = c(0, 1), name = "Row proportion") +
  scale_colour_identity() +
  labs(x = "cluster-mode annotation (cluster plurality)",
       y = "per-cell vote (final annotation)") +
  theme_qc(8) +
  theme(axis.text.y = element_text(size = 7),
        legend.position = "right")
ggsave(file.path(FIG_DIR, "qc_confusion_heatmap_clean.png"), p_conf_cl,
       width = 7.5, height = 6.5, dpi = 300, device = ragg::agg_png)
cat("Saved: qc_confusion_heatmap_clean.png\n")
