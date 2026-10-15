#!/usr/bin/env Rscript
# ============================================================
# 02_Normalize_Batch UMAP 批次矫正演示图 — 渲染部分（Windows 原生 R）
#
# 用途：Harmony 批次矫正效果示意图（before/after 双面板）。
#
# 输入：02_Normalize_Batch/data/umap_batch.csv（由 umap_batch_compute.R
#       在 WSL 生成；本脚本只做过滤与渲染，不碰 raw 数据）
# Device：ragg::agg_png（Windows 原生抗锯齿）
#
# 关键设计要点：
#   - 仅上皮谱系（Epithelial/Chief/Enteroendocrine/Parietal）
#   - 仅 A 侧（adjacent-normal）组织：同组织跨患者比较，斑块分离
#     就是纯批次效应，无 A↔B 瘤内生物学混淆；
#   - 患者 Pat01/Pat17/Pat21
#   - 原数据直出：不做任何逐细胞过滤或抽样；demo 子集仅按谱系 × 组织 × 患者三条声明性准则圈定
#   - 判读标准：before = 同色斑块（批次效应）；after = 颜色均匀
#     混合（批次结构消除）。本数据批次效应本身温和（同实验室、
#     同平台），after 面板保留少量患者结构是如实呈现
#
# 产出：02_Normalize_Batch/figures/umap_batch_before_after.png
#
# Usage (Windows native R 4.5.1)：
#   "E:\R-4.5.1\bin\Rscript.exe" 02_Normalize_Batch\scripts\umap_batch_figure.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# ── 路径 ─────────────────────────────────────────────────────
BASE     <- "02_Normalize_Batch"
DATA_DIR <- file.path(BASE, "data")
FIG_DIR  <- file.path(BASE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── demo 子集：上皮谱系 × A 侧组织 × 3 名示例患者 ────────────
d <- fread(file.path(DATA_DIR, "umap_batch.csv"))
cat(sprintf("Loaded %d cells, %d patients, %d annotations\n",
            nrow(d), length(unique(d$patient)), length(unique(d$annotation))))

EPI_TYPES <- c("Epithelial", "Chief", "Enteroendocrine", "Parietal")
SEL_PATIENTS <- c("Pat01", "Pat17", "Pat21")
d_sel <- d[annotation %in% EPI_TYPES &
           tissue == "A" &
           patient %in% SEL_PATIENTS]
cat(sprintf("Demo subset: %d cells\n", nrow(d_sel)))

# 三种最大区分、色盲安全的 Okabe-Ito 子集色
PAT_COLS <- setNames(c("#0072B2", "#D55E00", "#009E73"), SEL_PATIENTS)

# ── 共享样式 ────────────────────────────────────────────────
# 图例 = 每个面板右上角内侧（A 侧 = adjacent-normal，非"正常组织"）；
# 放大图例点（override.aes）并去掉 alpha，白底黑框保证覆盖在散点上可读
theme_umap <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.background = element_rect(colour = "black", fill = "white"),
        legend.key = element_blank(),
        legend.text = element_text(size = 10))

# ── Before 面板（pre-Harmony PCA 上的 UMAP）─────────────────
p_pre <- ggplot(d_sel, aes(umap1_pre, umap2_pre)) +
  geom_point(aes(colour = patient), size = 0.7, alpha = 0.65) +
  scale_colour_manual(values = PAT_COLS, name = NULL,
                      breaks = SEL_PATIENTS) +
  labs(x = "UMAP1", y = "UMAP2",
       title = "Before Harmony (epithelial, adjacent-normal tissue)") +
  theme_umap +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))

# ── After 面板（Harmony embedding 上的 UMAP）────────────────
p_post <- ggplot(d_sel, aes(umap1_post, umap2_post)) +
  geom_point(aes(colour = patient), size = 0.7, alpha = 0.65) +
  scale_colour_manual(values = PAT_COLS, name = NULL,
                      breaks = SEL_PATIENTS) +
  labs(x = "UMAP1", y = "UMAP2",
       title = "After Harmony (epithelial, adjacent-normal tissue)") +
  theme_umap +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))

# ── 拼接 + 保存（各面板自带图例，右上角内侧）────────────────
p_all <- p_pre | p_post
ggsave(file.path(FIG_DIR, "umap_batch_before_after.png"), p_all,
       width = 13, height = 6.4, dpi = 300, device = ragg::agg_png)
cat(sprintf("Saved: %s\n", file.path(FIG_DIR, "umap_batch_before_after.png")))
