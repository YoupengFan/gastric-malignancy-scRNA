#!/usr/bin/env Rscript
# ============================================================
# 02_Normalize_Batch PC 验证图 — 渲染部分（Windows 原生 R）
#
# 用途：Harmony 矫正前后对比（2 面板）——A 批次混合直方图，
#       B 批次信号 per PC。
#
# 输入：02_Normalize_Batch/data/pc_validation.csv（per-PC 指标，
#       由 validate_pcs_compute.R 在 WSL 生成，含 pre/post 两套
#       系列）+ fnn_mixing_percell.csv（每细胞 k=20 同样本邻居
#       比例 pre/post，batch mixing 直方图用）
# Device：ragg::agg_png（Windows 原生抗锯齿）
#
# 关键设计要点：
#   - 论文风格：无面板标题、无图下 caption——面板仅保留左上角
#     A/B 标记（plot_annotation tag），定量数值与解释文字交给
#     论文图注（数据仍可直接由 pc_validation.csv 计算）
#   - A 面板：batch mixing 背靠背水平直方图（茎叶图式）——
#     per-cell 同样本邻居比例在竖轴、细胞计数在横轴，pre 向左
#     （琥珀）、post 向右（蓝），与 B 面板 pre/post 顺序一致；
#     random 期望（1/52，GSE150290 52 个 sample = Harmony 的
#     batch key）深红横虚线。矫正前分布沿竖轴铺开（批次隔离），
#     矫正后向 random 线收拢——直接量化矫正效果本身。图例置于
#     面板内右上角（post 侧高 ratio 处留白）。
#   - B 面板：批次信号 eta2(sample ~ PC) per PC，pre/post 成对
#     对比 —— 矫正后应大幅下降（Harmony 的 batch key 即 sample，
#     52 个 10x run）
#   - 颜色：Okabe-Ito 子集 —— 琥珀黄 #E69F00（矫正前）、蓝
#     #0072B2（矫正后），色盲安全
#
# 产出：02_Normalize_Batch/figures/pc_validation.png
#
# Usage (Windows native R 4.5.1)：
#   "E:\R-4.5.1\bin\Rscript.exe" 02_Normalize_Batch\scripts\validate_pcs_figure.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# ── 路径 ─────────────────────────────────────────────────────
BASE     <- "."
DATA_DIR <- file.path(BASE, "02_Normalize_Batch/data")
FIG_DIR  <- file.path(BASE, "02_Normalize_Batch/figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 读计算端输出 ────────────────────────────────────────────
res <- fread(file.path(DATA_DIR, "pc_validation.csv"))
fnn <- fread(file.path(DATA_DIR, "fnn_mixing_percell.csv"))

N_SAMPLES <- 52   # GSE150290 52 个 sample = Harmony batch key

PRE  <- "#E69F00"   # 矫正前：琥珀黄（Okabe-Ito）
POST <- "#0072B2"   # 矫正后：蓝

theme_pc <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "inside")

# ── (A) batch mixing：同样本邻居比例背靠背水平直方图 ─────────
# ratio 在竖轴、计数在横轴；pre 向左（负值、琥珀）、post 向右
# （蓝），与 B 面板顺序一致。共享 bin（以两分布最大值为界）；
# random 线 1/52 深红横虚线 + 右端标注（post 柱末端外留白，
# 不压数据；标注不入党例，图例只留给 pre/post 两个系列）。
# 图例置于面板内右上角——post
# 高 ratio 处几乎无细胞，正是留白区
bins <- seq(0, max(fnn$same_ratio_pre, fnn$same_ratio_post) * 1.05,
            length.out = 51)
bw    <- bins[2] - bins[1]
h_pre  <- hist(fnn$same_ratio_pre,  breaks = bins, plot = FALSE)
h_post <- hist(fnn$same_ratio_post, breaks = bins, plot = FALSE)
df_mix <- rbind(
  data.frame(y = h_pre$mids,  n = -h_pre$counts,  grp = "Pre-Harmony"),
  data.frame(y = h_post$mids, n =  h_post$counts, grp = "Post-Harmony"))
df_mix$grp <- factor(df_mix$grp, levels = c("Pre-Harmony", "Post-Harmony"))

pA <- ggplot(df_mix, aes(x = n, y = y, fill = grp)) +
  geom_col(orientation = "y", width = bw, colour = NA, alpha = 0.85) +
  geom_hline(yintercept = 1 / N_SAMPLES, colour = "#8B0000",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = Inf, y = 1 / N_SAMPLES + 0.02, hjust = 1,
           vjust = 0, colour = "#8B0000", size = 3.2,
           label = sprintf("Random (1/%d)", N_SAMPLES)) +
  scale_fill_manual(values = c("Pre-Harmony" = PRE, "Post-Harmony" = POST),
                    name = NULL) +
  scale_x_continuous(labels = abs, name = "Number of cells") +
  scale_y_continuous(limits = c(0, max(bins)),
                     breaks = seq(0, floor(max(bins) * 10) / 10, 0.1),
                     name = "Same-sample ratio (k=20)") +
  theme_pc + theme(legend.position.inside = c(0.98, 0.98),
                   legend.justification = c(1, 1))

# ── (B) batch signal per PC：eta2(sample) pre/post ───────────
df_eta <- rbind(
  data.frame(pc = res$pc, val = res$eta2_sample_pre,  grp = "Pre-Harmony"),
  data.frame(pc = res$pc, val = res$eta2_sample_post, grp = "Post-Harmony"))
df_eta$grp <- factor(df_eta$grp, levels = c("Pre-Harmony", "Post-Harmony"))

pB <- ggplot(df_eta, aes(pc, val, fill = grp)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Pre-Harmony" = PRE, "Post-Harmony" = POST),
                    name = NULL) +
  scale_x_continuous(breaks = c(1, seq(5, 30, 5))) +
  labs(x = "PC", y = "\u03b7\u00b2(sample ~ PC)") +
  theme_pc + theme(legend.position.inside = c(0.99, 0.99),
                   legend.justification = c(1, 1))

# ── 拼接 + 保存（各面板图例内置 + 左上角 A/B 标记）──────────
# A/B 两面板各自带内置图例（A 在右上角留白处、B 同），不再收集
# 到共享置顶图例——图例就放在数据旁边，读图不用来回找；面板无
# 标题、图下无 caption，仅保留左上角 A/B 标记，定量数值与解释
# 文字交给论文图注
p_all <- pA | pB
p_all <- p_all + plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 14, face = "bold")))
ggsave(file.path(FIG_DIR, "pc_validation.png"), p_all,
       width = 13, height = 5.3, dpi = 300, device = ragg::agg_png)
cat(sprintf("Saved: %s\n", file.path(FIG_DIR, "pc_validation.png")))
