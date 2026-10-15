#!/usr/bin/env Rscript
# ============================================================
# 02_Normalize_Batch 批次 UMAP demo 选例诊断 — diag_batch_demo.R
#
# 用途：为 umap_batch_figure.R 的 before/after 示例患者选择提供
#       定量依据。5 组指标全部零重算读 umap_batch.csv：
#         1) 患者内 A/B 组织质心距离 pre/post（配对）
#         2) 样本级（52 个 10x run）离全局质心分离度 sep_pre
#         3) 上皮谱系子集患者斑块分离度 sep_pre/post + clump
#         4) 上皮 × 单组织（A、B 分开）患者分离度——纯批次
#            斑块，无 A↔B 组织混淆
#         5) A 侧患者 clump 排名（d_pre/s_pre）——Pat01/Pat17/
#            Pat21 入选 demo 图的依据
#
# 输入：02_Normalize_Batch/data/umap_batch.csv（umap_batch_compute.R
#       在 WSL 生成，列 = umap1/2_pre/post + 样本/患者/组织/注释）
# 产出：无文件——指标打印到控制台；患者选择结论已固化进
#       umap_batch_figure.R 的 SEL_PATIENTS 与 header 说明
#
# 关键设计要点：
#   - 零重算 + <1 min：只读 UMAP 坐标 CSV，data.table 即可
#   - sep_pre = 患者/样本质心离全局质心的中位距离 ÷ 全局散布
#     尺度；clump = 质心偏离度 ÷ 患者内散布（大 = 明显斑块）
#   - Pat13 A 侧上皮 sep_pre 第 1 但被排除：A 侧上皮 55.7% 为
#     嫌疑污染（全患者最高，Pat01/17/21 均 <4%），详见
#     03_Cluster_Annotate/scripts/
#     contamination_analysis/NOTES_suspected_contamination.md
#
# Usage (Windows native R 4.5.1)：
#   "E:\R-4.5.1\bin\Rscript.exe" 02_Normalize_Batch\scripts\diagnostics\diag_batch_demo.R
# ============================================================
suppressPackageStartupMessages(library(data.table))

BASE <- "02_Normalize_Batch"
d <- fread(file.path(BASE, "data", "umap_batch.csv"))
cat("== annotations ==\n"); print(d[, .N, by = annotation][order(-N)])
cat("\n== tissue values ==\n"); print(d[, .N, by = tissue])

# 全局散布尺度（pre/post 各自归一化的分母）
glob <- function(x1, y1, x2, y2) list(
  gx = median(x1), gy = median(y1),
  s_all = median(sqrt((x1 - median(x1))^2 + (y1 - median(y1))^2)),
  gxp = median(x2), gyp = median(y2),
  s_allp = median(sqrt((x2 - median(x2))^2 + (y2 - median(y2))^2)))
g <- glob(d$umap1_pre, d$umap2_pre, d$umap1_post, d$umap2_post)

# ── 1) 患者内 A/B 组织质心距离（配对）────────────────────
cat("\n== A/B centroid distance per patient (paired) ==\n")
tissues <- sort(unique(d$tissue))
if (length(tissues) == 2) {
  ab <- rbindlist(lapply(unique(d$patient), function(p) {
    sub <- d[patient == p]
    if (length(unique(sub$tissue)) < 2) return(NULL)
    a <- sub[tissue == tissues[1]]; b <- sub[tissue == tissues[2]]
    d_pre  <- sqrt((median(a$umap1_pre)  - median(b$umap1_pre))^2 +
                   (median(a$umap2_pre)  - median(b$umap2_pre))^2)
    d_post <- sqrt((median(a$umap1_post) - median(b$umap1_post))^2 +
                   (median(a$umap2_post) - median(b$umap2_post))^2)
    data.table(patient = p, nA = nrow(a), nB = nrow(b),
               d_pre = d_pre, d_post = d_post, ratio = d_post / d_pre)
  }))
  print(ab[order(-d_pre)])
  cat(sprintf("median d_pre = %.2f, median d_post = %.2f\n",
              median(ab$d_pre), median(ab$d_post)))
} else {
  cat(sprintf("tissues = %s (skip A/B pairing)\n",
              paste(tissues, collapse = ", ")))
}

# ── 2) 样本级（52 个 10x run）分离度 sep_pre ──────────────
cat("\n== sample-level sep_pre (52 samples, whole data) top 12 ==\n")
sm <- d[, .(n = .N,
            d_pre  = median(sqrt((umap1_pre  - g$gx)^2  + (umap2_pre  - g$gy)^2)),
            d_post = median(sqrt((umap1_post - g$gxp)^2 + (umap2_post - g$gyp)^2))),
        by = sample]
sm[, sep_pre := d_pre / g$s_all]
sm[, sep_post := d_post / g$s_allp]
print(sm[order(-sep_pre)][1:12])
cat(sprintf("sample sep_pre: min %.2f | med %.2f | max %.2f\n",
            min(sm$sep_pre), median(sm$sep_pre), max(sm$sep_pre)))

# ── 3) 上皮谱系子集，患者斑块 ─────────────────────────────
epi_types <- c("Chief", "Epithelial", "Enteroendocrine", "Parietal")
d_epi <- d[annotation %in% epi_types]
cat(sprintf("\n== epi subset: %d cells ==\n", nrow(d_epi)))
ge <- glob(d_epi$umap1_pre, d_epi$umap2_pre, d_epi$umap1_post, d_epi$umap2_post)
pe <- d_epi[, .(n = .N,
                d_pre  = median(sqrt((umap1_pre  - ge$gx)^2  + (umap2_pre  - ge$gy)^2)),
                s_pre  = median(sqrt((umap1_pre  - median(umap1_pre))^2 +
                                     (umap2_pre  - median(umap2_pre))^2)),
                d_post = median(sqrt((umap1_post - ge$gxp)^2 + (umap2_post - ge$gyp)^2))),
            by = patient]
pe[, sep_pre := d_pre / ge$s_all]
pe[, sep_post := d_post / ge$s_allp]
pe[, clump := d_pre / s_pre]
print(pe[order(-sep_pre)][1:12])
cat(sprintf("epi patient sep_pre: min %.2f | med %.2f | max %.2f\n",
            min(pe$sep_pre), median(pe$sep_pre), max(pe$sep_pre)))

# ── 4) 上皮 × 单组织：纯批次斑块（无组织混淆）────────────
cat("\n== epi + tissue (each tissue separately), patient sep_pre top 8 ==\n")
for (t in tissues) {
  dt <- d_epi[tissue == t]
  gt <- glob(dt$umap1_pre, dt$umap2_pre, dt$umap1_post, dt$umap2_post)
  pt <- dt[, .(n = .N,
               d_pre  = median(sqrt((umap1_pre  - gt$gx)^2  + (umap2_pre  - gt$gy)^2)),
               d_post = median(sqrt((umap1_post - gt$gxp)^2 + (umap2_post - gt$gyp)^2))),
           by = patient]
  pt[, sep_pre := d_pre / gt$s_all]
  pt[, sep_post := d_post / gt$s_allp]
  cat(sprintf("\n-- tissue %s (n=%d) --\n", t, nrow(dt)))
  print(pt[order(-sep_pre)][1:8])
  cat(sprintf("  sep_pre: min %.2f | med %.2f | max %.2f ;  sep_post med %.2f\n",
              min(pt$sep_pre), median(pt$sep_pre), max(pt$sep_pre),
              median(pt$sep_post)))
}

# ── 5) A 侧 clump 排名（per patient 的 d_pre/s_pre）──────
cat("\n== epi + tissue A: clump = d_pre/s_pre, all patients ==\n")
dtA <- d_epi[tissue == tissues[1]]
gtA <- glob(dtA$umap1_pre, dtA$umap2_pre, dtA$umap1_post, dtA$umap2_post)
pA <- dtA[, .(n = .N,
              d_pre = median(sqrt((umap1_pre  - gtA$gx)^2 + (umap2_pre  - gtA$gy)^2)),
              s_pre = median(sqrt((umap1_pre  - median(umap1_pre))^2 +
                                  (umap2_pre  - median(umap2_pre))^2))),
          by = patient]
pA[, clump := d_pre / s_pre]
print(pA[order(-clump)])
