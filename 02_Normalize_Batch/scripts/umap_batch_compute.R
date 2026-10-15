#!/usr/bin/env Rscript
# ============================================================
# 02_Normalize_Batch UMAP 批次可视化 — 计算部分（WSL）
#
# 用途：UMAP 矫正前后对比，可视化批次效应。
#   Before = pre-Harmony PCA 上的 UMAP（本脚本现算）
#   After  = Harmony embedding 上的 UMAP（03 的 RDS 中已有）
#   图（由 umap_batch_figure.R 在 Windows 原生 R 渲染）按患者
#   着色：矫正前患者成彩色斑块（批次效应）；矫正后颜色混合。
#
# 输入：03_Cluster_Annotate/data/gse150290_annotated.rds
# 产出：02_Normalize_Batch/data/umap_batch.csv
#       （barcode, sample, patient, tissue, annotation,
#        umap1_pre, umap2_pre, umap1_post, umap2_post）
#
# Usage (WSL):
#   export PATH="$HOME/softwares/bin:$PATH"
#   Rscript 02_Normalize_Batch/scripts/umap_batch_compute.R
#   # ~15 min：7.0 GB RDS 加载（~2-3 min）+ pca FindNeighbors +
#   #           RunUMAP on 160k cells x 50 dims（~5-10 min）
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(uwot)
})

BASE     <- "."
RDS_PATH <- file.path(BASE, "03_Cluster_Annotate/data/gse150290_annotated.rds")
OUT_DIR  <- file.path(BASE, "02_Normalize_Batch/data")
TMP_DIR  <- "/tmp/umap_batch_r"
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_PCS       <- 50
N_NEIGHBORS <- 15   # 与 03 聚类的 UMAP 一致

log_fh <- file(file.path(TMP_DIR, "umap_batch_compute.log"), open = "a")
log <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n")
  writeLines(line, log_fh)
  flush(log_fh)
}

# ── 读入 ─────────────────────────────────────────────────────
t0 <- Sys.time()
log("Loading annotated RDS (7.0 GB) ...")
obj <- readRDS(RDS_PATH)
log(sprintf("  %d cells x %d features [%.0fs]",
            ncol(obj), nrow(obj), difftime(Sys.time(), t0, units = "secs")))

# ── 矫正前 UMAP：pre-Harmony PCA ─────────────────────────────
# 直接用 uwot 在 50 维 pre-Harmony PCA 上计算（参数与 03 步的
# Seurat RunUMAP 相同：n_neighbors=15, min_dist=0.3）。Seurat 的
# RunUMAP 图路径需要 umap-learn 包，故绕过它直接调底层
# uwot 包。
t0 <- Sys.time()
log(sprintf("uwot::umap on pca 50-d (n_neighbors=%d, min_dist=0.3) ...",
            N_NEIGHBORS))
set.seed(42)
umap_pre <- umap(Embeddings(obj, "pca"),
                 n_neighbors = N_NEIGHBORS, min_dist = 0.3,
                 n_components = 2, verbose = TRUE)
colnames(umap_pre) <- c("UMAP_1", "UMAP_2")
log(sprintf("  umap_pre done [%.0fs]", difftime(Sys.time(), t0, units = "secs")))

# ── 导出两套 UMAP + meta ─────────────────────────────────────
t0 <- Sys.time()
log("Exporting umap_batch.csv ...")
umap_post <- Embeddings(obj, "umap")   # 基于 harmony，03 步已计算
out <- data.frame(
  barcode    = colnames(obj),
  sample     = obj$sample,
  patient    = obj$patient,
  tissue     = obj$tissue,
  annotation = obj$annotation,   # 逐细胞投票（03 v2.1 最终注释）
  umap1_pre  = round(umap_pre[, 1], 4),
  umap2_pre  = round(umap_pre[, 2], 4),
  umap1_post = round(umap_post[, 1], 4),
  umap2_post = round(umap_post[, 2], 4),
  stringsAsFactors = FALSE
)
fwrite(out, file.path(TMP_DIR, "umap_batch.csv"))
log(sprintf("  %d rows [%.0fs]", nrow(out), difftime(Sys.time(), t0, units = "secs")))

file.copy(file.path(TMP_DIR, "umap_batch.csv"),
          file.path(OUT_DIR, "umap_batch.csv"), overwrite = TRUE)
log("Copied outputs. Done.")
close(log_fh)
