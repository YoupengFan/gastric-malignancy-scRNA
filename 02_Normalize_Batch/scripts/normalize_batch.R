#!/usr/bin/env Rscript
# ============================================================
# 02_Normalize_Batch: log1p(CP10K) + CC/MT 回归 + Harmony
#
# 输入：01_QC/data/gse150290_qc_matrix.rds + gse150290_qc_meta.csv
#       + gse150290_doublet_scores.csv（01_QC stage 2，Scrublet）
# 产出：02_Normalize_Batch/data/gse150290_norm_harmony.rds（Seurat
#       对象：counts/data/scale.data layers、50 PCs + harmony
#       reduction、HVG flag）+ normalize_batch_params.json +
#       normalize_batch.log
#
# 关键设计要点：
#   - 流程顺序：归一化 → CC 评分 → HVG 3000（黑名单后精选）→
#     回归 S/G2M/pct_mt（残差，不缩放）→ 残差方差重排 top 2000 →
#     z-score（clip 10）→ PCA 50 → Harmony（key=sample；theta=2,
#     sigma=0.1, nclust=100）
#   - scale.factor = 10000（CP10K，Seurat 社区默认）
#   - Scrublet 异型 doublet 剔除（threshold 0.25）：异型 doublet 是
#     偏离类分析的混淆项（混合身份 = 伪偏离轴）；同型 doublet
#     下游相对无害（归一化）。刻意不做大细胞/高 UMI 硬过滤
#     （会误杀 Parietal/Chief/恶性大细胞）
#   - 技术噪音黑名单：^MT- / ^RPL / ^RPS / ^MRPL / ^MRPS / ^MALAT1$
#     （vst 排序后剔除）；CC 基因不进黑名单
#   - 与论文对照：Zhang 2021 在 scale 步骤回归 total UMI + tSNE；
#     我们用 CC/MT 回归 + Harmony（时代改进）
#
# Usage:
#   export PATH="$HOME/softwares/bin:$PATH"
#   Rscript 02_Normalize_Batch/scripts/normalize_batch.R
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(harmony)
})

# ── Paths ─────────────────────────────────────────────────────
QC_DIR   <- "01_QC/data"
OUT_DIR  <- "02_Normalize_Batch/data"
TMP_DIR  <- "/tmp/norm_r"
OUT_TMP  <- file.path(TMP_DIR, "out")

# ── 参数 ──────────────────────────────────────────────────────
N_HVG_RANK    <- 5000  # vst 排序深度（黑名单剔除后仍足够取 N_HVG_INITIAL）
N_HVG_INITIAL <- 3000  # 回归前的初始宽选 HVG
N_HVG_FINAL   <- 2000  # 残差方差重排后的最终 HVG
N_PCS         <- 50
SCALE_CLIP    <- 10    # sc.pp.scale(max_value=10)
# harmony >= 2.0：max.iter.harmony 改名 max_iter；epsilon.harmony 已
# 移除（由固定阈值 1e-4 的 early_stop 替代）
HARMONY_PARAMS <- list(theta = 2, sigma = 0.1, nclust = 100,
                       max_iter = 20, early_stop = TRUE)
HARMONY_KEY    <- "sample"

# 技术噪音黑名单（MT- / 核糖体 / MALAT1）；CC 基因不进黑名单
BLACKLIST_PATTERNS <- c("^MT-", "^RPL", "^RPS", "^MRPL", "^MRPS", "^MALAT1$")

dir.create(OUT_TMP, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
log_fh <- file(file.path(TMP_DIR, "normalize_batch.log"), open = "a")
log <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n")
  writeLines(line, log_fh)
  flush(log_fh)
}

# ── 1. 读入 QC 产出，构建 Seurat 对象 ─────────────────────────
t0 <- Sys.time()
log(sprintf("Loading QC matrix from %s ...", QC_DIR))
counts <- readRDS(file.path(QC_DIR, "gse150290_qc_matrix.rds"))
# 91 个重复 symbol（U1/7SK/Y_RNA/SNORA40/ZNF 家族等）——make.unique
# 供 Seurat 使用；行是互异的 ENSEMBL 位点，仅名字撞车
rownames(counts) <- make.unique(rownames(counts))
meta <- fread(file.path(QC_DIR, "gse150290_qc_meta.csv"))
setnames(meta, "pct_mt", "percent.mt")
meta <- as.data.frame(meta)
rownames(meta) <- meta$barcode
meta <- meta[colnames(counts), ]
log(sprintf("  %d cells x %d genes loaded [%.0fs]",
            ncol(counts), nrow(counts), difftime(Sys.time(), t0, units = "secs")))

# ── 1b. Doublet 过滤：01_QC 的逐样本 Scrublet 分数 ────────────
# 异型 doublet 是偏离类分析的混淆项（混合身份 = 伪偏离轴）；
# 不做大细胞/高 UMI 硬过滤（会误杀真大单细胞：Parietal/Chief/恶性）
db_csv <- file.path(QC_DIR, "gse150290_doublet_scores.csv")
if (file.exists(db_csv)) {
  db <- fread(db_csv)
  is_dbl <- db$predicted_doublet == TRUE
  keep_bc <- setdiff(colnames(counts), db$barcode[is_dbl])
  n_scored <- nrow(db)
  n_dbl <- sum(is_dbl)
  n_kept <- length(keep_bc)
  log(sprintf("  Scrublet doublets excluded: %d / %d scored (%.2f%%); %d cells remain",
              n_dbl, n_scored, 100 * n_dbl / n_scored, n_kept))
  counts <- counts[, keep_bc, drop = FALSE]
  meta <- meta[keep_bc, , drop = FALSE]
} else {
  log("  WARNING: gse150290_doublet_scores.csv not found; doublet filter SKIPPED")
}

obj <- CreateSeuratObject(counts = counts, meta.data = meta, min.cells = 0, min.features = 0)
rm(counts); gc()

# ── 2. Normalize: log1p(x * 10K / lib) — CP10K (Seurat default) ──
t0 <- Sys.time()
lib <- colSums(obj[["RNA"]]$counts)
log(sprintf("Library size: median = %.0f, mean = %.0f (informational; scale.factor is fixed 10K)",
            median(lib), mean(lib)))
log("Normalizing with scale.factor = 10000 (CP10K, Seurat default) ...")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
log(sprintf("  Normalized [%.0fs]", difftime(Sys.time(), t0, units = "secs")))

# ── 3. Cell cycle scoring (Tirosh 2015, Seurat canonical) ─────
t0 <- Sys.time()
log("CellCycleScoring (cc.genes) ...")
obj <- CellCycleScoring(obj, s.features = cc.genes$s.genes,
                        g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
log(sprintf("  S.Score range [%.2f, %.2f], G2M.Score range [%.2f, %.2f] [%.0fs]",
            min(obj$S.Score), max(obj$S.Score),
            min(obj$G2M.Score), max(obj$G2M.Score),
            difftime(Sys.time(), t0, units = "secs")))

# ── 4. HVG: vst on all genes, then drop blacklist, keep 3000 ──
t0 <- Sys.time()
log(sprintf("FindVariableFeatures (vst, nfeatures=%d, keep %d after blacklist) ...",
            N_HVG_RANK, N_HVG_INITIAL))
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = N_HVG_RANK)
ranked <- VariableFeatures(obj)
is_bl <- grepl(paste(BLACKLIST_PATTERNS, collapse = "|"), ranked)
n_bl_in_rank <- sum(is_bl)
hvg3000 <- ranked[!is_bl][seq_len(N_HVG_INITIAL)]
VariableFeatures(obj) <- hvg3000
log(sprintf("  %d blacklist genes in top-5000 removed; HVG = %d [%.0fs]",
            n_bl_in_rank, length(hvg3000), difftime(Sys.time(), t0, units = "secs")))

# ── 5. Regress S/G2M/pct_mt on HVG 3000 (residuals, NOT scaled) ──
t0 <- Sys.time()
log("ScaleData: regress S.Score + G2M.Score + percent.mt (center only, no scaling) ...")
obj <- ScaleData(obj, features = hvg3000,
                 vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"),
                 do.scale = FALSE, do.center = TRUE)
log(sprintf("  Regression done [%.0fs]", difftime(Sys.time(), t0, units = "secs")))

# ── 6. Re-rank by residual variance -> top 2000 ───────────────
t0 <- Sys.time()
log(sprintf("Re-ranking by residual variance: %d -> %d ...", N_HVG_INITIAL, N_HVG_FINAL))
sd <- GetAssayData(obj, layer = "scale.data")
rm_ <- rowMeans(sd)
rv <- rowSums((sd - rm_)^2) / (ncol(sd) - 1)
log(sprintf("  Residual variance range [%.4f, %.4f], median %.4f",
            min(rv), max(rv), median(rv)))
top2000 <- names(sort(rv, decreasing = TRUE))[seq_len(N_HVG_FINAL)]
VariableFeatures(obj) <- top2000
# 报告与回归前 top-2000 的重叠
orig_top <- hvg3000[seq_len(N_HVG_FINAL)]
log(sprintf("  Overlap with pre-regression top-%d: %d/%d (%.1f%%)",
            N_HVG_FINAL, sum(top2000 %in% orig_top), N_HVG_FINAL,
            sum(top2000 %in% orig_top) / N_HVG_FINAL * 100))
rm(sd, rm_, rv); gc()

# ── 7. z-score final 2000 HVG residuals + clip 10 ─────────────
t0 <- Sys.time()
log("z-score final 2000 HVG + clip ±10 ...")
sd <- GetAssayData(obj, layer = "scale.data")[top2000, , drop = FALSE]
rm_ <- rowMeans(sd)
rsd <- sqrt(rowSums((sd - rm_)^2) / (ncol(sd) - 1))
zs <- sweep(sd, 1, rm_, "-")
rm(sd, rm_); gc()
zs <- sweep(zs, 1, rsd, "/")
zs[zs > SCALE_CLIP] <- SCALE_CLIP
zs[zs < -SCALE_CLIP] <- -SCALE_CLIP
obj <- SetAssayData(obj, layer = "scale.data", new.data = zs)
rm(zs); gc()
log(sprintf("  Scaled + clipped [%.0fs]", difftime(Sys.time(), t0, units = "secs")))

# ── 8. PCA on final HVG ───────────────────────────────────────
t0 <- Sys.time()
log(sprintf("RunPCA (features=%d, npcs=%d) ...", length(top2000), N_PCS))
obj <- RunPCA(obj, features = top2000, npcs = N_PCS, verbose = FALSE)
pca_stdev <- obj[["pca"]]@stdev   # SeuratObject 的 Stdev() 签名已变；直接读 slot
var_ratio <- pca_stdev[1:N_PCS]^2 / sum(pca_stdev[1:N_PCS]^2)
cumvar <- cumsum(var_ratio)
log(sprintf("  PC1: %.4f (%.2f%%); PC1-5 cum: %.1f%%; PC1-10 cum: %.1f%%; PC1-30 cum: %.1f%%",
            var_ratio[1], var_ratio[1] * 100, cumvar[5] * 100, cumvar[10] * 100, cumvar[30] * 100))
log(sprintf("  PCA done [%.0fs]", difftime(Sys.time(), t0, units = "secs")))

# ── 9. Harmony on 50 PCs, key = sample ────────────────────────
t0 <- Sys.time()
log(sprintf("RunHarmony (group.by.vars=%s, theta=%.1f, sigma=%.1f, nclust=%d) ...",
            HARMONY_KEY, HARMONY_PARAMS$theta, HARMONY_PARAMS$sigma, HARMONY_PARAMS$nclust))
obj <- RunHarmony(obj, group.by.vars = HARMONY_KEY, reduction.use = "pca",
                  reduction.save = "harmony",
                  theta = HARMONY_PARAMS$theta, sigma = HARMONY_PARAMS$sigma,
                  nclust = HARMONY_PARAMS$nclust,
                  max_iter = HARMONY_PARAMS$max_iter,
                  early_stop = HARMONY_PARAMS$early_stop)
log(sprintf("  Harmony done [%.0fs]", difftime(Sys.time(), t0, units = "secs")))

# ── 10. Save ──────────────────────────────────────────────────
t0 <- Sys.time()
log("Saving Seurat object (RDS) ...")
saveRDS(obj, file.path(OUT_TMP, "gse150290_norm_harmony.rds"), compress = FALSE)
log(sprintf("  Saved [%.0fs], size %.1f GB",
            difftime(Sys.time(), t0, units = "secs"),
            file.info(file.path(OUT_TMP, "gse150290_norm_harmony.rds"))$size / 1e9))

params <- list(
  step = "02_Normalize_Batch",
  normalization = list(method = "LogNormalize", scale_factor = 10000,
                       scale_factor_type = "CP10K (Seurat default)"),
  cell_cycle = list(method = "Seurat CellCycleScoring", gene_set = "Tirosh 2015 cc.genes"),
  hvg = list(initial_n = N_HVG_INITIAL, final_n = N_HVG_FINAL,
             selection = "vst", blacklist = BLACKLIST_PATTERNS,
             re_rank = "residual variance after regression"),
  regression = list(vars = c("S.Score", "G2M.Score", "percent.mt"),
                    output = "centered residuals (not scaled)"),
  scaling = list(method = "z-score", clip = SCALE_CLIP),
  pca = list(n_pcs = N_PCS, features = length(top2000)),
  harmony = HARMONY_PARAMS,
  paper_comparison = "Zhang 2021 regressed total UMI at scale step + tSNE; we use CC/MT regression + Harmony (era consensus)"
)
writeLines(toJSON(params, pretty = TRUE, auto_unbox = TRUE),
           file.path(OUT_TMP, "normalize_batch_params.json"))

log("Copying outputs to repo dir ...")
file.copy(file.path(OUT_TMP, "gse150290_norm_harmony.rds"),
          file.path(OUT_DIR, "gse150290_norm_harmony.rds"), overwrite = TRUE)
file.copy(file.path(OUT_TMP, "normalize_batch_params.json"),
          file.path(OUT_DIR, "normalize_batch_params.json"), overwrite = TRUE)
file.copy(file.path(TMP_DIR, "normalize_batch.log"),
          file.path(OUT_DIR, "normalize_batch.log"), overwrite = TRUE)
log(sprintf("Done. Output: %s", file.path(OUT_DIR, "gse150290_norm_harmony.rds")))
close(log_fh)
